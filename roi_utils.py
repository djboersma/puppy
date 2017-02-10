# coding: utf-8
import dicom
import SimpleITK as sitk
import numpy as np
import matplotlib
import logging
import pylab
logger = logging.getLogger()

def list_roinames(ds):
    assert(len(ds.ROIContourSequence)==len(ds.StructureSetROISequence))
    return [str(ssroi.ROIName) for ssroi in ds.StructureSetROISequence]

def sum_of_angles(points,name="unspecified"):
    """
    Code to determine left/right handedness of contour.
    """
    phi=0.
    npoints = len(points)
    assert(npoints>2)
    assert(points.shape[1]==3)
    assert(len(set(points[:,2]))==1)
    pppoints=np.append(points[:,:2],points[:2,:2],axis=0)
    dpppoints=np.diff(pppoints,axis=0)
    dp0=dpppoints[:-1]
    nonzerodp0=(0<np.sum(dp0**2,axis=1))
    if not nonzerodp0.all():
        Ngood = np.sum(nonzerodp0) 
        if Ngood < 2:
            logger.warn("got a pathological contour: only {} out of {} have nonzero line segments".format(Ngood,npoints))
            return np.nan
        else:
            logger.warn("buggy contour: only {} out of {} have nonzero line segments, going to re-call this function with cleaned point set".format(Ngood,npoints))
            return sum_of_angles(points[nonzerodp0],name)
    # if the previous works as intended, then the following assert should always pass
    assert(nonzerodp0.all())
    dp1=dpppoints[1:]
    kross = dp0[:,0]*dp1[:,1] - dp0[:,1]*dp1[:,0]
    dots  = dp0[:,0]*dp1[:,0] + dp0[:,1]*dp1[:,1]
    norms = np.sqrt(np.sum(dp0**2,axis=1)*np.sum(dp1**2,axis=1))
    assert((norms>0).all())
    sinphi = kross/norms
    sinphi[sinphi>1]=1
    sinphi[sinphi<-1]=-1
    maskQ23=(dots<0)
    maskQ2=maskQ23*(sinphi>0)
    maskQ3=maskQ23*(sinphi<0)
    maskBAD=maskQ23*(sinphi==0)
    phi=np.arcsin(sinphi)
    phi[maskQ23]*=-1
    phi[maskQ2]+=np.pi
    phi[maskQ3]-=np.pi
    if maskBAD.any():
        logger.warn("{} contains {} points where the contour retreats 180 degrees on itself".format(name,np.sum(maskBAD)))
        return np.nan
    roundphi = int(np.round(np.sum(phi)*180/np.pi));
    if abs(roundphi) != 360:
        logger.warn("({}) weird sum of contour angles: {}, should be + or - 360 degrees".format(name,roundphi))
    return roundphi

class contour_layer(object):
    """
    List of 2D contours. The ones with positive orientation (sum of angles
    between successive contour segments is +360 degrees) will be used to for
    inclusion, the ones with negative orientation (sum of angles is
    -360 degrees) will be used for exclusion. All points of an exclusion
    contour should be included by an inclusion contour.
    """
    def __init__(self,points=None,ref=None,name="notset"):
        self.name = name
        self.ref = ref
        self.z = points[0,2]
        self.inclusion = []
        self.exclusion = []
        self.add_contour(points,ref)
    def add_contour(self,points,ref=None):
        assert(len(points)>2)
        assert(self.z == points[0,2])
        if self.ref is None:
            self.ref = ref
        elif not ref is None:
            assert(ref==self.ref)
        orientation = sum_of_angles(points)
        path = matplotlib.path.Path(points[:,:2])
        if np.around(orientation) == 360:
            self.inclusion.append(path)
        elif np.around(orientation) == -360:
            self.exclusion.append(path)
        else:
            logger.error("({}) got a very weird contour a sum of angles equal to {}; z={} ref={}".format(self.name,orientation,len(points),self.z))
        logger.debug("layer {} has {} inclusion path(s) and {} exclusion path(s)".format(self.ref,len(self.inclusion),len(self.exclusion)))
    def contains_points(self,xycoords):
        Ncoords = len(xycoords)
        assert(xycoords.shape == (Ncoords,2))
        flatmask = np.zeros(len(xycoords),dtype=bool)
        for q in self.inclusion:
            flatmask |= q.contains_points(xycoords)
        for p in self.exclusion:
            flatmask &= np.logical_not(p.contains_points(xycoords))
        return flatmask
    def check(self):
        assert(len(self.inclusion)>0) # really?
        for p in self.exclusion:
            ok = False
            for q in self.inclusion:
                if q.contains_path(p):
                    ok = True
                    logger.debug("({}) layer {} exclusion contour check OK".format(self.name,self.z))
                break
            if not ok:
                logger.critical("({}) exclusion contour at z={} not contained in any inclusion contour".format(self.name,self.z))
                raise RuntimeError("contour error")
        logger.debug("({}) layer {} check OK".format(self.name,self.z))

class region_of_interest(object):
    def __init__(self,ds,roi_id,verbose=False):
        roi_found = False
        assert(len(ds.ROIContourSequence)==len(ds.StructureSetROISequence))
        for roi,ssroi in zip(ds.ROIContourSequence,ds.StructureSetROISequence):
            if (str(roi_id)==str(roi.RefdROINumber)) or (str(roi_id)==str(ssroi.ROIName)):
                self.roinr = int(roi.RefdROINumber)
                self.roiname = str(ssroi.ROIName)
                roi_found = True
                break
            else: # debug
                logger.debug("{} != {}".format(roi_id,roi.RefdROINumber))
        if not roi_found:
            raise ValueError("ROI with id {} not found".format(roi_id))
        self.ncontours = len(roi.ContourSequence)
        self.npoints_total = sum([len(c.ContourData) for c in roi.ContourSequence])
        self.xmin=np.inf
        self.ymin=np.inf
        self.zmin=np.inf
        self.xmax=-np.inf
        self.ymax=-np.inf
        self.zmax=-np.inf
        # we are sort the contours by depth-coordinate
        self.contour_layers=[]
        self.zlist = []
        self.dz = 0.
        #self.contour_refs=[]
        for contour in roi.ContourSequence:
            ref = contour.ContourImageSequence[0].RefdSOPInstanceUID
            npoints = int(contour.NumberOfContourPoints)
            # check assumption on number of contour coordinates
            assert(len(contour.ContourData)==3*npoints)
            points = np.array([float(coord) for coord in contour.ContourData]).reshape(npoints,3)
            zvalues = set(points[:,2])
            # check assumption that all points are in the same xy plane (constant z)
            assert(len(zvalues)==1)
            zvalue = zvalues.pop()
            if zvalue in self.zlist:
                ic = self.zlist.index(zvalue)
                self.contour_layers[ic].add_contour(points,ref)
            else:
                self.contour_layers.append(contour_layer(points,ref))
                self.zlist.append(zvalue)
            self.xmin = min(self.xmin,np.min(points[:,0]))
            self.ymin = min(self.ymin,np.min(points[:,1]))
            self.zmin = min(self.zmin,zvalue)
            self.xmax = max(self.xmax,np.max(points[:,0]))
            self.ymax = max(self.ymax,np.max(points[:,1]))
            self.zmax = max(self.zmax,zvalue)
        if verbose:
            logger.info("roi {}={} has {} points on {} contours with z range [{},{}]".format(
                    self.roinr,self.roiname,self.npoints_total,self.ncontours,self.zmin,self.zmax))
        for layer in self.contour_layers:
            layer.check()
        dz = set(np.diff(self.zlist))
        if len(dz) == 1:
            self.dz = dz.pop()
        else:
            dz = set(np.diff(np.around(self.zlist,decimals=6)))
            if len(dz) == 1:
                self.dz = dz.pop()
            else:
                logger.warn("{} not one single z step: {}".format(self.roiname,", ".join([str(d) for d in dz])))
                self.dz = 0.
    def have_mask(self):
        return self.dz != 0.
    def get_mask(self,img):
        if not self.have_mask():
            logger.warn("Irregular z-values, masking not yet supported")
            return None
        dims=img.GetSize()
        if len(dims)!=3:
            logger.error("ERROR only 3d images supported")
            return None
        #logger.debug("create roi mask image object with dims={}".format(dims))
        roimask = sitk.Image(dims,sitk.sitkUInt8)
        roimask.CopyInformation(img)
        orig = roimask.GetOrigin()
        space = roimask.GetSpacing()
        contained = True
        for o,s,d,rmin,rmax in zip(orig,space,dims,[self.xmin,self.ymin,self.zmin],[self.xmin,self.ymin,self.zmin]):
            contained &= (int(np.round(rmin-o)/s) in range(d))
            contained &= (int(np.round(rmax-o)/s) in range(d))
        if not contained:
            logger.warn('DUIZEND BOMMEN EN GRANATEN')
        else:
            logger.debug('YAY: roi "{}" is contained in image'.format(self.roiname))
        #logger.debug("copied infor orig={} spacing={}".format(orig,space))
        zmin = orig[2] - 0.5*space[2]
        zmax = orig[2] + (dims[2]-0.5)*space[2]
        #logger.debug("zmin={} zmax={}".format(zmin,zmax))
        xpoints=np.linspace(orig[0],orig[0]+space[0]*dims[0],dims[0],False)
        ypoints=np.linspace(orig[1],orig[1]+space[1]*dims[1],dims[1],False)
        xymesh = np.meshgrid(xpoints,ypoints)
        xyflat = np.array([(x,y) for x,y in zip(xymesh[0].flat,xymesh[1].flat)])
        eps=0.001*np.abs(self.dz)
        #logger.debug("got point mesh")
        if zmin-eps>self.zmax or zmax+eps<self.zmin:
            logger.warn("WARNING: no overlap in z ranges")
            return roimask
        clayer0 = self.contour_layers[0]
        #logger.debug contour0pts.shape
        z0=clayer0.z
        #logger.debug("z0={}".format(z0))
        #logger.debug("going to loop over z planes in image")
        for iz in range(dims[2]):
            z=orig[2]+space[2]*iz
            icz = int(np.round((z-z0)/self.dz))
            if icz>=0 and icz<len(self.contour_layers):
                #logger.debug("masking z={} with contour {}".format(z,icz))
                #path = matplotlib.path.Path(self.contours[icz][:,0:2])
                flatmask = self.contour_layers[icz].contains_points(xyflat)
                logger.debug("got {} points inside".format(np.sum(flatmask)))
                for iflat,b in enumerate(flatmask):
                    if not b:
                        continue
                    ix = iflat % dims[0]
                    iy = iflat / dims[0]
                    roimask[ix,iy,iz]=1
        return roimask
    def get_dvh(self,img,nbins=100,dmin=None,dmax=None,debuglabel=None):
        logger.debug("starting dvh calculation")
        dims=img.GetSize()
        if len(dims)!=3:
            logger.error("ERROR only 3d images supported")
            return None
        logger.debug("got size = {}".format(dims))
        aimg = sitk.GetArrayFromImage(img)
        logger.debug("got array with shape {}".format(list(aimg.shape)))
        if dmin is None:
            dmin=np.min(aimg)
        if dmax is None:
            dmax=np.max(aimg)
        logger.debug("dmin={} dmax={}".format(dmin,dmax))
        itkmask=self.get_mask(img)
        logger.debug("got mask with size {}".format(itkmask.GetSize()))
        amask=(sitk.GetArrayFromImage(itkmask)>0)
        if debuglabel:
            print("START debugging "+debuglabel)
            amask0=amask*(aimg==0.0)
            ix0,iy0,iz0 = np.where(amask0)
            logger.debug('indices for zero dose voxels have {}/{}/{} different values (out of {}/{}/{})'.format(
                len(set(ix0[:])), len(set(iy0[:])), len(set(iz0[:])),
                len(ix0), len(iy0), len(iz0)))
            assert((aimg[ix0,iy0,iz0]==0.0).all())
            hix0,xedges=np.histogram(ix0,bins=np.arange(-0.5,aimg.shape[0]-0.49,1.0))
            hiy0,yedges=np.histogram(iy0,bins=np.arange(-0.5,aimg.shape[1]-0.49,1.0))
            hiz0,zedges=np.histogram(iz0,bins=np.arange(-0.5,aimg.shape[2]-0.49,1.0))
            hixy0,foo,bar=np.histogram2d(ix0,iy0,bins=(xedges,yedges))
            hiyz0,foo,bar=np.histogram2d(iy0,iz0,bins=(yedges,zedges))
            hizx0,foo,bar=np.histogram2d(iz0,ix0,bins=(zedges,xedges))
            ahixy0 = np.array(hixy0).T
            ahiyz0 = np.array(hiyz0).T
            ahizx0 = np.array(hizx0).T
            oldfig=pylab.gcf()
            fig=pylab.figure(num="dvhdebug_"+debuglabel,figsize=[10,15])
            pylab.suptitle('distribution of indices of voxels in "{}" with dose==0'.format(self.roiname))
            pylab.subplot(321)
            pylab.bar(xedges[:-1],hix0,lw=0,color='b')
            pylab.xlim(xedges[0],xedges[-1])
            pylab.xlabel('X index')
            pylab.subplot(322)
            pylab.pcolormesh(yedges,zedges,ahiyz0)
            #pylab.xlim(yedges[0],yedges[-1])
            #pylab.ylim(zedges[0],zedges[-1])
            pylab.xlabel('Y index')
            pylab.ylabel('Z index')
            pylab.colorbar()
            pylab.subplot(323)
            pylab.bar(yedges[:-1],hiy0,lw=0,color='b')
            pylab.xlabel('Y index')
            pylab.xlim(yedges[0],yedges[-1])
            pylab.subplot(324)
            pylab.pcolormesh(zedges,xedges,ahizx0)
            #pylab.xlim(zedges[0],zedges[-1])
            #pylab.ylim(xedges[0],xedges[-1])
            pylab.xlabel('Z index')
            pylab.ylabel('X index')
            pylab.colorbar()
            pylab.subplot(325)
            pylab.bar(zedges[:-1],hiz0,lw=0,color='b')
            pylab.xlim(zedges[0],zedges[-1])
            pylab.xlabel('Z index')
            pylab.subplot(326)
            pylab.pcolormesh(xedges,yedges,ahixy0)
            #pylab.xlim(xedges[0],xedges[-1])
            #pylab.ylim(yedges[0],yedges[-1])
            pylab.xlabel('X index')
            pylab.ylabel('Y index')
            pylab.colorbar()
            pylab.savefig('dvhdebug_{}.pdf'.format(debuglabel))
            pylab.savefig('dvhdebug_{}.png'.format(debuglabel))
            del fig
            if oldfig:
                pylab.figure(oldfig.number)
            amaskneg=amask*(aimg<0.0)
            logger.debug("got mask with {} selected voxels, {} of which have zero dose, {} are negative".format(np.sum(amask),np.sum(amask0),np.sum(amaskneg)))
            print("END debugging "+debuglabel)
        dhist,dedges = np.histogram(aimg[amask],bins=nbins,range=(dmin,dmax))
        logger.debug("got histogram with {} edges for {} bins".format(len(dedges),nbins))
        adhist=np.array(dhist,dtype=float)
        adedges=np.array(dedges,dtype=float)
        dsum=0.5*np.sum(adhist*adedges[:-1]+adhist*adedges[1:])
        dhistsum=np.sum(adhist)
        amasksum=np.sum(amask)
        adchist=np.cumsum(adhist)
        logger.debug("dhistsum={} amasksum={} adchist[-1]={}".format(dhistsum,amasksum,adchist[-1]))
        assert(amasksum==dhistsum)
        assert(amasksum==adchist[-1])
        logger.debug("survived assert")
        if dhistsum>0:
            logger.debug("getting dvh")
            dvh=-1.0*adchist/dhistsum+1.0
            logger.debug("got dvh with dsum={} dvh[0]={} adchist[0]={}".format(dsum,dvh[0],adchist[0]))
        else:
            logger.warn("dhistsum is zero or negative")
        return dvh, dedges, dhistsum, dsum
                


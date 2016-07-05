import dicom
import SimpleITK as sitk
import numpy as np
import matplotlib
import logging
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
    for i in range(npoints):
        # we might avoid an explicit for loop if numpy would have a "cyclic array"
        j=(i+1)%npoints
        k=(i+2)%npoints
        ji=points[j,:]-points[i,:]
        kj=points[k,:]-points[j,:]
        assert(ji[2]==0)
        assert(kj[2]==0)
        krossproduct=ji[0]*kj[1]-ji[1]*kj[0]
        norms=np.sqrt(np.sum(ji**2)*np.sum(kj**2))
        sinphi=krossproduct/norms
        if sinphi<=-1:
            phi+=180.
            logger.warn("OOPS full turnaround? sinphi={} krossproduct={} norms={}".format(sinphi,krossproduct,norms))
        elif sinphi<1:
            phi+=np.arcsin(sinphi)*180./np.pi
    if int(np.abs(np.round(phi))) != 360:
        logger.warn("contour {} has an anomalous sum of segment angles: {}".format(name,phi))
    return phi

class contour_layer(object):
    """
    List of 2D contours. The ones with positive orientation (sum of angles
    between successive contour segments is +360 degrees) will be used to for
    inclusion, the ones with negative orientation (sum of angles is
    -360 degrees) will be used for exclusion. All points of an exclusion
    contour should be included by an inclusion contour.
    """
    def __init__(self,points,ref):
        self.ref = ref
        self.z = points[0,2]
        self.add_contour(points)
    def add_contour(self,points,ref=None):
        assert(self.z == points[0,2])
        if not ref is None:
             assert(ref==self.ref)
        orientation = sum_of_angles(points)

class region_of_interest(object):
    def __init__(self,ds,roi_id,verbose=False):
        ok=False
        assert(len(ds.ROIContourSequence)==len(ds.StructureSetROISequence))
        for roi,ssroi in zip(ds.ROIContourSequence,ds.StructureSetROISequence):
            if (str(roi_id)==str(roi.RefdROINumber)) or (str(roi_id)==str(ssroi.ROIName)):
                self.roinr = int(roi.RefdROINumber)
                self.roiname = str(ssroi.ROIName)
                ok = True
                break
            else: # debug
                logger.debug("{} != {}".format(roi_id,roi.RefdROINumber))
        if not ok:
            raise ValueError("ROI with id {} not found".format(roi_id))
        self.ncontours = len(roi.ContourSequence)
        self.npoints_total = sum([len(c.ContourData) for c in roi.ContourSequence])
        zvalues = [float(c.ContourData[-1]) for c in roi.ContourSequence]
        self.zmin=min(zvalues)
        self.zmax=max(zvalues)
        if verbose:
            logger.info("roi {}={} has {} points on {} contours with z range [{},{}]".format(
                    self.roinr,self.roiname,self.npoints_total,self.ncontours,self.zmin,self.zmax))
        # we are sort the contours by depth-coordinate
        self.contour_layers=[]
        self.dz = 0.
        self.contour_refs=[]
        for i,contour in enumerate(roi.ContourSequence):
            npoints = int(contour.NumberOfContourPoints)
            assert(len(contour.ContourData)==3*npoints)
            points = np.array([float(coord) for coord in contour.ContourData]).reshape(npoints,3)
            zvalues = set(points[:,2])
            assert(len(zvalues)==1)
            self.contours.append(points)
            self.contour_refs.append(contour.ContourImageSequence[0].RefdSOPInstanceUID)
            if i==1:
                self.dz = points[0,2] - self.contours[0][0,2]
                if abs(self.dz)<1e-6:
                    self.dz = 0.
            elif self.dz>0. or self.dz<0.:
                checkdz = points[0,2] - self.contours[i-1][0,2]
                if abs(checkdz-self.dz)>1e-6:
                    logger.error("ERROR {} self.dz={} checkdz={} points[0,2]={} previous contourz={}".format(
                                self.roiname,self.dz,checkdz,points[0,2],self.contours[i-1][0,2]))
                    self.dz = 0.
    def have_mask(self):
        return self.dz != 0.
    def get_mask(self,img):
        if self.dz == 0.:
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
        #logger.debug("copied infor orig={} spacing={}".format(orig,space))
        zmin = orig[2] - 0.5*space[2]
        zmax = orig[2] + (dims[2]-0.5)*space[2]
        #logger.debug("zmin={} zmax={}".format(zmin,zmax))
        xpoints=np.linspace(orig[0],orig[0]+space[0]*dims[0],dims[0],False)
        ypoints=np.linspace(orig[1],orig[1]+space[1]*dims[1],dims[1],False)
        xymesh = np.meshgrid(xpoints,ypoints)
        eps=0.001*np.abs(self.dz)
        #logger.debug("got point mesh")
        if zmin-eps>self.zmax or zmax+eps<self.zmin:
            logger.warn("WARNING: no overlap in z ranges")
            return roimask
        contour0pts = self.contours[0]
        #logger.debug contour0pts.shape
        z0=contour0pts[0,2]
        #logger.debug("z0={}".format(z0))
        #logger.debug("going to loop over z planes in image")
        for iz in range(dims[2]):
            z=orig[2]+space[2]*iz
            icz = int(np.round((z-z0)/self.dz))
            if icz>=0 and icz<self.ncontours:
                #logger.debug("masking z={} with contour {}".format(z,icz))
                path = matplotlib.path.Path(self.contours[icz][:,0:2])
                flatmask = path.contains_points(np.array([(x,y) for x,y in zip(xymesh[0].flat,xymesh[1].flat)]))
                #logger.debug("got {} points inside".format(np.sum(flatmask)))
                for iflat,b in enumerate(flatmask):
                    if not b:
                        continue
                    ix = iflat % dims[0]
                    iy = iflat / dims[0]
                    roimask[ix,iy,iz]=1
        return roimask
    def get_dvh(self,img,nbins=100,dmin=None,dmax=None):
        logger.debug("starting dvh calculation")
        dims=img.GetSize()
        if len(dims)!=3:
            logger.error("ERROR only 3d images supported")
            return None
        logger.debug("got size = {}".format(dims))
        aimg = sitk.GetArrayFromImage(img)
        logger.debug("got array")
        if dmin is None:
            dmin=np.min(aimg)
        if dmax is None:
            dmax=np.max(aimg)
        logger.debug("dmin={} dmax={}".format(dmin,dmax))
        itkmask=self.get_mask(img)
        logger.debug("got mask with size {}".format(itkmask.GetSize()))
        amask=(sitk.GetArrayFromImage(itkmask)>0)
        logger.debug("got mask with {} selected voxels".format(np.sum(amask)))
        dhist,dedges = np.histogram(aimg[amask],bins=nbins,range=(dmin,dmax))
        logger.debug("got histogram with {} edges for {} bins".format(len(dedges),nbins))
        adhist=np.array(dhist,dtype=float)
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
            logger.debug("got dvh")
        else:
            logger.warn("dhistsum is zero or negative")
        return dvh, dedges, dhistsum
                


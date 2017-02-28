# coding: utf-8
import dicom
import SimpleITK as sitk
import numpy as np
import matplotlib # for useful Path class, not for plotting...
from debugging_crap import make_elaborate_debugging_plots
from bounding_box import bounding_box
import logging
logger = logging.getLogger()

def list_roinames(ds):
    """
    Return the names of the ROIs in a given dicom structure set as a list of strings.
    """
    assert(hasattr(ds,"StructureSetROISequence"))
    assert(hasattr(ds,"ROIContourSequence"))
    assert(len(ds.ROIContourSequence)==len(ds.StructureSetROISequence))
    return [str(ssroi.ROIName) for ssroi in ds.StructureSetROISequence]

def scrutinize_contour(points):
    """
    Test that `points` is a (n,3) array suitable for contour definition.
    """
    assert(hasattr(points,"shape"))  # is it an array?
    assert(len(points.shape)==2)     # a 2-dim array, I mean?
    assert(points.shape[0]>=3)       # we need at least 3 points for a contour
    assert(points.shape[1]==3)       # we deal with points in 3-d space ...
    assert(len(set(points[:,2]))==1) # ... but we assume they all have the same z-coordinate
    # maybe add some more tests

def sum_of_angles(points,name="unspecified",rounded=True,scrutinize=False):
    """
    Code to determine left/right handedness of contour.  We do this by
    computing the angles between each successive segment in the contour. By
    "segment" we mean the difference vector between successive contour points.
    The cross and dot products between two segments are proportional to sine
    and cosine of the angle between the segments, respectively.
    """
    phi=0.
    # first: check assumptions
    if scrutinize:
        scrutinize_contour(points)
    # we have more assumptions, but will check them later
    npoints = points.shape[0] 
    # for computation convenience we add the first segment
    # (first two points) to the end of the list of points
    pppoints=np.append(points[:,:2],points[:2,:2],axis=0)
    # compute difference vectors (segments)
    dpppoints=np.diff(pppoints,axis=0)
    assert(dpppoints.shape == (npoints+1,2))
    # array with all unique seqments
    dp0=dpppoints[:-1]
    # array with all unique seqments, first segment moved to the end
    dp1=dpppoints[1:]
    assert(dp0.shape == (npoints,2))
    # check that segments have nonzero length
    nonzerodp0=(0<np.sum(dp0**2,axis=1))
    if not nonzerodp0.all():
        Ngood = np.sum(nonzerodp0) 
        if Ngood < 3:
            logger.warn("got a pathological contour: only {} out of {} have nonzero line segments".format(Ngood,npoints))
            return np.nan
        else:
            logger.warn("BUGGY CONTOUR: only {} out of {} have nonzero line segments, going to re-call this function with cleaned point set".format(Ngood,npoints))
            # by applying the nonzero mask on the points vector,
            # we leave out the points that are identical to their next neighbor
            # TODO: check against too deep recursion level (not a big concern here)
            return sum_of_angles(points[nonzerodp0],name=name,rounded=rounded,scrutinize=True)
    # if the previous works as intended, then the following assert should always pass
    assert(nonzerodp0.all())
    # now do ordinary vector calculus: cross product, dot product and norms
    kross = dp0[:,0]*dp1[:,1] - dp0[:,1]*dp1[:,0]
    dots  = dp0[:,0]*dp1[:,0] + dp0[:,1]*dp1[:,1]
    norms = np.sqrt(np.sum(dp0**2,axis=1)*np.sum(dp1**2,axis=1))
    # this assert is maybe paranoid and superfluous
    assert((norms>0).all())
    sinphi = kross/norms
    # guard against anomalies due to rounding errors
    sinphi[sinphi>1]=1
    sinphi[sinphi<-1]=-1
    # which Quadrant are we in?
    # Q1: less or equal pi/2 to the left
    # Q2: more than pi/2 to the left
    # Q3: more than pi/2 to the right
    # Q4: less or equal pi/2 to the right
    maskQ23=(dots<0)
    maskQ2=maskQ23*(sinphi>0)
    maskQ3=maskQ23*(sinphi<0)
    maskBAD=maskQ23*(sinphi==0)
    phi=np.arcsin(sinphi)
    # arcsin returns phi in range -pi .. +pi
    # Q2 (phi>0): phi -> +pi - phi
    # Q3 (phi<0): phi -> -pi - phi
    phi[maskQ23]*=-1
    phi[maskQ2]+=np.pi
    phi[maskQ3]-=np.pi
    if maskBAD.any():
        logger.warn("{} contains {} points where the contour retreats 180 degrees on itself".format(name,np.sum(maskBAD)))
        logger.warn("this is fixable (remove one or two points) but I did not implement that fix yet.")
        # TODO: is a warning and returning NAN sufficient? Shouldn't we crash and burn here?
        return np.nan
    sum_phi_deg = np.sum(phi)*180/np.pi
    round_sum_phi_deg = int(np.round(sum_phi_deg));
    if round_sum_phi_deg == 360:
        logger.debug("({}) POSITIVE: inclusion contour".format(name))
    elif round_sum_phi_deg == -360:
        logger.debug("({}) NEGATIVE: exclusion contour".format(name))
    else:
        logger.warn("({}) weird sum of contour angles: {} degrees, should be + or - 360 degrees".format(name,sum_phi_deg))
    if rounded:
        return round_sum_phi_deg
    else:
        return sum_phi_deg

class contour_layer(object):
    """
    This is an auxiliary class for the `region_of_interest` class defined below.
    A `contour_layer` object describes the ROI at one particular z-value.
    It is basically a list of 2D contours. The ones with positive orientation
    (sum of angles between successive contour segments is +360 degrees) will be
    used to for inclusion, the ones with negative orientation (sum of angles is
    -360 degrees) will be used for exclusion. All points of an exclusion
    contour should be included by an inclusion contour.
    """
    def __init__(self,points=None,ref=None,name="notset"):
        self.name = name
        self.ref = ref
        self.inclusion = []
        self.exclusion = []
        if points is None:
            self.z = None
        else:
            self.z = points[0,2]
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
    def contains_point(self,point):
        assert(len(point) == 2)
        is_contained = False
        for q in self.inclusion:
            if q.contains_point(point):
                is_contained = True
                break
        for p in self.exclusion:
            if p.contains_point(point):
                is_contained = False
                break
        return is_contained
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
            logger.error("ROI with id {} not found; structure set contains: ".format(roi_id) + ", ".join(list_roinames(ds)))
            raise ValueError("ROI with id {} not found".format(roi_id))
        self.ncontours = len(roi.ContourSequence)
        self.npoints_total = sum([len(c.ContourData) for c in roi.ContourSequence])
        self.bb = bounding_box()
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
            self.bb.should_contain_all(points)
        if verbose:
            logger.info("roi {}={} has {} points on {} contours with z range [{},{}]".format(
                    self.roinr,self.roiname,self.npoints_total,self.ncontours,self.bb.zmin,self.bb.zmax))
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
    def __repr__(self):
        return "roi {} defined by contours in {} layers, {}".format(self.roiname,len(self.contour_layers), self.bb)
    def have_mask(self):
        return self.dz != 0.
    def get_mask(self,img,zrange=None):
        """
        For a given image, compute for every voxel whether it is inside the ROI or not.
        The `zrange` can be used to limit the z-range of the ROI.
        If specified, the `zrange` should be contained in the z-range of the given image.
        """
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
        #############################################################################################
        # check that the bounding box of this ROI is contained within the volume of the given image #
        #############################################################################################
        contained = True
        for o,s,d,rmin,rmax in zip(orig,space,dims,self.bb.mincorner(),self.bb.maxcorner()):
            contained &= (int(np.round(rmin-o)/s) in range(d))
            contained &= (int(np.round(rmax-o)/s) in range(d))
        if not contained:
            logger.warn('DUIZEND BOMMEN EN GRANATEN orig={} space={} dims={} bbroi={}'.format(orig,space,dims,self.bb))
        else:
            logger.debug('YAY: roi "{}" is contained in image'.format(self.roiname))
        #logger.debug("copied infor orig={} spacing={}".format(orig,space))
        # ITK: the "origin" has the coordinates of the *center* of the corner voxel
        # zmin and zmax are the z coordinates of the boundary of the volume
        zmin = orig[2] - 0.5*space[2]
        zmax = orig[2] + (dims[2]-0.5)*space[2]
        eps=0.001*np.abs(self.dz)
        #logger.debug("got point mesh")
        if zmin-eps>self.bb.zmax+self.dz or zmax+eps<self.bb.zmin-self.dz:
            logger.warn("WARNING: no overlap in z ranges")
            logger.warn("WARNING: img z range [{}-{}], roi z range [{}-{}]".format(zmin,zmax,self.zmin,self.zmax))
            return roimask
        if zrange is None:
            zrange=(zmin,zmax)
        else:
            assert(len(zrange)==2)
            assert(zrange[0]>=zmin)
            assert(zrange[1]<=zmax)
            if zrange[0]-eps>self.bb.zmax+self.dz or zrange[1]+eps<self.bb.zmin-self.dz:
                logger.warn("WARNING: no overlap in (restricted) z ranges")
                return roimask
        # logger.debug("zmin={} zmax={}".format(zmin,zmax))
        # xpoints and ypoints contain the x/y coordinates of the voxel centers
        xpoints=np.linspace(orig[0],orig[0]+space[0]*dims[0],dims[0],False)
        ypoints=np.linspace(orig[1],orig[1]+space[1]*dims[1],dims[1],False)
        xymesh = np.meshgrid(xpoints,ypoints)
        xyflat = np.array([(x,y) for x,y in zip(xymesh[0].flat,xymesh[1].flat)])
        clayer0 = self.contour_layers[0]
        #logger.debug contour0pts.shape
        z0 = clayer0.z
        #logger.debug("z0={}".format(z0))
        #logger.debug("going to loop over z planes in image")
        for iz in range(dims[2]):
            z = orig[2]+space[2]*iz # z coordinate in image/mask
            if z<zrange[0] or z>zrange[1]:
                continue
            icz = int(np.round((z-z0)/self.dz)) # layer index
            if icz>=0 and icz<len(self.contour_layers):
                logger.debug("INSIDE roi: z index mask/image iz={} (z={}) layer index icz={} (z={})".format(iz,z,icz,self.contour_layers[icz].z))
                flatmask = self.contour_layers[icz].contains_points(xyflat)
                logger.debug("got {} points inside".format(np.sum(flatmask)))
                for iflat,b in enumerate(flatmask):
                    if not b:
                        continue
                    ix = iflat % dims[0]
                    iy = iflat / dims[0]
                    roimask[ix,iy,iz]=1
                    x = orig[0]+space[0]*ix # x coordinate in image/mask
                    y = orig[1]+space[1]*iy # y coordinate in image/mask
                    assert(self.contour_layers[icz].contains_point(point=(x,y)))
            elif icz<0:
                logger.debug("BELOW roi: z index mask/image iz={} (z={}) layer index icz={} (z0={} dz={})".format(iz,z,icz,z0,self.dz))
            else:
                logger.debug("ABOVE roi: z index mask/image iz={} (z={}) layer index icz={} (z0={} dz={} nlayer={})".format(iz,z,icz,z0,self.dz,len(self.contour_layers)))
        return roimask
    def get_dvh(self,img,nbins=100,dmin=None,dmax=None,zrange=None,debuglabel=None):
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
        itkmask=self.get_mask(img,zrange)
        logger.debug("got mask with size {}".format(itkmask.GetSize()))
        amask=(sitk.GetArrayFromImage(itkmask)>0)
        if debuglabel:
            make_elaborate_debugging_plots(aimg,amask,debuglabel,self.roiname)
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

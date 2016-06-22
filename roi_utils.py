import dicom
import SimpleITK as sitk
import numpy as np

def list_roinames(ds):
    assert(len(ds.ROIContourSequence)==len(ds.StructureSetROISequence))
    return [str(ssroi.ROIName) for ssroi in ds.StructureSetROISequence]

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
            #else: # debug
            #    print("{} != {}".format(roi_id,roi.RefdROINumber))
        if not ok:
            raise ValueError("ROI with id {} not found".format(roi_id))
        self.ncontours = len(roi.ContourSequence)
        self.npoints_total = sum([len(c.ContourData) for c in roi.ContourSequence])
        zvalues = [float(c.ContourData[-1]) for c in roi.ContourSequence]
        self.zmin=min(zvalues)
        self.zmax=max(zvalues)
        if verbose:
            print("roi {}={} has {} points on {} contours with z range [{},{}]".format(
                    self.roinr,self.roiname,self.npoints_total,self.ncontours,self.zmin,self.zmax))
        self.contours=[]
        self.dz = 0.
        self.contour_refs=[]
        for i,contour in enumerate(roi.ContourSequence):
            npoints = int(contour.NumberOfContourPoints)
            assert(len(contour.ContourData)==3*npoints)
            points = numpy.array([float(coord) for coord in contour.ContourData]).reshape(npoints,3)
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
                    self.dz = 0.
    def get_mask(self,img):
        if self.dz == 0.:
            print("Irregular z-values, masking not yet supported")
            return None
        dims=img.GetSize()
        if len(dims)!=3:
            print("ERROR only 3d images supported")
            return None
        roimask = sitk.Image(dims,sitk.sitkUInt8)
        roimask.CopyInformation(img)
        orig = roimask.GetOrigin()
        space = roimask.GetSpacing()
        zmin = orig[2] - 0.5000001*space[2]
        zmax = orig[2] + (dims[2]-0.50000001)*space[2]
        xpoints=np.linspace(orig[0],orig[0]+space[0]*dims[0],dims[0],False)
        ypoints=np.linspace(orig[1],orig[1]+space[1]*dims[1],dims[1],False)
        xymesh = np.meshgrid(xpoints,ypoints)
        if zmin>self.zmax or zmax<self.zmin:
            print("WARNING: no overlap in z ranges")
            return roimask
        contour0pts = self.contours[0]
        #print contour0pts.shape
        z0=contour0pts[0,2]
        #print("z0={}".format(z0))
        for iz in range(dims[2]):
            z=orig[2]+space[2]*iz
            icz = int(np.round((z-z0)/self.dz))
            if icz>=0 and icz<self.ncontours:
                path = matplotlib.path.Path(self.contours[icz][:,0:2])
                flatmask = path.contains_points(np.array([(x,y) for x,y in zip(xymesh[0].flat,xymesh[1].flat)]))
                for iflat,b in enumerate(flatmask):
                    if not b:
                        continue
                    ix = iflat % dims[0]
                    iy = iflat / dims[0]
                    roimask[ix,iy,iz]=1
        return roimask
    def get_dvh(self,img,nbins=100,dmin=None,dmax=None):
        dims=img.GetSize()
        if len(dims)!=3:
            print("ERROR only 3d images supported")
            return None
        aimg = sitk.GetArrayFromImage(img)
        if dmin is None:
            dmin=np.min(aimg)
        if dmax is None:
            dmax=np.max(aimg)
        amask=(sitk.GetArrayFromImage(self.get_mask(img))>0)
        dhist,dedges = np.histogram(aimg[amask],bins=nbins,range=(dmin,dmax))
        adhist=np.array(dhist,dtype=float)
        dhistsum=np.sum(adhist)
        amasksum=np.sum(amask)
        adchist=np.cumsum(adhist)
        assert(amasksum==dhistsum)
        assert(amasksum==adchist[-1])
        if dhistsum>0:
            dvh=-1.0*adchist/dhistsum+1.0
        return dvh, dedges, dhistsum
                

    # CODE TO DETERMINE LEFT/RIGHT HANDEDNESS OF CONTOUR
    #        phi=0.
    #        for i in range(npoints):
    #            j=(i+1)%npoints
    #            k=(i+2)%npoints
    #            ji=points[j,:]-points[i,:]
    #            kj=points[k,:]-points[j,:]
    #            assert(ji[2]==0)
    #            assert(kj[2]==0)
    #            krossproduct=ji[0]*kj[1]-ji[1]*kj[0]
    #            norms=np.sqrt(np.sum(ji**2)*np.sum(kj**2))
    #            sinphi=krossproduct/norms
    #            if sinphi<-1:
    #                phi+=180.
    #                print("OOPS full turnaround? sinphi={} krossproduct={} norms={}".format(sinphi,krossproduct,norms))
    #            elif sinphi<1:
    #                phi+=np.arcsin(sinphi)*180./np.pi
    #        print("contour {} point 1=({}) point {}=({}) phi = {} degrees".format(i,points[0,:],npoints,points[-1,:],phi))

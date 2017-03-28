import numpy as np

class bounding_box(object):
    """
    Define ranges in which things are in 3D space.
    We could probably get rid of the repetitive code
    by using two 3d arrays or even a 6d one, but keeping
    xyz lingo is more geometrically intuitive.
    """
    def __init__(self,**kwargs):
        nkeys = len(kwargs.keys())
        if nkeys == 0:
            self.xmin=np.inf
            self.ymin=np.inf
            self.zmin=np.inf
            self.xmax=-np.inf
            self.ymax=-np.inf
            self.zmax=-np.inf
        elif nkeys > 1:
            raise RuntimeError("too many arguments ({}) to bounding box constructor: {}".format(nkeys,kwargs))
        elif "bb" in kwargs:
            bb = kwargs["bb"]
            self.xmin=bb.xmin
            self.ymin=bb.ymin
            self.zmin=bb.zmin
            self.xmax=bb.xmax
            self.ymax=bb.ymax
            self.zmax=bb.zmax
        elif "xyz" in kwargs:
            xyz = kwargs["xyz"]
            axyz = np.array(xyz,dtype=float)
            if not (axyz.shape == (3,2) or axyz.shape==(6,)):
                raise ValueError("wrong shape for xyz limits: {}".format(axyz.shape))
            xmin,xmax,ymin,ymax,zmin,zmax = axyz.flat[:]
            if not (xmin<=xmax):
                raise ValueError("xmin should be less than xmax, got xmin={} xmax={}".format(xmin,xmax)))
            if not (ymin<=ymax):
                raise ValueError("ymin should be less than ymax, got ymin={} ymax={}".format(ymin,ymax)))
            if not (zmin<=zmax):
                raise ValueError("zmin should be less than zmax, got zmin={} zmax={}".format(zmin,zmax)))
            self.xmin=float(xmin)
            self.ymin=float(ymin)
            self.zmin=float(zmin)
            self.xmax=float(xmax)
            self.ymax=float(ymax)
            self.zmax=float(zmax)
    def __repr__(self):
        return "bounding box [[{},{}],[{},{}],[{},{}]]".format(
                self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)
    def volume(self):
        return (self.xmax-self.xmin)*(self.ymax-self.ymin)*(self.zmax-self.zmin)
    def empty(self):
        vol = self.volume()
        return ((vol==0) or np.isinf(vol))
    def __eq__(self,rhs):
        lhsvol = self.volume()
        rhsvol = rhs.volume()
        if lhsvol != rhsvol:
            return False
        if lhsvol==0 and rhsvol==0:
            return True
        return ((self.mincorner()==rhs.mincorner()).all() and (self.maxcorner()==rhs.maxcorner()).all())
    def should_contain(self,point):
        self.xmin = min(self.xmin,np.min(point[0]))
        self.ymin = min(self.ymin,np.min(point[1]))
        self.zmin = min(self.zmin,np.min(point[2]))
        self.xmax = max(self.xmax,np.max(point[0]))
        self.ymax = max(self.ymax,np.max(point[1]))
        self.zmax = max(self.zmax,np.max(point[2]))
    def should_contain_all(self,points):
        self.xmin = min(self.xmin,np.min(points[:,0]))
        self.ymin = min(self.ymin,np.min(points[:,1]))
        self.zmin = min(self.zmin,np.min(points[:,2]))
        self.xmax = max(self.xmax,np.max(points[:,0]))
        self.ymax = max(self.ymax,np.max(points[:,1]))
        self.zmax = max(self.zmax,np.max(points[:,2]))
    def mincorner(self):
        return np.array([self.xmin,self.ymin,self.zmin])
    def maxcorner(self):
        return np.array([self.xmax,self.ymax,self.zmax])
    def contains(self,point,inner=False):
        assert(len(point>2))
        if inner:
            return ( point[0]>self.xmin and point[0]<self.xmax and
                     point[1]>self.ymin and point[1]<self.ymax and
                     point[2]>self.zmin and point[2]<self.zmax )
        else:
            return ( point[0]>=self.xmin and point[0]<=self.xmax and
                     point[1]>=self.ymin and point[1]<=self.ymax and
                     point[2]>=self.zmin and point[2]<=self.zmax )
    def encloses(self,bb,inner=False):
        return (self.contains(bb.mincorner(),inner) and self.contains(bb.maxcorner(),inner))

#######################################################################
# TESTING
#######################################################################
import unittest

class test_bounding_box(unittest.TestCase):
    def test_xyz_constructor(self):
        bb0 = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        self.assertEqual( bb0.xmin, 1 )
        self.assertEqual( bb0.xmax, 2 )
        self.assertEqual( bb0.ymin, 3 )
        self.assertEqual( bb0.ymax, 4 )
        self.assertEqual( bb0.zmin, 5 )
        self.assertEqual( bb0.zmax, 6 )
        bb1 = bounding_box(xyz=[1,2,3,4,5,6])
        self.assertEqual( bb1.xmin, 1 )
        self.assertEqual( bb1.xmax, 2 )
        self.assertEqual( bb1.ymin, 3 )
        self.assertEqual( bb1.ymax, 4 )
        self.assertEqual( bb1.zmin, 5 )
        self.assertEqual( bb1.zmax, 6 )
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[[1,2],[4,3],[5,6]])
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[[np.nan,2],[4,3],[5,6]])
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[[1,2],[3],[5,6,7]])
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[[1,2],[3],[5,6]])
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[1,2,3,4,5,6,7])
        with self.assertRaises(ValueError):
            bbxyz_wrong = bounding_box(xyz=[1,2,3,4,5])
    def test_max_one_key_constructor(self):
        with self.assertRaises(RuntimeError):
            bb0 = bounding_box()
            bbxyz = bounding_box(xyz=[[1,2],[3,4],[5,6]],bb=bb0)
    def test_bb_constructor(self):
        bbxyz = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        bbxyz2 = bounding_box(bb=bbxyz)
    def test_default_constructor(self):
        bb0 = bounding_box()
        for b in [bb.xmin,bb.ymin,bb.zmin]:
            self.assertTrue( np.isinf(b) and b>0)
        for b in [bb.xmax,bb.ymax,bb.zmax]:
            self.assertTrue( np.isinf(b) and b<0)
    def test_equal(self):
        bb0 = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        bb1 = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        bb2 = bounding_box(xyz=[[1,2],[3,5],[5,6]])
        assertEqual(bb0,bb1)
        assertFalse(bb0,bb2)
        bb0 = bounding_box()
        bb1 = bounding_box()
        assertEqual(bb0,bb1)
        assertFalse(bb0,bb2)
        bb0 = bounding_box()
        bb1 = bounding_box(bb0)
        assertEqual(bb0,bb1)
        assertFalse(bb1,bb2)
        bb1 = bounding_box(xyz=[[1,2],[3,4],[5,5]])
        assertEqual(bb0,bb1) # both bounding boxes are empty
        assertFalse(bb1,bb2)
    def test_corners(self):
        bbxyz = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        assertEqual(bbxyz.mincorner(), np.array([1.,3.,5.]))
        assertEqual(bbxyz.maxcorner(), np.array([2.,4.,6.]))
    def test_contains(self):
        bbxyz = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        xx,yy,zz = np.meshgrid(np.arange(0.5,2.6,1.0), np.arange(2.5,4.6,1.0), np.arange(4.5,6.6,1.0) )
        for i,point in enumerate(zip(xx.flat,yy.flat,zz.flat)):
            if i==13:
                assertTrue(bbxyz.contains(point))
            else:
                assertFalse(bbxyz.contains(point))
    def test_grow(self):
        bbxyz = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        for point in np.array(range(30)).reshape(10,3):
            bbxyz.should_contain(point)
            assertTrue(bbxyz.contains(point))
        assertFalse(bbxyz.contains((0,0,0))
        assertFalse(bbxyz.contains((0,0,30))
        assertFalse(bbxyz.contains((-1,10,10))
        bbxyz2 = bounding_box(xyz=[[1,2],[3,4],[5,6]])
        bbxyz2.should_contain_all( np.array(range(30)).reshape(10,3) )
        assertEqual(bbxyz,bbxyz2)

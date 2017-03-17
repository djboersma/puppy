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
    #    assert(xmin<=xmax)
    #    assert(ymin<=ymax)
    #    assert(zmin<=zmax)
    #    self.xmin=float(xmin)
    #    self.ymin=float(ymin)
    #    self.zmin=float(zmin)
    #    self.xmax=float(xmax)
    #    self.ymax=float(ymax)
    #    self.zmax=float(zmax)
    def __repr__(self):
        return "bounding box [[{},{}],[{},{}],[{},{}]]".format(
                self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)
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
    def contains(self,point):
        assert(len(point>2))
        return ( point[0]>=self.xmin and point[0]<=self.xmax and
                 point[1]>=self.ymin and point[1]<=self.ymax and
                 point[2]>=self.zmin and point[2]<=self.zmax )
    def encloses(self,bb):
        return ( bb.xmin>=self.xmin and bb.xmax<=self.xmax and
                 bb.ymin>=self.ymin and bb.ymax<=self.ymax and
                 bb.zmin>=self.zmin and bb.zmax<=self.zmax )


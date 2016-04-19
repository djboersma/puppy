import numpy as np

class three_point_parabola(object):
    """
    This utility class can be used to quickly get
    a better approximation of the location of an extremum
    in a series of (x[i],y[i]) points (sorted such that x[i] are
    monotonically increasing). The idea is that you find
    a point x[i],y[i] such that y[i-1] and y[i+1] are both
    either higher or lower than y[i], 
    """
    def __init__(self,x3,y3):
        if len(x3) != 3:
            raise ValueError("x3 arg should have length 3")
        if len(y3) != 3:
            raise ValueError("x3 arg should have length 3")
        if x3[0]>=x3[1] or x3[1]>=x3[2]:
            raise ValueError("x3 should be monotonically increasing")
        # maybe some more checks are needed
        sx10=x3[1]+x3[0]
        sx21=x3[2]+x3[1]
        dx10=x3[1]-x3[0]
        dx21=x3[2]-x3[1]
        dx20=x3[2]-x3[0]
        dx012=dx10*dx21*dx20
        dy10=y3[1]-y3[0]
        dy21=y3[2]-y3[1]
        self.a=(dy21*dx10-dy10*dx21)/dx012
        self.b=(dy10*sx21*dx21-dy21*sx10*dx10)/dx012
        self.c=np.mean(y3-self.a*x3**2-self.b*x3)
    def get_extremum(self):
        if np.abs(self.a)>0:
            xext=-0.5*self.b/self.a
            fext=self.c-0.25*self.b**2/self.a
        else:
            xext=None
            fext=None
        return (xext,fext)
    def has_maximum(self):
        return (self.a<0)
    def has_minimum(self):
        return (self.a>0)
    def is_linear(self):
        return (a==0.) # maybe better: abs(a)<epsilon
    def evaluate(self,x):
        xext,fext=self.get_extremum()
        return self.a*(x-xext)**2+fext

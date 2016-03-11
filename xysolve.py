def xysolve(x,y):
    """
    Simple function to fit a straight line y=ax+b through a bunch of (x,y)
    points, using least squares method (analytical solution) y=a*x+b.
    Returns a, b and D=N*Sxx-Sx**2. If D is zero, then a and b will be None.
    (D is the determinant of linear equations for minimizing the sum of
    squares. N=len(x), and Sx and Sxx are the sums of x and x**2, respectively).
    """
    assert(len(x)==len(y))
    N=len(x)
    Sxx=np.sum(x**2)
    Sx=np.sum(x)
    Sy=np.sum(y)
    Sxy=np.sum(x*y)
    D=N*Sxx-Sx**2
    if np.abs(D)>0:
        a=(N*Sxy-Sx*Sy)/D
        b=(Sxx*Sy-Sx*Sxy)/D
    else:
        print("WARNING: equation singular")
        a=None
        b=None
    return (a,b,D)

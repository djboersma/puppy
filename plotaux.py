import numpy as np

def roundmin(val,margin,roundval):
    """
    Utility function to compute the largest integer multiple of `roundval` that is
    at least `margin` less than the minimum of the `val` array.
    """
    return roundval*np.floor((np.min(val)-margin)/float(roundval))
def roundmax(val,margin,roundval):
    """
    Utility function to compute the smallest integer multiple of `roundval` that is
    at least `margin` larger than the maximum of the `val` array.
    """
    return roundval*np.ceil((np.max(val)+margin)/float(roundval))

def roundlim(val,margin,roundval):
    return [roundmin(val,margin,roundval),roundmax(val,margin,roundval)]

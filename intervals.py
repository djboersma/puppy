#import logging
#logger = logging.getLogger()

def get_interval_overlap_fraction(t0,t1,s0,s1):
    """
    Auxiliary function to compute what fraction of [s0,s1] has an overlap with [t0,t1]
    """
    assert(s1>s0)
    assert(t1>=t0)
    if t1<s0:
        return 0.
    if t0>s1:
        return 0.
    if t1>=s1 and t0<=s0:
        return 1.
    return float(min(t1,s1)-max(t0,s0))/float(s1-s0)

def get_interval_overlap(t0,t1,s0,s1):
    """
    Auxiliary function to compute the size of the overlap of [s0,s1] with [t0,t1]
    """
    assert(s1>s0)
    assert(t1>=t0)
    if t1<s0:
        return 0.
    if t0>s1:
        return 0.
    if t1>=s1 and t0<=s0:
        return float(s1-s0)
    return float(min(t1,s1)-max(t0,s0))

#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np
import os

def find_bb_indices(img,threshold=0):
    """
        Find the indices that define the sub-image of an image that exclude
        the volume of the image that is completely below the threshold value.
    """
    sizes = img.GetSize();
    n=len(sizes)
    assert(n==3)
    aimg = sitk.GetArrayFromImage(img).swapaxes(0,2)
    done = False
    ie = np.array(sizes,dtype=int)
    ib = np.zeros(n,dtype=int)
    while not done:
        done = True
        while ib[0]<ie[0] and (aimg[ib[0],:,:]<threshold).all():
            ib[0] += 1
            done = False
        while ib[0]<ie[0] and (aimg[ie[0]-1,:,:]<threshold).all():
            ie[0] -= 1
            done = False
        while ib[1]<ie[1] and (aimg[:,ib[1],:]<threshold).all():
            ib[1] += 1
            done = False
        while ib[1]<ie[1] and (aimg[:,ie[1]-1,:]<threshold).all():
            ie[1] -= 1
            done = False
        while ib[2]<ie[2] and (aimg[:,:,ib[2]]<threshold).all():
            ib[2] += 1
            done = False
        while ib[2]<ie[2] and (aimg[:,:,ie[2]-1]<threshold).all():
            ie[2] -= 1
            done = False
        if not done:
            print "again"
    return (ib,ie)


if __name__ == '__main__':
    import argparse, os
    parser = argparse.ArgumentParser(description='Python script to copy a 3D MHD image file with air layers removed.')
    parser.add_argument('-I','--inputpath',type=str,dest='INPUTPATH',help='input mhd file name')
    parser.add_argument('-O','--outputpath',type=str,dest='OUTPUTPATH',help='output mhd file name')
    parser.add_argument('-t','--threshold',type=int,default=0,dest='THRESHOLD',help='minimum "interesting" voxel value to keep')
    myargs = parser.parse_args()
    #print("myargs = {}".format(myargs))
    assert(myargs.INPUTPATH != myargs.OUTPUTPATH)
    oldimg = sitk.ReadImage(myargs.INPUTPATH)
    ib,ie = find_bb_indices(oldimg,threshold=myargs.THRESHOLD)
    newimg = oldimg[ib[0]:ie[0],ib[1]:ie[1],ib[2]:ie[2]]
    print("origin old={} new={}".format(oldimg.GetOrigin(),newimg.GetOrigin()))
    print("size old={} new={}".format(oldimg.GetSize(),newimg.GetSize()))
    print("spacing old={} new={}".format(oldimg.GetSpacing(),newimg.GetSpacing()))
    sitk.WriteImage(newimg,myargs.OUTPUTPATH)

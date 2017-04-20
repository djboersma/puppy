#!/usr/bin/env python

import SimpleITK as sitk
import numpy as np

def zpad(img,nzbot=0,nztop=0):
    """
    Return an augmented copy of an image object: the central part is the same
    as the original, but we add nzbot (nztop) copies of the bottom (top) layer
    to the bottom (top).
    """
    if type(img) != sitk.SimpleITK.Image:
        raise TypeError("first arg of zpad should be a SimpleITK image, I got a {}".format(type(img)))
    if type(nzbot) != int or type(nztop) != int:
        raise TypeError("wrong type for nzbot({}) and/or nztop({}), should both be ints".format(type(nzbot),type(nztop)))
    if nzbot<0 or nztop<0 or nztop+nzbot==0:
        raise ValueError("nzbot and nztop should be nonnegative and at least one should be positive");
    origin = list(img.GetOrigin());
    spacing = list(img.GetSpacing());
    oldsize = list(img.GetSize());
    if len(oldsize) != 3:
        raise ValueError("this function should only be used with 3D images");
    if type(img[0,0,0]) == tuple:
        raise TypeError("for now, this function only works with scalar voxel values")
    nx,ny,nz = oldsize
    oldarray = sitk.GetArrayFromImage(img)
    assert(oldarray.shape == (nz,ny,nx) )
    newarray = oldarray.copy()
    newarray.resize( (nz+nzbot+nztop, ny, nx) )
    newarray[nzbot:nzbot+nz,:,:] = oldarray[0:nz,:,:]
    for iz in range(nzbot):
        newarray[iz,:,:] = oldarray[0,:,:]
    for iz in range(nztop):
        newarray[nzbot+nz+iz,:,:] = oldarray[-1,:,:]
    padded_img = sitk.GetImageFromArray(newarray)
    padded_img.SetSpacing(spacing)
    origin[2] -= nzbot * spacing[2]
    padded_img.SetOrigin(origin)
    return padded_img

#######################################################################
# TESTING
#######################################################################
import unittest

class test_zpad(unittest.TestCase):
    def setUp(self):
        self.a2 = np.random.randint(-1024,high=4096,size=(3,4))
        self.img2 = sitk.GetImageFromArray(self.a2)
        self.a3 = np.random.randint(-1024,high=4096,size=(3,4,5))
        #print("3d image has {} voxels with value zero".format(np.sum(self.a3==0)))
        self.img3 = sitk.GetImageFromArray(self.a3)
        self.a4 = np.random.randint(-1024,high=4096,size=(3,4,5,6))
        self.img4 = sitk.GetImageFromArray(self.a4)
    def test_wrong_types(self):
        with self.assertRaises(TypeError):
            zpad(self.a3,1,2); # not an image
        with self.assertRaises(TypeError):
            zpad(self.img3,1.,2); # float instead of int
        with self.assertRaises(TypeError):
            zpad(self.img3,1,2.); # float instead of int
        with self.assertRaises(TypeError):
            zpad(self.img3,1.,2.); # float instead of int
        with self.assertRaises(ValueError):
            zpad(self.img3,0,0); # at least one padding argument should be larger than 0
        with self.assertRaises(ValueError):
            zpad(self.img3,-5,10); # both ints should be positive
        with self.assertRaises(ValueError):
            zpad(self.img3,10,-5); # both ints should be positive
        with self.assertRaises(ValueError):
            zpad(self.img2,1,2); # image should be 3D
        with self.assertRaises(TypeError):
            zpad(self.img4,1,2); # image should be 3D and have scalar voxel type
    def test_bottom(self):
        for p in [1,5,15]:
            img3 = sitk.Image(self.img3)
            img3.MakeUnique()
            size3 = img3.GetSize()
            padded_img3 = zpad(img3,p,0)
            psize3 = padded_img3.GetSize()
            self.assertEqual(size3[0],psize3[0])
            self.assertEqual(size3[1],psize3[1])
            self.assertEqual(size3[2]+p,psize3[2])
            ai3 = sitk.GetArrayFromImage(img3).swapaxes(0,2)
            api3 = sitk.GetArrayFromImage(padded_img3).swapaxes(0,2)
            for iz in range(p):
                #print("testing bottom iz={} for p={}".format(iz,p))
                self.assertTrue((api3[:,:,iz]==ai3[:,:,0]).all())
    def test_top(self):
        for q in [1,5,15]:
            img3 = sitk.Image(self.img3)
            img3.MakeUnique()
            size3 = img3.GetSize()
            padded_img3 = zpad(img3,0,q)
            psize3 = padded_img3.GetSize()
            self.assertEqual(size3[0],psize3[0])
            self.assertEqual(size3[1],psize3[1])
            self.assertEqual(size3[2]+q,psize3[2])
            ai3 = sitk.GetArrayFromImage(img3).swapaxes(0,2)
            api3 = sitk.GetArrayFromImage(padded_img3).swapaxes(0,2)
            for iz in range(q):
                #print("testing top iz={} for q={}".format(iz,q))
                self.assertTrue((api3[:,:,-iz-1]==ai3[:,:,-1]).all())
    def test_both(self):
        for p in [0,1,5,15]:
            for q in [0,4,16]:
                if p==0 and q==0:
                    continue
                img3 = sitk.Image(self.img3)
                img3.MakeUnique()
                size3 = img3.GetSize()
                padded_img3 = zpad(img3,p,q)
                psize3 = padded_img3.GetSize()
                self.assertEqual(size3[0],psize3[0])
                self.assertEqual(size3[1],psize3[1])
                self.assertEqual(size3[2]+p+q,psize3[2])
                ai3 = sitk.GetArrayFromImage(img3).swapaxes(0,2)
                api3 = sitk.GetArrayFromImage(padded_img3).swapaxes(0,2)
                for iz in range(p):
                    #print("both: testing bottom iz={} for p={}".format(iz,q))
                    self.assertTrue((api3[:,:,iz]==ai3[:,:,0]).all())
                for iz in range(q):
                    #print("both: testing top iz={} for q={}".format(iz,q))
                    self.assertTrue((api3[:,:,-iz-1]==ai3[:,:,-1]).all())

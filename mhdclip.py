#!/usr/bin/env python

import os, sys
import numpy as np
from glob import glob


mhdlist=[]
zapdir=os.environ.get('PATIENT_DATA_COLLECTION',None)
for a in sys.argv[1:]:
    print "checking a=",a
    nglob=0
    for b in glob(a):
        print "checking b=",b
        if b[-4:]==".mhd":
            print "bingo"
            nglob+=1
            mhdlist.append(b)
    if nglob>0:
        continue
    if zapdir:
        print "checking zapdir/a=",os.path.join(zapdir,a)
        for b in glob(os.path.join(zapdir,a)):
            print "checking b=",b
            if b[-4:]==".mhd":
                print "bingo"
                nglob+=1
                mhdlist.append(b)

print mhdlist
#sys.exit(0)

import SimpleITK as sitk

def new_sitk_clamp(img):
    """
    Clamping (forcing pixel values within a range) can be done with
    sitk.Clamp(img,pixelvalue,lower,upper), but in SimpleITK version 0.8.* this
    function was buggy (output image always identical to input image).
    This wrapper function is used in case your SimpleITK installation seems new enough.
    """
    return sitk.Clamp(img,img.GetPixelIDValue(),-1024.,3000.)

def old_sitk_clamp(img):
    """
    Clamping (forcing pixel values within a range) can be done with
    sitk.Clamp(img,pixelvalue,lower,upper), but in SimpleITK version 0.8.* this
    function was buggy (output image always identical to input image).
    This wrapper function is used in case your SimpleITK installation seems too old.
    """
    aimg=sitk.GetArrayFromImage(img)
    nz,ny,nx=aimg.shape
    mx,my,mz=img.GetSize()
    print("nx={} ny={} nz={}".format(nx,ny,nz))
    print("mx={} my={} mz={}".format(mx,my,mz))
    assert( (mx,my,mz) == (nx,ny,nz) )
    mask=(aimg>3000)
    nover=np.sum(mask.flat)
    print("{} pixels have a value larger than 3000".format(nover))
    zi,yi,xi = np.indices((nz,ny,nx))
    for iz,iy,ix,pix0 in zip(zi[mask].flat,yi[mask].flat,xi[mask].flat,aimg[mask].flat):
        print("({},{},{}): going to force pix0={}={} down to 3000".format(ix,iy,iz,pix0,img.GetPixel(ix,iy,iz)))
        img.SetPixel(ix,iy,iz,3000)
    mask=(aimg<-1024)
    nunder=np.sum(mask.flat)
    print("{} pixels have a value less than -1024".format(nunder))
    for ix,iy,iz,pix0 in zip(xi[mask].flat,yi[mask].flat,zi[mask].flat,aimg[mask].flat):
        print("({},{},{}): going to force pix0={} up to -1024".format(ix,iy,iz,pix0))
        img.SetPixel(ix,iy,iz,-1024)
    return img

if sitk.Version_MajorVersion() == 0 and sitk.Version_MinorVersion() <= 8:
    print("WARNING: using manual clamping because your SimpleITK version {} ".format(sitk.Version_VersionString())
          +"seems to be an old one (older than July 2015) with a bug in the Clamp function.")
    sitk_clamp = old_sitk_clamp
else:
    print("YAY: using native sitk clamping function because " +
          "your SimpleITK version seems to be new enough: {}.".format(sitk.Version_VersionString()))
    sitk_clamp = new_sitk_clamp

for mhd in mhdlist:
    img=sitk.ReadImage(mhd)
    cimg = sitk_clamp(img)
    cmhd=mhd.replace(".mhd",".clipped.mhd");
    sitk.WriteImage(cimg,cmhd)

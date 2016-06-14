#!/usr/bin/env python

import os, sys
import numpy as np
from glob import glob


mhdlist=[]
zapdir='/Users/montecarlo/Analysis/DICOM/Hodgkin/zaptestdata'
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
for mhd in mhdlist:
    img=sitk.ReadImage(mhd)
    aimg=sitk.GetArrayFromImage(img)
    nz,ny,nx=aimg.shape
    mx,my,mz=img.GetSize()
    print("nx={} ny={} nz={}".format(nx,ny,nz))
    print("mx={} my={} mz={}".format(mx,my,mz))
    assert( (mx,my,mz) == (nx,ny,nz) )
    mask=(aimg>3000)
    zi,yi,xi = np.indices((nz,ny,nx))
    for iz,iy,ix,pix0 in zip(zi[mask].flat,yi[mask].flat,xi[mask].flat,aimg[mask].flat):
        print("({},{},{}): going to force pix0={}={} down to 3000".format(ix,iy,iz,pix0,img.GetPixel(ix,iy,iz)))
        img.SetPixel(ix,iy,iz,3000)
    mask=(aimg<-1024)
    for ix,iy,iz,pix0 in zip(xi[mask].flat,yi[mask].flat,zi[mask].flat,aimg[mask].flat):
        print("({},{},{}): going to force pix0={} up to -1024".format(ix,iy,iz,pix0))
        img.SetPixel(ix,iy,iz,-1024)
    cmhd=mhd.replace(".mhd",".clipped.mhd");
    print("writing clamped file",cmhd);
    sitk.WriteImage(img,cmhd)

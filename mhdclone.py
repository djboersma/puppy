#!/usr/bin/env python

import os, sys
import numpy as np
from glob import glob

nclones=10

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
for mhd in mhdlist:
    img=sitk.ReadImage(mhd)
    print("Going to make 10 clones of {}".format(nclones,mhd))
    for iclone in range(nclones):
        cmhd=mhd.replace(".mhd",".{}.mhd".format(iclone));
        print("writing clone {}...".format(cmhd))
        sitk.WriteImage(img,cmhd)

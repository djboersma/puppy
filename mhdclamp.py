#!/usr/bin/env python

import os, sys
import numpy as np
from glob import glob
import SimpleITK as sitk
import logging
logger = logging.getLogger()

def get_mhd_list(filepats):
    mhdlist=[]
    for a in sys.argv[1:]:
        logger.debug("checking pattern={}".format(a))
        nglob=0
        for b in glob(a):
            logger.debug("checking path={}".format(b))
            if b[-4:]==".mhd":
                logger.debug("bingo")
                nglob+=1
                mhdlist.append(b)
        if nglob>0:
            continue
        zapdir=os.environ.get('PATIENT_DATA_COLLECTION',None)
        if zapdir:
            fpat = os.path.join(zapdir,a)
            logger.debug("checking pattern=zapdir/a={}".format(fpat))
            for b in glob(fpat):
                logger.debug("checking b={}".format(b))
                if b[-4:]==".mhd":
                    logger.debug("bingo")
                    nglob+=1
                    mhdlist.append(b)
        if nglob == 0:
            logger.warn("NO MHD FILES FOUND corresponding to a={}".format(a))
    logger.debug("MHD LIST = " + ",".join(mhdlist))
    return mhdlist

def itk_clamp(mhd):
    """
    Clamping (forcing pixel values within a range) can be done with
    sitk.Clamp(img,pixelvalue,lower,upper), but in SimpleITK version 0.8.* this
    function was buggy (output image always identical to input image).
    This wrapper function is used in case your SimpleITK installation seems new enough.
    """
    img=sitk.ReadImage(mhd)
    cmhd=mhd.replace(".mhd",".clamped.mhd");
    img_clamped=sitk.Clamp(img,img.GetPixelIDValue(),-1024,3000)
    sitk.WriteImage(img_clamped,cmhd)

def np_clamp(mhd):
    """
    Clamping (forcing pixel values within a range) can be done with
    sitk.Clamp(img,pixelvalue,lower,upper), but in SimpleITK version 0.8.* this
    function was buggy (output image always identical to input image).
    This wrapper function is used in case your SimpleITK installation seems too old.
    """
    img=sitk.ReadImage(mhd)
    aimg=sitk.GetArrayFromImage(img)
    nz,ny,nx=aimg.shape
    mx,my,mz=img.GetSize()
    logger.debug("nx={} ny={} nz={}".format(nx,ny,nz))
    logger.debug("mx={} my={} mz={}".format(mx,my,mz))
    assert( (mx,my,mz) == (nx,ny,nz) )
    mask=(aimg>3000)
    zi,yi,xi = np.indices((nz,ny,nx))
    for iz,iy,ix,pix0 in zip(zi[mask].flat,yi[mask].flat,xi[mask].flat,aimg[mask].flat):
        logger.debug("({},{},{}): going to force pix0={}={} down to 3000".format(ix,iy,iz,pix0,img.GetPixel(ix,iy,iz)))
        img.SetPixel(ix,iy,iz,3000)
    mask=(aimg<-1024)
    for ix,iy,iz,pix0 in zip(xi[mask].flat,yi[mask].flat,zi[mask].flat,aimg[mask].flat):
        logger.debug("({},{},{}): going to force pix0={} up to -1024".format(ix,iy,iz,pix0))
        img.SetPixel(ix,iy,iz,-1024)
    cmhd=mhd.replace(".mhd",".clipped3.mhd");
    logger.debug("writing clamped file",cmhd);
    sitk.WriteImage(img,cmhd)

def clamp_mhds(fpats):
    if 100*sitk.Version_MajorVersion() + sitk.Version_MinorVersion()>=9:
        logger.debug("Congratulations, your SimpleITK installation is new enough: {}.".format(sitk.Version_VersionString()))
        clamp = itk_clamp
    else:
        logger.warn("Your SimpleITK installation is old ({}), with a buggy Clamp implementation.".format(sitk.Version_VersionString()))
        clamp = np_clamp
    for mhd in get_mhd_list(fpats):
        clamp(mhd)

def get_opts():
    import argparse
    parser = argparse.ArgumentParser("foobar")
    parser.add_argument("-v","--verbose",action='store_true',help="print debugging information?")
    parser.add_argument("filepattern",nargs='+',help="list of mhd file patterns")
    args=parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    logger.debug("args={}".format(args))
    return args

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    args=get_opts()
    clamp_mhds(args.filepattern)

#!/usr/bin/env python

"""
development snippet: how to use DicomIO and ROIs from a structure set
to get a bounding box corresponding to the "external" structure, and then
crop the image files to get mean lean input for the simulations...
"""

import os
import SimpleITK as sitk
import numpy as np
import logging
import time
import argparse
logger = logging.getLogger()
try:
    from RtStruct import RtStruct
except ImportError:
    logger.error("in order to use this script you'll need the DicomIO project in your python path")
    raise ImportError("in order to use this script you'll need the DicomIO project in your python path")

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-P','--plan',help='plan name')
    parser.add_argument('-l','--list',action='store_true',help='list plan names')
    parser.add_argument('-c','--clamp',action='store_true',help='do also clamp, not only crop')
    parser.add_argument("-v","--verbose",action='store_true',help="print debugging information?")
    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(level=logging.DEBUG)
    if not hasattr(args,'plan'):
        args.plan = 'all'
    logger.debug("args={}".format(args))
    return args

def stringpos(pos):
    return "({0:6.1f},{1:6.1f},{2:6.1f})".format(pos[0],pos[1],pos[2])

def boxvolume(pos1,pos2):
    volume=np.prod(np.array(pos2)-np.array(pos1))
    return volume

def string2pos(pos1,pos2):
    return "[{0},{1}]".format(stringpos(pos1),stringpos(pos2))


def get_dictionary():
    logger.debug("going to get DICOM dir")
    dicomdir = os.path.join(os.environ['DICOM_DIRECTORY'],'DICOM_4D_STUDIE')
    logger.debug("going to get MHD dir")
    patdir = os.environ['PATIENT_DATA_COLLECTION']
    if not os.path.isdir(dicomdir):
        raise IOError('cannot find DICOM dir {}'.format(dicomdir))
    logger.debug("filling dict")
    the_dict = dict()
    the_dict['40_P1'] = dict(rs = os.path.join(dicomdir,'40_P1','RS.Test_C_hodgkin-uas-040.dcm'),
                             mhd = os.path.join(patdir,'pat40_lunga00.mhd') )
    the_dict['44'] = dict(rs = os.path.join(dicomdir,'44','STRUCTURESET.dcm'),
                          mhd = os.path.join(patdir,'pat44_lunga00.mhd') )
    the_dict['51_P1'] = dict(rs = os.path.join(dicomdir,'51_P1','RS.Test_C_hodgkin-051.dcm'),
                             mhd = os.path.join(patdir,'pat51_lunga00.mhd') )
    the_dict['53_SHORT'] = dict(rs = os.path.join(dicomdir,'53_SHORT','STRUCTURESET.dcm'),
                                mhd = os.path.join(patdir,'pat53_lunga00.mhd') )
    the_dict['53_LONG'] = dict(rs = os.path.join(dicomdir,'53_LONG','STRUCTURESET.dcm'),
                               mhd = os.path.join(patdir,'pat53_lunga00.mhd') )
    logger.debug("checking dict")
    for plan,data in the_dict.items():
        logger.debug("plan="+plan+"rs="+data['rs'])
        logger.debug("plan="+plan+"mhd="+data['mhd'])
        assert(os.path.exists(data['rs']))
        assert(os.path.exists(data['mhd']))
    logger.debug("returning dict")
    return the_dict

def get_structure_set(p=None,clamp=False,verbose=False,roiname = '4D External', altroiname='External_4D'):
    logger.debug("going to get dictionary")
    the_dict = get_dictionary()
    if type(p)==list:
        planlist=p
    elif p is None or p=="all":
        planlist=the_dict.keys()
    elif type(p)==str:
        planlist=[p]
    else:
        raise RuntimeError("don't know what to do with {}".format(p))
    for plan in planlist:
        rsfile = the_dict[plan]['rs']
        mhdfile = the_dict[plan]['mhd']
        logger.debug("going to get RS file {}".format(rsfile))
        ss = RtStruct(rsfile)
        logger.debug("going to parse RS file {}".format(rsfile))
        before = time.time()
        ss.parse()
        after = time.time()
        logger.debug("Parsing took {tparse} seconds (wall time), found {Nroi} ROIs, {ok} {roiname}.".format(
            tparse=after-before,
            Nroi=len(ss.all_roi_names()),
            ok=("including" if roiname in ss.all_roi_names() else "NOT INCLUDING"),
            roiname=roiname))
        if roiname in ss.all_roi_names():
            bb = ss.roi_bounding_box(roiname)
        elif altroiname in ss.all_roi_names():
            bb = ss.roi_bounding_box(altroiname)
        else:
            logger.warn("OOPS: {ss} contains neither {roi} nor {altroi}".format(ss=rsfile,roi=roiname,altroi=altroiname))
            logger.warn("OOPS: available roi names are: {all}".format(all="; ".join(ss.all_roi_names())))
            continue
        volbb = boxvolume(bb[0],bb[1])
        mhd = sitk.ReadImage(mhdfile)
        origin = np.array(mhd.GetOrigin())
        anti_origin = np.array(mhd.GetOrigin())+np.array(mhd.GetSize())*np.array(mhd.GetSpacing())
        volmhd = boxvolume(origin,anti_origin)
        logger.info("plan={plan:8s} BB={bb} mhd={mhd} volBB/volMHD={Vratio:4.1f}%".format(
            plan = plan,
            bb = string2pos(bb[0],bb[1]),
            mhd = string2pos(origin,anti_origin),
            Vratio = 100.*volbb/volmhd))
        spacing = np.array(mhd.GetSpacing())
        ibb0 = np.round((np.array(bb[0])-origin)/spacing).astype(int)
        ibb1 = np.round((np.array(bb[1])-origin)/spacing).astype(int)
        sizes = np.array(mhd.GetSize())
        if (ibb0<0).any() or (ibb0>=sizes).any() or (ibb1<0).any() or (ibb1>=sizes).any():
            logging.error('some BB indexes out of range, ibb0={} ibb1={} sizes={}'.format(ibb0,ibb1,sizes))
        else:
            logging.info('crop: ibb0={} ibb1={} sizes={}'.format(ibb0,ibb1,sizes))
            oneoneone = np.ones(3,dtype=int)
            wantsize=ibb1-ibb0+oneoneone
            # now we add one voxel margin wherever we are not already on the edge
            ilow=(ibb0>0).astype(int)
            wantsize+=ilow
            ibb0-=ilow
            ihigh=(ibb1+oneoneone<sizes).astype(int)
            wantsize+=ihigh
            for phase in range(0,100,10):
                phase_str = "{0:02d}".format(phase)
                inmhdfile = mhdfile.replace("00",phase_str)
                outmhdfile = "cropped_"+os.path.basename(inmhdfile)
                # ITK filter names are bonkers! The 'Crop' filter does not actually crop...
                cropmhd = sitk.RegionOfInterest( sitk.ReadImage(inmhdfile), wantsize, ibb0)
                # check check check
                csizes = np.array(cropmhd.GetSize())
                if phase == 0:
                    logging.info('crop: old={old}, new={new}, want={want}, kept {pct:.2f} % of voxels'.format(
                        old=stringpos(sizes),
                        new=stringpos(csizes),
                        want=stringpos(wantsize),
                        pct=np.prod(csizes)*100./np.prod(sizes)))
                # save output
                if clamp:
                    cropclampmhd = sitk.Clamp(cropmhd,cropmhd.GetPixelIDValue(),-1024,3000)
                    sitk.WriteImage(cropclampmhd,outmhdfile.replace('.mhd','.clamped.mhd'))
                else:
                    sitk.WriteImage(cropmhd,outmhdfile)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    args = get_options()
    if args.list:
       d = get_dictionary()
       logger.info("plans in dictionary:\n" + "\n".join(d.keys()))
    else:
       get_structure_set(args.plan,args.clamp,args.verbose)





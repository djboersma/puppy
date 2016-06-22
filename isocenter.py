#!/usr/bin/env python

import dicom
import sys
import logging
logger = logging.getLogger()

def get_isocenter(planfile,field=0):
    plan = dicom.read_file(planfile)
    print("plan has {} beams".format(len(plan.IonBeamSequence)))
    beam0 = plan.IonBeamSequence[field]
    icp0 = beam0.IonControlPointSequence[0]
    iso = icp0.IsocenterPosition
    return (float(iso[0]),float(iso[1]),float(iso[2]))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) == 2:
        arg1=str(sys.argv[1])
        logging.debug("sys.argv[1]={}".format(arg1))
        logging.info("iso={} mm".format(get_isocenter(arg1)))
    if len(sys.argv) == 3:
        arg1=str(sys.argv[1])
        arg2=int(sys.argv[2])
        print("sys.argv[1]={} sys.argv[2]={}".format(arg1,arg2))
        print("iso={} mm".format(get_isocenter(arg1,arg2)))
    else:
        raise RuntimeError("I need exactly one argument (a DICOM planfile)")

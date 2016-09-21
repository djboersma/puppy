#!/usr/bin/env python

import os
import dicom
import argparse
import logging
from eclipse_proton_plan import eclipse_proton_plan
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description='Python script to extract useful information out of DICOM plan file.')
parser.add_argument('-P','--plan',dest='PLAN',help='plan file', required=True)
parser.add_argument('-F','--field',dest='FIELD',default='',
                    help='If a name or number is specified, use only a single field out of a multi-field plan (default: use all fields)')
args = parser.parse_args()
filepaths = [args.PLAN]
filepath = None
if 'DICOM_DIRECTORY' in os.environ:
    filepaths.append( os.path.join(os.environ['DICOM_DIRECTORY'],args.PLAN) )
for fp in filepaths:
    if os.path.exists(fp):
        filepath = fp
if filepath is None:
    raise IOError("cannot find plan file {}".format(args.PLAN))

logger.debug('the plan file is {}'.format(filepath))
if args.FIELD:
    logger.info('using only field number {}'.format(args.FIELD))
else:
    logger.debug('using all fields')

plan = eclipse_proton_plan(filepath)
logger.info("{}".format(plan))

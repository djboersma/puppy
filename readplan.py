#!/usr/bin/env python

import os
import dicom
import argparse
import logging
from eclipse_proton_plan import eclipse_proton_plan
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)
myargs = None

def get_options():
    global myargs
    parser = argparse.ArgumentParser(description='Python script to extract useful information out of DICOM plan file.')
    parser.add_argument('-P','--path',dest='PATH',help='path of plan file')
    parser.add_argument('-l','--list',dest='LIST',default=False,action='store_true')
    parser.add_argument('-F','--field',dest='FIELD',default='',nargs=1,
                        help='If a name or number is specified, use only a single field out of a multi-field plan (default: use all fields)')
    myargs = parser.parse_args()
    #return args


def find_and_show_plan_file():
    global myargs
    filepaths = [myargs.PATH]
    filepath = None
    if 'DICOM_DIRECTORY' in os.environ:
        filepaths.append( os.path.join(os.environ['DICOM_DIRECTORY'],myargs.PATH) )
    for fp in filepaths:
        if os.path.exists(fp):
            filepath = fp
            break
    if filepath is None:
        raise IOError("cannot find plan file {}".format(myargs.PATH))

    logger.debug('the plan file is {}'.format(filepath))
    if myargs.FIELD:
        logger.info('using only field number {}'.format(myargs.FIELD))
    else:
        logger.debug('using all fields')
    plan = eclipse_proton_plan(filepath,True)
    logger.info("{}".format(plan))

def main():
    global myargs
    get_options()
    logger.debug('path={}'.format(myargs.PATH))
    logger.debug('field={}'.format(myargs.FIELD))
    logger.debug('list={}'.format(myargs.LIST))
    if myargs.LIST:
        find_and_show_plan_file()

main()

#!/usr/bin/env python

import numpy as np
import os,sys
from roi_utils import region_of_interest, get_intersection_volume, list_roinames
import dicom
import logging
logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

logger.info("PATIENT 44")
ss44 = dicom.read_file(os.path.join(os.environ['DICOM_DIRECTORY'],'DICOM_4D_STUDIE','44','STRUCTURESET.dcm'))
#print("\n".join(list_roinames(ss)))
itv44 = region_of_interest(ss44,'ITV 0+20+50+80+M')
ctv44 = region_of_interest(ss44,'CTV 50')
logger.info('itv volume is {}'.format(itv44.get_volume()))
logger.info('ctv volume is {}'.format(ctv44.get_volume()))
logger.info('itv ctv intersection volume is {}'.format(get_intersection_volume([itv44,ctv44])))

logger.info("PATIENT 51")
ss51 = dicom.read_file(os.path.join(os.environ['DICOM_DIRECTORY'],'DICOM_4D_STUDIE','51_P1','RS.Test_C_hodgkin-051.dcm'))
ctv51 = region_of_interest(ss51,'4D CTV')
logger.info('ctv volume is {}'.format(ctv51.get_volume()))

logger.info("PATIENT 53")
ss53 = dicom.read_file(os.path.join(os.environ['DICOM_DIRECTORY'],'DICOM_4D_STUDIE','53_SHORT','STRUCTURESET.dcm'))
ctv53 = region_of_interest(ss51,'CTV50')
logger.info('ctv volume is {}'.format(ctv53.get_volume()))

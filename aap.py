#!/usr/bin/env python

import numpy as np
import os,sys
from roi_utils import region_of_interest, get_intersection_volume, list_roinames
import dicom

ss = dicom.read_file(os.path.join(os.environ['DICOM_DIRECTORY'],'DICOM_4D_STUDIE','44','STRUCTURESET.dcm'))

print("\n".join(list_roinames(ss)))

itv = region_of_interest(ss,'ITV 0+20+50+80+M')
ctv = region_of_interest(ss,'CTV 50')

print('itv volume is {}'.format(itv.get_volume()))
print('ctv volume is {}'.format(ctv.get_volume()))
print('itv ctv intersection volume is {}'.format(get_intersection_volume([itv,ctv])))

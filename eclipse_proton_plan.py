import dicom
import logging
import numpy as np

logger = logging.getLogger()

def get_if_available(ds,name,defaultvalue="NOT AVAILABLE"):
    """
    Retrieve data corresponding to a tag.
    If tag is not available in the given data set, use the default value.
    """
    tag = dicom.datadict.tag_for_name(name)
    if ds.has_key(tag):
        value = ds[tag].value
    else:
        value = defaultvalue
    logger.debug("{}={}".format(name,value))
    return value

def is_close(x,y,eps=1e-6):
    sumabs = np.abs(x)+np.abs(y)
    absdif = np.abs(x-y)
    ok = (sumabs == 0) or (absdif < eps*0.5*sumabs)
    return ok

class eclipse_proton_layer(object):
    def __init__(self,ctrlpnt,cumsumchk=[],verbose=False):
        if verbose:
            logger.debug('{}. ion beam sequence element with type {}'.format(i,type(ion_beam)))
            for k in ion_beam.keys():
                if dicom.datadict.dictionary_has_tag(k):
                    kw = dicom.datadict.keyword_for_tag(k)
                else:
                    kw = "(UNKNOWN)"
                logger.debug('k={} keyword={}'.format(k,kw))
        nspot =int(ctrlpnt.NumberOfScanSpotPositions)
        assert(ctrlpnt.NominalBeamEnergyUnit == 'MEV')
        if nspot == 1:
            self.w = np.array([float(ctrlpnt.ScanSpotMetersetWeights)])
        else:
            self.w = np.array([float(w) for w in ctrlpnt.ScanSpotMetersetWeights])
        assert( nspot == len( self.w ) )
        assert( nspot*2 == len( ctrlpnt.ScanSpotPositionMap ) )
        self.cpindex = int(ctrlpnt.CPIndex)
        self.msw = float(ctrlpnt.CumulativeMetersetWeight)
        if cumsumchk:
            logger.debug("msw={0:14.8g} cumsumchk={1:14.8g} diff={2:14.8g}".format( self.msw, cumsumchk[0], self.msw - cumsumchk[0]))
            assert( is_close(self.msw,cumsumchk[0]) )
        self.energy = float(ctrlpnt.NominalBeamEnergy)
        self.npainting = int(ctrlpnt.NumberOfPaintings)
        xy = np.array([float(pos) for pos in ctrlpnt.ScanSpotPositionMap]).reshape(nspot,2)
        self.x = np.array(xy[:,0])
        self.y = np.array(xy[:,1])
        self.spot_id = str(ctrlpnt.ScanSpotTuneID)
        self.dx = float(ctrlpnt.ScanningSpotSize[0])
        self.dy = float(ctrlpnt.ScanningSpotSize[1])
        #assert( self.msw>0 )
        if (self.w<=0).any():
            logger.debug("%s".format(self))
            logger.error("layer number {} has {} spots with zero weight:".format(self.cpindex,np.sum(self.w<=0)))
            assert(np.sum(self.w<=0) == len(self.w))
        assert( (self.w>0).all() ^ (self.w==0.0).all() )
        wsum=np.sum(self.w)
        logger.debug("layer number {} has {} spots, of which {} with zero weight, cumsum={}, sum(w)={}".format(self.cpindex,len(self.w),np.sum(self.w<=0),self.msw,wsum))
        # assert( (wsum-self.msw) < 1e-6*(wsum+self.msw) )
        assert( self.npainting == 1 )
        cumsumchk[0]+=wsum
    def __repr__(self):
        txt = "layer number {}, energy={} MeV, msw={}, dx={} dy={}".format(self.cpindex, self.energy, self.msw, self.dx, self.dy)
        for x,y,w in zip(self.x,self.y,self.w):
            txt += "\nx={}, y={} msw={}".format(x,y,w)
        return txt

class eclipse_proton_field(object):
    def __init__(self,ionbeam,verbose=False):
        if verbose:
            logger.debug('{}. ion beam sequence element with type {}'.format(i,type(ion_beam)))
            for k in ion_beam.keys():
                if dicom.datadict.dictionary_has_tag(k):
                    kw = dicom.datadict.keyword_for_tag(k)
                else:
                    kw = "(UNKNOWN)"
                logger.debug('k={} keyword={}'.format(k,kw))
        if not ionbeam.RadiationType == 'PROTON':
            raise RuntimeError('No protons!?? What the hell?')
        self.name = ionbeam.BeamName
        self.number = int(ionbeam.BeamNumber)
        self.msw = ionbeam.FinalCumulativeMetersetWeight
        self.dose_unit = ionbeam.PrimaryDosimeterUnit
        self.gantry_angle = ionbeam.IonControlPointSequence[0].GantryAngle
        self.nlayers = len(ionbeam.IonControlPointSequence)
        assert(self.nlayers%2 == 0)
        cumsumchk=[0.]
        self.layers = [ eclipse_proton_layer(ionbeam.IonControlPointSequence[2*i],cumsumchk) for i in range(self.nlayers/2) ]
        assert(len(self.layers)*2 == self.nlayers)
        assert( is_close(self.msw,cumsumchk[0]) )
    def numpify(self):
        biglist = []
        for l in self.layers:
            energy = l.energy
            dx = l.dx
            dy = l.dy
            for x,y,w in zip(l.x,l.y,l.w):
                biglist.append([x,y,w,energy,dx,dy])
        names = "x","y","w","energy","dx","dy"
        return names,np.array(biglist)

class eclipse_proton_plan(object):
    def __init__(self,filename):
        self.new_file(filename)
    def new_file(self,filename):
        self.filename = filename
        self.plan = dicom.read_file(filename)
        self.basic_checks()
        self.RTPlanDate = get_if_available(self.plan,'RTPlanDate')
        self.RTPlanLabel = get_if_available(self.plan,'RTPlanLabel')
        self.RTPlanDescription = get_if_available(self.plan,'RTPlanDescription')
        self.PatientBirthDate = get_if_available(self.plan,'PatientBirthDate')
        self.PatientID = get_if_available(self.plan,'PatientID')
        self.PatientSex = get_if_available(self.plan,'PatientSex')
        self.PatientName = get_if_available(self.plan,'PatientName')
        self.fields = [ eclipse_proton_field(ion_beam) for ion_beam in self.plan.IonBeamSequence ]
    def basic_checks(self):
        if not self.plan.has_key( dicom.datadict.tag_for_name('IonBeamSequence') ):
            logger.error("This DICOM file '{}' does not seem to be a PLAN file, because it does not contain an Ion Beam Sequence.".format(self.filename))
            raise RuntimeError("this DICOM file does not seem to contain Ion Beams.")
        #logger.debug('plan label={} descr={} date={}'.format(plan[rtplan_label_tag].value,plan[rtplan_descr_tag].value,plan[rtplan_date_tag].value))
        ion_beams = self.plan.IonBeamSequence
        logger.debug('the ion beams data element has value type {}'.format(type(ion_beams)))
        logger.debug('the ion beam sequence has {} ion beams'.format(len(ion_beams)))
        if len(ion_beams) == 0:
            logger.error("Ion beams sequence is actually empty.")
            raise RuntimeError("No beams, no glory!")
    def __repr__(self):
        txt = "Plan: '{}' created on {}, description: '{}'\n".format(self.RTPlanLabel,self.RTPlanDate,self.RTPlanDescription)
        txt += "Patient: ID={}, sex={}, name={}, birth date = {}'".format(self.PatientID,self.PatientSex,self.PatientName,self.PatientBirthDate)
        for field in self.fields:
            txt += "\nField '{}' has {} layers, angle={}, first/last layer E={}/{} MeV, dose={} {}, {} spots.".format(
                    field.name,len(field.layers), field.gantry_angle,
                    field.layers[0].energy,field.layers[-1].energy,
                    field.msw,field.dose_unit,
                    sum([len(layer.w) for layer in field.layers]))
        return txt

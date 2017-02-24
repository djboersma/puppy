import pylab
import numpy as np
import logging
logger = logging.getLogger()

def make_elaborate_debugging_plots(aimg,amask,debuglabel,roiname):
    logger.debug("START debugging "+debuglabel)
    amask0=amask*(aimg==0.0)
    iz0,iy0,ix0 = np.where(amask0)
    iz,iy,ix = np.where(amask)
    logger.debug('indices for zero dose voxels have {}/{}/{} different values (out of {}/{}/{})'.format(
        len(set(ix0[:])), len(set(iy0[:])), len(set(iz0[:])),
        len(ix0), len(iy0), len(iz0)))
    assert((aimg[iz0,iy0,ix0]==0.0).all())
    hix0,xedges=np.histogram(ix0,bins=np.arange(-0.5,aimg.shape[2]-0.49,1.0))
    hiy0,yedges=np.histogram(iy0,bins=np.arange(-0.5,aimg.shape[1]-0.49,1.0))
    hiz0,zedges=np.histogram(iz0,bins=np.arange(-0.5,aimg.shape[0]-0.49,1.0))
    izmax = np.array(hiz0).argmax()
    hixy0,foo,bar=np.histogram2d(ix0,iy0,bins=(xedges,yedges))
    hiyz0,foo,bar=np.histogram2d(iy0,iz0,bins=(yedges,zedges))
    hizx0,foo,bar=np.histogram2d(iz0,ix0,bins=(zedges,xedges))
    rxmin=xedges[np.min(ix)]
    rxmax=xedges[np.max(ix)+1]
    rymin=yedges[np.min(iy)]
    rymax=yedges[np.max(iy)+1]
    rzmin=zedges[np.min(iz)]
    rzmax=zedges[np.max(iz)+1]
    ahixy0 = np.array(hixy0).T
    ahiyz0 = np.array(hiyz0).T
    ahizx0 = np.array(hizx0).T
    oldfig=pylab.gcf()
    fig=pylab.figure(num="dvhdebug_"+debuglabel,figsize=[10,15])
    pylab.suptitle('distribution of indices of voxels in "{}" with dose==0'.format(roiname))
    pylab.subplot(321)
    pylab.bar(xedges[:-1],hix0,lw=0,color='b')
    pylab.xlim(xedges[0],xedges[-1])
    ymin,ymax=pylab.ylim()
    pylab.plot([rxmin]*2,[ymin,ymax],'r-')
    pylab.plot([rxmax]*2,[ymin,ymax],'r-')
    pylab.xlabel('X index')
    pylab.subplot(322)
    pylab.pcolormesh(yedges,zedges,ahiyz0)
    pylab.plot([rymin,rymax,rymax,rymin,rymin],[rzmin,rzmin,rzmax,rzmax,rzmin],'y-')
    pylab.xlim(yedges[0],yedges[-1])
    pylab.ylim(zedges[0],zedges[-1])
    pylab.xlabel('Y index')
    pylab.ylabel('Z index')
    pylab.colorbar()
    pylab.subplot(323)
    pylab.bar(yedges[:-1],hiy0,lw=0,color='b')
    ymin,ymax=pylab.ylim()
    pylab.plot([rymin]*2,[ymin,ymax],'r-')
    pylab.plot([rymax]*2,[ymin,ymax],'r-')
    pylab.xlabel('Y index')
    pylab.xlim(yedges[0],yedges[-1])
    pylab.subplot(324)
    pylab.pcolormesh(zedges,xedges,ahizx0)
    pylab.plot([rzmin,rzmax,rzmax,rzmin,rzmin],[rxmin,rxmin,rxmax,rxmax,rxmin],'y-')
    pylab.xlim(zedges[0],zedges[-1])
    pylab.ylim(xedges[0],xedges[-1])
    pylab.xlabel('Z index')
    pylab.ylabel('X index')
    pylab.colorbar()
    pylab.subplot(325)
    pylab.bar(zedges[:-1],hiz0,lw=0,color='b')
    pylab.xlim(zedges[0],zedges[-1])
    ymin,ymax=pylab.ylim()
    pylab.plot([rzmin]*2,[ymin,ymax],'r-')
    pylab.plot([rzmax]*2,[ymin,ymax],'r-')
    pylab.xlabel('Z index')
    pylab.subplot(326)
    pylab.pcolormesh(xedges,yedges,ahixy0)
    pylab.plot([rxmin,rxmax,rxmax,rxmin,rxmin],[rymin,rymin,rymax,rymax,rymin],'y-')
    pylab.xlim(xedges[0],xedges[-1])
    pylab.ylim(yedges[0],yedges[-1])
    pylab.xlabel('X index')
    pylab.ylabel('Y index')
    pylab.colorbar()
    pylab.savefig('dvhdebug_{}.pdf'.format(debuglabel))
    pylab.savefig('dvhdebug_{}.png'.format(debuglabel))
    del fig
    for iz in range(amask.shape[0]):
        nmask = np.sum(amask[iz,:,:])
        if nmask == 0:
            logger.debug("iz={} nmask={} dosesum={}".format(iz,nmask,np.sum(aimg[iz,:,:])))
            continue
        debugname="overkill_{}_{}".format(debuglabel,iz)
        fig=pylab.figure(num="overkill",figsize=[10,10])
        fig.clf()
        pylab.subplot(221)
        pylab.pcolormesh(aimg[izmax,:,:])
        pylab.title("dose at iz={}".format(iz))
        pylab.colorbar()
        pylab.subplot(222)
        pylab.pcolormesh(amask[izmax,:,:])
        pylab.title("mask at iz={}".format(iz))
        pylab.colorbar()
        pylab.subplot(223)
        pylab.pcolormesh(aimg[izmax,:,:]*amask[izmax,:,:])
        pylab.title("dose within mask at iz={}".format(iz))
        pylab.colorbar()
        mask0 = amask[iz,:,:]*(aimg[iz,:,:]==0)
        nmask0 = np.sum(mask0)
        if nmask0 > 0:
            pylab.subplot(224)
            pylab.pcolormesh(aimg[izmax,:,:]*amask[izmax,:,:])
            pylab.title("zero dose voxels within mask at iz={}".format(iz))
            pylab.colorbar()
        pylab.savefig('{}.pdf'.format(debugname))
        pylab.savefig('{}.png'.format(debugname))
        del fig
    if oldfig:
        fig = pylab.figure(oldfig.number)
    amaskneg=amask*(aimg<0.0)
    logger.debug("got mask with {} selected voxels, {} of which have zero dose, {} are negative".format(np.sum(amask),np.sum(amask0),np.sum(amaskneg)))
    print("END debugging "+debuglabel)



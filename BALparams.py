#!/usr/bin/env python
'''
--------------------------------------------------------------------------------
Author: Jesse A. Rogerson, rogerson@yorku.ca
Additional Credits: Catherine J. Grier
                    Patrick B. Hall,
                    Daniel E. Vandenberk

This code calculates the BALnicity index of a quasar given some
user-supplied spectrum.

In the literature, the definition of the BALnicity index has gone through
many iterations. This code is meant to be a one-stop-shop for all definitons.
It allows the user to either
    a) choose from a set of predefined indexes
    b) manually decide values of some or all of the
       parameters to define their own BALnicity Index

This code has predefined:
1. BI - The original
The definition comes from Weymann et al. 1991, ApJ, 373, 23
balidx = integrate(3000,25000,(1-f(v)/0.9)*C)
f(v) = normalized flux as function of velocity from CIV center
C = 1 when 1-f(v)/0.9 is continuously > 0 for > 2000km/s,
  = 0 otherwise -- so C is just a mask

The integration starts at 3000km/s from the CIV center to avoid "associated"
systems.  It ends at 25000 to avoid the SiIV region.  This is equivalent to
masking out the < 3000km/s and > 25000km/s regions.
NB: this definition EXCLUDES all absorption prior to the
    2000km/s cut off.

2. BI_0 - defined:
A version of original BI that extends all the way to zero velocity on the
lower limit of the integral
All other rules are the same as BI
Gibson et al., 2008, ApJ, 675, 985

3.AIT - Trump et al., 2006, ApJS, 165, 1

4. AI450
The 'Absorption Index' originally defined in
--Hall et al. 2002, ApJS, 141, 267
But it was re-defined in:
--Paris et al. 2012, A&A, 548, 66

5.DI
The 'Detection Index' defined in
Paris et al. 2012, A&A, 548, 66
It is the same as the original BI but it INCLUDES all
absorption over the whole trough, not just after the
vmin value.

FOR HELP:
----------
$> ./BALparams.py -h

usage: BALparams.py [-h] [-index INDEX] [-zerr ZERR] [-lam LAM] [-vlo VLO]
                    [-vhi VHI] [-v V] [-f F] [-inc INC] [-out OUT] [-pop POP]
                    [-writeout WRITEOUT] [-smooth SMOOTH]
                    file zem

positional arguments:
  file          /path/to/deredshifted/normalized/ASCII/spectrum
  zem           redshift of target

optional arguments:
  -h, --help          show this help message and exit
  -index INDEX        The BALnicity Index to be measure (BI, BI_0, AIT, AI450)
  -zerr ZERR          error in the redshift of the target (default to 0.0)
  -lam LAM            rest wavelength of BAL ion (default to CIV)
  -vlo VLO            low velocity cut-off (default to 3000.0 km/s)
  -vhi VHI            high velocity cut-off (default to 25000.0 km/s)
  -v V                minimum continguously under f (default to 2000.0)
  -f F                amount to stay under for v (default to 0.9)
  -inc INC            Include absorption before v_min (y/n)? (default is no)
  -out OUT            The name of the output plot (default BALplot.eps)
  -pop POP            Do you want to have the plot pop-up (y/n)?(default is no)
  -writeout WRITEOUT  The name of output file (default is BAL_BI.dat
  -smooth SMOOTH      The amount of smoothing to be applied to the spectra
                      (default is no smoothing)

To RUN on DEFAULT mode:
----------
The only required inputs from the user are:
file: an ASCII spectrum of a quasar. It should be normalized
      to avoid bad pixels. It should be in the observed frame.
      (it will be deredshifted in-program)
      the first three columns must be:
            #wavelength flux flux_err
zem:  The Redshift of the quasar
(these are the 'positional arguments' from the 'help' above)

With these, the user can execute:

$> ./BALparams.py file zem

Upon executing without any other parameters specified,
BALparams.py will default to the ORIGINAL BI from
Weymann et al. 1991, ApJ, 373, 23

EXAMPLES of OTHER ways to RUN:
----------
a) There are a few pre-set BALnicity indexes that are
available in this script: BI, BI_0, AI450, DI

If the user specifies one of these via the '-index' parameter
on the commandline, it will automatically set the code to
calculate that index.

$> ./BALparams.py file zem -index AI450

b) Fully manual
The user may choose to define any (or all) of the
'optional arguments' in the 'help' section. However,
in order to not defautlt to the original BI values, the
user MUST specify '-index man,' (see below)

e.g.,

$> ./BALparams.py file zem -index man -vlo 0 -vhi 35000
--or--
$> ./BALparams.py file zem -inded man -vhi 35000 -zerr 0.056 -v 1000 -inc y

And any 'optional parameters' you do not set will remain
the default BI values.

NOTES:
----------

a) If you set '-index' to be one of the available pre-set indexes, it will
override any additional optional parameters you set.

b) The code does not validate your optional parameters. So if you define a
'-vlo' that is higher in than the defined '-vhi,' the code will NOT yell at you.

c) The code will print out the configuration of the variables each time it
runs, so you will know what you set.

HISTORY
--------------------------------------------------------------------------------
2015-05-30 - JAR - created
2015-06-01 - JAR - major changes made
2015-06-10 - JAR - condensed everything into a single BALnicity function,
                   which can adapt to whatever measurement you want
2015-06-11 - JAR - changed name to BALparams.py (from testBI.py)
2015-06-24 - JAR - fixed bug regarding index selection.
                   in order to set parameters individually users must now set:
                   '-index man'
2015-06-25 - JAR - The plotting function now colours in the parts of the
                   spectrum that meets the BALnicity conditions with a blue
                   shading.
                 - The above required removing the step where we isolate only
                   the velocity regime the users was interested in. As a result
                   the code now includes the entire spectrum array for the
                   whole calculation process. This DOES slow down the code a
                   bit. Will be the focus of improvements later
                 - added a Min Velocity of trough, just because.
2015-06-26 - JAR - was crashing if BALnicity was equal=0. added an error-catch
                   to make sure it doesn't crash when this happens
2015-07-21 - CJG - added a file output flag '-writeout'
                 - added a smoothing flag '-smooth'. boxcar smoothing if
                   requested by user.
--------------------------------------------------------------------------------
'''
import numpy as np
from sys import argv
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import copy as cp
import argparse

#global variables
lightspeed=299792.458 #km/s
civ_0=1548.202 #Ang
siv_0=1396.76 #Ang

#for validation later
yes=set(['yes','y','YES','Y',True,1])
no=set(['no','n','NO','N',False,0])
indexes=set(['BI','BI_0','AIT','AI450','DI','man'])

def BALnicity(**kwargs):
    '''
    Calculating the BALnicity Index
    '''
    #Step1 - is the BALnicity requested allowed?
    #      - crash if the submitted measurement type isn't proper
    measure=kwargs['index']
    if measure not in indexes:
        print "I don't know what to measure because '",measure,"'"
        print "is not part of my programming. Crashing."
        sys.exit()
    #      - read in spectrum given on command line
    spectrum=np.genfromtxt(kwargs['file'],usecols=(0,1,2))
    #      - set individual kwargs to variable names (easier)
    zem=kwargs['zem']
    zerr=kwargs['zerr']
    lam_0=kwargs['lam']
    vlolimit=kwargs['vlo']
    vhilimit=kwargs['vhi']
    vmin=kwargs['v']
    flim=kwargs['f']
    inc=kwargs['inc']
    smooth = kwargs['smooth']

    #Step2 - if they specified a BI, then override their options
    if measure=='BI':
        lam_0,vlolimit,vhilimit,vmin,flim,inc=civ_0,3000.0,25000.0,2000.0,0.9,'n'
    elif measure=='BI_0':
        lam_0,vlolimit,vhilimit,vmin,flim,inc=civ_0,0.0,25000.0,2000.0,0.9,'n'
    elif measure=='AIT':
        lam_0,vlolimit,vhilimit,vmin,flim,inc=civ_0,0.0,29000.0,1000.0,0.9,'y'
    elif measure=='AI450':
        lam_0,vlolimit,vhilimit,vmin,flim,inc=civ_0,0.0,25000.0,450.0,0.9,'y'
    elif measure =='DI':
        lam_0,vlolimit,vhilimit,vmin,flim,inc=civ_0,3000.0,25000.0,2000.0,0.9,'y'

    print '--------------------------------------------------------'
    print 'BALnicity measure:',measure
    print 'spectrum:',kwargs['file']
    print 'z=',zem,'+/-',zerr
    print 'BAL ion wavelength (Ang):',lam_0
    print 'BAL region:',vlolimit,'< v (km/s) <',vhilimit
    print 'must be under f(v)=',flim,'contigously for:',vmin,'km/s'
    if inc in yes:
        print 'Including absorption from before',vmin
    else:
        print 'Only including absorption AFTER', vmin
    print '-'
    #Step3 - switch to rest-frame
    spectrum[:,0]=spectrum[:,0]/(1.+zem)

    #Step4 - Convert spectrum to velocity from CIV
    spectrum=lam2vel(spectrum,civ_0)

    #Step5 - pull out ONLY the relevant velocity region
    #      - (set by vlolimit,vhilimit)
    #      - (decided it better to work with 1D lists from here on)
    #lam,flux,flux_err,vbal,dvbal=velCut(spectrum,vlolimit,vhilimit)
    #
    #Chagned this step to:
    #Step5 - individualize the spectrum into lists.
    lam=cp.deepcopy(spectrum[:,0])
    if smooth > 0:
        f2 = cp.deepcopy(spectrum[:,1])
        flux=SmoothBoxCar(f2, smooth)
        flux_err=cp.deepcopy(spectrum[:,2])
    if smooth == 0:
        flux=cp.deepcopy(spectrum[:,1])
        flux_err=cp.deepcopy(spectrum[:,2])
    vbal=cp.deepcopy(spectrum[:,3])
    dvbal=cp.deepcopy(spectrum[:,4])

    #Step6 - Calculate differential "equivalent width" (1-f(v)/0.9)
    #      - NB: you COULD change 0.9 by changing flim...
    vw=[1.0-(val/flim) for val in flux]
    ve=[(err/flim)**2 for err in flux_err]
    #Step7 - Calculate the Cvalues
    #     7a) First, the C values are calculated in the same manner as the BI.
    C=Cvalues(vbal,vlolimit,vhilimit,vmin,vw)
    #     7b) This is only done for AI450 at the moment
    #      7b Second, C' values are calculated by redoing this from high v to low.
    #         using the SAME C array, and writing over/changing it as needed
    #      - C' is 1 or 0 -- 1 if absorption < $frac*continuum contiguously
    #        for > vmin km/s; C' is set to 0 for the first vmin km/s.
    if measure=='AI450' or kwargs['inc'] in yes:
        C=Cpr(C,vbal,vlolimit,vhilimit,vmin,vw)

    #Step8 - determine begin/end lambda via C list
    #      - (this is for plotting purposes mostly)
    start,finish=lam_0,0
    finish = next((lam[i] for i,c in enumerate(C) if c==1),None)
    start = next((lam[len(lam)-1-i] for i,c in enumerate(C[::-1]) if c==1),None)
    #Step9 - Print results
    print 'Results:'
    if 1 in C:
        BI=0
        errBI=0
        for i,v in enumerate(vbal):
            BI+=vw[i]*C[i]*dvbal[i]
            errBI+=ve[i]*C[i]*dvbal[i]
        print 'BALnicity index = ',BI,'+/-',math.sqrt(errBI),' (km/s)'
        #
        # Round to nearest 100 km/s
        BIround = (100.*(int(BI / 100.) + int(2*(BI % 100.)/100.)))
        print 'BALnicity index = ',BIround,' (rounded to nearest 100km/s)'
        #
        # Find Vmax
        _v=[v for i,v in enumerate(vbal) if C[i] >0]
        vmax,vmin=max(_v),min(_v)
        verr=zerr*lightspeed
        print 'Max Velocity = ',vmax,'+/-',verr
        #
        # Round Vmax to nearest 150km/s
        vmaxround = (150.*(int(vmax / 150.) + int(2*(vmax % 150.)/150.)))
        print 'Max Vmax = ',vmaxround,'(rounded to nearest 150 km/s)'
        #
        print 'Min Velocity = ',vmin,'+/-',verr
        #
        # Round Vmax to nearest 150km/s
        vminround = (150.*(int(vmin / 150.) + int(2*(vmin % 150.)/150.)))
        print 'Max Vmax = ',vminround,'(rounded to nearest 150 km/s)'
        #
        # Find \chi^2_{trough}
        rms=math.sqrt(np.mean(np.array([flux_err[i] for i,v in enumerate(C) if v==1])**2))
        arr=[((1-flux[i])/rms)**2 for i,v in enumerate(C) if v==1]
        N=len(arr)
        chi2=(1./N)*math.fsum(arr)
        print 'The reduced Chi^2 = ',chi2
        print '(Definied in Paris et al. 2012, A&A, 548, 66)'
        #Step10 - plot results
        plotSpec(spectrum,vlolimit,vhilimit,start,finish,zem,flim,kwargs['out'],kwargs['pop'],C)
    else:
        print 'BALnicity index = 0.0 +/- 0.0'
        print 'Must not be a BAL?'
        #Setting all values to 0 if BI = 0 so it can still be output to a file
        #without crashing if there's no BAL.
        BI, errBI, BIround, vmax, vmin, verr, vmaxround, vmin, vminround, chi2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    #Write out results to a file.
    #NB: Right now this just appends the results, so if the file already exists,
    #it won't write over it-- it'll just add to it. Could be problematic.
    outfile = open(kwargs['writeout'], 'a')
    s = '%8.5f      %5.5f      %5.1f     %8.3f     %5.3f      %5.3f      %8.3f     %5.3f      %5.3f       %5.3f  \n'%(BI, math.sqrt(errBI), BIround, vmax, verr, vmaxround, vmin, verr, vminround, chi2)
    outfile.write(s)
    outfile.close()
    print '--------------------------------------------------------'
    return
    #--------------------------------------------------------

def SmoothBoxCar(x,N):
    '''
    BoxCar smoothing function, optional
    '''
    boxcar=np.ones(N)
    return convolve(x, boxcar/boxcar.sum())

def velCut(spectrum,vlolimit,vhilimit):
    '''This function pulls out the relevant velocity range
    returns 5 individual lists
    '''
    lam_orig=cp.deepcopy(spectrum[:,0])
    flux_orig=cp.deepcopy(spectrum[:,1])
    flux_err_orig=cp.deepcopy(spectrum[:,2])
    vbal_orig=cp.deepcopy(spectrum[:,3])
    dvbal_orig=cp.deepcopy(spectrum[:,4])
    lam,flux,flux_err,vbal,dvbal=[],[],[],[],[]
    #note: list.insert(0,val) will put each value at the beginning of the list
    #      this is useful to have it sorted highest/lowest
    for i,value in enumerate(vbal_orig):
        if value > vlolimit and value < vhilimit:
            lam.insert(0,lam_orig[i])
            flux.insert(0,flux_orig[i])
            flux_err.insert(0,flux_err_orig[i])
            vbal.insert(0,value)
            dvbal.insert(0,dvbal_orig[i])
    #return lam,flux,flux_err,vbal,dvbal
    return lam_orig,flux_orig,flux_err_orig,vbal_orig,dvbal_orig
    #--------------------------------------------------------

def Cpr(C_temp,vbal,vlolimit,vhilimit,vmin,vw):
    '''
    #--------------------------------------------------------
    #C' values are calculated by redoing the Cvalues from high v to low.
    #using the SAME C array, and writing over/changing it as needed
    #C' is 1 or 0 -- 1 if absorption < $frac*continuum contiguously
    #for > vmin km/s; C' is set to 0 for the first vmin km/s.
    #--------------------------------------------------------
    '''
    Cprime=C_temp
    i=len(vbal)-1
    for v in reversed(vbal):
        if v>vlolimit and v<(vhilimit-vmin):
            __vbal=np.float64(vhilimit)-vbal
            _vbal=__vbal[i]-__vbal
            _vw=[vw[j] for j,temp in enumerate(_vbal) if temp>=0 and temp<=vmin]
            if min(_vw)>0:
                Cprime[i]=1
            _vw=[vw[j] for j,temp in enumerate(_vbal) if temp>=(-0.5*vmin) and temp<=(0.5*vmin)]
            if min(_vw)>0:
                Cprime[i]=1
        i-=1
    return Cprime

def Cvalues(vbal,vlolimit,vhilimit,vmin,vw):
    '''
    #--------------------------------------------------------
    #Calculating the values of 'C'
    #C = 1 when 1-f(v)/0.9 is continuously > 0 for > v_min
    #  = 0 otherwise -- so C is just a mask
    #for the original BI sets v_min = 2000.0 km/s
    #however vmin has been changed in various publications
    #this function will take whatever vmin you give it
    #--------------------------------------------------------
    '''
    C=(np.zeros(len(vbal),int)).tolist() #starting with all C's =0
    for i,v in enumerate(vbal):
        if v>(vlolimit+vmin) and v<vhilimit:
            #v acts as the reference point, sliding along vbal
            _vbal=v-vbal #this works because im using lists (not numpy)
            #if it's 0<v<2000.0, then pull out the 'equivalent width' values
            _vw=[vw[j] for j,temp in enumerate(_vbal) if temp>=0 and temp<=vmin]
            #if all of the values of _vw are > 1, that means for 2000 km/s
            #before this current v, all values of f(v) were below 0.9
            #therefore set C to 1
            if min(_vw)>0:
                C[i]=1
    return C
    #--------------------------------------------------------

def lam2vel(spec,rest=civ_0):
    '''
    #--------------------------------------------------------
    #Function takes a spectrum (lam,flux,fluxerr) and returns
    #a new spectrum array with two appended columns:
    #velocity, and dv
    #velocity - measured realtive to some rest lambda
    #dv       - the pixel width in velocity space
    #
    #If a rest wavelength of the BAL ion is not provided
    #it assumes rest=civ_0
    #
    #Acknowledgements:
    #this code is based off Kate Grier
    #
    #See Hall et al. 2002,ApJS,141,267 for equation reference:
    #v = c*dz/(1+z) = c*(z_CIV - wl_BAL/wl_rest + 1)/(1+z_CIV)
    #--------------------------------------------------------
    '''
    vbal=np.zeros(len(spec[:,0]))
    dvbal=np.zeros(len(spec[:,0]))
    for i,w in enumerate(spec[:,0]):
        minlambda=(rest-rest)/rest
        maxlambda=(w-rest)/rest
        Rmin=1./(1+minlambda)
        Rmax=1./(1+maxlambda)
        vel_min=-lightspeed*(Rmin**2-1)/(Rmin**2+1)
        vel_max=-lightspeed*(Rmax**2-1)/(Rmax**2+1)
        vbal[i]=(vel_min-vel_max)
        #calculating dv infinitesimal (just bin-width)
        if i>0 and i<len(spec[:,0]):
            dvbal[i]=math.fabs(vbal[i]-vbal[i-1])
    #in order to concatenate must reshape
    spec=np.concatenate((spec,np.reshape(vbal,(len(vbal),1)),
                    np.reshape(dvbal,(len(dvbal),1))),axis=1)
    #return the original spectrum, with two more columns
    return spec
    #--------------------------------------------------------

def plotSpec(spec,vlolimit,vhilimit,s,f,zem,flim,filename,p,C):
    '''Plotting the spectrum along with the window of BALnicity'''
    #xlimits=[s-500,f+500]
    xlimits=[1200,1600]
    ylimits=[0,2]
    fig = plt.figure()
    ax1=fig.add_subplot(111)
    plt.rc('text',usetex=True)
    plt.rc('font',family='sans-serif')

    #split into two x-axes (for rest vs. observed)
    ax1.plot(spec[:,0],spec[:,1],'k',linewidth=1.0)
    #plotting the continuum at 1.0
    ax1.plot((1000,2000),(1,1),'k')
    #ax1.annotate('cont',xy=(((f-s)/2)+s,1.1),xytext=(((f-s)/2)+s,1.05))
    #plotting the 0.9 limit
    ax1.plot((1000,2000),(flim,flim),'k--')
    #ax1.annotate('0.9',xy=(s-20,0.8),xytext=(s-20,0.8))

    #vertical lines on BAL limits
    ax1.plot((s,s),(0,4),'r--',linewidth=1.0)
    ax1.annotate('stop',xy=(s-8,0.4),color='r',xytext=(s-8,0.4))
    ax1.plot((f,f),(0,4),'r--',linewidth=1.0)
    ax1.annotate('start',xy=(f+5,0.4),color='r',xytext=(f+5,0.4))
    ax1.plot((vlolimit,vlolimit),(0,4),'--')
    ax1.plot((vhilimit,vhilimit),(0,4),'--')
    #plt.plot((civ_0,civ_0),(0,4),'k--')
    C=np.array(C)
    #ax1.fill_between(spec[:,0],spec[:,1],0.9, where=(C==1))
    ax1.fill_between(spec[:,0],spec[:,1],0.9, where=(C==1))

    ax1.set_xlim(xlimits[0],xlimits[1])
    ax1.set_ylim(ylimits[0],ylimits[1])
    ax1.set_ylabel('Normalized Flux Density')
    ax1.set_xlabel('Rest-frame Wavelength \AA')
    ax1.set_yticks((0.0,0.5,1.0,1.5,2.0,flim))

    #2nd axis
    ax2=ax1.twiny()
    ax2.set_xlim(xlimits[0]*(1+zem),xlimits[1]*(1+zem))
    ax2.set_xticks([4400,4600,4800,5000,5200,5400])
    ax2.set_xlabel('Observed Wavelength (\AA)')
    ax2.xaxis.set_minor_locator(MultipleLocator(50))

    plt.savefig(filename)
    print 'Saved a figure for you called:',filename
    if p in yes:
        print '...Pushing figure to X11, popping up...'
        plt.show()
    return
    #--------------------------------------------------------

#-------------------------------------------------------------------------------
#read in from command line
parser=argparse.ArgumentParser()
parser.add_argument('file', type=str, help='/path/to/deredshifted/normalized/ASCII/spectrum')
parser.add_argument('zem', type=float, help='redshift of target')
parser.add_argument('-index', type=str, default='BI', help='The BALnicity Index to be measure (BI, BI_0, AIT, AI450)')
parser.add_argument('-zerr', type=float, default=0.0, help='error in the redshift of the target (default to 0.0)')
parser.add_argument('-lam', type=float, default=civ_0, help='rest wavelength of BAL ion (default to CIV)')
parser.add_argument('-vlo', type=float, default=3000.0, help='low velocity cut-off (default to 3000.0 km/s)')
parser.add_argument('-vhi', type=float, default=25000.0, help='high velocity cut-off (default to 25000.0 km/s)')
parser.add_argument('-v', type=float, default=2000.0, help='minimum continguously under f (default to 2000.0)')
parser.add_argument('-f', type=float, default=0.9, help='amount to stay under for v (default to 0.9)')
parser.add_argument('-inc', type=str, default='n', help='Include absorption before v_min (y/n)? (default is no) ')
parser.add_argument('-out', type=str, default='BALplot.eps', help='The name of the output plot (default BALplot.eps)')
parser.add_argument('-pop', type=str, default='n', help='Do you want to have the plot pop-up (y/n)?(default is no)')
parser.add_argument('-writeout', type=str, default='BAL_BI.dat', help='The name of output file (default is BAL_BI.dat')
parser.add_argument('-smooth', type=int, default = 0, help='The amount of smoothing to be applied to the spectra (default is no smoothing)')
kwargs=vars(parser.parse_args())

BALnicity(**kwargs)

#!/usr/bin/env python
'''
--------------------------------------------------------------------------------
Author: Jesse A. Rogerson, jesserogerson.com, rogerson@yorku.ca
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
2015-06-12 - JAR - original commit to github.com/jesserogerson
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
2015-09-02 - JAR - fixed error found in def Cpr(): (reported by CJG)
2015-09-03 - CJG - added 'multiple trough' functionality
                 - now reports how many troughs there are, what their individual
                   BALnicities are. Still reports the total as well.
2015-09-03 - JAR - moved the actual BAL calc to a new function called 'BALcalc'
                 - moved smoothing into the lam2vel() function
2015-09-04 - JAR - modified 'multiple troughs' functionality
                 - created new function troughRegions()
                   it accepts a C array of 1's and 0's and tells you where
                   the troughs are.
2015-10-06 - JAR - added a spectra idenfifier (gem/sdss/boss) to the writeout
                 - augmented current label in writeout
                 - added lambda ranges for BALs, helsp for EW later
2015-10-07 - JAR - fixed plotting issue. Now program creates individual BAL
                   plots for each 'sdss','boss', or 'gem'
2015-10-22 - JAR - forced the CIV global to 1545... will catch associated civ
                 - added a new function EWcalc(). It calculates EW, centroid
                   velocity, and average depth, of any trough identified.
                 - removed the CIV global change to 1545. don't like it.
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
import string
import os
from scipy.ndimage.filters import convolve

#global variables
lightspeed=299792.458 #km/s
civ_0=1548.202 #Ang
#civ_0=1545.000 #Ang - this will allow us to catch 'associated' systems
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
    #find which survey spectrum is from
    typ=kwargs['file'].split('.')[1]
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
    objName=os.path.splitext(kwargs['file'])[0]

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
        print 'Only including absorption AFTER',vmin
    print '-'

    #Step3 - switch to rest-frame
    spectrum[:,0]=spectrum[:,0]/(1.+zem)

    #Step4 - Convert spectrum to velocity from CIV, smooth if requested
    spectrum=lam2vel(spectrum,smooth,civ_0)

    #Step5 - pull out ONLY the relevant velocity region
    #      - (set by vlolimit,vhilimit)
    #      - (decided it better to work with 1D lists from here on)

    #Step5 - individualize the spectrum into lists.
    lam=cp.deepcopy(spectrum[:,0])
    flux=cp.deepcopy(spectrum[:,1])
    flux_err=cp.deepcopy(spectrum[:,2])
    vbal=cp.deepcopy(spectrum[:,3])
    dvbal=cp.deepcopy(spectrum[:,4])

    #Step6 - Calculate differential "equivalent width" (1-f(v)/0.9)
    vw=[1.0-(val/flim) for val in flux] #allowable range: 0 < vw < 1.0
    ve=[(err/flim)**2 for err in flux_err]

    #Step7 - Calculate the Cvalues
    C,regions=Cvalues(vbal,vlolimit,vhilimit,vmin,vw) #first pass
    #redo the C array if we're including stuff BEFORE vmin cutoff
    if measure=='AI450' or kwargs['inc'] in yes:
        C,regions=Cpr(C,vbal,vlolimit,vhilimit,vmin,vw)

    #Step8 - Determine how many troughs there are!
    #      - and make a trough dictionary
    numTroughs= len(regions) #how many troughs?
    troughLimits=[vbal[reg].tolist() for reg in regions] #finds the values in vbal of those differences
    print 'Number of BAL troughs found: ', numTroughs
    troughDict={} #dictionary helps keep of track of stuff
    letterList=list(string.ascii_uppercase)
    count=0
    for t in troughLimits:
        troughDict[letterList[count]]=t
        count+=1

    #Step9 - Print Results
    if numTroughs>0:
        for t in troughDict:
            # Determine begin/end velocity (this is for plotting purposes mostly)
            start,finish=troughDict[t][0], troughDict[t][1]
            #Pull out the region of this specific trough
            region = np.where(np.logical_and(vbal <= start, vbal >= finish))
            vbalNew = vbal[region]
            CNew = np.ones(np.size(region))
            vwNew = np.array(vw)[region]
            dvbalNew = np.array(dvbal)[region]
            veNew = np.array(ve)[region]
            lamNew = np.array(lam)[region]
            fluxNew = np.array(flux)[region]
            flux_errNew = np.array(flux_err)[region]

            #use BALcalc() to find all values of interest
            label=objName+' '+t
            BI,errBI,BIround,vmax,verr,vmaxround,vmin,verr,vminround,lamMin,lamMax,chi2=BALcalc(label,vbalNew,CNew,vwNew,dvbalNew,veNew,lamNew,fluxNew,flux_errNew,zerr)
            label=objName[4:]
            EW,sigmaEW, dmax7, sigmaDmax7, v_cent, dBAL, sigmadBAL=EWcalc(label,vmax,vmin,lamMin,lamMax,lam,flux,flux_err)
            print '*Writing results for trough',t,'to file'
            #Write out results to a file.
            #NB: Right now this just appends the results, so if the file already exists,
            #it won't write over it-- it'll just add to it. Could be problematic.
            outfile = open(kwargs['writeout'], 'a')
            s = '%s     %s      %s      %8.5f      %5.5f      %5.1f     %8.3f     %5.3f      %5.3f      %8.3f     %5.3f      %5.3f       %8.3f      %8.3f  %5.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f     \n'%(label,typ,t, BI, errBI, BIround, vmax, verr, vmaxround, vmin, verr, vminround, lamMin,lamMax,chi2,EW,sigmaEW, dmax7, sigmaDmax7, v_cent, dBAL, sigmadBAL)
            outfile.write(s)
            outfile.close()

        #do BALcalc() on the WHOLE vel range as well
        label=objName[4:]
        BI,errBI,BIround,vmax,verr,vmaxround,vmin,verr,vminround,lamMin,lamMax,chi2=BALcalc(label,vbal,C,vw,dvbal,ve,lam,flux,flux_err,zerr)
        print '*Writing results for TOTAL BAL to file'
        outfile = open('TOTAL_'+kwargs['writeout'], 'a')
        s = '%s      %s     %i      %8.5f      %5.5f      %5.1f     %8.3f     %5.3f      %5.3f      %8.3f     %5.3f      %5.3f      %8.3f      %8.3f     %5.3f  \n'%(label,typ,numTroughs, BI, errBI, BIround, vmax, verr, vmaxround, vmin, verr, vminround, lamMin, lamMax, chi2)
        outfile.write(s)
        outfile.close()
        #Step10 - plot results
        plotSpec(spectrum,vlolimit,vhilimit,zem,flim,(kwargs['out'].split('.')[0]+'_'+typ+'.'+kwargs['out'].split('.')[1]),kwargs['pop'],C,troughDict)

    else:#the case where no troughs were found
        print 'BALnicity index = 0.0 +/- 0.0'
        print 'Must not be a BAL?'
        #Setting all values to 0 if BI = 0 so it can still be output to a file
        #without crashing if there's no BAL.
        BI,errBI,BIround,vmax,vmin,verr,vmaxround,vmin,vminround,lamMin,lamMax,chi2,EW,sigmaEW, dmax7, sigmaDmax7, v_cent, dBAL, sigmadBAL = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0
        #Write out results to a file.
        #NB: Right now this just appends the results, so if the file already exists,
        #it won't write over it-- it'll just add to it. Could be problematic.
        label=objName[4:]
        outfile = open(kwargs['writeout'], 'a')
        s = '%s      %s     %i      %8.5f      %5.5f      %5.1f     %8.3f     %5.3f      %5.3f      %8.3f     %5.3f      %5.3f      %8.3f   %8.3f    %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f  \n'%(label,typ,numTroughs, BI, errBI, BIround, vmax, verr, vmaxround, vmin, verr, vminround, lamMin, lamMax, chi2,EW,sigmaEW, dmax7, sigmaDmax7, v_cent, dBAL, sigmadBAL)
        outfile.write(s)
        outfile.close()
    print '--------------------------------------------------------'
    return
    #--------------------------------------------------------

def EWcalc(name,vmax,vmin,lamMin,lamMax,lam,flux,flux_err):
    '''
    calculating various parameters: EW, average depth, etc.
    '''
    #BI,errBI,BIround,vmax,verr,vmaxround,vmin,verr,vminround,lamMin,lamMax,chi2
    #what are the RLFs? will be used for
    limits=[lamMin,lamMax]
    parmFile='norm'+name+'.parm'
    lowest=10.
    with open(parmFile,'r') as f:
        s=f.readlines()[-2]
    f.close()
    line=s.strip().split('=')[1]
    temp=map(float,line.split(','))
    RLF=[[temp[0],temp[1]]]
    for t in range(2,len(temp)-1,2):
        RLF.insert(0,[temp[t],temp[t+1]])
    #identify the indicies that reflect the given RLF windows
    w=0
    w=np.array([],dtype=int)
    for bounds in RLF:
        temp_w=[index for index,value in enumerate(lam) if value > bounds[0] and value < bounds[1]]
        w=np.concatenate((w,temp_w))
        temp_w=[]
    print '------------------------------------------------------------'
    print 'Measureing Equivalent Width'
    print 'Trough Range:'+str(limits)
    #print 'Corresponding velocity range: '+str(round(high_v,0))+' < v < '+str(round(low_v,0))
    print 'Corresponding velocity range: '+str(round(vmax,0))+' < v < '+str(round(vmin,0))
    print 'Corresponding velocity width: '+str(round(vmin-vmax,0))
    print 'Step1: Measure the uncertainty in the w continuum flux'
    print '       in the adjacent normalization windows (\DeltaFc)'
    print 'Normalization windows:',RLF
    temp_f=0.
    temp_ferr=0.
    number=len(w)
    #Calculating the Mean and Uncertainty in the Mean
    #using the Normalization windows
    for i in w:
        temp_f=temp_f+flux[i]#add up flux in each bin
        temp_ferr=temp_ferr+(flux_err[i]**2)#add up uncertainy squared
        number=number+1
    meanFc=temp_f/number#calculate the mean
    deltaFc=math.sqrt(temp_ferr/(number-1))#RMS noise of the distribution
    deltaFc=deltaFc/math.sqrt(number)#the uncertainty of the mean
    print '\Delta F_c = '+str(round(deltaFc,2))
    print 'Step2: Measure EW and uncertainy on EW'
    temp=0.
    tempA=0.
    tempB=0.
    sigma_low=0.
    low_avg=7
    N_v=0.
    vCent=0.
    vCentWeight=0.
    dBAL=0.
    dBALdev=0.
    pixWidth=0
    cont=1.0
    #Calculating values such as EW, lowest point, N_v etc.
    #the actual measurement loop
    for i in range(len(lam)):
        if lam[i]>=limits[0] and lam[i]<=limits[1]:
            #Measure of dBAL, average depth from normalized continuum
            dBAL+=(1-flux[i])
            pixWidth+=1
            #Measure centroid velocity, weighted, cumulative
            lambdatest=(lam[i]-civ_0)/civ_0
            Rtest=1./(1+lambdatest)
            veltest=lightspeed*(Rtest**2-1)/(Rtest**2+1)
            vCent+=veltest*(1-flux[i])
            vCentWeight+=1-flux[i]
            veltest,lambdatest,Rtest=0,0,0
            #Measure N_v (for column density)
            #N_v=N_v+(-1.*math.log(flux[i],math.e))
            #EWmeasure
            temp=temp+(1.0-(flux[i]/cont))*(lam[i]-lam[i-1])
            #delta EW measure
            tempA=tempA+((flux[i]/cont)*(lam[i]-lam[i-1])) #continuum
            tempB=tempB+((lam[i]-lam[i-1])*(flux_err[i]/cont))**2 #fluxerror
            #pulling out the lowest value and RMS
            if flux[i] <= lowest:
                #test_newlowest=(flux[i-3]+flux[i-2]+flux[i-1]+flux[i]+flux[i+1]+flux[i+2]+flux[i+3])/7.
                #print test_newlowest
                newlowest=np.mean(flux[i-low_avg:i+low_avg+1])
            if newlowest <= lowest:
                lowest=cp.copy(newlowest) ###
                sigma_low=math.sqrt(np.sum((flux_err[i-low_avg:i+low_avg+1])**2)/(low_avg-1))
                #sigma_low=math.sqrt((flux_err[i-3]**2+flux_err[i-2]**2+flux_err[i-1]**2
                #                     +flux_err[i]**2+flux_err[i+1]**2+flux_err[i+2]**2+flux_err[i+3]**2)/6)
                sigma_low=sigma_low/math.sqrt(low_avg)
    for i in range(len(lam)):
        if lam[i]>=limits[0] and lam[i]<=limits[1]:
            dBALdev+=((1-flux[i])-(dBAL/pixWidth))**2
    dBALdev=math.sqrt(dBALdev/pixWidth) #Standard Deviation of dBAL
    EW=temp
    tempA=(tempA*(deltaFc/cont))**2
    deltaEW=math.sqrt(tempA+tempB)
    #print 'N_CVI ='+str(round((3.7679e+14/(1549.055*0.286))*N_v))
    #print 'Equivalent Width = '+str(round(EW,3))
    #print 'Uncertainty continuum = '+str(round(tempA,3))
    #print 'Uncertainty flux = '+str(round(tempB,3))
    #print 'Total propagated uncertainty = '+str(round(deltaEW,3))
    print ''
    print 'Maximum Depth (dmax7, avg over 7 pixels):'+str(round(1-lowest,3))+' +/- '+str(round(sigma_low,3))
    print 'The Centroid Velocity is:'+str(round((vCent/vCentWeight),3))
    print 'Average BAL trough depth (dBAL):'+str(round(dBAL/pixWidth,3))+' +/- '+str(round((dBALdev/math.sqrt(pixWidth)),3))
    print 'RETURN: '+str(round(EW,3))+' +/- '+str(round(deltaEW,3))
    print '------------------------------------------------------------'
    #return EW,sigmaEW, dmax7, sigmaDmax7, v_cent, dBAL, sigmadBAL
    return EW,deltaEW,(1-lowest),sigma_low,(vCent/vCentWeight),(dBAL/pixWidth),(dBALdev/math.sqrt(pixWidth))
    #--------------------------------------------------------

def BALcalc(name,vbal,C,vw,dvbal,ve,lam,flux,flux_err,zerr):
    '''
    This is where the actual BI and other associated
    values are calculated.
    '''
    print '-----'
    print 'Here are the results for trough ',name
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
    # Find Vmax,Vmin
    _v=[v for i,v in enumerate(vbal) if C[i] >0]
    vmax,vmin=max(_v),min(_v)
    verr=zerr*lightspeed
    print 'Max Velocity = ',vmax,'+/-',verr
    vmaxround = (150.*(int(vmax / 150.) + int(2*(vmax % 150.)/150.)))
    print 'Min Velocity = ',vmin,'+/-',verr
    vminround = (150.*(int(vmin / 150.) + int(2*(vmin % 150.)/150.)))
    # Find lamMax, lamMin
    _l=[l for i,l in enumerate(lam) if C[i] > 0]
    lamMin,lamMax=min(_l),max(_l)
    #
    # Find \chi^2_{trough}
    rms=math.sqrt(np.mean(np.array([flux_err[i] for i,v in enumerate(C) if v==1])**2))
    arr=[((1-flux[i])/rms)**2 for i,v in enumerate(C) if v==1]
    N=len(arr)
    chi2=(1./N)*math.fsum(arr)
    print 'The reduced Chi^2 = ',chi2
    print '(Definied in Paris et al. 2012, A&A, 548, 66)'
    print '-----'
    return BI, math.sqrt(errBI), BIround, vmax, verr, vmaxround, vmin, verr, vminround, lamMin, lamMax, chi2
    #--------------------------------------------------------

def smoothBoxCar(x,N):
    '''
    BoxCar smoothing function, optional
    '''
    boxcar=np.ones(N)
    return convolve(x, boxcar/boxcar.sum())
    #--------------------------------------------------------

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
    #works on the premise:
    #if the value of the Carray was 1 previously, and the next
    #continuum value is below the limit, then the next value in
    #the Carray should be a 1 one. If it isn't... change it.
    Cprime=cp.deepcopy(C_temp)
    for i,v in enumerate(vbal):
        if v>vlolimit and v<vhilimit:
            if Cprime[i-1]==1 and vw[i]>0.0 and vw[i]<1.0:
                Cprime[i]=1
    return Cprime,troughRegions(Cprime)
    #--------------------------------------------------------

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
    return C,troughRegions(C)
    #--------------------------------------------------------

def troughRegions(C):
    '''
    Finds the sections of the C array that are 1's, returns
    a 2D-list, where each entry is the start/stop of a trough
    '''
    start,stop=0,0
    s=False
    regions=[]
    for i,c in enumerate(C):
        if c==1 and s==False:
            start=i
            s=True
        elif c==0 and s==True:
            stop=i-1
            s=False
            regions.append([start,stop])
        if s==True and i==len(C)-1:
            stop=i
            regions.append([start,stop])
    return regions
    #--------------------------------------------------------

def lam2vel(spec,smooth,rest=civ_0):
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
    #smooth if asked for
    if smooth > 0:
        spec[:,1]=smoothBoxCar(spec[:,1], smooth)
    #return the original spectrum, with two more columns
    return spec
    #--------------------------------------------------------

def plotSpec(spec,vlolimit,vhilimit,zem,flim,filename,p,C,troughDict):
    '''Plotting the spectrum along with the window of BALnicity'''
    print filename
    sys.exit()
    plt.rc('text',usetex=True)
    plt.rc('font',family='sans-serif')
    fig = plt.figure()
    ax1=fig.add_subplot(1,1,1)
    xlimits=[1200,1600]
    ylimits=[0,2]

    #top axis - or axis 1
    ax1.set_xlim(xlimits[0],xlimits[1])
    ax1.set_ylim(ylimits[0],ylimits[1])
    ax1.set_ylabel('Normalized Flux Density')
    ax1.set_xlabel('Rest-frame Wavelength \AA')
    #ax1.set_yticks((0.0,0.5,1.0,1.5,2.0,flim))

    ax1.plot(spec[:,0],spec[:,1],'k',linewidth=1.0)
    #plotting the continuum at 1.0
    ax1.plot((1000,2000),(1,1),'k')
    #plotting the 0.9 limit
    ax1.plot((1000,2000),(flim,flim),'k--')
    C=np.array(C)
    ax1.fill_between(spec[:,0],spec[:,1],flim, where=(C==1))
    ax1.annotate(str(filename),xy=(1275,(ylimits[1]*0.1)))

    #calculating analogus x boundary
    minlambda=(xlimits[0]-civ_0)/civ_0
    maxlambda=(xlimits[1]-civ_0)/civ_0
    Rmin=1./(1+minlambda)
    Rmax=1./(1+maxlambda)
    vel_min=lightspeed*(Rmin**2-1)/(Rmin**2+1)
    vel_max=lightspeed*(Rmax**2-1)/(Rmax**2+1)

    #2nd axis

    ax2=ax1.twiny()
    ax2.set_xlim(vel_min,vel_max)
    ax2.set_ylim(ylimits[0],ylimits[1])
    ax2.set_xlabel('velocity (km/s)')

    #ERROR - this doesn't seem to be plotting properly
    #        another problem for another day
    for t in troughDict:
        ax2.annotate(t, xy=(troughDict[t][0],1.5),xytext=(troughDict[t][0],1.5))
    #    #print troughDict[t]
    #    ax2.plot([np.array(troughDict[t][0]),np.array(troughDict[t][0])], ylimits, color = 'r', marker = 'o')
    #    ax2.plot([np.array(troughDict[t][1]),np.array(troughDict[t][1])], ylimits, color = 'r', marker = 'o')

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
parser.add_argument('file', type=str, help='/path/to/normalized/ASCII/spectrum')
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

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import os
import glob
import ROOT
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/Users/kath/icecube2019/pytools')  ## So it can find things without links
import livetime
import hdfchain
import ic86_geometry


# In[2]:


## Load the data
mfiles11 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2011/*.h5')
mfiles12 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2012/*.h5')
mfiles13 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2013/*.h5')
mfiles14 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2014/*.h5')
mfiles15 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2015/*.h5')
mfiles16 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2016/*.h5')
mfiles17 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2017/*.h5')
mfiles18 = glob.glob('/Users/kath/data/testruns_L3reprocess/IC86.2018/*.h5')
mfiles11.sort()
mfiles12.sort()
mfiles13.sort()
mfiles14.sort()
mfiles15.sort()
mfiles16.sort()
mfiles17.sort()
mfiles18.sort()

'''
mfiles11 = mfiles11[0:10]    #<---- COMMENT THIS IN FOR TESTING!
mfiles12 = mfiles12[0:10]   
mfiles13 = mfiles13[0:10]   
mfiles14 = mfiles14[0:10]   
mfiles15 = mfiles15[0:10]   
mfiles16 = mfiles16[0:10]   
mfiles17 = mfiles17[0:10]   
mfiles18 = mfiles18[0:10]   
'''

t11 = hdfchain.HDFChain(mfiles11)
t12 = hdfchain.HDFChain(mfiles12)
t13 = hdfchain.HDFChain(mfiles13)
t14 = hdfchain.HDFChain(mfiles14)
t15 = hdfchain.HDFChain(mfiles15)
t16 = hdfchain.HDFChain(mfiles16)
t17 = hdfchain.HDFChain(mfiles17)
t18 = hdfchain.HDFChain(mfiles18)

years = [2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]
#years = [2012, 2013, 2014, 2015, 2016, 2017, 2018]

trees = {2011:t11, 2012:t12, 2013:t13, 2014:t14, 2015:t15, 2016:t16, 2017:t17, 2018:t18}
colors = {2011:'black', 2012:'red', 2013:'orange', 2014:'yellow', 2015:'green', 2016:'blue', 2017:'violet', 2018:'magenta'}


# In[3]:


## Livetime, Nevents
nev = {}
lt = {}
for y in years:
    t = trees[y]
    
    lt[y] = livetime.livetime_from_deltaT(trees[y],10,0)[0]
    print("Livetime of %d = %f"%(y, lt[y]))

## Fit success and other common things we'll need over and over
lap_status = {}
prescale = {}
stdfilter = {}
iffilter = {}
st2filter = {}
for y in years:
    t = trees[y]
    lap_status[y] = t.get_node("/Laputop").col('fit_status')[:]
    if y>2011:
        prescale[y] = t.get_node("/IceTop_EventPrescale").col('value')[:]
    else:
        prescale[y] = np.ones_like(lap_status[y])  ##<--- because 2011 is broken
    stdfilter[y] = t.get_node("/IceTop_StandardFilter").col('value')[:]
    iffilter[y] = t.get_node("/IceTop_InFillFilter").col('value')[:]
    if y>2015:
        st2filter[y] = t.get_node("/IceTop_TwoStationFilter").col('value')[:]
    else:
        st2filter[y] = np.zeros_like(lap_status[y])




    


# In[4]:


## PULSE TIMING: IceTopHLCSeedRTPulses
for y in years:
    print("Processing ", str(y))
    t = trees[y]
    recocut = (lap_status[y]==0)*(stdfilter[y]==1)
    
    pname = "IceTopHLCSeedRTPulses"
    #pname = "CleanedHLCTankPulses"


    ptime = t.get_node("/"+pname).col('time')[:]
    print(ptime)
    w = np.ones_like(ptime)/lt[y]

    ## Three different time scales (wide, medium, zoom-in)
    plt.figure(1)
    plt.hist(ptime,label=y, weights=w,histtype='step',color=colors[y],bins=100,range=(0,1e7),log=True)
    plt.xlabel(pname+" time (ns)")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_time1.png")

    plt.figure(2)
    plt.hist(ptime,label=y, weights=w,histtype='step',color=colors[y],bins=100,range=(0,1e5),log=True)
    plt.xlabel(pname+" time (ns)")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_time2.png")

    plt.figure(3)
    plt.hist(ptime,label=y, weights=w,histtype='step',color=colors[y],bins=100,range=(9000,15000),log=True)
    plt.xlabel(pname+" time (ns)")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_time3.png")

    ## NOTE: All years have times that go past 1e6...


# In[5]:


## PULSES CHARGE: IceTopHLCSeedRTPulses
for y in years:
    print("Processing ", str(y))
    t = trees[y]
    recocut = (lap_status[y]==0)*(stdfilter[y]==1)

    pname = "IceTopHLCSeedRTPulses"
    #pname = "CleanedHLCTankPulses"

    pcharge = t.get_node("/"+pname).col('charge')[:]
    print(pcharge)
    w = np.ones_like(pcharge)/lt[y]

    ## How many of them are NAN?
    #plt.figure(1)
    nans = np.isnan(pcharge)
    nnans = len(pcharge[nans])
    print("%d: NaNs are %d out of %d (%f percent)"%(y, nnans, len(nans), nnans/len(nans)*100))
    
    ## Turn the NAN's into an actual number, so that they show up on the plot
    pcharge[nans] = 0.001
    
    plt.hist(np.log10(pcharge),label=y, weights=w,histtype='step',color=colors[y],bins=100,range=(-3.0,4.0),log=True)
    plt.xlabel(pname+" charge (VEM)")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_charge.png")




    


# In[6]:


## PULSES WIDTH: IceTopHLCSeedRTPulses
for y in years:
    print("Processing ", str(y))
    t = trees[y]
    recocut = (lap_status[y]==0)*(stdfilter[y]==1)

    pname = "IceTopHLCSeedRTPulses"
    #pname = "CleanedHLCTankPulses"

    pwidth = t.get_node("/"+pname).col('width')[:]
    print(pwidth)
    w = np.ones_like(pwidth)/lt[y]
    plt.hist(pwidth,label=y, weights=w,histtype='step',color=colors[y],bins=100,range=(0,600),log=True)
    plt.xlabel(pname+" width (ns)")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_width.png")


    ## NOTE: looks like 2017 or 2018 has a few widths > 1000



# In[7]:


## PULSES STR/OM/PMT: IceTopHLCSeedRTPulses
for y in years:
    print("Processing ", str(y))
    t = trees[y]
    recocut = (lap_status[y]==0)*(stdfilter[y]==1)

    pname = "IceTopHLCSeedRTPulses"
    #pname = "CleanedHLCTankPulses"

    pstr = t.get_node("/"+pname).col('string')[:]
    pom = t.get_node("/"+pname).col('om')[:]
    ppmt = t.get_node("/"+pname).col('pmt')[:]

    plt.figure(1)
    w = np.ones_like(pstr)/lt[y]
    plt.hist(pstr,label=y, weights=w,histtype='step',color=colors[y],bins=82,range=(0,82),log=True)
    plt.xlabel(pname+" String number")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_stringnr.png")


    plt.figure(2)
    plt.hist(pom,label=y, weights=w,histtype='step',color=colors[y],bins=4,range=(61,65),log=True)
    plt.xlabel(pname+" OM number")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_omnr.png")


    #plt.figure(3)
    #plt.hist(pom,label=y, weights=w,histtype='step',color=colors[y],log=True)
    #plt.xlabel(pname+" PMT number")
    #plt.ylabel("Nevents (norm by livetime)")
    #plt.legend(loc='upper right')



# In[8]:


## PULSES INFO: IceTopHLCSeedRTPulses
for y in years:
#for y in [2016,2017,2018]:
    print("Processing ", str(y))
    t = trees[y]
    recocut = (lap_status[y]==0)*(stdfilter[y]==1)
    
    cut = np.ones_like(stdfilter[y],dtype=bool)   ## Just plot everything
    ## What cut to apply?
    #cut = st2filter[y]==0
    
    pname = "IceTopHLCSeedRTPulses"
    #pname = "CleanedHLCTankPulses"

    pnchan = t.get_node("/"+pname+"_info").col('nchan')[:]
    pnchan1 = t.get_node("/"+pname+"_info").col('nchan_1hit')[:]
    pnstr = t.get_node("/"+pname+"_info").col('nstrings')[:]
    pnhit = t.get_node("/"+pname+"_info").col('nhit')[:]
    pqtot = t.get_node("/"+pname+"_info").col('tot_charge')[:]
    pt1 = t.get_node("/"+pname+"_info").col('first_time')[:]
    plength = t.get_node("/"+pname+"_info").col('length')[:]

    #cut = pnchan>0   ## Non-empty pulse-series only
    
    ## How much to zoom in?
    nstrmax = 10
    ntankmax = nstrmax*2
    
    plt.figure(1)
    w = np.ones_like(pnchan)/lt[y]
    plt.hist(pnchan[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=ntankmax+1,range=(0,ntankmax+1),log=True)
    plt.xlabel(pname+" Nch")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_Nch.png")

    plt.figure(2)
    plt.hist(pnchan1[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=ntankmax+1,range=(0,ntankmax+1),log=True)
    plt.xlabel(pname+" Nch 1st hit")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_Nch1.png")

    plt.figure(3)
    plt.hist(pnstr[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=nstrmax+1,range=(0,nstrmax+1),log=True)
    plt.xlabel(pname+" Nstr")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_Nstr.png")

    plt.figure(4)
    plt.hist(pnhit[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=ntankmax+1,range=(0,ntankmax+1),log=True)
    plt.xlabel(pname+" Nhit")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_Nhit.png")

    plt.figure(5)
    ## Some of these are "NAN" too!  Because one of the pulses is.
    nans = np.isnan(pqtot)
    nnans = len(pqtot[nans])
    print("%d: Qtot NaNs are %d out of %d (%f percent)"%(y, nnans, len(nans), nnans/len(nans)*100))
    ## Turn the NAN's into an actual number, so that they show up on the plot
    pqtot[nans] = 0.1
    plt.hist(np.log10(pqtot[cut]),label=y, weights=w[cut],histtype='step',color=colors[y],bins=100,range=(-1.0,4.5),log=True)
    plt.xlabel(pname+" Qtotal")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_Qtot.png")

    plt.figure(6)
    plt.hist(pt1[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=100,range=(0,1e5),log=True)
    plt.xlabel(pname+" First time")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_1sttime.png")

    plt.figure(7)
    plt.hist(plength[cut],label=y, weights=w[cut],histtype='step',color=colors[y],bins=100,range=(0,4000),log=True)
    plt.xlabel(pname+" Length between 1st and last hit")
    plt.ylabel("Nevents (norm by livetime)")
    plt.legend(loc='upper right')
    plt.savefig("plots/yoy_"+pname+"_length.png")


 


# In[9]:


t11.close()
t12.close()
t13.close()
t14.close()
t15.close()
t16.close()
t17.close()
t18.close()


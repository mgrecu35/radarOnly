fnameDPR='../DPRData/Atlantic/2A.GPM.DPR.V8-20180723.20180830-S183511-E200744.025593.V06A.HDF5'
fnameDPR='../monthly/SEAsia/2A.GPM.DPR.V8-20180723.20180602-S021248-E034522.024198.V06A.HDF5'
from netCDF4 import Dataset
import numpy as np
import pickle


import glob

def readSfcRainS(fname):
    fh=Dataset(fname,'r')
    sfcRain=fh['NS/SLV/precipRateNearSurface'][:,:]
    return sfcRain
#rainL=[]
#for f in fs:
#    rain=readSfcRainS(f)
#    rainL.append(rain.sum())
#rainL=np.array(rainL)
#inds=np.argsort(rainL)
import pickle
#pickle.dump(inds,open('sortedInd.pklz','wb'))

indx=range(327)
#stop
def readSfcRainCmb(fname,latmin,latmax):
    fh=Dataset(fname,'r')
    lat=fh['NS/Latitude'][:,24]
    a=np.nonzero((lat-latmin)*(lat-latmax)<0)
    sfcRain=fh['MS/surfPrecipTotRate'][a[0],:]
    #pRate=fh['NS/precipTotRate'][a[0],:,:]
    binNodes=fh['MS/phaseBinNodes'][a[0],:]
    lat=fh['MS/Latitude'][a[0],:]
    fh.close()
    return sfcRain,binNodes,a, lat

def readSfcRainCmbX(fname,latmin,latmax):
    fh=Dataset(fname,'r')
    lat=fh['FS/Latitude'][:,24]
    a=np.nonzero((lat-latmin)*(lat-latmax)<0)
    sfcRain=fh['FS/surfPrecipTotRate'][a[0],12:37]
    #pRate=fh['NS/precipTotRate'][a[0],:,:]
    binNodes=fh['FS/phaseBinNodes'][a[0],12:37]
    lat=fh['FS/Latitude'][a[0],12:37]
    fh.close()
    return sfcRain,binNodes,a, lat

def readSfcRain(fname,a):
    fh=Dataset(fname,'r')
    sfcRain=fh['NS/SLV/precipRateNearSurface'][a[0],12:37]
    lat=fh['NS/Latitude'][a[0],12:37]
    lon=fh['NS/Longitude'][a[0],12:37]
    sfcType=fh['NS/PRE/landSurfaceType'][a[0],12:37]
    h0=fh['NS/VER/heightZeroDeg'][a[0],12:37]
    bz=fh['NS/VER/binZeroDeg'][a[0],12:37]
    bcf=fh['NS/PRE/binClutterFreeBottom'][a[0],12:37]
    btop=fh['NS/PRE/binStormTop'][a[0],12:37]
    pType=(fh['NS/CSF/typePrecip'][a[0],12:37]/1e7).astype(int)
    binBB=fh['NS/CSF/binBBPeak'][a[0],12:37]
    bsfc=fh['NS/PRE/binRealSurface'][a[0],12:37]
    return lat,lon,sfcRain,sfcType,h0,pType,binBB,bsfc
import glob
#fs=[]
#for i in range(1,31):
#    fs1=glob.glob("/gpmdata/2018/08/%2.2i/Xradar/2A.GPM.DPR*"%i)
#    fs.extend(sorted(fs1))
#print(len(fs))

rainL=[]

import matplotlib.pyplot as plt
import matplotlib

sfcCMB_st_L=np.zeros((25),float)
sfcCMB6_st_L=np.zeros((25),float)
sfcCMBX_st_L=np.zeros((25),float)
sfcDPR_st_L=np.zeros((25),float)
countR_st_L=np.zeros((25),float)

sfcCMB_cv_L=np.zeros((25),float)
sfcCMB6_cv_L=np.zeros((25),float)
sfcCMBX_cv_L=np.zeros((25),float)
sfcDPR_cv_L=np.zeros((25),float)
countR_cv_L=np.zeros((25),float)

sfcCMB_st_O=np.zeros((25),float)
sfcCMB6_st_O=np.zeros((25),float)
sfcCMBX_st_O=np.zeros((25),float)
sfcDPR_st_O=np.zeros((25),float)
countR_st_O=np.zeros((25),float)

sfcCMB_cv_O=np.zeros((25),float)
sfcCMB6_cv_O=np.zeros((25),float)
sfcCMBX_cv_O=np.zeros((25),float)
sfcDPR_cv_O=np.zeros((25),float)
countR_cv_O=np.zeros((25),float)

fs=sorted(glob.glob("V6X/2B*"))
latmin=10
latmax=45

for f in fs[:60]:
    sfcRainCMB,nodes,a,latCMB=readSfcRainCmb(f,latmin,latmax)
    print(f)
    sdate=f.split('.')[-3]
    year=sdate[0:4]
    mm=sdate[4:6]
    day=sdate[6:8]
    orb=f.split('.')[-2]
    print(f,orb)
    fdpr=glob.glob('/gpmdata/%s/%s/%s/radar/2A.GPM.DPR.V8*%s*'%(year,mm,day,orb))
    fcmb=glob.glob('/gpmdata/%s/%s/%s/radar/2B.GPM.DPR*COR*%s*'%(year,mm,day,orb))
    fcmb1=glob.glob('/itedata/ITE745/%s/%s/%s/radar/2B.GPM.DPR*COR*%s*'%(year,mm,day,orb))
    print(fdpr[0])
    lat,lon,sfcRain,sfcType,h0,pType,binBB,bsfc=readSfcRain(fdpr[0],a)
    sfcRainCMBv6,nodes,a,latCMB=readSfcRainCmb(fcmb[0],latmin,latmax)
    sfcRainCMBX,nodesX,aX,latCMBX=readSfcRainCmbX(fcmb1[0],latmin,latmax)
    #stop
    a=np.nonzero(pType>0)
    b=np.nonzero(sfcType[a]!=0)
    if len(b[0])>0:
        if sfcRainCMBv6[a][b].max()>50:
            ind=np.argmax(sfcRainCMBv6[a][b])
            print('nodes=',nodes[a[0][b[0][ind]],a[1][b[0][ind]],:])
            print(sfcRainCMBv6[a][b].max(),sfcRainCMB[a][b].max(),sfcRainCMBX[a][b].max())

    for i,j in zip(a[0],a[1]):
        if sfcRain[i,j]==sfcRain[i,j]:
            if pType[i,j]==1:
                if sfcType[i,j]!=0:
                    sfcCMB_st_L[j]+=sfcRainCMB[i,j]
                    sfcCMB6_st_L[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_st_L[j]+=sfcRainCMBX[i,j]
                    sfcDPR_st_L[j]+=sfcRain[i,j]
                    countR_st_L[j]+=1
                else:
                    sfcCMB_st_O[j]+=sfcRainCMB[i,j]
                    sfcCMB6_st_O[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_st_O[j]+=sfcRainCMBX[i,j]
                    sfcDPR_st_O[j]+=sfcRain[i,j]
                    countR_st_O[j]+=1
            if pType[i,j]==2:
                if sfcType[i,j]!=0:
                    sfcCMB_cv_L[j]+=sfcRainCMB[i,j]
                    sfcCMB6_cv_L[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_cv_L[j]+=sfcRainCMBX[i,j]
                    sfcDPR_cv_L[j]+=sfcRain[i,j]
                    countR_cv_L[j]+=1
                else:
                    sfcCMB_cv_O[j]+=sfcRainCMB[i,j]
                    sfcCMB6_cv_O[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_cv_O[j]+=sfcRainCMBX[i,j]
                    sfcDPR_cv_O[j]+=sfcRain[i,j]
                    countR_cv_O[j]+=1
plt.figure(figsize=(8,6))
plt.suptitle('Midlatitude Summer')                    
plt.subplot(221)
plt.plot(sfcCMB_st_L/countR_st_L)
plt.plot(sfcCMB6_st_L/countR_st_L)
plt.plot(sfcCMBX_st_L/countR_st_L)
plt.plot(sfcDPR_st_L/countR_st_L)
plt.legend(['CV7','CV6','CX_V6','DV6'])
plt.title('Stratiform Land')
plt.subplot(223)
plt.plot(sfcCMB_st_O/countR_st_O)
plt.plot(sfcCMB6_st_O/countR_st_O)
plt.plot(sfcCMBX_st_O/countR_st_O)
plt.plot(sfcDPR_st_O/countR_st_O)
plt.title('Stratiform Ocean')
plt.subplot(222)
plt.plot(sfcCMB_cv_L/countR_cv_L)
plt.plot(sfcCMB6_cv_L/countR_cv_L)
plt.plot(sfcCMBX_cv_L/countR_cv_L)
plt.plot(sfcDPR_cv_L/countR_cv_L)
plt.title('Convective Land')
plt.subplot(224)
plt.plot(sfcCMB_cv_O/countR_cv_O)
plt.plot(sfcCMB6_cv_O/countR_cv_O)
plt.plot(sfcCMBX_cv_O/countR_cv_O)
plt.plot(sfcDPR_cv_O/countR_cv_O)
plt.title('Convective Ocean')
plt.tight_layout()
plt.savefig('condRainMidLat_v1.png')

#plt.plot(0.5*(np.array(countBB)+np.array(countnoBB)))
#plt.xlabel('Ray')
#plt.ylabel('Counts')
#plt.legend(['BB','noBB','Average'])


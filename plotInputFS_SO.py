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
def readSfcRainCmb(fname,latmin,latmax,a):
    fh=Dataset(fname,'r')
    lat=fh['NS/Latitude'][:,24]
    a0=np.nonzero((lat-latmin)*(lat-latmax)<0)
    pType=(fh['NS/Input/precipitationType'][a[0],:]/1e7).astype(int)
    sfcRain=fh['NS/surfPrecipTotRate'][a[0],:]
    #pRate=fh['NS/precipTotRate'][a[0],:,:]
    binNodes=fh['NS/phaseBinNodes'][a[0],:]
    lat=fh['NS/Latitude'][a[0],:]
    #print(fh['DiagGroup'])
    tb=fh['DiagGroup/dset1'][a[0],:,:]
    fh.close()
    return sfcRain,binNodes,a, lat, pType, tb

def readSfcRainCmbX(fname,latmin,latmax):
    fh=Dataset(fname,'r')
    lat=fh['FS/Latitude'][:,24]
    a=np.nonzero((lat-latmin)*(lat-latmax)<0)
    sfcRain=fh['FS/surfPrecipTotRate'][a[0],:]
    pType=(fh['FS/Input/precipitationType'][a[0],:]/1e7).astype(int)
    #pRate=fh['NS/precipTotRate'][a[0],:,:]
    binNodes=fh['FS/phaseBinNodes'][a[0],:]
    lat=fh['FS/Latitude'][a[0],:]

    fh.close()
    return sfcRain,binNodes,a, lat, pType

def readSfcRain(fname,a,latmin,latmax):
    fh=Dataset(fname,'r')
    lat=fh['FS/Latitude'][:,24]
    a=np.nonzero((lat-latmin)*(lat-latmax)<0)
    sfcRain=fh['NS/SLV/precipRateNearSurface'][a[0],:]
    lat=fh['NS/Latitude'][a[0],:]
    lon=fh['NS/Longitude'][a[0],:]
    sfcType=fh['NS/PRE/landSurfaceType'][a[0],:]
    h0=fh['NS/VER/heightZeroDeg'][a[0],:]
    bz=fh['NS/VER/binZeroDeg'][a[0],:]
    bcf=fh['NS/PRE/binClutterFreeBottom'][a[0],:]
    btop=fh['NS/PRE/binStormTop'][a[0],:]
    pType=(fh['NS/CSF/typePrecip'][a[0],:]/1e7).astype(int)
    binBB=fh['NS/CSF/binBBPeak'][a[0],:]
    bsfc=fh['NS/PRE/binRealSurface'][a[0],:]
    return lat,lon,sfcRain,sfcType,h0,pType,binBB,bsfc

def readSfcRainFS(fname,a,latmin,latmax):
    fh=Dataset(fname,'r')
    lat=fh['FS/Latitude'][:,24]
    a=np.nonzero((lat-latmin)*(lat-latmax)<0)
    sfcRain=fh['FS/SLV/precipRateNearSurface'][a[0],:]
    lat=fh['FS/Latitude'][a[0],:]
    lon=fh['FS/Longitude'][a[0],:]
    sfcType=fh['FS/PRE/landSurfaceType'][a[0],:]
    h0=fh['FS/VER/heightZeroDeg'][a[0],:]
    bz=fh['FS/VER/binZeroDeg'][a[0],:]
    bcf=fh['FS/PRE/binClutterFreeBottom'][a[0],:]
    btop=fh['FS/PRE/binStormTop'][a[0],:]
    pType=(fh['FS/CSF/typePrecip'][a[0],:]/1e7).astype(int)
    binBB=fh['FS/CSF/binBBPeak'][a[0],:]
    bsfc=fh['FS/PRE/binRealSurface'][a[0],:]
    zm=fh['FS/PRE/zFactorMeasured'][a[0],:,:,:]
    snowIce=fh['FS/PRE/snowIceCover'][a[0],:]
    return lat,lon,sfcRain,sfcType,h0,pType,binBB,bsfc,zm,bz,bcf,snowIce
import glob
#fs=[]
#for i in range(1,31):
#    fs1=glob.glob("/gpmdata/2018/08/%2.2i/Xradar/2A.GPM.DPR*"%i)
#    fs.extend(sorted(fs1))
#print(len(fs))

rainL=[]

import matplotlib.pyplot as plt
import matplotlib

sfcCMB_st_L=np.zeros((49),float)
sfcCMB6_st_L=np.zeros((49),float)
sfcCMBX_st_L=np.zeros((49),float)
sfcDPR_st_L=np.zeros((49),float)
countR_st_L=np.zeros((49),float)

sfcCMB_cv_L=np.zeros((49),float)
sfcCMB6_cv_L=np.zeros((49),float)
sfcCMBX_cv_L=np.zeros((49),float)
sfcDPR_cv_L=np.zeros((49),float)
countR_cv_L=np.zeros((49),float)

sfcCMB_st_O=np.zeros((49),float)
sfcCMB6_st_O=np.zeros((49),float)
sfcCMBX_st_O=np.zeros((49),float)
sfcDPR_st_O=np.zeros((49),float)
countR_st_O=np.zeros((49),float)

sfcCMB_cv_O=np.zeros((49),float)
sfcCMB6_cv_O=np.zeros((49),float)
sfcCMBX_cv_O=np.zeros((49),float)
sfcDPR_cv_O=np.zeros((49),float)
countR_cv_O=np.zeros((49),float)

fs=sorted(glob.glob("V6X/2B*"))
latmin=-60
latmax=-40
countMiss=np.zeros((4),float)

tbL=[]
zmL=[]
bzdL=[]
bcfL=[]
coordL=[]
bsfcL=[]
sfcRainL=[]
jL=[]

tbL2=[]
bzdL2=[]
bcfL2=[]
coordL2=[]
bsfcL2=[]
sfcRainL2=[]
jL2=[]

for f in fs[:300]:
    sfcRainCMB,nodes,ac,latCMBX,pTypeCMB=readSfcRainCmbX(f,latmin,latmax)
    print(f)
    sdate=f.split('.')[-3]
    year=sdate[0:4]
    mm=sdate[4:6]
    day=sdate[6:8]
    orb=f.split('.')[-2]
    print(f,orb)
    fdpr=glob.glob('/itedata/ITE749/%s/%s/%s/radar/2A.GPM.DPR.V9*%s*'%(year,mm,day,orb))
    fcmb=glob.glob('/gpmdata/%s/%s/%s/radar/2B.GPM.DPR*COR*%s*'%(year,mm,day,orb))
    fcmb1=glob.glob('/itedata/ITE745/%s/%s/%s/radar/2B.GPM.DPR*COR*%s*'%(year,mm,day,orb))
    print(fdpr[0])
    #latmin=latCMBX[:,24].min()-0.0
    #latmax=latCMBX[:,24].max()+0.0
    lat,lon,sfcRain,sfcType,h0,pType,binBB,bsfc,zm,bzd,\
        bcf,snowIce=readSfcRainFS(fdpr[0],ac,latmin,latmax)

    print(lat[0]-latCMBX[0])
    #stop
    sfcRainCMBv6,nodes,a,\
        latCMB,pTypeCMB,tb=readSfcRainCmb(fcmb[0],latmin,latmax,ac)
    print(latCMBX[0]-latCMB[0])
    sfcRainCMBX,nodesX,aX,latCMBX2,pTypeCMBX=readSfcRainCmbX(fcmb1[0],latmin,latmax)
    #
    #stop
    a=np.nonzero(pType>0)
    b=np.nonzero(sfcType[a]!=0)
    if len(b[0])>0:
        if sfcRainCMBv6[a][b].max()>50:
            ind=np.argmax(sfcRainCMBv6[a][b])
            print('nodes=',nodes[a[0][b[0][ind]],a[1][b[0][ind]],:])
            print(sfcRainCMBv6[a][b].max(),sfcRainCMB[a][b].max(),sfcRainCMBX[a][b].max())

    for i,j in zip(a[0],a[1]):
        if i==sfcRain.shape[0]-1:
            continue
        if sfcRain[i,j]==sfcRain[i,j]:
            if pType[i,j]==1:
                if sfcType[i,j]!=0:
                    sfcCMB_st_L[j]+=sfcRainCMB[i,j]
                    sfcCMB6_st_L[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_st_L[j]+=sfcRainCMBX[i,j]
                    sfcDPR_st_L[j]+=sfcRain[i,j]
                    countR_st_L[j]+=1
                    if pTypeCMB[i,j]<1:
                        countMiss[0]+=1
                else:
                    sfcCMB_st_O[j]+=sfcRainCMB[i,j]
                    sfcCMB6_st_O[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_st_O[j]+=sfcRainCMBX[i,j]
                    sfcDPR_st_O[j]+=sfcRain[i,j]
                    countR_st_O[j]+=1
                    if snowIce[i,j]==0:
                        bsfcL.append(bsfc[i,j,0])
                        bzdL.append(bzd[i,j])
                        zmL.append(zm[i,j,100:,:])
                        tbL.append(tb[i,j,:])
                        coordL.append([lon[i,j],lat[i,j]])
                        sfcRainL.append(sfcRainCMB[i,j])
                        bcfL.append(bcf[i,j])
                        jL.append([j,1])

                    if pTypeCMB[i,j]<1:
                        countMiss[1]+=1
            if pType[i,j]==2:
                if sfcType[i,j]!=0:
                    sfcCMB_cv_L[j]+=sfcRainCMB[i,j]
                    sfcCMB6_cv_L[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_cv_L[j]+=sfcRainCMBX[i,j]
                    sfcDPR_cv_L[j]+=sfcRain[i,j]
                    countR_cv_L[j]+=1
                    if pTypeCMB[i,j]<1:
                        countMiss[2]+=1
                else:
                    sfcCMB_cv_O[j]+=sfcRainCMB[i,j]
                    sfcCMB6_cv_O[j]+=sfcRainCMBv6[i,j]
                    sfcCMBX_cv_O[j]+=sfcRainCMBX[i,j]
                    sfcDPR_cv_O[j]+=sfcRain[i,j]
                    countR_cv_O[j]+=1
                    if snowIce[i,j]==0:
                        bsfcL.append(bsfc[i,j,0])
                        bzdL.append(bzd[i,j])
                        zmL.append(zm[i,j,100:,:])
                        tbL.append(tb[i,j,:])
                        sfcRainL.append(sfcRainCMB[i,j])
                        coordL.append([lon[i,j],lat[i,j]])
                        bcfL.append(bcf[i,j])
                        jL.append([j,2])
                    if pTypeCMB[i,j]<1:
                        countMiss[3]+=1
    
    a=np.nonzero(pType==0)
    b=np.nonzero(sfcType[a]!=0)
    for i,j in zip(a[0],a[1]):
        if i==sfcRain.shape[0]-1:
            continue
        if snowIce[i,j]==0:
            bsfcL2.append(bsfc[i,j,0])
            bzdL2.append(bzd[i,j])
            #zmL.append(zm[i,j,100:,:])
            tbL2.append(tb[i,j,:])
            coordL2.append([lon[i,j],lat[i,j]])
            sfcRainL2.append(sfcRainCMB[i,j])
            bcfL2.append(bcf[i,j])
            jL2.append([j,0])
    print('Missed Ocean =',countMiss[1]/countR_st_O.sum(),countMiss[3]/countR_cv_O.sum(),\
              (countMiss[1]+countMiss[3])/(countR_st_O.sum()+countR_cv_O.sum()))
    print('Missed Land =',countMiss[0]/countR_st_L.sum(),countMiss[2]/countR_cv_L.sum(),\
              (countMiss[0]+countMiss[2])/(countR_st_L.sum()+countR_cv_L.sum()))
    bL=(sum(sfcCMB_cv_L)+ sum(sfcCMB_st_L))/\
              (sum(sfcCMB6_cv_L)+ sum(sfcCMB6_st_L))
    bLx=(sum(sfcCMBX_cv_L)+ sum(sfcCMBX_st_L))/\
              (sum(sfcCMB6_cv_L)+ sum(sfcCMB6_st_L))
    bO=(sum(sfcCMB_cv_O)+ sum(sfcCMB_st_O))/\
              (sum(sfcCMB6_cv_O)+ sum(sfcCMB6_st_O))
    bOx=(sum(sfcCMBX_cv_O)+ sum(sfcCMBX_st_O))/\
              (sum(sfcCMB6_cv_O)+ sum(sfcCMB6_st_O))
    print('syst_diff=> Land=',bL,bLx,' Ocean=',bO,bOx)
plt.figure(figsize=(12,8))
plt.suptitle('SO Winter')                    
plt.subplot(221)
plt.plot(sfcCMB_st_L/countR_st_L)
plt.plot(sfcCMB6_st_L/countR_st_L)
plt.plot(sfcCMBX_st_L/countR_st_L)
plt.plot(sfcDPR_st_L/countR_st_L)
plt.legend(['CV7','CV6','ITE745','DITE749'],ncols=2)
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
plt.savefig('condRainSO_cmp2.png')

#plt.plot(0.5*(np.array(countBB)+np.array(countnoBB)))
#plt.xlabel('Ray')
#plt.ylabel('Counts')
#plt.legend(['BB','noBB','Average'])

# 1.1383305980412695 1.2234138453398562  Ocean= 0.9886373783085047 1.0046354827886301
import xarray as xr
tbx=xr.DataArray(tbL,dims=['ns','nt'])
zmx=xr.DataArray(zmL,dims=['ns','nz','n2'])
bsfcx=xr.DataArray(bsfcL,dims=['ns'])
bcfx=xr.DataArray(bcfL,dims=['ns'])
bzdx=xr.DataArray(bzdL,dims=['ns'])
jx=xr.DataArray(jL,dims=['ns','n2'])
coordx=xr.DataArray(coordL,dims=['ns','n2'])
sfcRainX=xr.DataArray(sfcRainL,dims=['ns'])
d=xr.Dataset({"Tb":tbx,"Zm":zmx,"bsfc":bsfcx,"bzd":bzdx,"j":jx,"sfcRain":sfcRainX,\
                  "bcf":bcfx,"coord":coordx})
d.to_netcdf("SO_traininDataSet_2.nc")

tbx2=xr.DataArray(tbL2,dims=['ns','nt'])
bsfcx2=xr.DataArray(bsfcL2,dims=['ns'])
bcfx2=xr.DataArray(bcfL2,dims=['ns'])
bzdx2=xr.DataArray(bzdL2,dims=['ns'])
jx2=xr.DataArray(jL2,dims=['ns','n2'])
coordx2=xr.DataArray(coordL2,dims=['ns','n2'])
sfcRainX2=xr.DataArray(sfcRainL2,dims=['ns'])
d2=xr.Dataset({"Tb":tbx2,"bsfc":bsfcx2,"bzd":bzdx2,"j":jx2,"sfcRain":sfcRainX2,\
                  "bcf":bcfx2,"coord":coordx2})
d2.to_netcdf("SO_traininDataSet_noPrecip.nc")

pickle.dump([sfcCMB_cv_O,sfcCMB_st_O,countR_cv_O,countR_st_O],\
                open('SO_ratesAndCounts.pklz','wb'))

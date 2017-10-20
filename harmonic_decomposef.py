#!/usr/bin/env python3
"""
Created on Wed Sep 14 07:24:26 2016
@author: dipankar saikia # deep07004@gmail.com
"""
import sys
import numpy as np
from scipy import linalg as linalg
import matplotlib.pylab as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d


def plot_it (trinfo,dataR,dataT,twb,twe):
    t_end=25
    zoom=100
    t_start=trinfo[0,0]
    ntwb=np.int((twb+np.abs(t_start))/trinfo[0,1])
    ntwe=np.int((twe+np.abs(t_start))/trinfo[0,1])
    nn=int(np.ceil((t_end-t_start)/trinfo[0,1]))
    tt_axis=np.linspace(t_start,t_end,nn)
    fig1 = plt.figure(figsize=(12,12))
    plt.subplot(1,2,1)
    plt.title('Radial')
    plt.ylim([-10,399])
    for i in range(0,nrf):
        plt.plot(tt_axis,dataR[0:nn,i]*zoom+trinfo[i,3],'k-')
        a=np.argmax(dataR[ntwb:ntwe,i])*trinfo[0,1]+twb
        b=np.max(dataR[ntwb:ntwe,i])*zoom+trinfo[i,3]
        plt.plot(a,b,'rx')
    ax=fig1.gca()
    ax.add_patch(patches.Rectangle((twb, -50), (twe-twb), 500,alpha=0.20))
    plt.subplot(1,2,2)
    plt.title('Transverse')
    plt.ylim([-10,399])
    for i in range(0,nrf):
        plt.plot(tt_axis,dataT[0:nn,i]*zoom+trinfo[i,3],'k-')
    ax=fig1.gca()
    ax.add_patch(patches.Rectangle((twb, -50), (twe-twb), 500,alpha=0.20))
  #  plt.show()

maxamp=0.05

LCS=True # If data is in left handed R-T-Z system where T -> +90 from R
# e.g Park and Levin fomulation is in right jhanded R-T-Z system
# Raysum is in left handed coordinate system 

if LCS:
    csm=-1.0
else:
    csm=1.0

try:
    fp=open(sys.argv[1],'r')
except:
    sys.exit()

fpin=fp.readlines()
fp.close()
nl=len(fpin)
trinfo = np.zeros([1000,4])
j=-1
k=0
for i in range(0,len(fpin)):
    l=fpin[i].split()
    if len(l)==4:
        j=j+1
        trinfo[j,0]=np.float(l[0])
        trinfo[j,1]=np.float(l[1])
        trinfo[j,2]=np.int(l[2])
        trinfo[j,3]=np.float(l[3])
        k=0
        if (i==0):
            npts=np.int(l[2])
            dataR = np.zeros([npts,1000])
            dataT = np.zeros([npts,1000])
    else :
        dataR[k,j]=np.float(l[0])
        dataT[k,j]=np.float(l[1])*csm
        k=k+1

nrf=j+1
ts=(np.argmax(dataR[:,0])*trinfo[0,1])+trinfo[0,0]-3.0
print (np.argmax(dataR[:,0]))
plot_it (trinfo,dataR,dataT,ts,ts+6)


tt=np.linspace(0,trinfo[0,1]*trinfo[0,2],trinfo[0,2])+trinfo[0,0]

# convert the data to frequency domain

dataRF_in = np.fft.rfft(dataR,axis=0)
dataTF_in = np.fft.rfft(dataT,axis=0)
oldnf=len(dataRF_in[:,0])
print ('Frequency no' , oldnf)
oldff=np.linspace(0,oldnf,oldnf)
nf=npts
ff=np.linspace(0,oldnf,nf)
dataRF = np.zeros([nf,nrf],dtype=complex)
dataTF = np.zeros([nf,nrf],dtype=complex)

for i in range(0,nrf):
    print(i)
    func1 = interp1d(oldff,dataRF_in[:,i].real, kind='slinear')
    func2 = interp1d(oldff,dataTF_in[:,i].real, kind='slinear')
    dataRF[:,i].real = func1(ff)
    dataTF[:,i].real = func2(ff)
    func1 = interp1d(oldff,dataRF_in[:,i].imag, kind='slinear')
    func2 = interp1d(oldff,dataTF_in[:,i].imag, kind='slinear')
    dataRF[:,i].imag = func1(ff)
    dataTF[:,i].imag = func2(ff)
print ('Frequency no' , nf)

# Modelled radial component
op=np.zeros([2*nrf,10])
Ac=np.array(np.zeros(nf),dtype=complex)
Bc=np.array(np.zeros(nf),dtype=complex)
Cc=np.array(np.zeros(nf),dtype=complex)
Dc=np.array(np.zeros(nf),dtype=complex)
Ec=np.array(np.zeros(nf),dtype=complex)
Fc=np.array(np.zeros(nf),dtype=complex)
Gc=np.array(np.zeros(nf),dtype=complex)
Hc=np.array(np.zeros(nf),dtype=complex)
Ic=np.array(np.zeros(nf),dtype=complex)
Jc=np.array(np.zeros(nf),dtype=complex)
data=np.zeros(2*nrf)
bootstrap=100
A=np.zeros([npts,bootstrap])
B=np.zeros([npts,bootstrap])
C=np.zeros([npts,bootstrap])
D=np.zeros([npts,bootstrap])
E=np.zeros([npts,bootstrap])
F=np.zeros([npts,bootstrap])
G=np.zeros([npts,bootstrap])
H=np.zeros([npts,bootstrap])
I=np.zeros([npts,bootstrap])
J=np.zeros([npts,bootstrap])
for bn in range(0,bootstrap):
    print('Bootstrap no.: ',bn)
    listbn=np.random.randint(0,nrf,nrf)
    i=0    
    for k in listbn:
        j=i+nrf
        theta = trinfo[k,3]*np.pi/180.0 
        op[i,0] = 3.0
        op[i,1] = np.cos(theta)
        op[i,2] = np.sin(theta)
        op[i,3] = np.cos(2.0*theta)
        op[i,4] = np.sin(2.0*theta)
        op[i,5] = 0.0
        op[i,6] = np.cos(theta)
        op[i,7] = np.sin(theta)
        op[i,8] = np.cos(2.0*theta)
        op[i,9] = np.sin(2.0*theta)
        op[j,0] = 0.0
        op[j,1] = np.sin(theta)
        op[j,2] = -1.0*np.cos(theta)
        op[j,3] = np.sin(2.0*theta)
        op[j,4] = -1.0*np.cos(2.0*theta)
        op[j,5] = 3.0
        op[j,6] = -1.0*np.sin(theta)
        op[j,7] = np.cos(theta)
        op[j,8] = -1.0*np.sin(2.0*theta)
        op[j,9] = np.cos(2.0*theta)
        i+=1
    
    invop=linalg.pinv2(op,check_finite=True)
    
    for i in range(0,nf):
        data[0:nrf]=dataRF.real[i,listbn]
        data[nrf:2*nrf]=dataTF.real[i,listbn]
        [Ac.real[i],Bc.real[i],Cc.real[i],Dc.real[i],Ec.real[i],Fc.real[i],Gc.real[i],Hc.real[i],Ic.real[i],Jc.real[i]] = np.matmul(invop,data)
        data[0:nrf]=dataRF.imag[i,listbn] 
        data[nrf:2*nrf]=dataTF.imag[i,listbn]
        [Ac.imag[i],Bc.imag[i],Cc.imag[i],Dc.imag[i],Ec.imag[i],Fc.imag[i],Gc.imag[i],Hc.imag[i],Ic.imag[i],Jc.imag[i]] = np.matmul(invop,data)
    # convert to time domain
    A[:,bn]=np.fft.irfft(Ac)[0:npts]
    B[:,bn]=np.fft.irfft(Bc)[0:npts]
    C[:,bn]=np.fft.irfft(Cc)[0:npts]
    D[:,bn]=np.fft.irfft(Dc)[0:npts]
    E[:,bn]=np.fft.irfft(Ec)[0:npts]
    F[:,bn]=np.fft.irfft(Fc)[0:npts]
    G[:,bn]=np.fft.irfft(Gc)[0:npts]
    H[:,bn]=np.fft.irfft(Hc)[0:npts]
    I[:,bn]=np.fft.irfft(Ic)[0:npts]
    J[:,bn]=np.fft.irfft(Jc)[0:npts]
    
AF,BF,CF,DF,EF,FF,GF,HF,IF,JF = np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3])

AF[:,1] = np.mean(A,axis=1)
AF[:,0] = AF[:,1]-np.std(A,axis=1)
AF[:,2] = AF[:,1]+np.std(A,axis=1)
BF[:,1] = np.mean(B,axis=1)
BF[:,0] = BF[:,1]-np.std(B,axis=1)
BF[:,2] = BF[:,1]+np.std(B,axis=1)
CF[:,1] = np.mean(C,axis=1)
CF[:,0] = CF[:,1]-np.std(C,axis=1)
CF[:,2] = CF[:,1]+np.std(C,axis=1)
DF[:,1] = np.mean(D,axis=1)
DF[:,0] = DF[:,1]-np.std(D,axis=1)
DF[:,2] = DF[:,1]+np.std(D,axis=1)
EF[:,1] = np.mean(E,axis=1)
EF[:,0] = EF[:,1]-np.std(E,axis=1)
EF[:,2] = EF[:,1]+np.std(E,axis=1)
FF[:,1] = np.mean(F,axis=1)
FF[:,0] = FF[:,1]-np.std(F,axis=1)
FF[:,2] = FF[:,1]+np.std(F,axis=1)
GF[:,1] = np.mean(G,axis=1)
GF[:,0] = GF[:,1]-np.std(G,axis=1)
GF[:,2] = GF[:,1]+np.std(G,axis=1)
HF[:,1] = np.mean(H,axis=1)
HF[:,0] = HF[:,1]-np.std(H,axis=1)
HF[:,2] = HF[:,1]+np.std(H,axis=1)
IF[:,1] = np.mean(I,axis=1)
IF[:,0] = IF[:,1]-np.std(I,axis=1)
IF[:,2] = IF[:,1]+np.std(I,axis=1)
JF[:,1] = np.mean(J,axis=1)
JF[:,0] = JF[:,1]-np.std(J,axis=1)
JF[:,2] = JF[:,1]+np.std(J,axis=1)


fig = plt.figure(figsize=(8.3,11.7),frameon=False)
xlimit=[-2,15]
fig.subplots_adjust(hspace=0.5)
plt.subplot(5,2,1)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('Constanat (K=0)')
plt.plot(tt[0:npts],AF[0:npts])
plt.subplot(5,2,3)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('cos (K=1)')
plt.plot(tt[0:npts],BF[0:npts])
plt.subplot(5,2,5)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('sin (K=1)')
plt.plot(tt[0:npts],CF[0:npts])
plt.subplot(5,2,7)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('cos (K=2)')
plt.plot(tt[0:npts],DF[0:npts])
plt.subplot(5,2,9)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('sin (K=2)')
plt.plot(tt[0:npts],EF[0:npts])

# Unmodelled

plt.subplot(5,2,2)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('Constanat (K=0)')
plt.plot(tt[0:npts],FF[0:npts])
plt.subplot(5,2,4)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('cos (K=1)')
plt.plot(tt[0:npts],GF[0:npts])
plt.subplot(5,2,6)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('sin (K=1)')
plt.plot(tt[0:npts],HF[0:npts])
plt.subplot(5,2,8)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('cos (K=2)')
plt.plot(tt[0:npts],IF[0:npts])
plt.subplot(5,2,10)
plt.ylim([-maxamp,maxamp])
plt.xlim(xlimit)
plt.title('sin (K=2)')
plt.plot(tt[0:npts],JF[0:npts])
###
fp=open('modeled.dat','w')
for i in range(0,npts):
    out=str(tt[i])+'\t'+str(AF[i,1])+'\t'+str(BF[i,1])+'\t'+str(CF[i,1])+'\t'+str(DF[i,1])+'\t'+str(EF[i,1])+'\n'
    fp.write(out)
fp.close()
fp=open('unmodeled.dat','w')
for i in range(0,npts):
    out=str(tt[i])+'\t'+str(FF[i,1])+'\t'+str(GF[i,1])+'\t'+str(HF[i,1])+'\t'+str(IF[i,1])+'\t'+str(JF[i,1])+'\n'
    fp.write(out)
fp.close()
### Analysis
t1=0.0
t2=4.0
n1=int(np.ceil((np.abs(trinfo[0,0])+t1)/trinfo[0,1]))
n2=int(np.ceil((np.abs(trinfo[0,0])+t2)/trinfo[0,1]))
print (n1,n2)
eta=np.arctan2(np.max(np.abs(BF[n1:n2,1])),np.max(np.abs(CF[n1:n2,1])))*180.0/np.pi
print('Strike :',eta)
eta=np.arctan2(np.max(np.abs(DF[n1:n2,1])),np.max(np.abs(EF[n1:n2,1])))*180.0/np.pi
print('Strike :',eta)
plt.savefig('plot.pdf')
plt.show()

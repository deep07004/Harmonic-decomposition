#!/usr/bin/env python3
"""
Created on Wed Sep 14 07:24:26 2016

@author: dipankar (deep07004@gmail.com)
"""
import sys
import numpy as np
from scipy import linalg as linalg
import matplotlib.pylab as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d


def plot_it (trinfo,dataR,dataT,twb,twe):
    t_end=twe
    zoom=100
    t_start=trinfo[0,0]
    ntwb=np.int((twb+np.abs(t_start))/trinfo[0,1])
    ntwe=np.int((twe+np.abs(t_start))/trinfo[0,1])
    nn=int(np.ceil((t_end-t_start)/trinfo[0,1]))
    time=np.linspace(t_start,t_end,nn)
    fig = plt.figure(figsize=(4.5,6.0))
    plt.clf()
    ax1 = fig.add_axes([0.075, 0.075, 0.375, 0.875])
    ax2 = fig.add_axes([0.55, 0.075, 0.375, 0.875])
    ax1.set_title('Radial')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Back-azimuth (deg)')
    ax1.set_xlim(twb, twe)
    ax1.set_ylim(-5, 370)
    ax1.set_yticks(())
    ax2.set_title('Transverse')
    ax2.set_xlabel('Time (s)')
    #ax2.set_ylabel('Back-azimuth (deg)')
    ax2.set_xlim(twb, twe)
    ax2.set_ylim(-5, 370)
    for i in range(0,nrf):
        y = trinfo[i,3]
        data = dataR[0:nn,i]
        maxval = 60
        ax1.plot(time,y+data*maxval,'k-',linewidth=0.5)
        ax1.fill_between(time, y, y+data*maxval, where=data+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
        ax1.fill_between(time, y, y+data*maxval, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
        data = dataT[0:nn,i]
        ax2.plot(time,y+data*maxval,'k-',linewidth=0.5)
        ax2.fill_between(time, y, y+data*maxval, where=data+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
        ax2.fill_between(time, y, y+data*maxval, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
    fig.savefig('baz-section.pdf',dpi=300)

# plotting limits
xlimit=[-2,10]  #  Time
bootstrap=200
zoom=5

LCS=False # If data is in left handed R-T-Z system where T -> +90 from R
# e.g Park and Levin fomulation is in right handed R-T-Z system
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
print (np.argmax(dataR[:,0]))
plot_it (trinfo,dataR,dataT,xlimit[0],xlimit[1])


tt=np.linspace(0,trinfo[0,1]*trinfo[0,2],np.int(trinfo[0,2]))+trinfo[0,0]

# convert the data to frequency domain

dataRF_in = np.fft.rfft(dataR,axis=0)
dataTF_in = np.fft.rfft(dataT,axis=0)
oldnf=len(dataRF_in[:,0])
oldff=np.linspace(0,oldnf,oldnf)
nf = npts*2
print ('Frequency no' , oldnf, nf)
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

#print ('Frequency no' , nf)

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
    print('Bootstrap no.: ',bn+1)
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
        op[j,5] = 0.3
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
    A[:,bn]=np.fft.irfft(Ac*3.0)[0:npts]
    B[:,bn]=np.fft.irfft(Bc)[0:npts]
    C[:,bn]=np.fft.irfft(Cc)[0:npts]
    D[:,bn]=np.fft.irfft(Dc)[0:npts]
    E[:,bn]=np.fft.irfft(Ec)[0:npts]
    F[:,bn]=np.fft.irfft(Fc*0.3)[0:npts]
    G[:,bn]=np.fft.irfft(Gc)[0:npts]
    H[:,bn]=np.fft.irfft(Hc)[0:npts]
    I[:,bn]=np.fft.irfft(Ic)[0:npts]
    J[:,bn]=np.fft.irfft(Jc)[0:npts]
    
AF,BF,CF,DF,EF,FF,GF,HF,IF,JF = np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3]),np.zeros([npts,3]),np.zeros([npts,3]),\
                                np.zeros([npts,3])

AF[:,1] = np.mean(A,axis=1)
AF[:,0] = AF[:,1]-1.98*np.std(A,axis=1)
AF[:,2] = AF[:,1]+1.98*np.std(A,axis=1)
BF[:,1] = np.mean(B,axis=1)
BF[:,0] = BF[:,1]-1.98*np.std(B,axis=1)
BF[:,2] = BF[:,1]+1.98*np.std(B,axis=1)
CF[:,1] = np.mean(C,axis=1)
CF[:,0] = CF[:,1]-1.98*np.std(C,axis=1)
CF[:,2] = CF[:,1]+1.98*np.std(C,axis=1)
DF[:,1] = np.mean(D,axis=1)
DF[:,0] = DF[:,1]-1.98*np.std(D,axis=1)
DF[:,2] = DF[:,1]+1.98*np.std(D,axis=1)
EF[:,1] = np.mean(E,axis=1)
EF[:,0] = EF[:,1]-1.98*np.std(E,axis=1)
EF[:,2] = EF[:,1]+1.98*np.std(E,axis=1)
FF[:,1] = np.mean(F,axis=1)
FF[:,0] = FF[:,1]-1.98*np.std(F,axis=1)
FF[:,2] = FF[:,1]+1.98*np.std(F,axis=1)
GF[:,1] = np.mean(G,axis=1)
GF[:,0] = GF[:,1]-1.98*np.std(G,axis=1)
GF[:,2] = GF[:,1]+1.98*np.std(G,axis=1)
HF[:,1] = np.mean(H,axis=1)
HF[:,0] = HF[:,1]-1.98*np.std(H,axis=1)
HF[:,2] = HF[:,1]+1.98*np.std(H,axis=1)
IF[:,1] = np.mean(I,axis=1)
IF[:,0] = IF[:,1]-1.98*np.std(I,axis=1)
IF[:,2] = IF[:,1]+1.98*np.std(I,axis=1)
JF[:,1] = np.mean(J,axis=1)
JF[:,0] = JF[:,1]-1.98*np.std(J,axis=1)
JF[:,2] = JF[:,1]+1.98*np.std(J,axis=1)

## Plotting the harmonic components
maxval = zoom
fig = plt.figure(figsize=[6.,6.])
plt.clf()
ax1 = fig.add_axes([0.025, 0.3, 0.45, 0.65])
ax2 = fig.add_axes([0.525, 0.3, 0.45, 0.65])
ax3 = fig.add_axes([0.025, 0.075, 0.45, 0.15])
ax4 = fig.add_axes([0.525, 0.075, 0.45, 0.15])
ax1.set_title('Anisotropy/Dip')
#ax1.set_xlabel('Time (s)')
ax1.set_ylabel('')
ax1.set_yticks(())
ax1.set_xlim(xlimit[0], xlimit[1])
ax1.set_ylim(-0.5, 4.75)
ax2.set_title('Unmodelled')
#ax2.set_xlabel('Time (s)')
ax2.set_ylabel('')
ax2.set_yticks(())
ax2.set_xlim(xlimit[0], xlimit[1])
ax2.set_ylim(-0.5, 4.75)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('')
ax3.set_yticks(())
ax3.set_xlim(xlimit[0], xlimit[1])
ax3.set_ylim(-0.5, 1.75)
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('')
ax4.set_yticks(())
ax4.set_xlim(xlimit[0], xlimit[1])
ax4.set_ylim(-0.5, 1.75)
i=4
for f in [AF, BF, CF, DF, EF]:
    data = f[0:npts,1] * maxval
    datal = f[0:npts,0] * maxval
    datau = f[0:npts,2] * maxval
    ax1.plot(tt[0:npts],i+data,'k-',linewidth=0.5)
    ax1.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 >= 0.,facecolor='cyan', linewidth=0)
    ax1.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 <= 0.,facecolor='cyan', linewidth=0)
    ax1.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax1.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
    i = i-1
m=(np.abs(xlimit[1])/12.0)**0.9
ax1.annotate('constant',xy=(xlimit[1]-3.4*m,4.12))
ax1.annotate('cos($\phi$)',xy=(xlimit[1]-2.585*m,3.12))
ax1.annotate('sin($\phi$)',xy=(xlimit[1]-2.5*m,2.12))
ax1.annotate('cos(2$\phi$)',xy=(xlimit[1]-2.985*m,1.12))
ax1.annotate('sin(2$\phi$)',xy=(xlimit[1]-2.85*m,0.12))
ax1.annotate('(N-S)',xy=(xlimit[1]-2.25*m,2.75))
ax1.annotate('(E-W)',xy=(xlimit[1]-2.3*m,1.75))
ax1.annotate('(N-S)',xy=(xlimit[1]-2.25*m,0.75))
ax1.annotate('(E-W)',xy=(xlimit[1]-2.3*m,-0.25))
i=4
for f in [FF, GF, HF, IF, JF]:
    data = f[0:npts,1] * maxval
    datal = f[0:npts,0] * maxval
    datau = f[0:npts,2] * maxval
    ax2.plot(tt[0:npts],i+data,'k-',linewidth=0.5)
    ax2.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 >= 0.,facecolor='cyan', linewidth=0)
    ax2.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 <= 0.,facecolor='cyan', linewidth=0)
    ax2.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax2.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
    i = i-1
ax2.text(xlimit[0]-1.45*m,3.925,'H0',fontsize=14)
ax2.text(xlimit[0]-1.45*m,2.925,'H1',fontsize=14)
ax2.text(xlimit[0]-1.45*m,1.925,'H2',fontsize=14)
ax2.text(xlimit[0]-1.45*m,0.925,'H3',fontsize=14)
ax2.text(xlimit[0]-1.45*m,-0.090,'H4',fontsize=14)
ax2.annotate('constant',xy=(xlimit[1]-3.4*m,4.15))
ax2.annotate('cos($\phi$)',xy=(xlimit[1]-2.585*m,3.12))
ax2.annotate('sin($\phi$)',xy=(xlimit[1]-2.5*m,2.12))
ax2.annotate('cos(2$\phi$)',xy=(xlimit[1]-2.985*m,1.12))
ax2.annotate('sin(2$\phi$)',xy=(xlimit[1]-2.85*m,0.12))
ax2.annotate('(N-S)',xy=(xlimit[1]-2.25*m,2.75))
ax2.annotate('(E-W)',xy=(xlimit[1]-2.3*m,1.75))
ax2.annotate('(N-S)',xy=(xlimit[1]-2.25*m,0.75))
ax2.annotate('(E-W)',xy=(xlimit[1]-2.3*m,-0.25))

data = np.sqrt((BF[0:npts,1])**2+ (CF[0:npts,1])**2) * maxval
ax3.fill_between(tt[0:npts], 1, 1+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
ax3.plot(tt[0:npts], 1+data,'k-',linewidth=0.5)
data = np.sqrt((DF[0:npts,1])**2+ (EF[0:npts,1])**2) * maxval
ax3.fill_between(tt[0:npts], 0, data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
ax3.plot(tt[0:npts], data,'k-',linewidth=0.5)
ax3.annotate('$\sqrt{H1^{2}+H2^{2}}$',xy=(xlimit[1]-4.5*m,1.15))
ax3.annotate('$\sqrt{H3^{2}+H4^{2}}$',xy=(xlimit[1]-4.5*m,0.15))
data = np.sqrt((GF[0:npts,1])**2+ (HF[0:npts,1])**2) * maxval
ax4.fill_between(tt[0:npts], 1, 1+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
ax4.plot(tt[0:npts], 1+data,'k-',linewidth=0.5)
data = np.sqrt((IF[0:npts,1])**2+ (JF[0:npts,1])**2) * maxval
ax4.fill_between(tt[0:npts], 0, data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
ax4.plot(tt[0:npts], data,'k-',linewidth=0.5)
ax4.annotate('$\sqrt{H1^{2}+H2^{2}}$',xy=(xlimit[1]-4.5*m,1.15))
ax4.annotate('$\sqrt{H3^{2}+H4^{2}}$',xy=(xlimit[1]-4.5*m,0.15))
fig.savefig('aniso-dip.pdf',dpi=300)

###
fp=open('modelled.dat','w')
fp1=open('unmodelled.dat','w')
for i in range(0,npts):
    out = "{0:0.3F} {1:0.4F} {2:0.4F} {3:0.4F} {4:0.4F} {5:0.4F} {6:0.4F} {7:0.4F} {8:0.4F} {9:0.4F} {10:0.4F} {11:0.4F} {12:0.4F} {13:0.4F} {14:0.4F} {15:0.4F}\n"\
    .format(tt[i],AF[i,0],AF[i,1],AF[i,2],BF[i,0],BF[i,1],BF[i,2],CF[i,0],CF[i,1],CF[i,2],DF[i,0]\
    ,DF[i,1],DF[i,2],EF[i,0],EF[i,1],EF[i,2])
    fp.write(out)
    out1 = "{0:0.3F} {1:0.4F} {2:0.4F} {3:0.4F} {4:0.4F} {5:0.4F} {6:0.4F} {7:0.4F} {8:0.4F} {9:0.4F} {10:0.4F} {11:0.4F} {12:0.4F} {13:0.4F} {14:0.4F} {15:0.4F}\n"\
    .format(tt[i],FF[i,0],FF[i,1],FF[i,2],GF[i,0],GF[i,1],GF[i,2],HF[i,0],HF[i,1],HF[i,2],IF[i,0]\
    ,IF[i,1],IF[i,2],JF[i,0],JF[i,1],JF[i,2])
    fp1.write(out1)
fp.close()
fp1.close()

#plt.show()

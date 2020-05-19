#!/usr/bin/env python3
import sys, os, argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches

## Parse the arguments ##
parser = argparse.ArgumentParser()
parser.add_argument("-n",type=int,default=1,help="No of anisotropic/dipping layers to analyse. Default is 1.")
parser.add_argument("-f",type=str,default="",help="Prefix to the output files of harmonic decomposition. Default is none, which\
    means files to be analysed are modelled.dat and unmodelled.dat.")
parser.add_argument("--horizontal",action="store_true",help="Force to look for horizonal anisotropic layer.")
parser.add_argument("--dip",action="store_true",help="Force to look for dipping isotropic layer.")
parser.add_argument("-z",type=float,default=10.0,help="Zoom level of the plot. Default is 10.")
parser.add_argument("-w",type=float,default=2.0,help="Width of the highlighted region. Default is 2s.")
args = parser.parse_args()
## Argument parsing complete ##
xlimit=[-2,12]  #  Time

modelled_file = args.f+'modelled.dat'
unmodelled_file = args.f+'unmodelled.dat'

if os.path.exists(modelled_file):
    modelled = open(modelled_file).readlines()
else:
    print('File',modelled_file,'doesnot exist')
    sys.exit()
if os.path.exists(unmodelled_file):
    unmodelled = open(unmodelled_file).readlines()
else:
    print('File',unmodelled_file,'doesnot exist')
    sys.exit()
npts = len(modelled)
A = np.zeros([npts,15])
B = np.zeros([npts,15])
tt = np.zeros(npts)
for i in range(0,npts):
    f = modelled[i].split()
    ff = unmodelled[i].split()
    tt[i] = np.float(f[0])
    for j in range(1,16):
        A[i,j-1] = np.float(f[j])
        B[i,j-1] = np.float(ff[j]) 
window = args.w
dt = np.int(1.0/(tt[1]-tt[0]))
m_vec1 = np.sqrt(A[:,4]**2+A[:,7]**2)
m_vec2 = np.sqrt(A[:,10]**2+A[:,13]**2)
um_vec1 = np.sqrt(B[:,4]**2+B[:,7]**2)
um_vec2 = np.sqrt(B[:,10]**2+B[:,13]**2)
if args.horizontal:
    tmp_vec = np.copy(m_vec2)
else:
    tmp_vec = np.copy(m_vec1)
for nlay in range(1,args.n+1):
    print('\v############ Layer',nlay,'#############')
    if nlay>1:
        tmp_vec[n-dt:n+dt]=0.0
    sortindex = np.argsort(tmp_vec)
    sortindex = sortindex[::-1]
    n = sortindex[0]
    t1 = tt[n]
    print('\vDepth to the layer: ',t1,'s')
    ## Estimate alpha ##
    if args.dip:
        alpha = 90.*np.pi/180. - np.arctan2(A[n,4],A[n,7])
        theta = alpha*180.0/np.pi
        print("The direction of dip is {0:0.2f}".format(theta))
    else:
        if args.horizontal:
            alpha = (90.*np.pi/180.- np.arctan2(A[n,10],A[n,13]))/2.0
            theta = alpha*180.0/np.pi
            if theta > 90:
                theta = theta - 180.
            print("\tSymmetry axis direction {0:0.2f} deg".format(theta))
        else:
            #alpha = 90.*np.pi/180.-np.arctan(A[n,4]/A[n,7])
            alpha = 90.*np.pi/180. - np.arctan2(A[n,4],A[n,7])
            R1 = A[n,10] * np.cos(2*alpha)
            R2 = A[n,13] * np.sin(2*alpha)
            print("\v{0:0.4f} {1:0.4f} {2:0.4f}".format(R1,R2,alpha*180./np.pi))
            print('\vIf the signal is from top of the anisotropic layer:')
            if R1 >0 and R2>0:
                theta = 180.+alpha*180.0/np.pi
                if theta > 180:
                    theta = theta - 360 
                print("\tFast axis symmetry with symmetry axis direction {0:0.2f} deg".format(theta))
            else:
                theta = alpha*180.0/np.pi
                if theta > 180:
                    theta = theta - 360
                print("\tSlow axis symmetry with symmetry axis direction {0:0.2f} deg".format(theta))
            
            print('\vIf the signal is from bootom of the anisotropic layer:')
            if R1 >0 and R2>0:
                theta = 180.+alpha*180.0/np.pi
                if theta > 180:
                    theta = theta - 360 
                print("\tSlow axis symmetry with symmetry axis direction {0:0.2f} deg".format(theta))
            else:
                theta = alpha*180.0/np.pi
                if theta > 180:
                    theta = theta - 360
                print("\tFast axis symmetry with symmetry axis direction {0:0.2f} deg".format(theta))

    ## Plotting ##
    maxval = args.z
    fig = plt.figure(figsize=[6,6.0])
    plt.clf()
    ax1 = fig.add_axes([0.025, 0.075, 0.45, 0.875])
    ax2 = fig.add_axes([0.525, 0.075, 0.45, 0.875])
    ax1.set_title('Anisotropy/Dip')
    ax1.set_ylabel('')
    ax1.set_yticks(())
    ax1.set_xlim(xlimit[0], xlimit[1])
    ax1.set_ylim(-0.25, 6.75)
    ax1.set_xlabel('Time (s)')
    ax2.set_title('Unmodelled')
    ax2.set_ylabel('')
    ax2.set_yticks(())
    ax2.set_xlim(xlimit[0], xlimit[1])
    ax2.set_ylim(-0.25, 6.75)
    ax2.set_xlabel('Time (s)')

    i=6
    for j in range(1,15,3):
        data = A[:,j] * maxval
        datal = A[:,j-1] * maxval
        datau = A[:,j+1] * maxval
        ax1.plot(tt[0:npts],i+data,'k-',linewidth=0.5)
        ax1.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 >= 0.,facecolor='cyan', linewidth=0)
        ax1.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 <= 0.,facecolor='cyan', linewidth=0)
        ax1.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 >= 0.,facecolor='red', linewidth=0)
        ax1.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
        i = i-1
    data = m_vec1 * maxval
    ax1.fill_between(tt[0:npts], i, i+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax1.plot(tt[0:npts], i+data,'k-',linewidth=0.5)
    i = i-1
    data = m_vec2 * maxval
    ax1.fill_between(tt[0:npts], i, i+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax1.plot(tt[0:npts], i+data,'k-',linewidth=0.5)

    m=(np.abs(xlimit[1])/12.0)**0.9
    ax1.annotate('constant',xy=(xlimit[1]-3.4*m,6.12))
    ax1.annotate('cos($\phi$)',xy=(xlimit[1]-2.585*m,5.12))
    ax1.annotate('sin($\phi$)',xy=(xlimit[1]-2.5*m,4.12))
    ax1.annotate('cos(2$\phi$)',xy=(xlimit[1]-2.985*m,3.12))
    ax1.annotate('sin(2$\phi$)',xy=(xlimit[1]-2.85*m,2.12))
    ax1.annotate('(N-S)',xy=(xlimit[1]-2.25*m,4.75))
    ax1.annotate('(E-W)',xy=(xlimit[1]-2.3*m,3.75))
    ax1.annotate('(N-S)',xy=(xlimit[1]-2.25*m,2.75))
    ax1.annotate('(E-W)',xy=(xlimit[1]-2.3*m,1.75))
    ax1.annotate('$\sqrt{H1^{2}+H2^{2}}$',xy=(xlimit[1]-4.5*m,1.1))
    ax1.annotate('$\sqrt{H3^{2}+H4^{2}}$',xy=(xlimit[1]-4.5*m,0.1))
    ax1.add_patch(patches.Rectangle((t1-window/2.0,-1),window,10,alpha=0.25))
    i=6
    for j in range(1,15,3):
        data = B[:,j] * maxval
        datal = B[:,j-1] * maxval
        datau = B[:,j+1] * maxval
        ax2.plot(tt[0:npts],i+data,'k-',linewidth=0.5)
        ax2.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 >= 0.,facecolor='cyan', linewidth=0)
        ax2.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 <= 0.,facecolor='cyan', linewidth=0)
        ax2.fill_between(tt[0:npts], i, i+datal, where=datal+1e-6 >= 0.,facecolor='red', linewidth=0)
        ax2.fill_between(tt[0:npts], i, i+datau, where=datau+1e-6 <= 0.,facecolor='dodgerblue', linewidth=0)
        i = i-1
    data = um_vec1 * maxval
    ax2.fill_between(tt[0:npts], i, i+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax2.plot(tt[0:npts], i+data,'k-',linewidth=0.5)
    i = i-1
    data = um_vec2 * maxval
    ax2.fill_between(tt[0:npts], i, i+data, where=data+1e-6 >= 0.,facecolor='red', linewidth=0)
    ax2.plot(tt[0:npts], i+data,'k-',linewidth=0.5)
    ax2.text(xlimit[0]-1.45*m,5.925,'H0',fontsize=14)
    ax2.text(xlimit[0]-1.45*m,4.925,'H1',fontsize=14)
    ax2.text(xlimit[0]-1.45*m,3.925,'H2',fontsize=14)
    ax2.text(xlimit[0]-1.45*m,2.925,'H3',fontsize=14)
    ax2.text(xlimit[0]-1.45*m,1.925,'H4',fontsize=14)
    ax2.annotate('constant',xy=(xlimit[1]-3.4*m,6.15))
    ax2.annotate('cos($\phi$)',xy=(xlimit[1]-2.585*m,5.12))
    ax2.annotate('sin($\phi$)',xy=(xlimit[1]-2.5*m,4.12))
    ax2.annotate('cos(2$\phi$)',xy=(xlimit[1]-2.985*m,3.12))
    ax2.annotate('sin(2$\phi$)',xy=(xlimit[1]-2.85*m,2.12))
    ax2.annotate('(N-S)',xy=(xlimit[1]-2.25*m,4.75))
    ax2.annotate('(E-W)',xy=(xlimit[1]-2.3*m,3.75))
    ax2.annotate('(N-S)',xy=(xlimit[1]-2.25*m,2.75))
    ax2.annotate('(E-W)',xy=(xlimit[1]-2.3*m,1.75))
    ax2.annotate('$\sqrt{H1^{2}+H2^{2}}$',xy=(xlimit[1]-4.5*m,1.1))
    ax2.annotate('$\sqrt{H3^{2}+H4^{2}}$',xy=(xlimit[1]-4.5*m,0.1))
    plt.show()

## FFT Propagation Frensnell based mode
## Messy code 
## MJC 31-7-2017
######################################
## All units are in S.I.
######################################
## DONE
## Converted to work with optical wavelengths
## Added interpolate
## Added SuperGauss for input beam
## Removed bandwidth 
######################################
## TODO
## Powell lens
######################################

import numpy as np
import scipy.interpolate as scipy
import matplotlib.pyplot as plt
n = 2**9
xSize = 100e-3
xs=np.linspace(-xSize/2,xSize/2,n)
ys=xs
c=3e8
Xs,Ys=np.meshgrid(xs,ys)
Rs = np.sqrt(Xs**2+Ys**2)

sigx = 4e-3
sigy = 1.5e-3
offset = 0e-3
SGx = 4
SGy = 2
Ein = np.exp(-((Xs-offset)**2/(2*sigx**2))**SGx)*np.exp(-((Ys-offset)**2/(2*sigy**2))**SGy)

InterpolationRatio = 200

f = 5.6391e14
lam = c/f
k = 2*np.pi/lam
xinterp = xs/InterpolationRatio
yinterp = ys/InterpolationRatio
Xinterp,Yinterp = np.meshgrid(xinterp,yinterp)


focallength = 50e-3
dz = 0.5e-3
zn = 100
zrange=np.linspace(focallength-dz,focallength+dz,zn)

phaseThinLens = -((k*Rs**2)/(2*focallength));
Esx = Ein*np.exp(1j*phaseThinLens);

kxlim = (2*np.pi)/(xs[2]-xs[0])
kx = np.linspace(-kxlim,kxlim,n)
ky=kx
Eout = np.zeros([np.size(xs),np.size(ys),np.size(zrange)],dtype=complex)
zcounter =0
for z in zrange:
    zcounter = int(zcounter)
    xi = z*kx/k;
    yi = z*ky/k;
    Xi,Yi=np.meshgrid(xi,yi)
    Ri2 = Xi**2 + Yi**2
    Ri = np.sqrt(Ri2)
    dxs = xs[2]-xs[1]
    xSoffset = np.min(xs)
    ySoffset = np.min(ys)
    EiTransformx = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(Esx*np.exp(1j*k*Rs**2/(2*(z))))))
    Eix = -1j*(k/(2*np.pi))*(np.exp(1j*k*(z))/z)*(np.exp(1j*k*Ri2/(2*z)))*EiTransformx
    Eix = Eix*np.exp(-1j*k*focallength);
   # Erf = scipy.interpolate.RectBivariateSpline(Xinterp,Yinterp,Eix);
    Erf=Eix
    Eout[:,:,zcounter]=Erf
    plt.figure(2)
    plt.cla()
    plt.pcolormesh(xs,ys,np.abs(Eout[:,:,zcounter])**2)
    plt.title(str(z))
    plt.draw()
    plt.pause(0.001)
    plt.show()
    
    zcounter = zcounter +1 
    
Ex = Eout[int(np.size(xs)/2),:,:]
Ey = Eout[:,int(np.size(ys)/2),:]

plt.subplot(211)
plt.cla()
plt.pcolormesh(zrange,xi,np.abs(Ex)**2)
plt.subplot(212)
plt.cla()
plt.pcolormesh(zrange,yi,np.abs(Ey)**2)
import numpy as np
import os
from scipy.special import j1
import copy
'''
This module calculates the surface brightness of the Roman WFI PSF at a given angular offset from 
a given source. This code is optimized for large angular offsets from the source, and takes a different
approach than STPSF. Rather than numerically calculating the fourier transform of the pupil to create
a template, this code uses analytic expressions for the fourier transform of the aperture and mirror struts.
As such, the surface brightness can be calculated for any arbitrary angular offset efficiently, without having
to oversample the pupil image. The cost of this efficiency is that, due to not modeling the details of
the pupil, it is not as accurate on small scales. It is able to reproduce the diffraction spikes produced
by the mirror struts.

Calling the calculator is simple:

>p0=romanpsf.romanpsf(x_offset_arcsec,y_offset_arcsec,wave,bars=bardict)

where p0 is the normalized surface brightness of the psf at x and y offsets from the source.
The input x and y offsets are assumed to be in the ideal coordinate system (https://roman-docs.stsci.edu/data-handbook/wfi-data-levels-and-products/coordinate-systems) in arcsec
wave is the wavelength for which the PSF is calculated in meters
bardict is a dictionary describing the position-dependent parameters of the mirror struts as seen by the focal plane
additional useage to include the LOLO-diffuser struts is in development

defaultbarparam is a dictionary defined below that sets the default parameters of this dictionary (calculated
for a source at the center of SCA01).


In general, you want to recalculate these parameters for a given source position to account for the parallax
of the mirror struts at different points in the focal plane. The getbarparamtheta method allows for the
straightforward calculation of this. Useage is as follows

getbarparamtheta(thetax,thetay,bardict)

thetax and thetay describe the position of the source in observatory (v2, v3) coordinates (in arcsec)
the method will fill in the necessary parameters to bardict given the offset, assuming the true strut
positions (this method modifies bardict, so store a copy if necessary)

TODO:
account for elliptical distortion of primary, secondary
include diffraction spikes from LOLO-diffuser struts. They should be aligned with the mirror struts in a reference position, but are visible for an offset source
'''


i=complex(0,1)
pi=np.pi
arcsectorad=pi/180/60/60

defaulttruebars={}
defaulttruebars['1']={}
defaulttruebars['2']={}
defaulttruebars['3']={}
defaulttruebars['4']={}
defaulttruebars['5']={}
defaulttruebars['6']={}

defaulttruebars['1']['theta0']=1.25971588e+00
defaulttruebars['1']['thetax1']=-3.47664426e+00
defaulttruebars['1']['thetay1']=-3.22432189e+00
defaulttruebars['1']['thetax2']=-2.86434405e+02
defaulttruebars['1']['thetay2']=3.13288770e+02
defaulttruebars['1']['thetaxy']=-1.61687473e+02

defaulttruebars['2']['theta0']=2.91163892e+00
defaulttruebars['2']['thetax1']=3.30347730e+00
defaulttruebars['2']['thetay1']=6.21256328e+00
defaulttruebars['2']['thetax2']=2.45334483e+02
defaulttruebars['2']['thetay2']=-2.66565225e+0
defaulttruebars['2']['thetaxy']=1.38851363e+02

defaulttruebars['3']['theta0']=2.32505304e+00
defaulttruebars['3']['thetax1']=-2.65408765e-02
defaulttruebars['3']['thetay1']=-1.21373211e+00
defaulttruebars['3']['thetax2']=8.84447477e+01
defaulttruebars['3']['thetay2']=3.51039495e+01
defaulttruebars['3']['thetaxy']=6.76308874e+01

defaulttruebars['4']['theta0']=8.21877192e-01
defaulttruebars['4']['thetax1']=-1.24926163e+00
defaulttruebars['4']['thetay1']=-2.89823915e-01
defaulttruebars['4']['thetax2']=1.62371885e+01
defaulttruebars['4']['thetay2']=-9.04770843e+01
defaulttruebars['4']['thetaxy']=-5.48434331e-01

defaulttruebars['5']['theta0']=2.13975696e-01
defaulttruebars['5']['thetax1']=2.35089534e+00
defaulttruebars['5']['thetay1']=6.49051207e+00
defaulttruebars['5']['thetax2']=3.11491054e+02
defaulttruebars['5']['thetay2']=-1.78068099e+02
defaulttruebars['5']['thetaxy']=1.97766241e+02

defaulttruebars['6']['theta0']=1.89786968e+00
defaulttruebars['6']['thetax1']=6.34653878e-01
defaulttruebars['6']['thetay1']=-5.17084440e+00
defaulttruebars['6']['thetax2']=-3.48991711e+02
defaulttruebars['6']['thetay2']=1.80665508e+02
defaulttruebars['6']['thetaxy']=-2.24039326e+02

#center of chip1
thetax0=1312.9491452484797*arcsectorad
thetay0=-1040.7853726755036*arcsectorad

defaultbarparam={}
defaultbarparam['1']={}
defaultbarparam['2']={}
defaultbarparam['3']={}
defaultbarparam['4']={}
defaultbarparam['5']={}
defaultbarparam['6']={}

#effective length of optical path in m
defaultbarparam['zref']=18.813

#width of secondary mirror struts
defaultbarparam['1']['w']=7.46179521e-02
defaultbarparam['2']['w']=7.41656285e-02
defaultbarparam['3']['w']=7.41208867e-02
defaultbarparam['4']['w']=7.46782500e-02
defaultbarparam['5']['w']=7.44190186e-02
defaultbarparam['6']['w']=7.45386249e-02

#length of secondary mirror struts
defaultbarparam['1']['l']=9.79994893e-01
defaultbarparam['2']['l']=9.89084836e-01
defaultbarparam['3']['l']=9.76846302e-01
defaultbarparam['4']['l']=9.66738093e-01
defaultbarparam['5']['l']=9.74839045e-01
defaultbarparam['6']['l']=9.73334491e-01

#centroid of secondary mirror struts
#just results in a phase offset
defaultbarparam['1']['x0']=1.4335505999999998
defaultbarparam['1']['y0']=2.72694435

defaultbarparam['2']['x0']=3.08062875
defaultbarparam['2']['y0']=-0.05643900000000002

defaultbarparam['3']['x0']=1.6969326
defaultbarparam['3']['y0']=-2.49930705

defaultbarparam['4']['x0']=-1.70916105
defaultbarparam['4']['y0']=-2.4588590999999997

defaultbarparam['5']['x0']=-3.09379785
defaultbarparam['5']['y0']=-0.023516249999999933

defaultbarparam['6']['x0']=-1.6470781499999998
defaultbarparam['6']['y0']=2.74011345


#normalized amplitude of diffraction spikes
#investigation involving JWST PSFs finds that the diffraction spikes
#have a slightly lower amplitude than one would expect for their geometry
#this may be due to wierdness in the multi-plane diffraction or some other
#characteristic of the struts. For now, scaling the diffraction spikes by
#a certain amount seems to work.
#This technically messes with the normalization but not by a noticable
#amount because the struts are a small fraction of the area

defaultbarparam['1']['amp']=1
defaultbarparam['2']['amp']=1
defaultbarparam['3']['amp']=1
defaultbarparam['4']['amp']=1
defaultbarparam['5']['amp']=1
defaultbarparam['6']['amp']=1

def getbarparamtheta(thetaxi,thetayi,barparam,truebars=defaulttruebars):
    #get the apparent projected angle of the mirror struts
    #given the angular offset from the telescope axis
    
    #do quadratic x,y polynomial for theta instead of full geometric solution
    #the geometry must be just a bit too complicated

    dx=thetaxi-thetax0
    dy=thetayi-thetay0
    for j in range(6):
        bari=str(j+1)
        barparam[bari]['theta']=truebars[bari]['theta0']+truebars[bari]['thetax1']*dx+truebars[bari]['thetay1']*dy+truebars[bari]['thetax2']*dx**2+truebars[bari]['thetay2']*dy**2+truebars[bari]['thetaxy']*dx*dy

    
def getbardict(thetax,thetay,truebars=defaulttruebars):
    '''
    Thetas in radians.
    '''
    #same as above, but returns a new dict instead of modifying an existing one
    barparam=copy.deepcopy(defaultbarparam)
    
    dx=thetax-thetax0
    dy=thetay-thetay0

    for j in range(6):
        bari=str(j+1)
        barparam[bari]['theta']=truebars[bari]['theta0']+truebars[bari]['thetax1']*dx+truebars[bari]['thetay1']*dy+truebars[bari]['thetax2']*dx**2+truebars[bari]['thetay2']*dy**2+truebars[bari]['thetaxy']*dx*dy

    return barparam

getbarparamtheta(thetax0,thetay0,defaultbarparam,truebars=defaulttruebars)


def barpsf(x,y,k,w,l,x0,y0,rotation=0):
    #calculates the fourier transform pattern of a bar with a given center and rotation

    #rotate x and y to a coordinate system aligned with the bar
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)

    #the resulting Fourier transform is a sinc in x, and y, with phases shifted to the central x and y
    return np.sinc(k*w*xrot/2/np.pi)*np.sinc(k*l*yrot/2/np.pi)*np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)*l*w


def romanpsf(xoffset,yoffset,wave,barparam=defaultbarparam,x0=0,y0=0,flength=18.813,ainner=.783/2,aouter=2.553/2,pixsize=10E-6):
    #calculates the PSF of Roman's WFI
    
    #PSF is evaluated with x and y in meters
    #so convert angular offset to physical offset using the optical path length
    #none of the individual components need the optical path length, since it cancels out
    #however, the interference between the bar psfs and the mirror psfs does need it
    #since the bar psfs are phase-shifted relative to the mirror psfs
    x=xoffset*arcsectorad*flength
    y=yoffset*arcsectorad*flength

    #theta_perp for disk fourier transform
    theta=np.sqrt(x**2+y**2)
    
    k=2*np.pi/wave/flength

    #fourier transform of disk is J1 bessel function divided by k*a*theta_perp (https://adriftjustoffthecoast.wordpress.com/2013/06/06/2d-fourier-transform-of-the-unit-disk/)
    mirrorpsfouter=np.array(j1(k*aouter*theta)/(k*aouter*theta),dtype=complex)*2*np.pi*aouter**2
    mirrorpsfinner=np.array(j1(k*ainner*theta)/(k*ainner*theta),dtype=complex)*2*np.pi*ainner**2

    if type(theta)==np.ndarray:
        mirrorpsfouter[theta==0]=np.pi*aouter**2
        mirrorpsfinner[theta==0]=np.pi*ainner**2
    else:
        if theta==0:
            mirrorpsfouter=np.pi*aouter**2+0*i
            mirrorpsfinner=np.pi*ainner**2+0*i
    #subtract inner mirror from the outer
    final=mirrorpsfouter-mirrorpsfinner

    barpsfs=np.zeros_like(final)
    #subtract the bars
    bararea=[]
    for j in range(6):
        barpsfj=barparam[str(j+1)]['amp']*barpsf(x,y,k,barparam[str(j+1)]['w'],barparam[str(j+1)]['l'],barparam[str(j+1)]['x0'],barparam[str(j+1)]['y0'],rotation=barparam[str(j+1)]['theta'])
        barpsfs+=barpsfj
        #final-=barparam[str(j+1)]['amp']*barpsf(x,y,k,barparam[str(j+1)]['w'],barparam[str(j+1)]['l'],barparam[str(j+1)]['x0'],barparam[str(j+1)]['y0'],rotation=barparam[str(j+1)]['theta'])
        #keep track of the area of the bars for normalization
        bararea.append(barparam[str(j+1)]['w']*barparam[str(j+1)]['l'])
        
    final-=barpsfs
    #alow for a phase shift to a given x0 and y0 in case this is passed through more optics
    #but don't do it by default
    if not ((x0==0) and (y0==0)):
        final*=np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)

    #normalize the final pattern by k/2/pi/sqrt(area)*pixelsize
    pixscale=pixsize/flength/arcsectorad
    final*=k/2/np.pi/np.sqrt(np.pi*(aouter**2-ainner**2)-np.sum(bararea))*pixsize/pixscale
    
    return final

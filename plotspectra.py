#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 12:39:55 2024

@author: spec
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits
import pandas as pd



def get_spectra(fitsfile):
    '''
    function to obtain the values of the spectra
    '''
    img_data, img_header = fits.getdata(fitsfile, header=True)
    cdelta1 = img_header['CDELT1']
    crval1  = img_header['CRVAL1']
    npoints = img_header['NAXIS1']
    ll = np.arange(0,npoints)*cdelta1+crval1
    return ll, img_data

def plot_spectra(fitsfile):
    ll, img_data = get_spectra(fitsfile)
    # the following code is used in case we want to see only a part of the spectrum 
    '''
    for l in range(5000,6800,10):
        print(l,l+10)
        index = np.where(ll>l)[0][0]
        indexg = np.where(llg>l)[0][0]
        plt.plot(ll[index:index+1000],img_data[index:index+1000]/max(img_data), label='Star')
        plt.plot(llg[indexg:indexg+1000],img_datag[indexg:indexg+1000]/max(img_datag), label ='Sun')
        plt.legend()
        plt.xlabel('Wavelength ($\mathrm{\AA}$)')
        plt.ylabel('Normalized flux')
        plt.savefig('/home/spec/WORK/PEEC2024/plots/'+fitsfile.replace('/home/spec/WORK/PEEC2024/Spectra_rv/','')+'_'+str(l)+'.jpg')
        plt.show()
        #plt.plot(ll[index-5000:index+5000],img_data[index-5000:index+5000])
    '''
    # the following code shows the whole spectrum
    plt.plot(ll,img_data)
    plt.title(fitsfile.replace('/home/spec/WORK/PEEC2024/spectra_vsin2_4/',''))
    #plt.savefig('/home/spec/WORK/PEEC2024/plots/'+fitsfile.replace('/home/spec/WORK/PEEC2024/spectra_vsin2_4/','')+'.jpg')
    plt.show()
    return fitsfile

'list of spectra to correct:'
spectra_ = os.listdir('/home/spec/WORK/PEEC2024/spectra_vsin2_4')

directory="exp_stars_information.csv"
Table=pd.read_csv(directory)
llg,img_datag = get_spectra('/home/spec/Programs/ARES/sun_harps_ganymede.fits')
'plot'
for i in range (0,len(Table)):
    #if 'fits' in spectra_[i]:
    plot_spectra('/home/spec/WORK/PEEC2024/spectra_vsin2_4/'+Table["fits_name"][i].replace("_rv.fits",".fits"))  #spectra not corrected in rv
    plot_spectra('/home/spec/WORK/PEEC2024/Spectra_rv/'+Table["fits_name"][i]) #spectra corrected in rv
    plt.show()

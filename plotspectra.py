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

def get_spectra(fitsfile):
    img_data, img_header = fits.getdata(fitsfile, header=True)
    cdelta1 = img_header['CDELT1']
    crval1  = img_header['CRVAL1']
    npoints = img_header['NAXIS1']
    ll = np.arange(0,npoints)*cdelta1+crval1
    return ll, img_data

def plot_spectra(fitsfile):
    ll, img_data = get_spectra(fitsfile)
    plt.plot(ll,img_data)
    plt.title(fitsfile)
    plt.show()
    return fitsfile

'list of spectra to correct:'
spectra_ = os.listdir('/home/spec/WORK/PEEC2024/spectra_vsin2_4')


'correction'
#for i in range (0,len(spectra_)):
 #   if 'fits' in spectra_[i]:
plot_spectra('/home/spec/WORK/PEEC2024/spectra_vsin2_4/'+'CoRoT-32_SOPHIE_HE_2020.fits')
plt.show()
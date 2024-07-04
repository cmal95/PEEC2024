#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:15:13 2024

@author: spec
"""

import numpy as np
import pandas as pd
from astropy.io import fits
import os
from PyAstronomy import pyasl
import vsini_code

PATH_TO_SAVE='/home/spec/WORK/PEEC2024/Spectra_rv/'

def get_spectra(fitsfile):
    img_data, img_header = fits.getdata(fitsfile, header=True)
    cdelta1 = img_header['CDELT1']
    crval1  = img_header['CRVAL1']
    npoints = img_header['NAXIS1']
    ll = np.arange(0,npoints)*cdelta1+crval1
    return ll, img_data, cdelta1


def treat_rv(name_file,l1,l2,ref):
    
    """
    This function does the correction of the spectrum in radial velocity.
    Parameters:
    -name_file: The name of the file to correct the spectrum.
    -l1: Value of the initial wavelenght
    -l2: value of the final wavelenght
    -ref: reference spectrum
    """
    ll_r, flux_r, snr = get_spectra(name_file)
    ll_r0, flux_r0 , snr0= get_spectra(ref)   #Using a rest frame solar spectrum defined on top
    ini = np.where(ll_r>l1)[0][0]
    inf = np.where(ll_r>l2)[0][0]
    ini2 = np.where(ll_r0>l1)[0][0]
    inf2 = np.where(ll_r0>l2)[0][0]
    #ll_r, flux_r = ll_r[ini:inf],flux_r[ini:inf]
    #ll_r0, flux_r0 = ll_r0[ini2:inf2], flux_r0[ini2:inf2]
    rv, cc = pyasl.crosscorrRV(ll_r[ini:inf], flux_r[ini:inf], ll_r0[ini2:inf2], flux_r0[ini2:inf2], -100., 100., 0.1, skipedge=300)
    maxind = np.argmax(cc)
    rv = rv[maxind]
    print("REF0: Cross-correlation function is maximized at dRV = ", rv, " km/s")
    
    c=299792.458 ### light velocity
    
    File_op=fits.open(name_file)
    Data=File_op[0].data
    Header=File_op[0].header

    nome0=os.path.basename(name_file.split(".fits")[0])

    File_op.close()
    wave_start = Header['CRVAL1']
    wave_delta = Header['CDELT1']
    wave = np.arange(Data.size) * wave_delta + wave_start

    wave_corr=wave/(1+rv/c)

    wave_int=np.arange(wave_corr[0],wave_corr[-1],wave_delta)

    Flux_int=np.interp(wave_int,wave_corr,Data)

    Header["CDELT1"]=wave_delta
    Header["CRVAL1"]=wave_int[0]
    Header["COMMENT"]= "Corrigido na velocidade radial."
    Header["RV"] = rv
    New_name=nome0 + "_rv.fits"
    print(nome0)

    fits.writeto(PATH_TO_SAVE + New_name,Flux_int,Header,overwrite=True)
    
'reference spectrum:'
ref = '/home/spec/Programs/ARES/sun_harps_ganymede.fits'

'list of spectra to correct:'
spectra_ = os.listdir('/home/spec/WORK/PEEC2024/spectra_vsin2_4')

directory="/home/spec/WORK/vsinimoog24/exp_stars_information.csv"
Table=pd.read_csv(directory)

'correction'
#for i in range (0,len(Table)):
    #if 'fits' in spectra_[i]:
treat_rv('/home/spec/WORK/vsinimoog24/spectra/'+Table["fits_name"][0].replace("_rv.fits",".fits"),5000,6800, ref)
#vsini_code.main()
 
'''for l in range(5000,6800,10):
    treat_rv('/home/spec/WORK/PEEC2024/spectra_vsin2_4/'+Table["fits_name"][0].replace("_rv.fits",".fits"),l,l+10, ref)
    vsini_code.main(p  )
    print(l,l+10)
'''



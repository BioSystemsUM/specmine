#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 15:04:11 2017

@author: scardoso
"""

def read_varian_spec_raw(spectrum_directory, fid_filename, procpar_filename, zero_filling=True, apodization=True):
    from nmrglue.fileio import varian, convert
    from nmrglue.process import pipe_proc
    
    #Read varian files:
    res=varian.read(dir=spectrum_directory, fid_file=fid_filename, procpar_file=procpar_filename)
    varian_data=res[1]
    varian_dic=res[0]    
    
    #Get main parameters for ppm scale:
    sw=float(varian_dic["procpar"]["sw"]["values"][0])
    obs=float(varian_dic["procpar"]["sfrq"]["values"][0])
    car=float(varian_dic["procpar"]["reffrq"]["values"][0])
    
    #Convert varian to pipe:
    universal_varian_dic=varian.guess_udic(varian_dic, varian_data)
    universal_varian_dic[0]["sw"]=sw
    universal_varian_dic[0]["obs"]=obs
    universal_varian_dic[0]["label"]=varian_dic["procpar"]["tn"]["values"][0]
    universal_varian_dic[0]["car"]=car
    C = convert.converter() 
    C.from_varian(varian_dic, varian_data, universal_varian_dic) 
    pipe_dic, varian_pipe_data = C.to_pipe()
    pipe_dic["FDF1SW"]=sw
    pipe_dic["FDF1OBS"]=obs
    pipe_dic["FDF1LABEL"]=universal_varian_dic[0]["label"]
    pipe_dic["FDF1CAR"]=car
    
    #If zero filling:
    if zero_filling:
        pipe_dic, varian_data=pipe_proc.zf(pipe_dic, varian_data)
    
    #If apodization:
    if apodization:
        lb_d=float(varian_dic["procpar"]["lb"]["values"][0])
        pipe_dic, varian_data=pipe_proc.em(pipe_dic, varian_data, lb=lb_d)
    
    #Fourier Transform:
    pipe_dic, varian_data=pipe_proc.ft(pipe_dic, varian_data, auto=True, inv=True)    
    
    #Phase correction:
    p_zero=float(varian_dic["procpar"]["rp"]["values"][0])
    p_one=float(varian_dic["procpar"]["lp"]["values"][0])
    pipe_dic, varian_data=pipe_proc.ps(pipe_dic, varian_data, p0=p_zero, p1=p_one, inv=True)

    #Baseline Correction:
    pipe_dic, varian_data=pipe_proc.cbf(pipe_dic, varian_data)

    #Remove imaginary numbers:    
    dic, data = pipe_proc.di(pipe_dic, varian_data)
    
    #Calculate ppm scale:
    ppm=list()
    ppm_f=(obs*1000000 + sw/2 - (car*1000000)) / car
    ppm_width=sw/car
    ppm_i=ppm_f-ppm_width
    n=int(len(data))
    ppm_step=ppm_width/n
    ppm.append(ppm_i)
    for i in range(1, n):
        ppm.append(ppm[i-1]+ppm_step)

    return (ppm, abs(data))

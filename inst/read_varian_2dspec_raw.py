# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 18:33:50 2020

@author: Bruno
"""

def read_varian_spec2d_raw(spectrum_directory,fid_filename,procpar_filename,zero_filling=True,apodization=True):
    from nmrglue.fileio import varian, convert
    from nmrglue.process import pipe_proc

    # Read varian files
    res = varian.read(dir=spectrum_directory,fid_file=fid_filename,procpar_file=procpar_filename)
    varian_data = res[1]
    varian_dic = res[0]

    # Get main parameters for ppm scale:
    sw = float(varian_dic["procpar"]["sw"]["values"][0]) #direct dimension
    sw1 = float(varian_dic["procpar"]["sw1"]["values"][0]) #indirect dimension

    obs = float(varian_dic["procpar"]["sfrq"]["values"][0]) #direct dimension
    obs1 = float(varian_dic["procpar"]["dfrq"]["values"][0]) #indirect dimension

    car = float(varian_dic["procpar"]["reffrq"]["values"][0]) #direct dimension
    car1 = float(varian_dic["procpar"]["reffrq1"]["values"][0]) #indirect dimension


    # Convert varian to pipe
    universal_varian_dic = varian.guess_udic(varian_dic,varian_data)
    
    ## direct dimension
    universal_varian_dic[0]["sw"] = sw
    universal_varian_dic[0]["obs"] = obs
    universal_varian_dic[0]["label"] = varian_dic["procpar"]["tn"]["values"][0]
    universal_varian_dic[0]["car"] = car
    
    ## indirect dimension
    universal_varian_dic[1]["sw"] = sw1
    universal_varian_dic[1]["obs"] = obs1
    universal_varian_dic[1]["label"] = varian_dic["procpar"]["dn"]["values"][0]
    universal_varian_dic[1]["car"] = car1


    C = convert.converter()
    C.from_varian(varian_dic,varian_data,universal_varian_dic)
    pipe_dic, varian_pipe_data = C.to_pipe()

    ## direct dimension
    pipe_dic["FDF2SW"] = sw
    pipe_dic["FDF2OBS"] = obs
    pipe_dic["FDF2LABEL"] = universal_varian_dic[0]["label"]
    pipe_dic["FDF2CAR"] = car
    
    ## indirect dimension
    pipe_dic["FDF1SW"] = sw1
    pipe_dic["FDF1OBS"] = obs1
    pipe_dic["FDF1LABEL"] = universal_varian_dic[1]["label"]
    pipe_dic["FDF1CAR"] = car1
    
    # process the direct dimension
    ## zero filling
    if zero_filling:
        pipe_dic, varian_pipe_data = pipe_proc.zf(pipe_dic,varian_pipe_data)
    ## apodization
    if apodization:
        lb_d = float(varian_dic["procpar"]["lb"]["values"][0])
        pipe_dic, varian_pipe_data = pipe_proc.em(pipe_dic,varian_pipe_data, lb = lb_d)
        # for lorentz-to-gauss
        # pipe_dic, varian_pipe_data = pipe_proc.gm(pipe_dic, varian_pipe_data)
    ## Fourier transform
    pipe_dic, varian_pipe_data = pipe_proc.ft(pipe_dic,varian_pipe_data,auto=True)
    ## Phase Correction
    p_zero = float(varian_dic["procpar"]["rp"]["values"][0])
    p_one = float(varian_dic["procpar"]["lp"]["values"][0])
    pipe_dic,varian_pipe_data = pipe_proc.ps(pipe_dic,varian_pipe_data,p0=p_zero,p1=p_one)
    ## Remove imaginary numbers:
    pipe_dic, varian_pipe_data = pipe_proc.di(pipe_dic,varian_pipe_data)
    
    # process the indirect dimension
    pipe_dic, varian_pipe_data = pipe_proc.tp(pipe_dic,varian_pipe_data)
    ## zero filling
    if zero_filling:
        pipe_dic, varian_pipe_data = pipe_proc.zf(pipe_dic,varian_pipe_data)
    ## apodization
    if apodization:
        lb_d1 = float(varian_dic["procpar"]["lb1"]["values"][0])
        pipe_dic, varian_pipe_data = pipe_proc.em(pipe_dic,varian_pipe_data,lb = lb_d1)
        # for lorentz-to-gauss
        # pipe_dic, varian_pipe_data = pipe_proc.gm(pipe_dic, varian_pipe_data)
    ## Fourier transform
    pipe_dic,varian_pipe_data = pipe_proc.ft(pipe_dic,varian_pipe_data,auto=True)
    ## Phase Correction
    p_zero1 = float(varian_dic["procpar"]["rp1"]["values"][0])
    p_one1 = float(varian_dic["procpar"]["lp1"]["values"][0])
    pipe_dic, varian_pipe_data = pipe_proc.ps(pipe_dic,varian_pipe_data,p0=p_zero1,p1=p_one1)
    ## Remove imaginary numbers:
    pipe_dic, varian_pipe_data = pipe_proc.di(pipe_dic,varian_pipe_data)
    
    dic, data = pipe_proc.tp(pipe_dic,varian_pipe_data)

    # Calculate both ppm scales:
    ## direct dimension
    ppm = list()
    ppm_f = (obs*1000000 + sw/2 - (car*1000000)) / car
    ppm_width = sw/car
    ppm_i = ppm_f - ppm_width
    n = int(data.shape[1])
    ppm_step = ppm_width/n
    ppm.append(ppm_i)
    for i in range(1,n):
        ppm.append(ppm[i-1]+ppm_step)
    new_ppm = ppm[::-1]
    
    ## indirect dimension
    # handling homonuclear cases
    if car1 == car and sw == sw1:
        ppm1 = list()
        ppm1.append(ppm_i)
        n1 = int(data.shape[0])
        ppm_step1 = ppm_width/n1
        for i in range(1, n1):
            ppm1.append(ppm1[i - 1] + ppm_step1)
        new_ppm1 = ppm1[::-1]
    else:
        ppm1 = list()
        ppm_f1 = (obs1*1000000 + sw1/2 - (car1*1000000)) / car1
        ppm_width1 = sw1/car1
        ppm_i1 = ppm_f1 - ppm_width1
        n1 = int(data.shape[0])
        ppm_step1 = ppm_width1/n1
        ppm1.append(ppm_i1)
        for i in range(1,n1):
            ppm1.append(ppm1[i-1]+ppm_step1)
        new_ppm1 = ppm1[::-1]

    ppms = [new_ppm, new_ppm1]
    
    return (ppms,abs(data))
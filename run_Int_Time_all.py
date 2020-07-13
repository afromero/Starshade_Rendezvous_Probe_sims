#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 14:28:38 2018

@author: romerowo
"""

from pylab import *
import os
import time

HIP_list = [
            8102,  # tau Ceti
            37279, # Procyon A
            108870,# epsilon Indi A
            32349, # Sirius A
            19849, #omicron 2 Eridani
            97649, # Altair
            99240, # delta Pavonis
            15510, #82 Eridani
            96100, #sigma Draconis
            2021,  # beta Hyi
            61317,  #beta CVn
            22449, #1 Ori
            113368, #Fomalhaut
            17378, #Delta Eridani
            78072,  #Gamma Serpentis
            1599,   #Zeta Tucanae
            27072, # Gamma Leporis
            64394,  #Beta Comae Berenices
            105858, #Gamma Pavonis
            14632,  #24 Capricorni
            67927,  #Eta Bootis
            91262,  #Vega
            57757,  #Beta Virginis
            86974, #Mu Herculis
            16537, # epsilon Eridani
            37826]  #Pollux
            

outdir = 'pub_runs'
count = 0
for hip in HIP_list:
    print(hip)
    com = './Integration_Time_full_spec.py -hip %d -outdir %s > %s/star_%d.txt'%(hip, outdir, outdir, hip)
    print(com)
    os.system(com)
    time.sleep(10.0)
    count+=1
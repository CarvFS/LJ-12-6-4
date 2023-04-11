'''
Script to run calculations using a set of ions:
 - use 1* while setting 'entropicDecomp = 0' or 2* for 'entropicDecomp = 1'
'''

import subprocess
import os
import pandas as pd
import scipy.signal
import numpy as np

# Functions to read data files, get first maximum and minimum positions

def read_vv(filename):
    return pd.read_fwf(filename, 
                     skiprows=3, 
                     infer_nrows=200).drop('#', axis=1)

def first_max(df, separation = False):
    maxima = df.drop('SEPARATION', axis=1).apply(lambda x: scipy.signal.argrelextrema(x.values, np.greater))
    if separation:
        maxima = maxima.apply(lambda x: x[0][0])
        df = df.drop('SEPARATION', axis=1).apply(lambda x: 
                      df.SEPARATION.iloc[maxima[x.name]])
        return df
    else:
        return maxima.apply(lambda x: x[0][0])

def first_min(df, separation = False):
    minima = df.drop('SEPARATION', axis=1).apply(
        lambda x: scipy.signal.argrelextrema(x.values, np.less))
    if separation:
        minima = minima.apply(lambda x: x[0][0])
        df = df.drop('SEPARATION', axis=1).apply(lambda x: 
                      df.SEPARATION.iloc[minima[x.name]])
        return df
    else:
        return minima.apply(lambda x: x[0][0])

###############################################################################################################

ions = ['Na+','MG','Cl-']

r = open("results","w+")

for i in ions:
    # write make_mdl script
    f = open("make_mdl_test.sh","w+")
    f.write("python $HOME/ambermd/amber.git/AmberTools/src/rismtools/lj-12-6-4/MDL_script/MDL_Generator.py "
    "-lrc leaprc.water.spce "
    "-frc frcmod.ions1lm_1264_spce frcmod.ions234lm_1264_spce "
    "-seq WAT %s "
    "-fn spce-ion "
    "-target %s \n" % (i,i))
    f.close()

    # run make_mdl script
    subprocess.call(['sh', './make_mdl_test.sh'])

    # run files to get rism1d results
    subprocess.call(['sh', './scriptLJ.sh'])

    # read data files

    df = read_vv('pse3.gvv')
    nvv = read_vv('pse3.nvv')

    # find ion-O distance

    iod = nvv.loc[first_max(df)["%s:O" % i]]["SEPARATION"]

    # find coordination number

    cn = nvv.loc[first_min(df)["%s:O" % i]]["%s:O" % i]

    # calculate hydration free energy

    # 1*
    headers = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 7, nrows=0).columns[1:]
    therm = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 8, nrows = 4, header = None, names = headers)
    hfe = therm.loc['Excess_chemical_potential_PR'][["%s" % i]]["%s" % i]

    # 2*
    #headers = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 8, nrows=0).columns[1:]
    #therm = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 9, nrows = 8, header = None, names = headers)
    #dh = therm.loc['Solvation_energy'][["%s" % i]]["%s" % i]
    #mtds =therm.loc['-Temperature*solvation_entropy_PR'][["%s" % i]]["%s" % i]
    #hfe = dh + mtds
    
    # write results
    r.write("%s: HFE = %s, IOD = %s, CN = %s \n" % (i,str(hfe),str(iod),str(cn)))

    # clean files to run another ion (should also work by cleaning only the .sav files)
    os.system('rm kh.* pse2.* pse3.*')

    
r.close()
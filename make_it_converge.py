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

############################ monovalent ions
# CS is not converging, yet. Add NH4+, H+'s and H3O+
#ions = ['LI','Na+','K+','RB','TL','CU1','AG','F','Cl-','BR','IOD'] # 12-6-4 ions***lm

#ions = ['LI','Na+','K+','RB','AG','F','Cl-','BR'] # 12-6 ions***lm

############################# divalent ions
# Be, FE2, HG, CA and BA are not converging, yet
#ions = ['CU', 'NI', 'ZN', 'CO', 'Cr', 'MG', 'V2+', 'MN', 'CD', 'Sn', 'SR'] # 12-6-4 ions***lm

############################ trivalent ions
# AL, FE, CR, IN, Tl, Y, PR, Nd, SM, EU3, GD3, TB, Dy, Er, Tm and LU are not converging, yet
# ions = ['LA', 'CE']

ions = ['AL']

r = open("make_it_converge.out", "w+")
#r = open("results-LJ-12-6-4-tri-LM.out","w+")
#r = open("results-LJ-12-6-dival.out","w+")
r.write("        HFE,      IOD,      CN \n")
for i in ions:
    # write make_mdl script for LJ 12-6-4
    
    f = open("make_mdl_test.sh","w+")
    f.write("python $HOME/ambermd/amber.git/AmberTools/src/rismtools/lj-12-6-4/MDL_script/MDL_Generator.py "
    "-lrc leaprc.water.spce "
    "-frc frcmod.ions1lm_1264_spce frcmod.ions234lm_1264_spce "
    "-seq WAT %s "
    "-fn spce-ion "
    "-target %s \n" % (i,i))
    f.close()
    
    # write make_mdl script for LJ 12-6

   # f = open("make_mdl_test.sh","w+")
   # f.write("python $HOME/ambermd/amber.git/AmberTools/src/rismtools/lj-12-6-4/MDL_script/MDL_Generator.py "
   # "-lrc leaprc.water.spce "
   # "#-frc frcmod.ionsjc_spce "
   # "-frc frcmod.ions1lm_126_spce frcmod.ions234lm_126_spce "
   # "-seq WAT %s "
   # "-fn spce-ion " % i) 4
   # f.close()

    # run make_mdl script
    subprocess.call(['sh', './make_mdl_test.sh'])

    # clean files kh.*, pse2.*, pse3.*
    os.system('rm kh.* pse2.* pse3.*')

    closure_list = ['kh','pse2','pse3']
    for indx,c in enumerate(closure_list):
        
        if indx == 0:
            D = [40.00,50.00,55.26]
        else:
            D = [55.26]
        
        #D=[40.00,50.00,55.26]
        
        for indx2,dens in enumerate(D):
            #print("%d" % temp)
            #input("press to continue")
            f = open("Run.mir.%s" % c, "w+")
            f.write("#!/bin/sh \n")
            f.write("\n")
            if indx>0:
                f.write("prev_closure=\"%s\" \n" % closure_list[indx-1])
            f.write("closure=\"%s\" \n" % c)
            f.write("cat > $closure.inp << EOF \n")
            f.write("&PARAMETERS \n")
            f.write("outlist = 'uxghcbtensq', THEORY='DRISM', closure='${closure^^}', \n")
            f.write("exchem_SC=1,exchem_SM=1 \n")
            f.write("extra_precision=0 \n")
            f.write("entropicDecomp=0 \n")
            f.write("!grid \n")
            f.write("NR=16384, DR=0.025, rout=100, kout=30., \n")
            f.write("!MDIIS \n")
            f.write("mdiis_nvec=50, mdiis_del=0.3, tolerance=1.e-9, \n")
            f.write("!iter \n")
            f.write("ksave=-1, progress=1, maxstep=10000,\n")
            f.write("!EStat \n")
            f.write("SMEAR=1, ADBCOR=0.5, \n")
            f.write("!bulk solvent properties \n")
            f.write("temperature=298,DIEps=78.497, \n")
            f.write("NSP=1 \n")
            f.write("/ \n")
            f.write("&SPECIES \n")
            f.write("ndens=2 \n")
            f.write("DENSITY=%f,0.0d0 \n" % dens)
            f.write("model=\"spce-ion.mdl\" \n")
            f.write("/ \n")
            f.write("EOF\n")
            f.write("\n")
            if indx>0:
                f.write("cp $prev_closure.sav $closure.sav \n")
            f.write("$AMBERHOME/bin/rism1d $closure")
            f.close()
    
            # run files to get rism1d results
            #subprocess.call(['sh', './scriptLJ.sh'])
            os.system('./Run.mir.%s' % c)
            
'''
    # get hydration free energy
    
    # 1*
    headers = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 7, nrows=0).columns[1:]
    therm = pd.read_csv('pse3.therm', delim_whitespace=True, skiprows = 8, nrows = 4, header = None, names = headers)
    hfe = therm.loc['Excess_chemical_potential_PR'][["%s" % i]]["%s" % i]

    # read data files
    
    if i == 'IOD':
        i = 'I'
    elif i == 'CU1':
        i = 'CU'

    df = read_vv('pse3.gvv')
    nvv = read_vv('pse3.nvv')

    # find ion-O distance

    iod = nvv.loc[first_max(df)["%s:O" % i]]["SEPARATION"]

    # find coordination number

    cn = nvv.loc[first_min(df)["%s:O" % i]]["%s:O" % i]
    
    # write results
    r.write("%s: %s, %s, %s \n" % (i,str(hfe),str(iod),str(cn)))

    # clean files to run another ion (should also work by cleaning only the .sav files)
    #os.system('rm kh.* pse2.* pse3.*')

'''    
r.close()

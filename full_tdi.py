# full_tdi.py
'''Script t  consolidate and analyze classical MonteCarlo and DFT thermodynamic integration 
using Pandas dataframes'''

from parse_vasp import *
import os
import subprocess
import pandas as pd
import numpy as np
import socket
import os

from subprocess import Popen, PIPE
import time

host = socket.gethostname()

runDir="/u/smwahl/scr/workspace"
saveDir="/u/smwahl/dat/"

# analysis scripts
runinfo="/u/smwahl/scripts/vasp_automation/runinfo"
eos="/u/smwahl/scripts/eos1"
cmc_result="/u/smwahl/code/python/vasp_automation/cmc_result"
cmc_einstein="/u/smwahl/code/python/vasp_automation/cmc_einstein"
lambda_cci="/u/smwahl/scripts/vasp_automation/lambda_cci"

# Run directories
cmcStr = "lFeSiO160 lSiO161 lSiO162 lFeSiO161 lFeSiO162 lSiO163 lFeSiO163 lSiO164 lFeSiO164 lSiO165 lFeSiO165 lSiO166 lSiO167 lFeSiO166"

cmcEinsteinStr = ""
cmcEinsteinRuns = cmcEinsteinStr.split()
cmcPPRuns = cmcStr.split()
cmcRuns = cmcEinsteinRuns + cmcPPRuns
cmcDirs = [ os.path.join(runDir,name) for name in cmcRuns ]

dftStr = "lFeSiO125 lFeSiO126 lFeSiO127 lFeSiO128 lFeSiO129 lSiO126 lSiO127 lSiO128 lSiO129 lSiO130 lSiO131 lSiO132 lSiO133 lSiO134 lSiO135 lFeSiO130 lFeSiO131 lFeSiO132 lFeSiO133 lFeSiO134 lFeSiO135 lFeSiO136 lFeSiO137 lFeSiO138 lFeSiO139 lSiO136 lSiO137 lSiO138 lSiO139 lSiO140 lFeSiO140 lFeSiO141 lFeSiO142 lFeSiO143 lFeSiO144 lSiO141 lSiO142 lSiO143 lSiO144 lSiO145 lFeSiO145 lFeSiO146 lFeSiO147 lFeSiO148 lFeSiO149 lSiO146 lSiO147 lSiO148 lSiO149 lSiO150 lFeSiO150 lFeSiO151 lFeSiO152 lFeSiO153 lFeSiO154 lSiO151 lSiO152 lSiO153 lSiO154 lSiO155 lSiO156 lSiO157 lSiO158 lSiO159 lSiO160 lFeSiO155 lFeSiO156 lFeSiO157 lFeSiO158 lFeSiO159"
dftRuns = dftStr.split()
dftDirs = [ os.path.join(runDir,name) for name in dftRuns ]

target_pressures = [ 50., 100., 400. ]

# Get information about cmc runs

(stdout, stderr) = Popen([runinfo]+cmcRuns, stdout=PIPE).communicate()
print stdout

cmcInfo = []
for line in stdout.split('\n'):
    sline = line.split()
    try:
        assert sline[3] == 'cmc'
#        print sline
        id = sline[0]
        system = sline[1].replace('[','')
        temp = float(sline[5].replace('T=','').replace(',','') )
        volume = float(sline[6].replace('Vi=','') )
#        print [id,system,temp,volume]
        cmcInfo.append([id,system,temp,volume,host])
    except:
        pass

    
cmc = pd.DataFrame(cmcInfo,columns=['id','system','temp','volume','hostname'])
cmc['dir'] = cmcDirs
cmc.set_index('id',inplace=True,drop=False)


# Get information from about dft runs
(stdout, stderr) = Popen([runinfo]+dftRuns, stdout=PIPE).communicate()
print stdout

dftInfo = []
for line in stdout.split('\n'):
    sline = line.split()
    try:
        assert sline[3] == 'cci' or sline[3] == 'einstein' # check for correct format
        id = sline[0]
        system = sline[1].replace('[','')
        temp = float(sline[4].replace('T=','').replace(',','') )
        volume = float(sline[5].replace('V=','').replace(',','') ) 
        lam = float(sline[12].replace('lamd=','') )
        k_spring = sline[11].replace('K=','').replace(',$','')
        dftInfo.append([id,system,temp,volume,k_spring,lam,host])
    except:
        pass


dft = pd.DataFrame(dftInfo,columns=['id','system','temp','volume','k_spring','lambda','hostname'])
dft['dir'] = dftDirs

# Generate cmc results, saving to a file

timestr = time.strftime("%Y%m%d-%H%M%S")

cmc_result_file = "/u/smwahl/dat/cmc_" + timestr + ".dat"
print cmc_result_file

(stdout1, stderr) = Popen([cmc_einstein]+cmcEinsteinRuns, stdout=PIPE).communicate()
(stdout2, stderr) = Popen([cmc_result]+cmcPPRuns, stdout=PIPE).communicate()
print stdout1
print stdout2

f = open(cmc_result_file,'w')
f.write(stdout1)
f.write(stdout2)
f.close()

nrow = cmc.shape[0]

cmcResultInfo = []
for line in stdout2.split('\n')[-(nrow+1):-1]:
    sline = line.split()
    print sline
    try:
        id = sline[0]
        result_line = line
        F_an = float(sline[4])
        F_cmc = float(sline[6])
        F_class = float(sline[8])
        cmcResultInfo.append([id,result_line,F_an,F_cmc,F_class,cmc_result_file])
    except:
        pass

cmc_tmp = pd.DataFrame(cmcResultInfo,columns=['id','result','F_an','F_cmc','F_class','result_file'])
cmc_tmp.set_index('id',inplace=True)

cmc = cmc.join(cmc_tmp)

# Find each separate integration

group_vol = dft.groupby('volume')

tdi_list = []
for group in group_vol:
    df = group[1]
    runs = df.id.tolist()
    l1 = df[ df['lambda'] == 1.0 ] 
    a =l1.loc[:,['id','system','temp','volume']]
    a['dft_tdi_runs'] = [ runs]
    tdi_list.append(a)

tdi = pd.concat(tdi_list)
tdi.columns = ['dft_id', 'system','temp','volume','dft_to_cl_ids']

# runs to run eos script on
eos_runs = tdi.dft_id.tolist()

# Find thermodynamic quantities of the DFT system and infer target pressure
# given a list of possibilities
target_pressures = [ 50., 100., 400.]
(stdout, stderr) = Popen([eos]+eos_runs, stdout=PIPE).communicate()
print stdout

eosInfo = []
for line in stdout.split('\n'):
    sline = line.split()
    print sline
    try:
        assert len(sline) > 10 # ignore empty lines
        id = sline[0]
        system = sline[1]
        temp = float(sline[2])
        time = ( float(sline[6]),float(sline[7]) )
        vol = float(sline[9])
        p = (float(sline[12]), float(sline[13]) )
        pv = (float(sline[16]), float(sline[17]) )
        u = (float(sline[20]), float(sline[21]) )
        h = (float(sline[24]), float(sline[25]) )
        ptar = target_pressures[ np.abs(np.array(target_pressures) - p[0]).argmin()]
        eosInfo.append([id,system,ptar,temp,time,vol,p,pv,u,h])
        print id,system,ptar,temp,time,vol,p,pv,u,h

    except:
        pass

dft_eos = pd.DataFrame(eosInfo,columns=['id','system','P_target','T', 'run_time', \
        'V','P', 'PV', 'U','H'])

# Index both tables by the lambda=1 run id
tdi.set_index('dft_id',inplace=True,drop=False)
dft_eos.set_index('id',inplace=True,drop=False)

tdi['P_target'] = dft_eos['P_target']

tdi = tdi.sort(['system','P_target','temp'])
dft_eos = dft_eos.sort(['system','P_target','T'])



# link tdi entries to their respectice cmc runs
# This is best done by matching temperature and volume, however the output of volume from
# various scripts have different number of significant figures

     
tmp_tdi = tdi.set_index(tdi['volume'].apply(np.round,args=[5]))
tmp_cmc = cmc.set_index(cmc['volume'].apply(np.round,args=[5]))
tmp_tdi['cmc_id'] = tmp_cmc['id']

tdi = tmp_tdi.set_index('dft_id',drop=False)

# Calculate the thermodynamic integration

tdi['num_lambda'] = tdi.dft_to_cl_ids.apply(len)

#timestr = time.strftime("%Y%m%d-%H%M%S")


tmp_cmc_file = "/u/smwahl/dat/tmp_cmc.dat"
tdi_result_file = "/u/smwahl/dat/tdi_" + timestr + ".dat"
print tdi_result_file

f = open(tdi_result_file,'w')

result_list = []

for idx, row in tdi.iterrows():

    # store cmc results in a temporary file
    tmpf = open(tmp_cmc_file,'w')
    cmcrun = row['cmc_id']
    tmpf.write(cmc.loc[cmcrun]['result'])
    tmpf.close()

    # run lambda_cci
    nlam = row['num_lambda']
    runs = row['dft_to_cl_ids']
    Pkbar = row['P_target']*10. # convert pressure from GPa to kbar

    print ' '.join([lambda_cci,tmp_cmc_file,str(nlam),str(Pkbar)]+runs)

    (stdout, stderr) = Popen([lambda_cci,tmp_cmc_file,str(nlam),str(Pkbar)]+runs, stdout=PIPE).communicate()
    print stdout

    # store results
    outlines = stdout.split('\n')
    dVlines = [line.split() for line in outlines[3:3+nlam] ]
    lams = [ float(x[1]) for x in dVlines]
    dVcell = [ ( float(x[6]), float(x[7]) ) for x in dVlines ]
    dV_var = [ float(x[10]) for x in dVlines ]
    dV_var_eV = [ float(x[13]) for x in dVlines ]

    for line in outlines:
        sline = line.split()
        if len(sline) > 0:
            if sline[0] == 'Quadratic':
                 F_dft_cl = ( float(sline[8]), float(sline[10]) )     
            if sline[0] == 'F_DFT=':
                F_dft = (float(sline[1]), float(sline[2] ) ) 
                P_targetV = (float(sline[4]), float(sline[5]) ) 
                U_dft =  (float(sline[7]), float(sline[8]) ) 
                G_dft = (float(sline[10]), float(sline[11]) ) 
                break

    result_list.append([idx, lams, dVcell, dV_var, dV_var_eV, F_dft_cl, \
        F_dft,P_targetV,U_dft,G_dft,tdi_result_file ])
    print result_list[-1]

    f.write(stdout)

f.close()


print result_list
tdi_tmp = pd.DataFrame(result_list,columns=['id', 'lambdas', 'dVcell', 'dV_var', \
        'dV_var_eV', 'F_cl_dft', 'F_dft','P_targetV','U_dft','G_dft','result_file'] )
tdi_tmp.set_index('id',inplace=True)


# may wish to check before joining

tdi = tdi.join(tdi_tmp)


tabDir = os.path.join(runDir,'tables')

# load dataFrames
tdi_old = pd.load('tdi.df')
cmc_old = pd.load('cmc.df')
dft_old = pd.load('dft.df')
dft_eos_old = pd.load('dft_eos.df')

# combine new cmc and dft tables with existing ones
cmc_comb = cmc_old.append(cmc)
dft_comb = dft_old.append(dft)
tdi_comb = tdi_old.append(tdi)
dft_eos_comb = dft_eos_old.append(dft_eos)

# save DataFrames
tdi.save(tabDir+'tdi_'+timestr+'.df')
cmc.save(tabDir+'cmc'+timestr+'.df')
dft.save(tabDir+'dft'+timestr+'.df')
dft_eos.save(tabDir+'dft_eos'+timestr+'.df')

tdi_comb.save(tabDir+'tdi_all_'+timestr+'.df')
cmc_comb.save(tabDir+'cmc_all_'+timestr+'.df')
dft_comb.save(tabDir+'dft_all_'+timestr+'.df')
dft_eos_comb.save(tabDir+'dft_eos_all_'+timestr+'.df')

tdi_comb.save('tdi.df')
cmc_comb.save('cmc.df')
dft_comb.save('dft.df')
dft_eos_comb.save('dft_eos.df')

# load dataFrames
tdi = pd.load('tdi.df')
cmc = pd.load('cmc.df')
dft = pd.load('dft.df')
dft_eos = pd.load('dft_eos.df')



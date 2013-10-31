# full_tdi.py
'''Script t  consolidate and analyze classical MonteCarlo and DFT thermodynamic integration 
using Pandas dataframes'''

from parse_vasp import *
import os
import subprocess
import pandas as pd
import numpy as np
import socket

host = socket.gethostname()

runDir="/u/smwahl/scr/workspace"
saveDir="/u/smwahl/dat/"

# analysis scripts
runinfo="/u/smwahl/scripts/vasp_automation/runinfo"
eos="/u/smwahl/scripts/eos1"
cmc_result="/u/smwahl/code/python/vasp_automation/cmc_result"
lambda_cci="/u/smwahl/scripts/vasp_automation/lambda_cci"

cmcStr = "lFe299 lFe300 lFe301 lFe302 lFe303 lFe304 lFe305 lFe306 lFe307 lFeMgO398 lFeMgO399 lFeMgO400 lFeMgO401 lFeMgO402 lFeMgO403 lFeMgO404 lFeMgO405 lFeMgO406 lMgO267 lMgO268 lMgO269 lMgO270 lMgO271 lMgO272 lMgO273 lMgO274"
cmcRuns = cmcStr.split(' ')
cmcDirs = [ os.path.join(runDir,name) for name in cmcRuns ]

dftStr = "lFe254 lFe255 lFe256 lFe257 lFe258 lFe259 lFe260 lFe261 lFe262 lFe263 lFe264 lFe265 lFe266 lFe267 lFe268 lFe269 lFe270 lFe271 lFe272 lFe273 lFe274 lFe275 lFe276 lFe277 lFe278 lFe279 lFe280 lFe281 lFe282 lFe283 lFe284 lFe285 lFe286 lFe287 lFe288 lFe289 lFe290 lFe291 lFe292 lFe293 lFe294 lFe295 lFe296 lFe297 lFe298 lFeMgO353 lFeMgO354 lFeMgO355 lFeMgO356 lFeMgO357 lFeMgO358 lFeMgO359 lFeMgO360 lFeMgO361 lFeMgO362 lFeMgO363 lFeMgO364 lFeMgO365 lFeMgO366 lFeMgO367 lFeMgO368 lFeMgO369 lFeMgO370 lFeMgO371 lFeMgO372 lFeMgO373 lFeMgO374 lFeMgO375 lFeMgO376 lFeMgO377 lFeMgO378 lFeMgO379 lFeMgO380 lFeMgO381 lFeMgO382 lFeMgO383 lFeMgO384 lFeMgO385 lFeMgO386 lFeMgO387 lFeMgO388 lFeMgO389 lFeMgO390 lFeMgO391 lFeMgO392 lFeMgO393 lFeMgO394 lFeMgO395 lFeMgO396 lFeMgO397 lMgO227 lMgO228 lMgO229 lMgO230 lMgO231 lMgO232 lMgO233 lMgO234 lMgO235 lMgO236 lMgO237 lMgO238 lMgO239 lMgO240 lMgO241 lMgO242 lMgO243 lMgO244 lMgO245 lMgO246 lMgO247 lMgO248 lMgO249 lMgO250 lMgO251 lMgO257 lMgO258 lMgO259 lMgO260 lMgO261 lMgO262 lMgO263 lMgO264 lMgO265 lMgO266"

dftRuns = dftStr.split(' ')
dftDirs = [ os.path.join(runDir,name) for name in dftRuns ]

target_pressures = [ 50., 100., 400. ]


from subprocess import Popen, PIPE

(stdout, stderr) = Popen([runinfo]+cmcRuns, stdout=PIPE).communicate()
print stdout

cmcInfo = []
for line in stdout.split('\n'):
    sline = line.split(' ')
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


# Get information from runinfo
(stdout, stderr) = Popen([runinfo]+dftRuns, stdout=PIPE).communicate()
print stdout

dftInfo = []
for line in stdout.split('\n'):
    sline = line.split(' ')
    try:
        assert sline[3] == 'cci,' # check for correct format
        id = sline[0]
        system = sline[1].replace('[','')
        temp = float(sline[4].replace('T=','').replace(',','') )
        volume = float(sline[5].replace('V=','').replace(',','') ) 
        lam = float(sline[12].replace('lamd=','') )
        dftInfo.append([id,system,temp,volume,lam,host])
    except:
        pass


dft = pd.DataFrame(dftInfo,columns=['id','system','temp','volume','lambda','hostname'])
dft['dir'] = dftDirs

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

# Generate cmc results, saving to a file

import time
timestr = time.strftime("%Y%m%d-%H%M%S")

cmc_result_file = "/u/smwahl/dat/cmc_" + timestr + ".dat"
print cmc_result_file

(stdout, stderr) = Popen([cmc_result]+cmcRuns, stdout=PIPE).communicate()
print stdout

f = open(cmc_result_file,'w')
f.write(stdout)
f.close()

nrow = cmc.shape[0]

cmcResultInfo = []
for line in stdout.split('\n')[-(nrow+1):-1]:
    sline = line.split()
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
cmc = cmc_tmp.set_index('id',inplace=True)

# link tdi entries to their respectice cmc runs
# This is best done by matching temperature and volume, however the output of volume from
# various scripts have different number of significant figures

     
tmp_tdi = tdi.set_index(tdi['volume'].apply(np.round,args=[5]))
tmp_cmc = cmc.set_index(cmc['volume'].apply(np.round,args=[5]))
tmp_tdi['cmc_id'] = tmp_cmc['id']

tdi = tmp_tdi.set_index('dft_id')

# Calculate the thermodynamic integration


# save database file

#pd.write_frame

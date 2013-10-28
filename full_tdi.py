# full_tdi.py
'''Script t  consolidate and analyze classical MonteCarlo and DFT thermodynamic integration'''

from parse_vasp import *
import os
import subprocess
import pandas as pd
import numpy as np

runDir="/u/smwahl/scr/workspace"
saveDir="/u/smwahl/dat/"

# scripts
runinfo="/u/smwahl/scripts/vasp_automation/runinfo"
eos="/u/smwahl/scripts/eos1"

cmcStr = "lFe299 lFe300 lFe301 lFe302 lFe303 lFe304 lFe305 lFe306 lFe307 lFeMgO398 lFeMgO399 lFeMgO400 lFeMgO401 lFeMgO402 lFeMgO403 lFeMgO404 lFeMgO405 lFeMgO406 lMgO267 lMgO268 lMgO269 lMgO270 lMgO271 lMgO272 lMgO273 lMgO274"
cmcRuns = cmcStr.split(' ')
cmcDirs = [ os.path.join(runDir,name) for name in cmcRuns ]

dftStr = "lFe254 lFe255 lFe256 lFe257 lFe258 lFe259 lFe260 lFe261 lFe262 lFe263 lFe264 lFe265 lFe266 lFe267 lFe268 lFe269 lFe270 lFe271 lFe272 lFe273 lFe274 lFe275 lFe276 lFe277 lFe278 lFe279 lFe280 lFe281 lFe282 lFe283 lFe284 lFe285 lFe286 lFe287 lFe288 lFe289 lFe290 lFe291 lFe292 lFe293 lFe294 lFe295 lFe296 lFe297 lFe298 lFeMgO353 lFeMgO354 lFeMgO355 lFeMgO356 lFeMgO357 lFeMgO358 lFeMgO359 lFeMgO360 lFeMgO361 lFeMgO362 lFeMgO363 lFeMgO364 lFeMgO365 lFeMgO366 lFeMgO367 lFeMgO368 lFeMgO369 lFeMgO370 lFeMgO371 lFeMgO372 lFeMgO373 lFeMgO374 lFeMgO375 lFeMgO376 lFeMgO377 lFeMgO378 lFeMgO379 lFeMgO380 lFeMgO381 lFeMgO382 lFeMgO383 lFeMgO384 lFeMgO385 lFeMgO386 lFeMgO387 lFeMgO388 lFeMgO389 lFeMgO390 lFeMgO391 lFeMgO392 lFeMgO393 lFeMgO394 lFeMgO395 lFeMgO396 lFeMgO397 lMgO227 lMgO228 lMgO229 lMgO230 lMgO231 lMgO232 lMgO233 lMgO234 lMgO235 lMgO236 lMgO237 lMgO238 lMgO239 lMgO240 lMgO241 lMgO242 lMgO243 lMgO244 lMgO245 lMgO246 lMgO247 lMgO248 lMgO249 lMgO250 lMgO251 lMgO257 lMgO258 lMgO259 lMgO260 lMgO261 lMgO262 lMgO263 lMgO264 lMgO265 lMgO266"

dftRuns = dftStr.split(' ')
dftDirs = [ os.path.join(runDir,name) for name in dftRuns ]

target_pressures = [ 50., 100., 400. ]

supbrocess

from subprocess import Popen, PIPE

(stdout, stderr) = Popen([runinfo]+cmcRuns, stdout=PIPE).communicate()
print stdout

cmcInfo = []
for line in stdout.split('\n')[:-1]:
    sline = line.split(' ')
    try:
        id = sline[0]
        system = sline[1].replace('[','')
        temp = float(sline[5].replace('T=','').replace(',','') )
        volume = float(sline[6].replace('Vi=','') )
        cmcInfo.append([id,system,temp,volume])
    except:
        pass

cmc['dir'] = cmcDirs
cmc = pd.DataFram(cmcInfo)



(stdout, stderr) = Popen([runinfo]+dftRuns, stdout=PIPE).communicate()
print stdout
dftInfo = []
for line in stdout.split('\n')[:-1]:
    sline = line.split(' ')
    try:
        id = sline[0]
        system = sline[1].replace('[','')
        temp = float(sline[5].replace('T=','').replace(',','') )
        volume = float(sline[6].replace('Vi=','') )
        dftInfo.append([id,system,temp,volume])
    except:
        pass

dft['dir'] = dftDirs
dft = pd.DataFram(dftInfo)

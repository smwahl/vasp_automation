# full_tdi.py
'''Script consolidate and analyze classical MonteCarlo and DFT thermodynamic integration 
using Pandas dataframes'''

#cd ~/code/python/vasp_automation

# import modules
from parse_vasp import *
import os
import subprocess
import pandas as pd
import numpy as np
import socket
from subprocess import Popen, PIPE
import time
import glob


# helper functions
def countInFile(query,paths):
    '''Count the number of occurences of a string in a file or list of files'''

    pathl = []
    if isinstance(paths,str):
        pathl.append(paths)
    elif isinstance(paths,list):
        pathl = paths
    else:
        raise Exception('Not a string or list of strings')

    count = 0
    for path in pathl:
        with open(path) as f:
            for line in f:
                if query in line:
                    count +=1
    return count

def updateData(old, new, idxcol=None):
    '''Defines the default way in which a single master dataframe is updated with new
    data.

    This appends new data as rows and drops any previous rows with the same index.

    By default the comparisons are done by index, if 'idxcol' is specified, then 
    duplicates are checked for in that column rather than in the indices 

    The motivation is having identical MD runs which have been run longer
    than the previously have. We therefore plan to save updated dataframes with the
    full runs. For comparing run lengths we would then compare between the updated 
    dataframe and a saved version of the the old dataframe with the shorter run
    lengths.'''

    tmp = pd.concat([old,new])
    if idxcol is None:
        tmp['tmp_idx'] = tmp.index.values
        tmp.drop_duplicates(cols='tmp_idx',take_last=True,inplace=True)
        comb = tmp.drop('tmp_idx',axis=1)
    else:
        comb = tmp.drop_duplicates(cols=idxcol,take_last=True)

    return comb

# Functions for generating dataframes
def makeCMC(cmcPPruns=[],cmcEinsteinRuns=[]):
    '''Generate a pandas dataframe storing information about cmc runs, given lists of
    ids for liquid runs using pairpotentials or solid runs using Einstein potentials.'''

    global host, addDate, kpoints, functional, runDir, runinfo, cmc_result, \
            cmc_einstein, timestr
    
    # combine list of runs
    cmcRuns = cmcEinsteinRuns + cmcPPRuns
    cmcDirs = [ os.path.join(runDir,name) for name in cmcRuns ]

    # Get information about cmc runs
    (stdout, stderr) = Popen([runinfo]+cmcRuns, stdout=PIPE).communicate()
    print stdout

    cmcInfo = []
    for line in stdout.split('\n'):
        sline = line.split()
        print sline
        try:
            assert sline[3] == 'cmc'
    #        print sline
            id = sline[0]
            system = sline[1].replace('[','')
            temp = float(sline[5].replace('T=','').replace(',','') )
            volume = float(sline[6].replace('Vi=','') )
            print [id,system,temp,volume]
            path = os.path.join(runDir,id)
            try:
                date =  os.path.getmtime(os.path.join(runDir,id,'run.scr'))
            except:
                date = None
            cmcInfo.append([id,system,temp,volume,host,path,date,addDate,kpoints, \
                functional])
        except:
            pass

    print cmcInfo

    cmc = pd.DataFrame(cmcInfo,columns=['id','system','temp','volume','hostname', \
            'path','date','add_date','kpoints','functional'])
    cmc['dir'] = cmcDirs
    cmc.set_index('id',inplace=True,drop=False)


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

    # add the cmc results to the dataframe
    cmc_tmp = pd.DataFrame(cmcResultInfo,columns=['id','result','F_an','F_cmc',\
            'F_class','result_file'])
    cmc_tmp.set_index('id',inplace=True)

    cmc = cmc.join(cmc_tmp)
    return cmc 


def makeDFT(dftRuns):
    '''Generate a pandas dataframe storing data about a collection of dft runs, given 
        a list of ids for those runs.'''

    global host, addDate, kpoints, functional, runDir, saveDir, runinfo, \
            lambda_cci, timestr

    dftDirs = [ os.path.join(runDir,name) for name in dftRuns ]
    # Get information from about dft runs

    (stdout, stderr) = Popen([runinfo]+dftRuns, stdout=PIPE).communicate()
    print stdout

    dftInfo = []
    for line in stdout.split('\n'):
        sline = line.split()
        try:
            assert sline[3] == 'cci' or sline[3] == 'einstein' # check for correct format
            id = sline[0]
            print id
            system = sline[1].replace('[','')
            temp = float(sline[4].replace('T=','').replace(',','') )
            volume = float(sline[5].replace('V=','').replace(',','') ) 
            tstep = float(sline[7].replace('tstep=','').replace(',','') )
            cutoff = float(sline[8].replace('cutoff=','').replace(',','') )
            lam = float(sline[12].replace('lamd=','') )
            k_spring = sline[11].replace('K=','').replace(',$','')
            path = os.path.join(runDir,id)
            outcars = glob.glob(path + '/out*/OUTCAR') + [path+'/OUTCAR'] 
            try:
                nstep = countInFile('LOOP+',outcars)
                date =  os.path.getmtime(outcars[-1])
                runtime = nstep*tstep
            except:
                date = None

            dftInfo.append([id,system,temp,volume,k_spring,lam,host,tstep,cutoff,\
                    nstep,runtime,date,path,addDate,kpoints,functional])
        except:
            pass

    dft = pd.DataFrame(dftInfo,columns=['id','system','temp','volume','k_spring',\
            'lambda','hostname','tstep','cutoff','nstep','runtime','date','path',\
            'add_date','kpoints','functional'])
    dft['dir'] = dftDirs
    return dft


def makeTDI(dft):
    '''Generate a dataframe for each thermodynamic integration based on finding the dft
     tdi runs with matching volume. Returns two dataframes, one for the tdi and a 
     one with the corresponding eos results for the lambda = 1 case'''

    global host, addDate, kpoints, functional, runDir, saveDir, runinfo, lambda_cci

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
            eosInfo.append([id,system,ptar,temp,time,vol,p,pv,u,h,kpoints,functional])
            print id,system,ptar,temp,time,vol,p,pv,u,h

        except:
            pass

    dft_eos = pd.DataFrame(eosInfo,columns=['id','system','P_target','T', 'run_time', \
            'V','P', 'PV', 'U','H','kpoints','functional'])

    # Index both tables by the lambda=1 run id
    tdi.set_index('dft_id',inplace=True,drop=False)
    dft_eos.set_index('id',inplace=True,drop=False)

    tdi['P_target'] = dft_eos['P_target']

    tdi = tdi.sort(['system','P_target','temp'])
    dft_eos = dft_eos.sort(['system','P_target','T'])

    return tdi, dft_eos

def linkCMC(tdi,cmc):
    '''Link tdi entries to their respectice cmc runs. This is done by matching 
    temperature and volume, however the output of volume from
    various scripts have different number of significant figures'''

    global host, addDate, kpoints, functional, runDir, saveDir, eos, timestr
     
    tmp_tdi = tdi.set_index(tdi['volume'].apply(np.round,args=[5]))
    tmp_cmc = cmc.set_index(cmc['volume'].apply(np.round,args=[5]))
    tmp_tdi['cmc_id'] = tmp_cmc['id']

    tdi = tmp_tdi.set_index('dft_id',drop=False)

    # Calculate the thermodynamic integration

    tdi['num_lambda'] = tdi.dft_to_cl_ids.apply(len)

    timestr = time.strftime("%Y%m%d-%H%M%S")

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

        (stdout, stderr) = Popen([lambda_cci,tmp_cmc_file,str(nlam),str(Pkbar)] \
                +runs,stdout=PIPE).communicate()
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
            F_dft,P_targetV,U_dft,G_dft,tdi_result_file,kpoints,functional ])
        print result_list[-1]

        f.write(stdout)

    f.close()


    print result_list
    tdi_tmp = pd.DataFrame(result_list,columns=['id', 'lambdas', 'dVcell', 'dV_var', \
            'dV_var_eV', 'F_cl_dft', 'F_dft','P_targetV','U_dft','G_dft', \
            'result_file','kpoints','functional'] )
    tdi_tmp.set_index('id',inplace=True)

    # may wish to check before joining
    return tdi.join(tdi_tmp)


if __name__ == "__main__":

    # General informations ( these are treated as global variables in functions
    # adding data to the dataframes
    # Need to add this data separately
    host = socket.gethostname()
    addDate = time.time()
    kpoints = 'balderesci' 
    functional = 'PBE' # at the moment this is to distingush regular GGA 'PBE' and GGA+U 'PBE+U')



    # Generate cmc results, saving to a file
    timestr = time.strftime("%Y%m%d-%H%M%S")

    # directories for reading raw data and saving consolidated data
    runDir="/u/smwahl/scr/workspace"
    saveDir="/u/smwahl/dat/"
    tabDir = os.path.join(saveDir,'tables')

    # analysis scripts (passed as globals)
    runinfo="/u/smwahl/scripts/vasp_automation/runinfo"
    eos="/u/smwahl/scripts/eos1"
    cmc_result="/u/smwahl/code/python/vasp_automation/cmc_result"
    cmc_einstein="/u/smwahl/code/python/vasp_automation/cmc_einstein"
    lambda_cci="/u/smwahl/scripts/vasp_automation/lambda_cci"

    # global host, addDate, kpoints, functional, runDir, saveDir, runinfo, eos, cmc_result, cmc_einstein, lambda_cci

    # Run directories (note separate arrays for liquid and solid cmc runs)
#    cmcStr = "lFeMgO509 lFeMgO510 lFeMgO511 lFeMgO512"
    cmcStr="lFe384 lMgO315 lFeMgO544 lFe385 lMgO316 lFeMgO545 lFe386 lMgO317 lFeMgO546"
    cmcPPRuns = cmcStr.split()
    cmcEinsteinStr = ""
    cmcEinsteinRuns = cmcEinsteinStr.split()
    #cmcRuns = cmcEinsteinRuns + cmcPPRuns
    #cmcDirs = [ os.path.join(runDir,name) for name in cmcRuns ]

#    dftStr = "lFeMgO489 lFeMgO490 lFeMgO491 lFeMgO492 lFeMgO493 lFeMgO494 lFeMgO495 lFeMgO496 lFeMgO497 lFeMgO498 lFeMgO499 lFeMgO500 lFeMgO501 lFeMgO502 lFeMgO503 lFeMgO504 lFeMgO505 lFeMgO506 lFeMgO507 lFeMgO508"
    dftStr = "lFe363 lFe364 lFe365 lFe366 lFe367 lMgO294 lMgO295 lMgO296 lMgO297 lMgO298 lFeMgO523 lFeMgO524 lFeMgO525 lFeMgO526 lFeMgO527 lFe368 lFe369 lFe370 lFe371 lFe372 lMgO299 lMgO300 lMgO301 lMgO302 lMgO303 lFeMgO528 lFeMgO529 lFeMgO530 lFeMgO531 lFeMgO532 lFe373 lFe374 lFe375 lFe376 lFe377 lMgO304 lMgO305 lMgO306 lMgO307 lMgO308 lFeMgO533 lFeMgO534 lFeMgO535 lFeMgO536 lFeMgO537"
    dftRuns = dftStr.split()
    dftDirs = [ os.path.join(runDir,name) for name in dftRuns ]

    #target_pressures = [ 50., 100., 400. ]

    # load old dataFrames
#    tab_num = '20140110-142950'
#    tab_num = '20140219-164724'
#    tab_num = '20140219-182020'
#    tab_num = '20140219-182020'
    tab_num = '20140327-144321'
    tdi_old = pd.load(tabDir + '/' + 'tdi_all_' + tab_num + '.df')
    cmc_old = pd.load(tabDir + '/' + 'cmc_all_' + tab_num + '.df')
    dft_old = pd.load(tabDir + '/' + 'dft_all_' + tab_num + '.df')
    dft_eos_old = pd.load(tabDir + '/' + 'dft_eos_all_' + tab_num + '.df')

    # generate the dataframes
    dft = makeDFT(dftRuns)
    if len(cmcPPRuns + cmcEinsteinRuns) > 0:
        cmc = makeCMC(cmcPPRuns,cmcEinsteinRuns)
    tdi, dft_eos = makeTDI(dft)

    # if necessary link dft runs to old cmc table
    # (e.g.) if only updating values for existing DFT runs
    if len(cmcPPRuns + cmcEinsteinRuns) == 0:
        cmc = cmc_old

    tdi = linkCMC(tdi,cmc)

    #link corresponding cmc runs to tdi (this could probably instead be handled with a merge)

    # combine new cmc and dft tables with existing ones
    if len(cmcPPRuns + cmcEinsteinRuns) > 0:
        cmc_comb = updateData(cmc_old,cmc)
    else:
        cmc_comb = cmc
    dft_comb = updateData(dft_old,dft,idxcol='id')
    tdi_comb = updateData(tdi_old,tdi)
    dft_eos_comb = updateData(dft_eos_old,dft_eos)

    # save DataFrames
    tdi.save(tabDir+'/tdi_'+timestr+'.df')
    cmc.save(tabDir+'/cmc_'+timestr+'.df')
    dft.save(tabDir+'/dft_'+timestr+'.df')
    dft_eos.save(tabDir+'/dft_eos_'+timestr+'.df')

    tdi_comb.save(tabDir+'/tdi_all_'+timestr+'.df')
    cmc_comb.save(tabDir+'/cmc_all_'+timestr+'.df')
    dft_comb.save(tabDir+'/dft_all_'+timestr+'.df')
    dft_eos_comb.save(tabDir+'/dft_eos_all_'+timestr+'.df')

    #tdi_comb.save('tdi.df')
    #cmc_comb.save('cmc.df')
    #dft_comb.save('dft.df')
    #dft_eos_comb.save('dft_eos.df')

    # load dataFrames

    # new data only
    #tdi = pd.load(tabDir+'/tdi_'+timestr+'.df')
    #cmc = pd.load(tabDir+'/cmc_'+timestr+'.df')
    #dft = pd.load(tabDir+'/dft_'+timestr+'.df')
    #dft_eos = pd.load(tabDir+'/dft_eos_'+timestr+'.df')

    #combined
    #tdi = pd.load('tdi.df')
    #cmc = pd.load('cmc.df')
    #dft = pd.load('dft.df')
    #dft_eos = pd.load('dft_eos.df')

    # print data location
    print 'Data directory: ' + saveDir
    print 'Table directory: '+ tabDir
    print 'Identifying time string: ' + timestr



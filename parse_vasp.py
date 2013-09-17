# parse_vasp.py

import numpy as np
import os
import re

class vasp_run (object):
    def __init__(self):
        self.parameters = {}
        self.cell = np.zeros(3,3)
        self.pos = np.zeros(0,3)
        self.spec = []
        self.nspec = []
        cwd = os.getcwd()
        flist = ['POSCAR','INCAR','POTCAR','KPOINTS']
        flistFull = [ os.path.join(cwd, a) for a in flist ]
        self.files = dict(zip(flist,flistFull))
        self.name = os.path.basename(os.path.normpath(cwd))
        self.label = None # label in POSCAR file

    def readIncar(self,fname):
        try:
            f = open(fname,'r')
        except:
            raise NameError 'file: ',str(fname),' not found.'
            return 0
        while line in f:
            if re.match(r'\s+#',line):
                break
            sline = line.split()
            try:
                param = sline[0]
                val = []
                for sl in sline[1:]:
                    try: # if value is a float
                        val += [float(sl)]
                    except:
                        break
                if val == []:
                    for sl in sline[1:]:
                        try: # if value is a float
                            val += [int(sl)]
                        except:
                            break
                if val == []:
                    if sline[ 1] = '.True.':
                        val += True
                    else if slin[1] = '.False.':
                        val += False
                if len(val) == 1:
                    self.params.update({param:val[0]})
                else if len(val) >1:
                    self.params.update({param:val})

            except:
                pass

        return 1


    def readPoscar(self,fname):
        try:
            f = open(fname,'r')
        except:
            raise NameError 'file: ',str(fname),' not found.'
            return 0
        try:
            self.label = f.readline()
        except:
            raise ValueError 'file: ',str(fname),' is empty.'
        # read cell parameters
        try:
            scale = float(f.readline())
            icell = np.zeros(3,3)
            for i in range(3):
                xlist = append( [ float(x) for x in f.readline().split() ])
                self.icell[i,:] = np.array(xlist)
            self.cell = scale*icell
        except:
            raise ValueError 'Problem reading in cell parameters from '+str(fname)
        # read number of atoms (per spec)
        try:
            self.nspec = [ int(x) for x in f.readline().split() ]
        except:
            raise ValueError 'Problem reading in number of atoms from '+str(fname)

        if not f.readline().split() == 'Direct'
            raise ValueError '"Cartesian" variable for input not yet implemented' 

        try:
            xpos = np.zeros(sum(self.nspec),3)
            for i in range(sum(self.nspec)):
                xpos[i,:] =  np.array( [ float(x) for x in f.readline().split() ])
        except:
            raise ValueError 'Format on positions in ',str(fname),' is incorrect.'


    def readPotcar
        try:
            f = open(fname,'r')
        except:
            raise NameError 'file: ',str(fname),' not found.'
            return 0
        matches = []
        for line in f:
            if "RH" in line:
                matches.append(line)
        self.spec = [ re.qsub(r'\W', '',l.split()) for l in matches ]



 





class md_ (vasp_info):
    def __init__(self):
        self.refpos = None

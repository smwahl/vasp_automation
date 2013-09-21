# super.py
''' Function to generate a supercell based on a VASP  position file (POSCAR,CONTCAR,
etc...). '''

import parse_vasp as vasp
import numpy as np
from numpy import linalg as la
import os

def supercell(superDim,directory=os.getcwd(),filename=None,outname=None):
    ''' Generate a supercell from a VASP position file, an ndarray containing integer
multiples of the 3 cell vectors to define the basis of the supercell. Prints results to a file'''

    # create a vasp_run object and read in CONTCAR
    run = vasp.vasp_run(directory)
    run.readIncar()
    run.readPotcar()

    if filename is None:
        fname = run.files['CONTCAR']
    else:
        fname = filename
    contPos = run.readPosFile(fname)

    # store values relavent to
    nspec = contPos['nspec']
    unitCell = contPos['cell']
    unitVol = np.dot(unitCell[0,:], np.cross(unitCell[1,:],unitCell[2,:]) )

    basis = contPos['pos']
    cartBasis = np.dot(basis,unitCell)
    
    superDimInv = la.inv(superDim)

    # generate parameters of the supercell described by superDim
    superCell = np.dot(superDim,unitCell)
    superLens = [ la.norm(superCell[i,:]) for i in range(3) ]
    superAngs = [ np.degrees(np.arccos(np.dot(superCell[i,:],superCell[j,:])
        /superLens[i]/superLens[j] )) for i,j in [ (1,2),(0,2),(0,1) ] ]
    superVol = np.dot(superCell[0,:], np.cross(superCell[1,:],superCell[2,:]) )

    print 'Supercell vectors:'
    print superCell
    print 'Vector lengths:'
    print 'a=',superLens[0],'b=',superLens[2],'c=',superLens[2]
    print 'Angles:'
    print 'alpha=',superAngs[0],'beta=',superAngs[2],'gamma=',superAngs[2]

    # Generate possible cartesian limits

    corners = np.array( [ [0, 0, 0 ],[ 0, 0, 1],[ 0, 1, 0],[ 1, 0, 0],[ 1, 1, 0]
        ,[ 1, 0, 1],[ 0, 1, 1],[1, 1, 1] ] )
    supCorners = np.dot(corners,superDim)
    max_ind = supCorners.max(0)
    min_ind = supCorners.min(0)


    positions = np.zeros([0,3])
    nspec_sup = []

    # split basis into portions
    basisBySpec = [ basis[sum(nspec[:i]):sum(nspec[:i+1]), : ]
            for i in range(len(nspec)) ]

    # iterate over different atom types
    for sbasis,Nbasis in zip(basisBySpec,nspec):

        # generate all possible positions within minimum and maximum indices
        # in fractional coordinates of primitive cell
        # This method should not lead to duplicate points !
        gen_pos = np.zeros([0,3])
        for i in range(min_ind[0],max_ind[0]):
            for j in range(min_ind[1],max_ind[1]):
                for k in range(min_ind[2],max_ind[2]):
                   gen_pos = np.vstack((gen_pos, sbasis + 
                       np.outer(np.ones(Nbasis), np.array([i,j,k]) ) ) )

        #translate positions into the supercell basis
        gen_pos_sup = np.dot(superDimInv.T, gen_pos.T).T
        Ngen,dim = gen_pos_sup.shape

        # eliminate positions that are outside of the new supercell
        count = 0
        for i in range(Ngen):
            pos = gen_pos_sup[i,:]
            # Check whether all indices are between 0 and 1 in supercell coordinates
            if pos.min() >= 0.0 and pos.max() <=1.0:
                positions = np.vstack((positions,pos))
                count += 1

        # store the number of atoms for species in the supercell
        nspec_sup.append(count) 

    print 'Number of Atoms:'
    print nspec_sup,sum(nspec_sup) 
#     print basis
#     print positions

    if outname is None:
        outfile = os.path.join(run.wd,''.join(['REFCAR',str(sum(nspec_sup))]) )
    else:
        outfile = outname
    run.writePosFile(outfile,superCell,nspec_sup,positions,label=
            'Supercell generated using super.py.' )
    
    # matlab code for calculating nearest mirror
    # shifts = [ 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1;
    #             -1 1 0; -1 0 1; 0 -1 1; -1 1 1;  1 -1 0; 1 0 -1; 0 1 -1; 1 -1 1;
    #              1 1 -1;  -1 -1 1;  -1 1 -1;  1 -1 -1];
    # shifted = shifts*vectors;
    # shift_len = sqrt( shifted(:,1).^2 + shifted(:,2).^2 + shifted(:,3).^2);
    # 
    # min_mirror = min(shift_len);

    # % vectors between opposing corners
    # oppos = [ shifted(7,:) ; shifted(6,:) - shifted(1,:);
    #           shifted(5,:)-shifted(2,:); shifted(4,:)-shifted(3,:)];
    # oppos_len = sqrt( oppos(:,1).^2 + oppos(:,2).^2 + oppos(:,3).^2);
    # max_oppos = max(oppos_len);


# Test runs functions from the commandline
if __name__ == '__main__':
    superDim = np.array([ [2,0,0],[0,2,0],[0,0,1] ] )
    supercell(superDim,'test_run')
    superDim = np.array([ [2,0,0],[0,1,-1],[0,2,1] ] )
    supercell(superDim,'test_run')

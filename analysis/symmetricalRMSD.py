"""
Calculate an RMSD that takes into account the symmetry of a protein.

Right now, this code is specific to rubisco, which has D4 symmetry.
It may generalize to other D4 proteins but it also may not, because we
reference specific subunits; it's also reliant on the interactor being
chain A, and then the rest being the rubisco large and then small subunits.
I (Serena Lotreck) have tried to refactor the code in such a way as to
make it as extensible as possilbe, but it's very much hard-coded specificity
at the moment!

Author: Josh Vermaas
"""
from vmd import Molecule, atomsel, molecule
import numpy as np
from scipy.spatial.transform import Rotation as R
from os.path import abspath


def get4x4rotation(i):
    """
    There are 7 new rotations we want to deal with beond the original
    orientation: the three rotations around "Z" and the two diagonal
    reflections.

    From ChatGPT:
        Here’s the full list of proper rotation operations for D₄
        (with rotation axis and angle):
        Operation Axis Angle
        R₀ [0, 0, 1] 0°
        R₁ [0, 0, 1] 90°
        R₂ [0, 0, 1] 180°
        R₃ [0, 0, 1] 270°
        R₄ [1, 0, 0] 180°
        R₅ [0, 1, 0] 180°
        R₆ [1, 1, 0] 180°
        R₇ [1, –1, 0] 180°

    parameters:
        i, int: which rotation we want

    returns:
        rr, np.array: the 4x4 rotation matrix 
    """
    # Define the rotation vector
    if i == 0:
        r = R.from_rotvec(90*np.array([0, 0, 1]), degrees=True)
    elif i == 1:
        r = R.from_rotvec(180*np.array([0, 0, 1]), degrees=True)
    elif i == 2:
        r = R.from_rotvec(270*np.array([0, 0, 1]), degrees=True)
    elif i == 3:
        r = R.from_rotvec(180*np.array([1,0,0]), degrees=True)
    elif i == 4:
        r = R.from_rotvec(180*np.array([0,1,0]), degrees=True)
    # We have to normalize the vector for the diagonal reflections so
    # we only rotate by 180 instead of more than that
    elif i == 5:
        r = R.from_rotvec(180*(1 / np.sqrt(2))*np.array([1,1,0]), degrees=True)
    elif i == 6:
        r = R.from_rotvec(180*(1 / np.sqrt(2))*np.array([1,-1,0]), degrees=True)

    # Make it into a 4x4 matrix
    rr = np.zeros((4,4), dtype=float)
    rr[:3,:3] = r.as_matrix()
    rr[3,3] = 1
    
    return rr


def alignProtein(mol, protein='rubisco'):
    """
    Align the protein to the X and Y axes.

    parameters:
        mol, Molecule: molecule object defining the protein
        protein, str, optional: name of the protein to align, default is rubisco

    returns: movesel: selection to move
    """
    # Align to X
    # Split the protein in half around the center to have something to do
    # alignment math on
    if protein == 'rubisco':
        AHalfProtein = atomsel("type CA and chain B C D E", frame=0)
        BHalfProtein = atomsel("type CA and chain F G H I", frame=0)
    else:
        pass # Would need implementation for other proteins
    vectortoalign = np.array(AHalfProtein.center()) - np.array(BHalfProtein.center()) # Align this vector to an axis.
    rotmatrix = R.align_vectors(np.array([[1,0,0]]), [vectortoalign])
    movesel = atomsel("all")
    movesel.frame=0
    
    # VMD wants a 4x4 matrix, not a 3x3 one.
    rr = np.zeros((4,4), dtype=float)
    rr[:3,:3] = rotmatrix[0].as_matrix()
    rr[3,3] = 1
    
    # Align the protein in all frames
    movesel.move(rr.T.flatten())

    # Now repeat for Y
    # Split in half in the other orientation
    if protein == 'rubisco':
        AHalfProtein = atomsel("type CA and chain B C H I", frame=0)
        BHalfProtein = atomsel("type CA and chain F G D E", frame=0)
    else:
        pass # Would need implementation for other proteins
    vectortoalign = np.array(AHalfProtein.center()) - np.array(BHalfProtein.center()) # Align this vector to an axis.
    rotmatrix = R.align_vectors(np.array([[0,1,0]]), [vectortoalign])
    movesel = atomsel("all")
    movesel.frame=0
    
    # VMD wants a 4x4 matrix, not a 3x3 one.
    rr = np.zeros((4,4), dtype=float)
    rr[:3,:3] = rotmatrix[0].as_matrix()
    rr[3,3] = 1

    # Align the protein in all frames
    movesel.move(rr.T.flatten())

    return movesel


def calculateSymmetricalRMSD(molDir, protein='rubisco'):
    """
    Perform alignments, rotations and calculate RMSD.

    parameters:
        molDir, str: directory path containing the 5 fold poses for the protein
        protein, str, optional: name of the protein to align, default is rubisco

    returns:
        rmsdmatrix, np.array: 
    """
    # Load molecule
    mol = Molecule.Molecule()
    for i in range(5):
        mol.load(f'{abspath(molDir)}/seed-1855_sample-{i}/model.cif', filetype='pdbx')

    # Align the proteins
    movesel = alignProtein(mol, protein='rubisco')

    # Center rubisco at the origin
    proteinsel = atomsel("type CA and not chain A") # "type CA" Specific to PyVMD
    proteinref = atomsel("type CA and not chain A", frame=0)
    for i in range(mol.numFrames()):
        proteinsel.frame = i
        movesel.frame = i
        # Line up rubisco
        movesel.move(proteinsel.fit(proteinref))
        # Center rubisco at the origin
        movesel.moveby(-1*np.array(proteinsel.center()))
    originalframes = mol.numFrames()
    proteinmove = atomsel("not chain A")
    for j in range(7): # Loop over different rotations, there are 8 total, 7 new ones.
        M = get4x4rotation(j)
        for i in range(originalframes): # Loop over duplicated frames
            molecule.dupframe(int(mol), i)
            proteinmove.frame = originalframes + i + 5 * j
            proteinmove.move(M.flatten())
    isel = atomsel("all")
    jsel = atomsel("all")
    rmsdmatrix = np.zeros((5,5), dtype=float)
    for i in range(5):
        isel.frame = i
        for j in range(i+1,5):
            minrmsd = 99999999999
            for k in range(8):
                jsel.frame = j + 5*k
                rmsd = isel.rmsdQCP(jsel)
                if rmsd < minrmsd:
                    minrmsd = rmsd
            rmsdmatrix[i][j] = minrmsd
    return rmsdmatrix















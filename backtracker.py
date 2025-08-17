from chemistry import Atom,Bond
from moleculeClass import Molecule
import copy
from tokenizer import*
from helperFunctions import*
from collections import deque



def getSkeletalStructuresDFS(molecule,moleculeList=None,visited=None):
    print("hi")
    """gets all skeletal structures using a DFS algorithm"""
    if moleculeList is None:
        moleculeList=[]
    if visited is None:
        visited=set()
    strMol=molecule.molToStr()
    if strMol in visited:
        return moleculeList
    visited.add(strMol)
    isAllBondedOrCircular=molecule.isBondedOrCircular()

    if isAllBondedOrCircular[0] and (not isAllBondedOrCircular[1]):
        if not (badCarbon(molecule) and len(moleculeList)>=1):
             moleculeList.append(molecule.cloneMolecule())
        return moleculeList
    length=len(molecule.atoms)
    for index in range(length-1):
        for index1 in range(index+1,length):
            atomOne=molecule.atoms[index]
            atomTwo=molecule.atoms[index1]
            canBond=atomOne.bondTo(atomTwo,"-")
            if canBond and molecule.isValid():
               if (hasOverOxyNitro(molecule) or
                     badHalogen(molecule)) and len(moleculeList)>=1:
                        atomOne.unBond(atomTwo)
                        continue
               getSkeletalStructuresDFS(molecule,moleculeList,visited)
            if canBond:
               atomOne.unBond(atomTwo)
    return moleculeList

    


    
    
def backtrackSkeletalStructure2(structure,moleculeList,
                               visited=None,bondList=None):
    """given a molecule that is all single bonded and has octets on all atoms
    add pi bonds until its complete structure inputted into this function 
    might have too many electrons"""
  #  prettyPrint(structure)
    if visited is None:
        visited=set()
    strMol=structure.molToStr()
    if strMol in visited:
        print("n")
        return moleculeList
    visited.add(strMol)
    if structure.isComplete():
             moleculeList.append(structure.cloneMolecule())
             return moleculeList
    bondList=structure.getBondList()

    #Atoms with higher bond capacity want to have more pi bonds
    for index in range(len(bondList)):
        bond=bondList[index]
       #cant add a pi bond to a triple bond
        if bond.type=="#":
            continue

        atomOne=bond.atomOne
        atomTwo=bond.atomTwo

        #if it has too little electrons return
        if structure.hasEnoughElectrons()==None:
            return moleculeList
        removeLpOne=atomOne.removeLonePair()
        removeLpTwo=atomTwo.removeLonePair()
        
    #if both lone pairs cant be removed, it means at least atom has no lone pairs
        if not (removeLpOne and removeLpTwo):
            continue
        bond.addPi()
        #after adding a pi, send the bond to the back of the list
        bondList.append(bondList.pop(index))
        backtrackSkeletalStructure2(structure,moleculeList,visited,bondList)

        bondList.insert(index,bondList[-1])
        bondList.pop(-1)
        bond.removePi()
        atomOne.addLonePair()
        atomTwo.addLonePair()
    
    return moleculeList
    


    



def backtrackSkeletalStructure(structure,moleculeList,
                                minScore=[9999999999],visited=set()):
    """given an all sigma bonded molecule with no lone pairs, adds pi bonds 
    and lp's until complete"""
    print('hi')
    strMol=structure.molToStr()
    if strMol in visited:
        return (moleculeList,minScore[0])
    visited.add(strMol)

    score=structure.getScore()
    #prune bad scores as they appear
    if score>=minScore[0]+2:
           # pass
            return (moleculeList,minScore[0])
    
    if structure.isComplete():
        moleculeList.append(structure.cloneMolecule())
        if score<minScore[0]:
            minScore[0]=score
        return (moleculeList,minScore[0])
    
    for atom in structure.atoms:
        atom.addLonePair()
        if structure.isValid():
            backtrackSkeletalStructure(structure,moleculeList,minScore,visited)
        atom.removeLonePair()

    bondList=structure.getBondList()
    for bond in bondList:
        if bond.type=="#":
            continue
        bond.addPi()
        if structure.isValid():
            backtrackSkeletalStructure(structure,moleculeList,minScore,visited)
        bond.removePi()
    return (moleculeList,minScore[0])
 

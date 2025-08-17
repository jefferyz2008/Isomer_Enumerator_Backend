from collections import deque
from chemistry import Atom,Bond
from moleculeClass import Molecule
import copy
from tokenizer import*
from helperFunctions import*

def getSkeletalStructuresBFS(molecule):
    #returns a list of all single bonded molecules
    moleculeList=[]
    visited=set()
    queue=deque()
    queue.append(molecule)
    visited.add(molecule.molToStr())
    while queue:
       # print(len(moleculeList))
        current=queue.popleft()
      #  prettyPrint(current)
        length=len(current.atoms)
        for index in range(length-1):
            for index1 in range(index+1,length):
                copyMol=current.cloneMolecule()
                copyAtom1=copyMol.atoms[index]
                copyAtom2=copyMol.atoms[index1]
                canBond=copyAtom1.bondTo(copyAtom2,"-")
                if not canBond:
                    continue

                strMol=copyMol.molToStr()
                if strMol in visited:
                   continue
                visited.add(strMol)
                #if oxygen has 3 single bonds
                #or nitrogren has 4 single bonds
                #or halogen has over 1 bond
                #abandon it if theres an alternative
                if (hasOverOxyNitro(copyMol) or
                     badHalogen(copyMol)) and len(moleculeList)>=1:
                        continue
                isBondedOrCircular=copyMol.isBondedOrCircular()
                #if all atoms are bonded and its not circular
                if (isBondedOrCircular[0] and (not isBondedOrCircular[1]) 
                    and copyMol.isValid()):
                    #checks for carbons with 1 bond
                    if badCarbon(copyMol) and len(moleculeList)>=1:
                        continue
                    moleculeList.append(copyMol.cloneMolecule())
                    continue
                if copyMol.isValid():
                    queue.append(copyMol)
    return moleculeList

def bfsSkeletalStructure(structure,minScore=None):
    moleculeList=[]
    #if it's complete return
    if minScore is None:
       minScore=float("inf")
    if structure.isComplete():
        moleculeList.append(structure)
        return (moleculeList,minScore)
    queue=deque()
    queue.append(structure.cloneMolecule())
    visited=set()
    while queue:
        current=queue.popleft()
    
        score=current.getScore()
        if score>minScore:
            continue
        strMol=current.molToStr()
        if strMol in visited:
            print(strMol)
            continue
        visited.add(strMol)
        if current.isComplete():
           moleculeList.append(current.cloneMolecule())
           score=current.getScore()
           if score<minScore:
               minScore=score
           continue
        allBonds=current.getBondList()
        for bond in allBonds:
            if not bond.type=="#":
                bond.addPi()
                if current.isValid():
                   queue.append(current.cloneMolecule())
                bond.removePi()

            atomOne=bond.atomOne
            atomOne.addLonePair()
            if current.isValid():
                queue.append(current.cloneMolecule())
            atomOne.removeLonePair()

            atomTwo=bond.atomTwo
            atomTwo.addLonePair()
            if current.isValid():
                queue.append(current.cloneMolecule())
            atomTwo.removeLonePair()
    return (moleculeList,minScore)
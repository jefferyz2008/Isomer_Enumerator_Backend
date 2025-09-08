from chemistry import Atom,Bond
from moleculeClass import Molecule
from tokenizer import*
import copy
from helperFunctions import*
import random

def getSkeletalStructuresDFS(molecule,moleculeList=None,visited=None):
    """gets all skeletal structures using a DFS algorithm"""
    if moleculeList is None:
        moleculeList=[]
    if visited is None:
        visited=set()
    strMol=molecule.molToStr()
    if strMol in visited:
        return moleculeList
    visited.add(strMol)
    if molecule.bondElectrons/2>=molecule.numAtoms-1:
       isAllBondedOrCircular=molecule.isBondedOrCircular()
    else:
        isAllBondedOrCircular=(False,False)
    #if it detects a cyclic molecule
    if isAllBondedOrCircular[1]:
        return moleculeList
    #the whole structure must be connected to be considered complete
    if isAllBondedOrCircular[0]:
        if not (badCarbon(molecule) and moleculeList):
             clone=molecule.cloneMolecule()
             clone.bondElectrons=molecule.currentElectrons
             moleculeList.append(clone)
        return moleculeList
    for index in range(molecule.numAtoms-1):
        atomOne=molecule.atoms[index]
        if (atomOne.electronDomains and moleculeList and 
            checkOverBonding(atomOne)):
            continue
        for index1 in range(index+1,molecule.numAtoms):
            atomTwo=molecule.atoms[index1]
            if (atomTwo.electronDomains and moleculeList
                 and checkOverBonding(atomTwo)):
                continue

            oldStrOne=atomOne.strAtom
            oldStrTwo=atomTwo.strAtom

            newBond=atomOne.sigmaBond(atomTwo)
            molecule.bondElectrons+=2


            #don't continue if bond unsuccesful
            if not newBond:
                continue
            else:
               molecule.currentElectrons+=2

            #eliminates peroxides if applicable
            if newBond and (newBond.sameTypeAtoms("O-O") and 
                moleculeList):
               atomOne.unBond(atomTwo,newBond)
               molecule.currentElectrons-=2
               atomOne.strAtom=oldStrOne
               atomTwo.strAtom=oldStrTwo
               continue

            #Both atoms not having an over octet means its valid
            if (newBond and atomOne.hasOctet()!=False 
                and atomTwo.hasOctet()!=False):
               atomOne.strAtom=atomOne.atomToStr()
               atomTwo.strAtom=atomTwo.atomToStr()
               getSkeletalStructuresDFS(molecule,moleculeList,visited)

            if newBond:
               atomOne.unBond(atomTwo,newBond)
               molecule.currentElectrons-=2
               molecule.bondElectrons-=2
               atomOne.strAtom=oldStrOne
               atomTwo.strAtom=oldStrTwo

               
    return moleculeList

    

def backtrackPiBondsLonePairs(structure,moleculeList,bondList,
                                minScore,visited):
    """given a skeletal structure, returns a list of complete molecules"""
    if structure.expandedOctet:
        return []
    if bondList is None:
        bondList=structure.getBondListNoH()

    #there can only be 1 perfect structure, so if there is, stop the search
    #if minScore[0]==0 and moleculeList:
      #  return (moleculeList,minScore)
    
    lowestScore=minScore[0]
    currentScore=structure.formalChargeSum
    if currentScore>lowestScore:
        return (moleculeList,minScore)
    
    strMol=structure.molToStr()
    if strMol in visited:
        return (moleculeList,minScore)
    visited.add(strMol)
    idealElectrons=structure.idealBondElectrons
    #if the structure has the ideal amount of bond electrons, 
    #it moves to lone pairs
    if structure.bondElectrons==idealElectrons:

        addedLonePairs=addLonePairs(structure,minScore[0])
        #this means that the score was too high while adding lone pairs
        if not addedLonePairs or (structure.formalChargeSum>lowestScore):
            #reverts to previous formal charge sum
            return (moleculeList,minScore)
        else:
            formalChargeSum=structure.formalChargeSum
            if formalChargeSum<=lowestScore:
                minScore[0]=formalChargeSum
                if addedLonePairs:
                    moleculeList.append(structure.cloneMolecule())
            structure.removeAllLonePairs()
            #we need to revert the formal Charge Sum back
            structure.formalChargeSum=currentScore
            return (moleculeList,minScore)
    
    octetDict=structure.octetDict
    for bond in sorted(bondList, key=lambda b:(b.numModifications,
                        -1*min(b.atomOne.prefBonds,b.atomTwo.prefBonds))):
        if bond.type=="≡":
            continue

        atomOne=bond.atomOne
        atomTwo=bond.atomTwo
        if octetDict[atomOne] or octetDict[atomTwo] or atomOne in {"H"}:
            continue
        oldStrOne=atomOne.strAtom
        oldStrTwo=atomTwo.strAtom

        bond.addPi()
        structure.bondElectrons+=2
        structure.currentElectrons+=2

        hasOctetOne=atomOne.hasOctet()
        hasOctetTwo=atomTwo.hasOctet()
        
        if hasOctetOne:
            octetDict[atomOne]=True
            formalChargeOne=abs(atomOne.getFormalCharge())
            structure.formalChargeSum+=formalChargeOne
        
        if hasOctetTwo:
            octetDict[atomTwo]=True
            formalChargeTwo=abs(atomTwo.getFormalCharge())
            structure.formalChargeSum+=formalChargeTwo
        
        atomOne.strAtom=atomOne.atomToStr()
        atomTwo.strAtom=atomTwo.atomToStr()
        backtrackPiBondsLonePairs(structure,moleculeList,bondList,
                                  minScore,visited)
        if hasOctetOne:
            octetDict[atomOne]=False
            structure.formalChargeSum-=formalChargeOne

        if hasOctetTwo:
            octetDict[atomTwo]=False
            structure.formalChargeSum-=formalChargeTwo
        
        bond.removePi()
        atomOne.strAtom=oldStrOne
        atomTwo.strAtom=oldStrTwo
        structure.bondElectrons-=2
        structure.currentElectrons-=2
    return (moleculeList,minScore)


def addLonePairs(structure,minScore): #takes in minscore for pruning
    """adds lone pairs on every atom of a molecule until its complete."""
    prevFormalChargeSum=structure.formalChargeSum
    octetDict=structure.octetDict
    for atom in structure.atoms:
        if atom.symbol in {"H","B","Al","Be"}:
            continue
        if octetDict[atom]:
            continue
        #adds until every atom has an octet
        atom.addUntilOctet()

        #if it already had an octet, don't add to formal charge sum
        atomFC=abs(atom.getFormalCharge())
        structure.formalChargeSum+=atomFC
        #reverts structure to its previous form.
        if structure.formalChargeSum>minScore:
            structure.formalChargeSum=prevFormalChargeSum
            structure.removeAllLonePairs()
            return False
    return True



#The following algorithms don't work as efficiently as this one

def backtrackSkeletalStructure(structure,moleculeList,bondList,atomList,
                                minScore=[float('inf')],visited=set()):
    """given an all sigma bonded molecule with no lone pairs, adds pi bonds 
    and lp's until complete"""
  #  prettyPrint(structure)
    if not bondList:
       #all bonds in the molecule except hydrogen, because we can't add pi
       #bonds to it
       bondList=structure.getBondListNoH()
    
    if not atomList:
         atomList=structure.atoms
         #we don't give H lone pairs
         atomList=[atom for atom in atomList if atom.symbol!="H"]

    #if we've found a perfect structure, we don't try anything else
    if moleculeList and minScore[0]==0 and not structure.expandedOctet:
        return (moleculeList,minScore[0])
    
    score=(structure.formalChargeSum if not structure.expandedOctet 
           else structure.sumFormalCharges())
    
    if score>minScore[0]:
       return (moleculeList,minScore[0])
    strMol=structure.molToStr()
    if strMol in visited:
        return (moleculeList,minScore[0])
    visited.add(strMol)

    completeYet=structure.isComplete()
    #in expanded octets, an atom having an octet doesn't dictate
    #that it's done adding electrons, so a different approach to formal 
    #charge pruning is needed.
    if completeYet:
        clone=structure.cloneMolecule()
        moleculeList.append(clone)
        if score<minScore[0]:
            minScore[0]=score
        return (moleculeList,minScore[0])
    octetDict=structure.octetDict
    for atom in atomList:
        #if it has an octet, no lone pairs can be added
        #except if its an expanded octet
        if octetDict[atom] and not atom.canExpandOctet:
            continue
        #saves string state for backtracking
        oldStr=atom.strAtom
        atom.addLonePair()

        structure.currentElectrons+=2
        atomHasOctet=atom.hasOctet()
        if atomHasOctet:
            octetDict[atom]=True
            atomFC=abs(atom.getFormalCharge())
            structure.formalChargeSum+=atomFC
        #hasOctet and hasEnoughElectrons 
        #will return None if it's octet isn't complete yet
        if atomHasOctet!=False and structure.hasEnoughElectrons()!=False:
            atom.strAtom=atom.atomToStr(True)
           #print("NEJFJNFIVEF") prettyPrint(structure)
            backtrackSkeletalStructure(structure,moleculeList,bondList
                                              ,atomList,minScore,visited)
        
        if atomHasOctet:
           #if we are removing it's octet, send it back to the last
           #formal charge sum
           octetDict[atom]=False
           structure.formalChargeSum-=atomFC
        atom.removeLonePair()
        structure.currentElectrons-=2
        atom.strAtom=oldStr
    for bond in sorted(bondList, key=lambda b:(b.numModifications,
                        -1*min(b.atomOne.prefBonds,b.atomTwo.prefBonds))):
        if bond.type=="≡":
            continue
        atomOne=bond.atomOne
        atomTwo=bond.atomTwo
        #we must make sure that none of the atoms have octets
        if ((octetDict[atomOne] and not atomOne.canExpandOctet) or 
            (octetDict[atomTwo] and not atomTwo.canExpandOctet)):
            print("Hi")
            continue
        #saves string states for backtracking
        oldStrOne=atomOne.strAtom
        oldStrTwo=atomTwo.strAtom
        print("JFNVFEV")
        assert(bond.addPi())
        prettyPrint(structure)
        structure.currentElectrons+=2
        structure.bondElectrons+=2
        atomOneHasOctet=atomOne.hasOctet()
        atomTwoHasOctet=atomTwo.hasOctet()
        if atomOneHasOctet:
            octetDict[atomOne]=True
            atomOneFC=abs(atomOne.getFormalCharge())
            structure.formalChargeSum+=atomOneFC
        if atomTwoHasOctet:
            octetDict[atomTwo]=True
            atomTwoFC=abs(atomTwo.getFormalCharge())
            structure.formalChargeSum+=atomTwoFC
        if (atomOneHasOctet!=False and atomTwoHasOctet!=False and 
            structure.hasEnoughElectrons()!=False):
            atomOne.strAtom=atomOne.atomToStr(True)
            atomTwo.strAtom=atomTwo.atomToStr(True)
            backtrackSkeletalStructure(structure,moleculeList,
                                            bondList,atomList,minScore,visited)

        """reverts molecule back to original state"""
        #if we removed an octet on any of the atoms, we subtract the formal 
        #charge from the sum
        if atomOneHasOctet:
            octetDict[atomOne]=False
            structure.formalChargeSum-=atomOneFC
        if atomTwoHasOctet:
            octetDict[atomTwo]=False
            structure.formalChargeSum-=atomTwoFC
        bond.removePi()
        structure.currentElectrons-=2
        structure.bondElectrons-=2
        atomOne.strAtom=oldStrOne
        atomTwo.strAtom=oldStrTwo
    
    return (moleculeList,minScore[0])



    





##################
def backtrackSkeletalStructure2(structure,moleculeList,
                               visited=None,bondList=None):
    #note this doesnt work yet
    """given a molecule that is all single bonded and has octets on all atoms
    add pi bonds until its complete structure inputted into this function 
    might have too many electrons"""
  #  prettyPrint(structure)
    if visited is None:
        visited=set()
    strMol=structure.molToStr()
    visited.add(strMol)
    if structure.isComplete():
             moleculeList.append(structure.cloneMolecule())
             return moleculeList
    bondList=structure.getBondList()

    #Atoms with higher bond capacity want to have more pi bonds
    for index in range(len(bondList)):
        bond=bondList[index]
       #cant add a pi bond to a triple bond
        if bond.type=="≡":
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

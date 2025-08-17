from chemistry import Bond, Atom
from helperFunctions import*
import copy
from collections import deque

class Molecule:
     def __init__(self,atoms,charge):
          #putting atoms with higher bonding capacity first
          # makes adding bonds easier
          atoms.sort(key=lambda atom: (atom.prefBonds,
                                        atom.atomicNumber), reverse=True)
          self.atoms=atoms
          self.charge=charge
          #number of electrons molecule should have
          self.numElectrons=sum(atom.valenceElectrons for
                  atom in self.atoms) - self.charge

     def sortAtoms(self):
         """sorts atoms for turning the molecule into a string"""
         self.atoms.sort(key=lambda atom:(len(atom.electronDomains),
                        atom.countBonds(),atom.getMolarMassSum(),atom.molarMass)
                        ,reverse=True)
     
     def getTotalElectrons(self):
         """returns the amount of electrons it does have"""
         totalElectrons=0
         for atom in self.atoms:
             for domain in atom.electronDomains:
                 ##each bond appears twice, so we divide by 2
                 if isinstance(domain,Bond):
                     totalElectrons+=domain.electrons/2
                 elif domain==":":
                     totalElectrons+=2
         return totalElectrons
     
     def hasEnoughElectrons(self):
         """returns True if the molecule has exactly enough electrons
            returns False if it has too many
            and returns None if it has too few"""
         totalElectrons=self.getTotalElectrons()
         numElectrons=self.numElectrons
         if totalElectrons>numElectrons:
             return False
         elif totalElectrons<numElectrons:
             return None
         return True

     def getMolarMass(self):
        """returns total molar mass of molecule"""
        return sum(atom.molarMass for atom in self.atoms)
     
     def isComplete(self):
         """returns True if the molecule is fully complete"""
         if not self.checkFormalCharges():
             return False
         for atom in self.atoms:
             if not atom.hasOctet():
                return False
         return self.hasEnoughElectrons()
    
     def checkFormalCharges(self):
         #sum of formal charges is the charge of molecule
         formalChargeSum=0
         for atom in self.atoms:
             formalCharge=atom.getFormalCharge()
             formalChargeSum+=formalCharge
         return formalChargeSum==self.charge
     

     def isValid(self):
         """returns true if molecule is valid"""
         for atom in self.atoms:
             if atom.hasOctet()==False:
                 return False
        # not too many electrons
         if self.hasEnoughElectrons()==False:
             return False
         return True
     
     def rearrange(self):
         """rearranges the electron domains of each atom"""
         for atom in self.atoms:
             atom.rearrange()
    
     def isBondedOrCircular(self):
        """returns a tuple with 2 bools
           the first represents whether the molecule is completely bonded
           the second represents if it's circular"""
        isCircular=False
        queue=deque()
        visited=set()
        atom=self.atoms[0]
        queue.append(atom)
        bondSet=set()
        while queue:
            current=queue.popleft()
            if not isCircular:
               if current in visited:
                   isCircular=True
                   break
            visited.add(current)

            for domain in current.electronDomains:
                atomOne=domain.atomOne
                atomTwo=domain.atomTwo
                if isinstance(domain,Bond):
                    if domain in bondSet:
                        continue
                    bondSet.add(domain)
                    if (not (current is atomOne)):
                        queue.append(atomOne)
                    else: 
                        queue.append(atomTwo)

        #want to see if every atom in the molecule is visited
        return (len(visited)==len(self.atoms),isCircular)
     
     def getBondList(self):
        """a list of all the bonds of the atom"""
        bondList=[]
        for atom in self.atoms:
             for domain in atom.electronDomains:
                 if isinstance(domain,Bond) and domain not in bondList:
                     bondList.append(domain)
        bondList.sort(key=lambda bond: 
                min(bond.atomOne.valenceElectrons,bond.atomTwo.valenceElectrons)
                ,reverse=True)
        return bondList
     
     def getBondListNoH(self):
         """returns all bonds except the ones that contain Hydrogen"""
         bondList=self.getBondList()
         return [bond for bond in bondList 
                if bond.atomOne.symbol != "H" and bond.atomTwo.symbol != "H"]
     
     def sumFormalCharges(self):
        """gets sum of all absolute value of formal charges in a molecule"""
        charges = [abs(atom.getFormalCharge()) for atom in self.atoms]
        return sum(charges)
     
     def getScore(self):
         """returns an int that says how stable the molecule is, lower
         scores preffered"""
         score=self.sumFormalCharges()
         bondList=self.getBondList()
         for bond in bondList:
             bondType=bond.type
             atomOne_=bond.atomOne
             atomTwo_=bond.atomTwo

             if bondType=="-":
                 #all these single bonds are unstable
                 #if atoms don't have octets, they can still omit these bonds
                 if not (atomOne_.hasOctet() and atomTwo_.hasOctet()):
                     continue
                 if bond.sameTypeAtoms("O-O"):
                     score+=10
                 if bond.sameTypeAtoms("N-N"):
                     score+=10
                 if bond.sameTypeAtoms("F-F"):
                     score+=20
                 if bond.sameTypeAtoms("O-F"):
                     score+=7
             #these triple bonds are unstable
             if bondType=="#":
                 
                 if atomOne_.symbol=="O" or atomTwo_.symbol=="O":
                     score+=15
                 if atomOne_.symbol=="F" or atomTwo_.symbol=="F":
                     score+=15
                    
            #these bonds with flourine are unstable
             if bond.sameTypeAtoms("C-F"):
                 if bondType=="#":
                     score+=30
                 if bondType=="=":
                     score+=10   
            #Oxygen doesn't like 3 bonds    
             if bond.atomOne.symbol=="O" and bond.atomOne.countBonds()>2:
                 score+=15
             if bond.atomTwo.symbol=="O" and bond.atomTwo.countBonds()>2:
                 score+=15

         return score
     
     def hasAtomThere(self,xPos,yPos):
         """returns true if an atom overlaps another atom"""
         for atom in self.atoms:
             if getDistance(atom.centerX,atom.centerY,xPos,yPos)<=10:
                 return True
         return False
         
     def allHasPositions(self):
         """checks if every atom is assigned to a centerX and centerY"""
         for atom in self.atoms:
             if not (atom.centerX and atom.centerY):
                 return False 
         return True
    
     def molToStr(self):
         """turns a molecule into a string unique to that molecule"""
         copyMol=self.cloneMolecule()
         copyMol.rearrange()
         copyMol.sortAtoms()
         returnStr=""
         for atom in copyMol.atoms:
             returnStr+=atom.symbol+"|"
             for domain in atom.electronDomains:
                 if isinstance(domain,Bond):
                     returnStr+=domain.type+domain.getOther(atom).symbol
                 else:
                     returnStr+=":"
         return returnStr

     def cloneMolecule(self):
        """returns a different object, but exact copy of all atoms and bonds"""
        oldToNew={}
        for atom in self.atoms:
            oldToNew[atom]=Atom(atom.symbol)
            for _ in range(atom.electronDomains.count(":")):
                oldToNew[atom].addLonePair()
        allBonds=self.getBondList()
        for bond in allBonds:
            atom1=oldToNew[bond.atomOne]
            atom2=oldToNew[bond.atomTwo]
            newBond=Bond(bond.type,atom1,atom2)
            atom1.electronDomains.append(newBond)
            atom2.electronDomains.append(newBond)
        atomList=[]
        for key in oldToNew:
            atomList.append(oldToNew[key])
        return Molecule(atomList,self.charge)
  
     def addUntilOctet(self):
         """for all atoms in the molecule, give them all octets"""
         for atom in self.atoms:
             atom.addUntilOctet()
         return True
     
     def assignPositions(self,width,height):#positionsList):
        """assign's the best positions to an atom for drawing"""
        queue=deque()
        bondSet=set()
        self.atoms[0].centerX=width/2
        self.atoms[0].centerY=height/2
        queue.append((None,self.atoms[0]))
        while queue:
            atomTuple=queue.popleft()
            current=atomTuple[1]
            predecessor=atomTuple[0]
            if predecessor:
                for dX,dY in [(-75,0),(75,0),(0,-75),(0,75)]:
                    if not self.hasAtomThere(predecessor.centerX+dX,
                                             predecessor.centerY+dY):
                        current.centerX=predecessor.centerX+dX
                        current.centerY=predecessor.centerY+dY
                        break
            for domain in current.electronDomains:
                if isinstance(domain,Bond) and domain not in bondSet:
                    queue.append((current,domain.getOther(current)))
                    bondSet.add(domain)
        for bond in self.getBondList():
            bond.startPos=(bond.atomOne.centerX,bond.atomOne.centerY)  
            bond.endPos=(bond.atomTwo.centerX,bond.atomTwo.centerY) 
        return True
     
     def hasSingleCenterAtom(self):
         """returns true if the molecule has a single central atom, used 
         for expanded octets"""
         molecule=self.cloneMolecule()
         molecule.atoms.sort(key=lambda atom: atom.valenceElectrons,
                             reverse=True)
         #the central atom has the most valence electrons
         centralAtom=molecule.atoms[0]
         visited=set()
         visited.add(centralAtom)
         #see if every atom can be explored through the central atom
         for domain in centralAtom.electronDomains:
             if isinstance(domain,Bond):
                 visited.add(domain.getOther(centralAtom))
         return len(visited)==len(self.atoms)
                 
         
     
     def getLonePairs(self):
         """assigns positions to lone pairs"""
         positions={}
         for atom in self.atoms:
             numLonePairs=atom.countDomains(":")
             if True:
                addedLPs=0
                for dX,dY in [(-75,0),(75,0),(0,-75),(0,75)]:
                    if not self.hasAtomThere(atom.centerX+dX,
                                             atom.centerY+dY):
                        if addedLPs>=numLonePairs:
                            print("hi")
                            break
                        positions[str(atom.centerX+dX/4)+","+
str(atom.centerY+dY/4)]="90" if ((dX,dY)==(-75,0) or (dX,dY)==(75,0)) else "180"
                        addedLPs+=1
                        print(addedLPs,numLonePairs)

         return positions
                 
             
             
         
         
     
     



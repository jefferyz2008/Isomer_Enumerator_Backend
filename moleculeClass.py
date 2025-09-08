from chemistry import Bond, Atom
from helperFunctions import*
import copy
from collections import deque
from specialDict import*
from queue import PriorityQueue
class Molecule:
     __slots__ = [
        "atoms", "charge", "formula", "numElectrons", "currentElectrons",
        "formalChargeSum", "bondElectrons", "idealBondElectrons",
        "octetDict", "numAtoms", "expandedOctet"
    ]
     def __init__(self,atoms,charge,formula):
          #putting atoms with higher bonding capacity first
          #to make backtracking more efficient
          atoms.sort(key=lambda atom: (atom.prefBonds,
                                      atom.atomicNumber), reverse=True)
          self.atoms=atoms
          self.charge=charge
          self.formula=formula
          #number of electrons molecule should have
          self.numElectrons=sum(atom.valenceElectrons for
                  atom in self.atoms) - self.charge
          
          #the number of electrons it does have
          self.currentElectrons=0

          #for formal charges
          self.formalChargeSum=0

          #number of electrons contained in bonds
          self.bondElectrons=0

          #the number of electrons that should be contained in bonds
          self.idealBondElectrons=(sum([atom.octetElectrons for atom in atoms])-
                                   self.numElectrons)

          #this is a dictionary with an atom as a key and a bool, true or false
          #representing whether the atom has an octet.
          self.octetDict=dictWithCounts()
          for atom in atoms:
              self.octetDict[atom]=True if (atom.hasOctet()==True) else False

          #controls expanded octet if it clearly doesnt have one
          self.numAtoms=len(atoms)
          if self.numAtoms<3 or self.numAtoms>8:
              for atom in atoms:
                  atom.canExpandOctet=False
          
          for atom in atoms:
              if self.count(atom.symbol)!=1:
                  atom.canExpandOctet=False
          self.expandedOctet=False
          for atom in atoms:
              if atom.canExpandOctet:
                  self.expandedOctet=True
                  
          
     def sortAtoms(self):
         """returns a list of sorted atoms for turning the molecule into a string"""
         atomList=self.atoms[:]
         atomList.sort(key=lambda atom:(atom.molarMass,len(atom.electronDomains),
                        atom.countBonds(),atom.getMolarMassSum(),atom.molarMass)
                        ,reverse=True)
         return atomList

     def count(self,symbol):
         assert(isinstance(symbol,str))
         """counts the number of a given type of atom in molecule"""
         return len([atom for atom in self.atoms if atom.symbol==symbol])
     
     
     def hasEnoughElectrons(self):
         """returns True if the molecule has exactly enough electrons
            returns False if it has too many
            and returns None if it has too few"""
         
         #expanded octet is an edge case:
         electrons=self.currentElectrons
         shouldHave=self.numElectrons
         return True if electrons==shouldHave else (False if 
                                        electrons>shouldHave else None)

     def getMolarMass(self):
        """returns total molar mass of molecule"""
        return sum(atom.molarMass for atom in self.atoms)
     
     def isComplete(self):
         """returns True if the molecule is fully complete"""
         octetDict=self.octetDict
         if octetDict.hasValue(False):
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
                min(bond.atomOne.prefBonds,bond.atomTwo.prefBonds)
                ,reverse=True)
        return bondList
     
     def getBondListNoH(self):
         """returns all bonds except the ones that contain Hydrogen"""
         bondList=self.getBondList()
         return [bond for bond in bondList 
                if bond.atomOne.symbol != "H" and bond.atomTwo.symbol != "H"]
     
     def sumFormalCharges(self):
        """gets sum of all absolute value of formal charges in a molecule
           used for expanded octet purposes"""
        charges = [abs(atom.getFormalCharge()) for atom in self.atoms if
                    (not (atom.canExpandOctet) and atom.hasOctet())]
        return sum(charges)
                       
     def getScore(self,isComplete=True):
                       #isComplete is there so we don't punish molecules that
                       #aren't complete yet
         """returns an int that says how stable the molecule is, lower
         scores preffered"""
         if self.expandedOctet and not (self.isComplete()):
             #if its expanded octet, we must recalculate formal charges
             #a different way
             score=self.sumFormalCharges()
         else:
             score=self.formalChargeSum
         bondList=self.getBondList()
         for bond in bondList:
             bondType=bond.type
             atomOne=bond.atomOne
             atomTwo=bond.atomTwo

             if bondType=="-" and isComplete:
                 #all these single bonds are unstable
                 if bond.sameTypeAtoms("O-O"):
                     score+=20
                 if bond.sameTypeAtoms("N-N"):
                     score+=10
                 if bond.sameTypeAtoms("F-F"):
                     score+=20
                 if bond.sameTypeAtoms("O-F"):
                     score+=7
                 if bond.sameTypeAtoms("O-N"):
                     score+=1/2
         if isComplete:
            score+= self.getFormalChargeSharing(self.sumFormalCharges())
          #  score+=self.punishFormalChargesElectroneg()
         return score
     
     
     
     def punishFormalChargesElectroneg(self):
        """
        Penalizes molecules for chemically unfavorable formal charges:
        - Negative FC on less electronegative atoms
        - Positive FC on highly electronegative atoms
        Returns an integer penalty score; higher = worse.
        """
        penalty = 0
        maxEN = 4.0  # Fluorine as reference max electronegativity

        for atom in self.atoms:
            fc = atom.getFormalCharge()
            en = atom.electroNegativity

            if fc < 0:
                # Negative FC on low EN atoms is bad
                # Low EN atoms have more penalty: (maxEN - EN)
                penalty += int((maxEN - en) * abs(fc) * 3)
            elif fc > 0:
                # Positive FC on high EN atoms is bad
                # High EN atoms have more penalty: EN itself
                penalty += int(en * abs(fc) * 3)

        return penalty

         
     def getFormalChargeSharing(self,sumCharges):
         """returns a number indicating how well spread the formal charge is"""
        
         #this is the absolute value of all formal charges
         sumCharges=self.sumFormalCharges()
         #number of atoms in the molecule that has a formal charge
         numChargeAtoms=len([atom for atom in  self.atoms if 
                                atom.getFormalCharge()])
         if not numChargeAtoms:
             return 0
         return sumCharges/numChargeAtoms
     
    
     def molToStr(self,polarityCheck=False):
         """turns a molecule into a string unique to that molecule"""

         atomList=self.sortAtoms()
         returnStr=""
         for atom in atomList:
             returnStr+=atom.strAtom
         return returnStr

     def cloneMolecule(self):
        """returns a different object, but exact copy of all atoms and bonds"""
        oldToNew={}
        for atom in self.atoms:
            oldToNew[atom]=Atom(atom.symbol)
            newAtom=oldToNew[atom]
            newAtom.canExpandOctet=atom.canExpandOctet
            for _ in range(atom.electronDomains.count(":")):
                newAtom.addLonePair()
            newAtom.currentElectrons=atom.currentElectrons
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
        newMolecule=Molecule(atomList,self.charge,self.formula)
        newMolecule.formalChargeSum=self.formalChargeSum
        newMolecule.expandedOctet=self.expandedOctet
        newMolecule.currentElectrons=self.currentElectrons
        newMolecule.bondElectrons=self.currentElectrons
        return newMolecule
  
     def addUntilOctet(self):
         """for all atoms in the molecule, give them all octets"""
         for atom in self.atoms:
             atom.addUntilOctet()
         return True
     
     def removeAllLonePairs(self):
         """removes all lone pairs from a molecule"""
         removedLonePairs=0
         for atom in self.atoms:
             while atom.removeLonePair():
                 removedLonePairs+=1
         self.currentElectrons-=removedLonePairs*2
         return True
             
     
     def getCenterAtom(self):
         """returns center Atom for expanded octets"""
         for atom in self.atoms:
             if self.count(atom.symbol)==1:
                 return atom
         return self.atoms[0]
    
     def hasSingleCenterAtom(self):
         """returns true if the molecule has a single central atom, used 
         for expanded octets"""
         molecule=self.cloneMolecule()
         molecule.atoms.sort(key=lambda atom: atom.valenceElectrons,
                             reverse=True)
         #the central atom has the most valence electrons
         centralAtom=self.getCenterAtom()
         if centralAtom.countDomains("bond")==1:
             return False
         visited=set()
         visited.add(centralAtom)
         #see if every atom can be explored through the central atom
         for domain in centralAtom.electronDomains:
             if isinstance(domain,Bond):
                 visited.add(domain.getOther(centralAtom))
         return len(visited)==len(self.atoms)
                 

     def split(self,atom,bond):
         """returns a 'molecule' reprsenting the frontier of the atom with 
         respect to the given bond. Used for determing polarity"""
         bondSet=set()
         bondSet.add(bond)
         queue=deque()
         queue.append(bond.getOther(atom))
         allAtoms=[]
         while queue:
             current=queue.popleft()
             allAtoms.append(current)
             for domain in current.electronDomains:
                 if isinstance(domain,Bond) and domain not in bondSet:
                     bondSet.add(domain)
                     queue.append(domain.getOther(current))
         newMolecule=Molecule(allAtoms,0,"")
         return newMolecule.molToStr(True)
    
     def getPolarity(self):
         """returns polar if molecule is polar nonpolar otherwise"""
         hasNotHydroCarbon=False
         #if it has lone pairs and it's not square planar or linear, then
         #it's nonpolar
         for atom in self.atoms:
             if (not hasNotHydroCarbon) and atom.symbol!="H" and atom.symbol!="C":
                 hasNotHydroCarbon=True
             if atom.countDomains("bond")==1:
                 continue
             vsepr=atom.getVSEPR()
             if atom.countDomains(":") and (not (vsepr[0]=="linear" 
                                                 or vsepr[0]=="Square Planar")):
                 return "polar"
         if not hasNotHydroCarbon:
             return "non-polar"
         
         #if theres an odd number of atoms or theres a central atom
         if len(self.atoms)%2==1 or self.hasSingleCenterAtom():
             for atom in self.atoms:
                 if atom.symbol=="H":
                     continue
                 if atom.countDomains("bond")==1:
                     continue
                 frontierSet=set()
                 for domain in atom.electronDomains:
                     if domain==":":
                         continue
                     frontierSet.add(self.split(atom,domain))
                 if len(frontierSet)==1:
              
                    return "non-polar"
             return "polar"
             
         else:
             for bond in self.getBondList():
                 if self.split(bond.atomOne,bond)==self.split(bond.atomTwo,bond):
                     return "non-polar"
             return "polar"
     def countPi(self):
         sum=0
         bondList=self.getBondList()
         for bond in bondList:
             if bond.type=="=":
                 sum+=1
             elif bond.type=="≡":
                 sum+=2
         return sum
     def countSigma(self):
         return len(self.getBondList())
     
     def updateAllSurroundingSets(self):
         for atom in self.atoms:
             atom.updateSurroundingSet()


         
         
         
##########ALL THE FOLLOWING FUNCTIONS ARE DEDICATED TO ASSIGNING ATOMS POSITIONS
#ON THE CANVAS

     def positionsToStr(self):
         #returns a string representing the positions of all the atoms
         returnStr=""
         self.rearrange()
         atomList=self.sortAtoms()
         for atom in atomList:
             returnStr+=(atom.atomToStr()+"|"+str(atom.centerX)+","+
                         str(atom.centerY))
         return returnStr
     
     def findLowAtom(self):
         lowest=float("inf")
         for atom in self.atoms:
             if atom.centerY<lowest:
                 lowest=atom.centerY
                 lowestAtom=atom
         return lowestAtom
     
     def findTopAtom(self):
        """Returns the atom with the smallest Y coordinate (top-most)"""
        topAtom = self.atoms[0]
        for atom in self.atoms:
            if atom.centerY < topAtom.centerY:
                topAtom = atom
        return topAtom

     def findLeftAtom(self):
        """Returns the atom with the smallest X coordinate (left-most)"""
        leftAtom = self.atoms[0]
        for atom in self.atoms:
            if atom.centerX < leftAtom.centerX:
                leftAtom = atom
        return leftAtom

     def findRightAtom(self):
        """Returns the atom with the largest X coordinate (right-most)"""
        rightAtom = self.atoms[0]
        for atom in self.atoms:
            if atom.centerX > rightAtom.centerX:
                rightAtom = atom
        return rightAtom
     
     def findCenter(self):
        """
        Returns the approximate center of the molecule using extremal atoms:
        - Top-most, left-most, and right-most atoms
        """
        topAtom = self.findTopAtom()
        leftAtom = self.findLeftAtom()
        rightAtom = self.findRightAtom()
        lowestAtom = self.findLowAtom()

        # Center X is midpoint between left-most and right-most atoms
        centerX = (leftAtom.centerX + rightAtom.centerX) / 2
        # Center Y is midpoint between top-most and lowest atom
        centerY = (topAtom.centerY + lowestAtom.centerY) / 2

        return (centerX, centerY)
            
     
     def findCenterAtom(self):
         """finds atom closest to the center"""
         centerX,centerY=self.findCenter()
         minDist=float('inf')
         centerAtom=self.atoms[0]
         for atom in self.atoms:
             distance=getDistance(centerX,centerY,atom.centerX,atom.centerY)
             if distance<minDist:
                 minDist=distance
                 centerAtom=atom
         return centerAtom
     
     
     def center(self):
         """centers molecule to center of canvas"""
         centerX,centerY=1200/2,600/2
         centerAtom=self.findCenterAtom()
         centerAtomX,centerAtomY=centerAtom.centerX,centerAtom.centerY
         xDist=centerX-centerAtomX
         yDist=centerY-centerAtomY
         for atom in self.atoms:
             atom.centerX+=xDist
            
         return True



         
                 
             

     def getPositionless(self):
         allPositionless=[]
         for atom in self.atoms:
             if not (atom.centerX and atom.centerY):
                 allPositionless.append(atom)
         return allPositionless
         
     def assignHorizontal(self,width,height):
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
                orderList=[(-75,0),(75,0)]
                for dX,dY in orderList:

                    if not self.hasAtomThere(predecessor.centerX+dX,
                                             predecessor.centerY+dY):
                        current.centerX=predecessor.centerX+dX
                        current.centerY=predecessor.centerY+dY
                        break
            #makes it so atoms with more bonds are positioned horizontally
            current.electronDomains.sort(key=lambda bond: len(
                self.split(current,bond)) if 
                isinstance(bond,Bond) else 0, reverse=True)
            for domain in current.electronDomains:
                if isinstance(domain,Bond) and domain not in bondSet:
                    queue.append((current,domain.getOther(current)))
                    bondSet.add(domain)
            
         
     def assignPositions(self,width,height):#positionsList):
        """assign's the best positions to an atom for drawing"""
        self.updateAllSurroundingSets()
        self.assignHorizontal(width,height)
        
        atomsToAssign=self.getPositionless()
        print(len(atomsToAssign))
        for atom in atomsToAssign:
            assignedPositions=False
            for atom1 in atom.surroundingSet:
                #if an atom around it has positions
                if atom1.centerX and atom1.centerY:
                    for dX,dY in [(0,-75),(0,75),(75,0),(-75,0)]:
                        if self.hasAtomThere(atom1.centerX+dX,atom1.centerY+dY):
                            continue
                        assignedPositions=True
                        atom.centerX=atom1.centerX+dX
                        atom.centerY=atom1.centerY+dY
                        break
            if not assignedPositions and atom.symbol!="H":
                atomsToAssign.append(atom)
        self.center()

        return True
     
     def getLonePairs(self):
         """assigns positions to lone pairs"""
         positions={}
         for atom in self.atoms:
             numLonePairs=atom.countDomains(":")
             print(numLonePairs)
             if True:
                addedLPs=0
                for dX,dY in [(-75,0),(75,0),(0,-75),(0,75)]:
                    alreadyThere=self.hasAtomThere(atom.centerX+dX,
                                             atom.centerY+dY)
                    if not alreadyThere or not atom.isBonded(alreadyThere):
                        if addedLPs>=numLonePairs:
                            break
                        positions[str(atom.centerX+dX/3)+","+
str(atom.centerY+dY/3)]=":-" if ((dX,dY)==(-75,0) or (dX,dY)==(75,0)) else ":|"
                        addedLPs+=1
         return positions
     
     def assignBonds(self):
        bonds={}
        for bond in self.getBondList():
            atomOne=bond.atomOne
            atomTwo=bond.atomTwo

            #type of the bond and - if it is horizontal or | if it is vertical
            typeAndOrientation=bond.type+("|" if 
                                atomOne.centerX==atomTwo.centerX else "-")
            
            bonds[str((atomOne.centerX+atomTwo.centerX)//2)+","
                +str((atomOne.centerY+atomTwo.centerY)//2)]=typeAndOrientation
        return bonds
     
     def hasAtomThere(self,xPos,yPos):
         """returns true if an atom overlaps another atom"""
         for atom in self.atoms:
             if getDistance(atom.centerX,atom.centerY,xPos,yPos)<=10:
                 return atom
         return False
     
     def allHasPositions(self):
         """checks if every atom is assigned to a centerX and centerY"""
         for atom in self.atoms:
             if not (atom.centerX and atom.centerY):
                 return False 
         return True
             
         
     
     



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
                if domain==":":
                    continue
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
     
     def maxBonds(self):
         """returns the number of bonds the atom with the most amount of bonds
         has"""
         maxBonds=-1*float('inf')
         for atom in self.atoms:
             if atom.countDomains("bond")>maxBonds:
                 maxBonds=atom.countDomains("bond")
         return maxBonds
             


     def split(self,atom,bond):
         """returns a 'molecule' reprsenting the frontier of the atom with 
         respect to the given bond. Used for determing polarity"""

        #in the direction of that bond
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
     
     def splitOtherDirection(self,atom,atom2):
         """gets frontier of molecule in the opposite direction of the atome"""
         if not atom2:
             return Molecule([],0,"")
         bondSet=set()
         queue=deque()
         for bond in atom.electronDomains:
             if bond==":":
                 continue
             if not bond.getOther(atom) is atom2:
                 bondSet.add(bond)
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
         return newMolecule
    
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
     def getRing(self):
         """gets all rings the molecule"""

         #it's not a ring
         if not self.isBondedOrCircular()[1]:
             return False
         



         
         
         
##########ALL THE FOLLOWING FUNCTIONS ARE DEDICATED TO ASSIGNING ATOMS POSITIONS
#ON THE CANVAS
     def getAverageDistance(self,x,y):
         distSum=0
         numAtoms=0
         for atom in self.atoms:
             if not (atom.centerX and atom.centerY):
                 continue
             distSum+=getDistance(x,y,atom.centerX,atom.centerY)
             numAtoms+=1
         return distSum/numAtoms
     
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
     
     
     def center(self,width,height):
         """centers molecule to center of canvas"""
         centerX,centerY=width/2,height/2
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
     
     def assignRing(self,width,height):
         
         pass
     
     def updatePredecessors(self):
         for atom in self.atoms:
             if atom.predecessor:
                 continue
             adjacentAtom=atom.getAdjacentAtom()
             if adjacentAtom:
                 atom.predecessor=adjacentAtom
         return True
              
            
     def recursiveAssignPositions(self, positionless=None,bondList=None):
        """recursively backtracks where atoms go on the canvas"""
        if positionless is None and bondList is None:
            bondList=self.getBondList()
            positionless = self.getPositionless()

        #returns true if everything has a position
        if not positionless: 
            return True
        
        #we do this because an atoms position is determined by the position 
        #of an atom next to it
        self.updatePredecessors()
        positionless.sort(key=lambda atom: (len(self.splitOtherDirection(atom,
                                                atom.predecessor).atoms)),
                                                reverse=True)
        atom=positionless.pop(0)
        predecessor=atom.predecessor
        if not predecessor:
            positionless.append(atom)
            if self.recursiveAssignPositions(positionless):
                    return True
            else:
                return False
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        #coordinate which results in optimal placement
        lengths=[65,130]
        for length in lengths:
            if length==130 and not (atom.countDomains("bond")>=3 or
                                     predecessor.countDomains("bond")>=3):
                continue
            if length==130 and atom.symbol=="H" or predecessor.symbol=="H":
                continue

            directions.sort(key=lambda direction: self.getAverageDistance
                            (predecessor.centerX+length*direction[0],
                                predecessor.centerY+length*direction[1]),
                                reverse=True)
            
            for direction in directions:
                dx = length * direction[0]
                dy = length * direction[1]
                newX = predecessor.centerX + dx
                newY = predecessor.centerY + dy
                newBond=atom.getBond(predecessor)
                if self.hasAtomThere(newX, newY):
                    continue
                atom.centerX=newX
                atom.centerY=newY
                shouldContinue=False
                if length==130:
                    newBond.isLong=True
                    for bond in bondList:
                        if (not bond.sameAtoms(newBond) and
                            bond.intersects(newBond)):
                            shouldContinue=True
                            break
                        if bond.intersectsAtom(atom):
                            shouldContinue=True
                            break
                for atom_ in self.atoms:
                    if shouldContinue:
                        break
                    if newBond.intersectsAtom(atom_):
                        shouldContinue=True
                        break
                    for bond_ in bondList:
                        if bond_.intersectsAtom(atom):
                            shouldContinue=True
                            break
                if not shouldContinue:
                    if self.recursiveAssignPositions(positionless,bondList):
                        return True 
                atom.centerX=0
                atom.centerY=0
                if length==130:
                    newBond.isLong=False
                
        #no valid placement found, put back for future attempts
        positionless.append(atom)
        return False
     




     def assignPositions(self,width,height):#positionsList):
        """assign's the best positions to an atom for drawing"""
        self.updateAllSurroundingSets()
        self.atoms[0].centerX=500
        self.atoms[0].centerY=300
        #self.assignHorizontal(width,height)
        self.recursiveAssignPositions()
        self.center(width,height)
        return
        
     def extendBonds(self):
         """makes longer bonds if the space isn't enough"""


     
     def getLonePairs(self):
         """assigns positions to lone pairs"""
         positions={}
         for atom in self.atoms:
             numLonePairs=atom.countDomains(":")
             if numLonePairs:
                addedLPs=0
                for dX,dY in [(-75,0),(75,0),(0,-75),(0,75)]:
                    alreadyThere=self.hasAtomThere(atom.centerX+dX,
                                             atom.centerY+dY)
                    if not alreadyThere or (not alreadyThere in 
                                            atom.surroundingSet):
                        if addedLPs>=numLonePairs:
                            break
                        positions[str(atom.centerX+dX/4)+","+
str(atom.centerY+dY/4)]=":-" if ((dX,dY)==(-75,0) or (dX,dY)==(75,0)) else ":|"
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
                +str((atomOne.centerY+atomTwo.centerY)//2)]=(typeAndOrientation+
                                                ("l" if bond.isLong else ""))
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
     


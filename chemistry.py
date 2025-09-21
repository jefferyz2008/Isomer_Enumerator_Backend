"""atom and bond classes and methods"""
from helperFunctions import*
from collections import deque
class Atom:
    
    __slots__ = [
        "symbol", "electronDomains", "mol", "molarMass", "valenceElectrons",
        "prefBonds", "atomicNumber", "canExpandOctet", "electroNegativity",
        "currentElectrons", "octetElectrons",
        "centerX", "centerY","strAtom","surroundingSet","alreadyChecked",
        "predecessor"
    ]
    def __init__(self,symbol):
        self.symbol=symbol

        #electron domains contain bonds and lone pairs
        self.electronDomains =[]

        #only these atoms are used in lewis structures
        allElements = {
    "H":  [1.01, 1, 1, 1, False, 2.20],
    "C":  [12.01, 4, 4, 6, False, 2.55],
    "N":  [14.01, 5, 3, 7, False, 3.04],
    "O":  [16.00, 6, 2, 8, False, 3.44],
    "S":  [32.07, 6, 2, 16, True, 2.58],
    "F":  [19.00, 7, 1, 9, False, 3.98],
    "Cl": [34.45, 7, 1, 17, True, 3.16],
    "Br": [79.70, 7, 1, 35, True, 2.96],
    "P":  [30.97, 5, 3, 15, True, 2.19],
    "I":  [126.90, 7, 1, 53, True, 2.66],
    "He": [4.00, 8, 4, 2, False, 0.0],  
    "Ne": [20.18, 8, 4, 10, False, 0.0],
    "Ar": [39.95, 8, 4, 18, True, 0.0],
    "Kr": [83.80, 8, 4, 36, True, 3.0],  # approximate
    "Xe": [131.29, 8, 4, 54, True, 2.6], # approximate
    "B":  [10.81, 3, 3, 5, False, 2.04],
    "Be": [9.01, 2, 2, 4, False, 1.57],
    "Si": [28.09, 4, 4, 14, True, 1.90],
    "Al": [26.98, 3, 3, 13, True, 1.61]
        }
               
        element=allElements[symbol]
        self.molarMass=element[0]
        self.valenceElectrons=element[1]
        self.prefBonds=element[2]
        self.atomicNumber=element[3]
        self.canExpandOctet=False and element[4]
        self.electroNegativity=element[5]
        
        #number of electrons attached to self
        self.currentElectrons=0

        #indicates number of electrons needed to complete an octet
        self.octetElectrons=(8 if symbol not in {"B","Al","H","Be"} else 
        (6 if symbol in {"B","Al"} else (4 if symbol=="Be" else 2)))

        self.strAtom=symbol+"|"

        self.surroundingSet=set()
        
        #only for indicating coordinates for drawing
        self.centerX=0
        self.centerY=0

        self.alreadyChecked=False
        

        #atom that this atom uses to get a position
        self.predecessor=None



    def getElectrons(self):
        """number of electrons attached to an atom"""
        electronSum=0
        for domain in self.electronDomains:
            if isinstance(domain,Bond):
                electronSum+=domain.electrons
            elif domain==":":
                electronSum+=2
        return electronSum

 
    def countDomains(self,item):
        """count how many of a given domain is attached to an atom"""
        count=0
        for domain in self.electronDomains:
            isBond=isinstance(domain,Bond)
            if isinstance(domain,str) and domain==":" and item==":":
                count+=1
                continue
            elif isBond and domain.type=="-" and item=="-":
                count+=1
                continue
            elif isBond and domain.type=="=" and item=="=":
                count+=1
                continue
            elif isBond and domain.type=="≡" and item=="≡":
                count+=1
                continue
            elif isBond and item=="bond":
                count+=1
                continue
        return count
    
    def countBonds(self):
        """counts number of bonds attached to an atom"""
        total=0
        for domain in self.electronDomains:
            if isinstance(domain,Bond):
                total+=domain.electrons/2
        return total
    
    def getMolarMassSum(self):
        """gets all the atoms connected to a given atom, and multiplies 
        their molar mass by the order of the bond between them"""
    #used for rearraning atoms, to get canonical strings
        return sum([bond.getOther(self).molarMass*bond.electrons/2
    +bond.electrons for bond in self.electronDomains if isinstance(bond,Bond)])
    

    def hasOctet(self):
        """returns true if an atom has an octet, False if its over an octet
        and None if the octet is incomplete"""
        electronSum=self.currentElectrons
        if self.canExpandOctet:
            if electronSum<8:
                return None
            else:
                return len(self.electronDomains)<=6 and electronSum<=16
        #all edge cases
        #hydrogen can only have 1 single bond
        if self.symbol=="H":
            if not self.electronDomains:
                return None
            if (self.countDomains("-")>1 or self.countDomains(":") or 
   self.countDomains("=") or self.countDomains("≡")):
                return False
            if self.countDomains("-")==1:
                return True
        #Be only takes 2 single bonds
        if self.symbol=="Be":
            if (self.countDomains(":") or self.countDomains("=") or 
                self.countDomains("≡") or self.countDomains("-")>2):
                return False
            if self.countDomains("-")==2:
                return True
            return None
        #Boron and Al each take 3 single bonds
        if self.symbol=="B" or self.symbol=="Al":
            if (self.countDomains(":") or self.countDomains("=") or 
                self.countDomains("≡") or self.countDomains("-")>3):
                return False
            if self.countDomains("-")==3:
                return True
            return None
        if electronSum==8:
            return True
        if electronSum<8:
            return None

        return False
    
    def addLonePair(self):
        """adds lone pair to atom"""
        self.currentElectrons+=2
        self.electronDomains.append(":")
        return True
    
    def removeLonePair(self):
        """removes lone pair from atom"""
        #checks if atom has a lone pair
        try:
            self.electronDomains.remove(":")
            self.currentElectrons-=2
            return True
        except:
            return False
       
    def sigmaBond(self,other):
        """single bonds self and other"""
        if self is other:
            #no bonding to same atoms
            return None
        """H can only have 1 single bond, so it should have nothing attached
        if it wants to bond"""
        selfSymbol=self.symbol
        otherSymbol=other.symbol
        if (selfSymbol=="H" and self.electronDomains) or (otherSymbol=="H" 
                                                and other.electronDomains):
            return None
        newBond=Bond("-",self,other)
        surroundingSet=self.surroundingSet
        otherSurroundingSet=other.surroundingSet
        if self in otherSurroundingSet or other in surroundingSet:
            return None
        
        self.electronDomains.append(newBond)
        other.electronDomains.append(newBond)
        other.currentElectrons+=2
        self.currentElectrons+=2
        #updates surroundings
        surroundingSet.add(other)
        otherSurroundingSet.add(self)
        return newBond
    
    def unBond(self,other,bond):
        """removes bond from both atoms"""
        self.currentElectrons-=2
        other.currentElectrons-=2
        self.electronDomains.remove(bond)
        other.electronDomains.remove(bond)
        self.surroundingSet.discard(other)
        other.surroundingSet.discard(self)
        return True
        
    def getFormalCharge(self):
        """returns formal charge of an atom"""
        valence=self.valenceElectrons
        for domain in self.electronDomains:
            valence-=(2 if domain==":" else domain.electrons/2)
        #noble gasses don't get formal charge
        if self.valenceElectrons==8 and valence:
            return 1000000
        return valence
    
    def atomToStr(self,lp=False,polarityCheck=False):
        """turns atom into a string for duplicate and polarity check"""
        self.rearrange(lp)
        returnStr=self.symbol+"|"
        for domain in self.electronDomains:
            if isinstance(domain,Bond):
                #for checking polarity, we don't care about bond type
                #because of resonance
                if polarityCheck:
                    returnStr+="-"+domain.getOther(self).symbol
                else:
                    returnStr+=domain.type+domain.getOther(self).symbol
        if lp:
            numLps=self.countDomains(":")
            for _ in range(numLps):
                returnStr+=":"
        
        return "("+returnStr+")"

                 
    def rearrange(self,lPs=False):
        """rearrages atom's bonds and lp's for hashing"""
        bonds = [b for b in self.electronDomains if not isinstance(b, str)]
        bonds.sort(key=lambda bond: ((bond.electrons,
                bond.getOther(self).atomicNumber)),reverse=True)
        newDomains = bonds
        if lPs: 
            newDomains+=[":"]*self.countDomains(":")
        self.electronDomains=newDomains
    
    def addUntilOctet(self):
        """adds lone pairs to an atom until it has octet. returns 
        how many lps it added"""
       #these atoms dont accept lps
        if self.symbol in {"Be","B","H","Al"}:
            return 0
        for numPairs in range((8-self.currentElectrons)//2):
            self.addLonePair()
        try:
           return numPairs*2
        except:
            return 0
        
    def isBonded(self,atom):
        """returns true if the 2 atoms are bonded"""
        for domain in self.electronDomains:
            if not isinstance(domain,Bond):
                continue
            if ((domain.atomOne is self and domain.atomTwo is atom) or 
                (domain.atomOne is atom and domain.atomTwo is self)):
                return True
        return False
    
    def updateSurroundingSet(self):
        """updates surroundings of all atoms"""
        returnSet=set()
        for bond in self.electronDomains:
            if not isinstance(bond,Bond):
                continue
            returnSet.add(bond.getOther(self))
        self.surroundingSet=returnSet
    
    def getAdjacentAtom(self):
        """adjacent atom with a position"""
        for domain in self.electronDomains:
            if (isinstance(domain,Bond) and domain.getOther(self).centerX and 
                domain.getOther(self).centerY):
                return domain.getOther(self)
        return None
    def nearestAssignedAtom(self):
        """nearest atom with a position"""
        queue=deque()
        queue.append(self)
        bondSet=set()
        while queue:
            current=queue.popleft()
            if current.centerX and current.centerY:
                return current
            for domain in current.electronDomains:
                if domain in bondSet or domain==":":
                    continue
                bondSet.add(domain)   
                queue.append(domain.getOther(current))
        return 
    

    
    def averageDistance(self,atomList):
        """returns distance between atom and all atoms in the list"""
        distSum=0
        numAtoms=0
        for atom in atomList:
            if not (atom.centerX and atom.centerY):
                continue
            distSum+=getDistance(self.centerX,self.centerY,
                                 atom.centerX,atom.centerY)
            numAtoms+=1
        return distSum/numAtoms if numAtoms else 0
    
    def isAssigned(self):
        """does it have a position"""
        return self.centerX and self.centerY
    
    def getBond(self,other):
        for domain in self.electronDomains:
            if domain!=":" and domain.getOther(self) is other:
                return domain
        return None
    





        
    
    #gets information like bond angles
    def getHybirdization(self):
        """returns hybirdization of self"""
        if len(self.electronDomains)==1:
            return "N/A"
        return "s"+"p"+str(len(self.electronDomains)-1)
    
    def getVSEPR(self):
        """returns vsepr geometry followed by bond angles """
        numBonds=self.countDomains("bond")
        if numBonds==1:
            return ["N/A","N/A"]
        numDomains=len(self.electronDomains)
        numLPs=self.countDomains(":")
        if numDomains==2:
            return ["linear","180°"]
        if numDomains==3:
            if numLPs:
                return ["bent","<120°"]
            else:
                return ["trigonal planar","120°"]
        if numDomains==4:
            if numLPs==2:
                return ["bent","<109°"]
            if numLPs==1:
                return ["trigonal pyrimidal","<109°"]
            else:
                return ["tetrahedral","109°"]
        if numDomains==5:
            if numLPs==3:
                return ["linear","180°"]
            if numLPs==2:
                return ["T-Shaped","<90°"]
            if numLPs==1:
                return ["SeeSaw","<120, <90°"]
            else:
                return ["trigonal bipyramidal","120°, 90°"]
        if numDomains==6:
            if numLPs==4:
                return ["linear","180°"]
            if numLPs==3:
                return ["T-Shaped","<90°"]
            if numLPs==2:
                return ["Square Planar","90°"]
            if numLPs==1:
                return ["square pyramidal","<90°"]
            else:
                return ["octahedral","90°"]
        return ["N/A","N/A"]
    
    
    def getRingDict(self):
         """gets a ring of atoms, starting with this atom"""
         predDict={}
         bondSet=set()
         queue=deque()
         visited=set()
         queue.append(self)
         while queue:
             current=queue.popleft()
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
                        predDict[atomOne]=current
                     else: 
                        queue.append(atomTwo)
                        predDict[atomTwo]=current
         return predDict
    

    def __repr__(self):
        return self.symbol
    def __hash__(self):
        return hash(id(self))
    def __str__(self):
        return self.symbol
    def __eq__(self,other):
        return self is other

class Bond:
     __slots__ = ["type", "electrons", "atomOne", "atomTwo", "numModifications"
                  ,"startX","endX","startY","endY","isLong"]
     def __init__(self,type,atomOne,atomTwo):
          self.type=type
          self.electrons=2 if type=="-" else (4 if type=="=" else 6)
          self.atomOne=atomOne
          self.atomTwo=atomTwo

          #for backtracking efficiency
          self.numModifications=0

          #for bonds coordinates
          self.startX=0
          self.startY=0
          self.endX=0
          self.endY=0

          self.isLong=False

     def sameAtoms(self,other):
         """Returns true if the atoms attached to each bond are the same atom"""
         selfAtomOne=self.atomOne
         selfAtomTwo=self.atomTwo
         otherAtomOne=other.atomOne
         otherAtomTwo=other.atomTwo
         return ((selfAtomOne is otherAtomOne 
            and selfAtomTwo is otherAtomTwo)  or (otherAtomTwo 
is selfAtomOne and selfAtomTwo is otherAtomOne))
     

     def sameTypeAtoms(self,other):
         """returns true if the atoms of both bonds are the same element"""
         selfAtomOne=self.atomOne.symbol
         selfAtomTwo=self.atomTwo.symbol
         return ((other[0]==selfAtomOne 
            and selfAtomTwo==other[2])  
            or (other[2]==selfAtomOne and 
selfAtomTwo==other[0]))
         
     
     def addPi(self):
         type=self.type
         if type=="≡":
             return False
         success=False
         atomOne=self.atomOne
         atomTwo=self.atomTwo
         self.numModifications+=1
         if type=="-":
             atomOne.currentElectrons+=2
             atomTwo.currentElectrons+=2
             self.type="="
             self.electrons=4
             success=True
         elif type=="=":
             atomOne.currentElectrons+=2
             atomTwo.currentElectrons+=2
             self.type="≡"
             self.electrons=6
             success=True
        
         return success
     
     def removePi(self):
         if self.type=="-":
             return False
         atomOne=self.atomOne
         atomTwo=self.atomTwo
         self.numModifications-=1
         if self.type=="≡":
             atomOne.currentElectrons-=2
             atomTwo.currentElectrons-=2
             self.type="="
             self.electrons=4
             return True
         elif self.type=="=":
             atomOne.currentElectrons-=2
             atomTwo.currentElectrons-=2
             self.type="-"
             self.electrons=2
             return True
         return False
     
     def getOther(self,atom):
         """returns the atom that is not the atom"""
         atomOne_=self.atomOne
         atomTwo_=self.atomTwo
         #gets the atom in the bond thats not the given atom
         if atom!=atomOne_ and atom!=atomTwo_:
             return False
         return atomOne_ if atom is atomTwo_ else atomTwo_
     
     def getFrontier(self,atom):
         """returns all atoms in the molecule accessible by the bond on the 
         side of the atom"""
         queue=deque()
         queue.append(atom)
         bondSet=set()
         bondSet.add(self)
         visited=set()
         while queue:
             current=queue.pop()
             visited.add(current)
             for domain in current.electronDomains:
                 if not isinstance(domain,Bond) or domain in bondSet:
                     continue
                 queue.append(domain.getOther(current))
         return visited
     
     
     def intersects(self,other):
         """returns true if bond intersects other bond"""
         atomOne=self.atomOne
         atomTwo=self.atomTwo
         otherAtomOne=other.atomOne
         otherAtomTwo=other.atomTwo
         if (atomOne is otherAtomOne or atomOne is otherAtomTwo 
             or atomTwo is otherAtomOne or atomTwo is otherAtomTwo):
             return False
         if (not atomOne.centerX and atomOne.centerY and atomTwo.centerX 
             and atomTwo.centerY):
             return False
         intersection=lineIntersection((atomOne.centerX,atomOne.centerY)
                        ,(atomTwo.centerX,atomTwo.centerY),
                        (otherAtomOne.centerX,otherAtomOne.centerY),
                        (otherAtomTwo.centerX,otherAtomTwo.centerY))
         if not intersection:
             return False
         if (intersection==(atomOne.centerX,atomOne.centerY) or 
             intersection==(atomTwo.centerX,atomTwo.centerY)):
             return False 
         return True
     
     def intersectsAtom(self,atom):
         atomOne=self.atomOne
         atomTwo=self.atomTwo
         if atom is atomOne or atom is atomTwo:
             return False
         centerX=atom.centerX
         centerY=atom.centerY
         return pointOnSegment((atomOne.centerX,atomOne.centerY),
                               (atomTwo.centerX,atomTwo.centerY),
                               (centerX,centerY))
         
         
         
     def __str__(self):
         return (self.atomOne.symbol+
                 self.type+self.atomTwo.symbol)
       
     def __repr__(self):
         return (self.atomOne.symbol+self.type+self.atomTwo.symbol)
     
     def __eq__(self,other):
         return other is self
     
     def __hash__(self):
        return hash(id(self))



     



     





     

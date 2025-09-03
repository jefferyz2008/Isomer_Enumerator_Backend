"""atom and bond classes and methods"""
class Atom:
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
        self.canExpandOctet=element[4]
        self.electroNegativity=element[5]
        
        #number of electrons attached to self
        self.currentElectrons=0

        #indicates number of electrons needed to complete an octet
        self.octetElectrons=(8 if symbol not in {"B","Al","H","Be"} else 
        (6 if symbol in {"B","Al"} else (4 if symbol=="Be" else 2)))

        #For backtracking efficiency
        self.numModifications=0

        #for formal charge
        self.formalCharge=0
        
        #only for indicating coordinates for drawing
        self.centerX=0
        self.centerY=0


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
        self.numModifications+=1
        self.currentElectrons+=2
        self.electronDomains.append(":")
        return True
    
    def removeLonePair(self):
        #checks if atom has a lone pair
        self.numModifications-=1
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
        if self.symbol=="H" and self.electronDomains:
            return None
        newBond=Bond("-",self,other)
        for domain in self.electronDomains:
           #two atoms can only have 1 bond between them
            if isinstance(domain,Bond) and domain.sameAtoms(newBond):
                return None
        self.electronDomains.append(newBond)
        other.electronDomains.append(newBond)
        other.currentElectrons+=2
        self.currentElectrons+=2
        return newBond
    
    def unBond(self,other,bond):
        self.currentElectrons-=2
        other.currentElectrons-=2
        self.electronDomains.remove(bond)
        other.electronDomains.remove(bond)
        return True
        
    def getFormalCharge(self):
        """returns formal charge of an atom"""
        valence=self.valenceElectrons
        for domain in self.electronDomains:
            valence-=(2 if domain==":" else domain.electrons/2)
        return valence
    
    def atomToStr(self,polarityCheck=False):
        returnStr=self.symbol+"|"
        for domain in self.electronDomains:
            if isinstance(domain,Bond):
                #for checking polarity, we don't care about bond type
                #because of resonance
                if polarityCheck:
                    returnStr+="-"+domain.getOther(self).symbol
                else:
                    returnStr+=domain.type+domain.getOther(self).symbol
            elif not polarityCheck:
                #for checking polarity, we don't care about lone pairs
                #because of resonance
                returnStr+=":"
        return returnStr

    
    def rearrange(self):
        """rearrages atom's bonds and lp's for hashing"""
        bonds = [b for b in self.electronDomains if not isinstance(b, str)]
        bonds.sort(key=lambda bond: ((bond.electrons,
                bond.getOther(self).atomicNumber)),reverse=True)
        lone_pairs = [lp for lp in self.electronDomains if isinstance(lp, str)]
        newDomains = bonds + lone_pairs
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
        
    
    #gets information like bond angles
    def getHybirdization(self):
        """returns hybirdization of self"""
        if len(self.electronDomains)==1:
            return "None"
        return "s"+"p"+str(len(self.electronDomains)-1)
    
    def getVSEPR(self):
        """returns vsepr geometry followed by bond angles """
        numDomains=len(self.electronDomains)
        numLPs=self.countDomains(":")
        if numDomains==2:
            return ["linear","180"]
        if numDomains==3:
            if numLPs:
                return ["bent","<120"]
            else:
                return ["trigonal planar","120"]
        if numDomains==4:
            if numLPs==2:
                return ["bent","<109"]
            if numLPs==1:
                return ["trigonal pyrimidal","<109"]
            else:
                return ["tetrahedral","109"]
        if numDomains==5:
            if numLPs==3:
                return ["linear","180"]
            if numLPs==2:
                return ["T-Shaped","<90"]
            if numLPs==1:
                return ["SeeSaw","<120, <90"]
            else:
                return ["trigonal bipyramidal","120, 90"]
        if numDomains==6:
            if numLPs==4:
                return ["linear","180"]
            if numLPs==3:
                return ["T-Shaped","<90"]
            if numLPs==2:
                return ["Square Planar","90"]
            if numLPs==1:
                return ["square pyramidal","<90"]
            else:
                return ["octahedral","90"]
        return ["None","None"]
    def __repr__(self):
        return self.symbol
    def __hash__(self):
        return hash(id(self))
    def __str__(self):
        return self.symbol
    def __eq__(self,other):
        return self is other

class Bond:
     def __init__(self,type,atomOne,atomTwo):
          self.type=type
          self.electrons=2 if type=="-" else (4 if type=="=" else 6)
          self.atomOne=atomOne
          self.atomTwo=atomTwo

          #for backtracking efficiency
          self.numModifications=0
     
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
         if self.type=="≡":
             return False
         atomOne=self.atomOne
         atomTwo=self.atomTwo
         self.numModifications+=1
         if self.type=="-":
             atomOne.currentElectrons+=2
             atomTwo.currentElectrons+=2
             self.type="="
             self.electrons=4
             return True
         if self.type=="=":
             atomOne.currentElectrons+=2
             atomTwo.currentElectrons+=2
             self.type="≡"
             self.electrons=6
             return True
         return False
     
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
         atomOne_=self.atomOne
         atomTwo_=self.atomTwo
         #gets the atom in the bond thats not the given atom
         if atom!=atomOne_ and atom!=atomTwo_:
             return False
         return atomOne_ if atom is atomTwo_ else atomTwo_
     

     def __str__(self):
         return (self.atomOne.symbol+
                 self.type+self.atomTwo.symbol)
       
     def __repr__(self):
         return (self.atomOne.symbol+self.type+self.atomTwo.symbol)
     
     def __eq__(self,other):
         return other is self
     
     def __hash__(self):
        return hash(id(self))





     





     

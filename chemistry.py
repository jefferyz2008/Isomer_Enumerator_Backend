import copy
class Atom:
    def __init__(self,symbol):
        self.symbol=symbol

        #electron domains contain bonds and lone pairs
        self.electronDomains =[]


        #only these atoms are used in lewis structures
        allElements = {
    "H":  [1.01, 1, 1, 1, False],
    "C":  [12.01, 4, 4, 6, False],
    "N":  [14.01, 5, 3, 7, False],
    "O":  [16.00, 6, 2, 8, False],
    "S":  [32.07, 6, 2, 16, True],
    "F":  [19.00, 7, 1, 9, False],
    "Cl": [34.45, 7, 1, 17, True],
    "Br": [79.70, 7, 1, 35, True],
    "P":  [30.97, 5, 3, 15, True],
    "I":  [126.90, 7, 1, 53, True],
    "He": [4.00, 8, 4, 2, False],
    "Ne": [20.18, 8, 4, 10, False],
    "Ar": [39.95, 8, 4, 18, True], 
    "Kr": [83.80, 8, 4, 36, True],
    "Xe": [131.29, 8, 4, 54, True],
    "B":  [10.81, 3, 3, 5, False],
    "Be": [9.01, 2, 2, 4, False],
    "Si": [28.09, 2, 4, 14, True],
    "Al": [26.98, 3, 3, 13, True]
}
        self.molarMass=allElements[symbol][0]
        self.valenceElectrons=allElements[symbol][1]
        self.prefBonds=allElements[symbol][2]
        self.atomicNumber=allElements[symbol][3]
        self.canExpandOctet=allElements[symbol][4]
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
            if isinstance(domain,str) and domain==":" and item==":":
                count+=1
            elif isinstance(domain,Bond) and domain.type=="-" and item=="-":
                count+=1
            elif isinstance(domain,Bond) and domain.type=="=" and item=="=":
                count+=1
            elif isinstance(domain,Bond) and domain.type=="#" and item=="#":
                count+=1
            elif isinstance(domain,Bond) and item=="bond":
                count+=1
        return count
    
    def countBonds(self):
        """counts number of bonds attached to an atom"""
        total=0
        for domain in self.electronDomains:
            if isinstance(domain,Bond):
                total+=domain.electrons/2
        return total
    
    def getMolarMassSum(self):
        """gets all the atoms connected to a given atom, and multiplies their molar
        mass by the order of the bond between them"""

    #used for rearraning atoms, to hash the molecule
        return sum([bond.getOther(self).molarMass*bond.electrons/2
            for bond in self.electronDomains if isinstance(bond,Bond)])
    
    def hasOctet(self):
        """returns true if an atom has an octet, False if its over an octet
        and None if the octet is incomplete"""
        electronSum=self.getElectrons()
        #hydrogen can only have 1 single bond
        if self.symbol=="H":
            if (self.countDomains("-")>1 or self.countDomains(":") or 
   self.countDomains("=") or self.countDomains("#")):
                return False
            if self.countDomains("-")==1:
                return True
            return None
        #Be only takes 2 single bonds
        if self.symbol=="Be":
            if (self.countDomains(":") or self.countDomains("=") or 
                self.countDomains("#") or self.countDomains("-")>2):
                return False
            if self.countDomains("-")==2:
                return True
            return None
        #Boron and Al each take 3 single bonds
        if self.symbol=="B" or self.symbol=="Al":
            if (self.countDomains(":") or self.countDomains("=") or 
                self.countDomains("#") or self.countDomains("-")>3):
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
        self.electronDomains.append(":")
        return True
    
    def removeLonePair(self):
        #checks if atom has a lone pair
        if ":" in self.electronDomains:
            self.electronDomains.remove(":")
            return True
        return False
    
    def bondTo(self,other,type):
        """creates a new bond object between self and other"""
        if self is other:
            #no bonding to same atoms
            return False
        #no bonding non-atoms
        if not isinstance(other,Atom):
            return False
        """H can only have 1 single bond, so it should have nothing attached
        if it wants to bond"""
        if self.symbol=="H":
            if self.electronDomains:
                return False
        newBond=Bond(type,self,other)

        for domain in self.electronDomains:
           #two atoms can only have 1 bond between them
            if isinstance(domain,Bond) and domain.sameAtoms(newBond):
                return False

        self.electronDomains.append(newBond)
        other.electronDomains.append(newBond)
        return True
    
    def unBond(self,other):
        for domain in self.electronDomains:
            if (isinstance(domain,Bond) and
                ((domain.atomOne is self and domain.atomTwo is other) or 
                 (domain.atomOne is other and domain.atomTwo is self))):
                     self.electronDomains.remove(domain)
                     other.electronDomains.remove(domain)
                     return True
        return False
    
    def removeLonePair(self):
        for domain in self.electronDomains:
            if domain==":":
                self.electronDomains.remove(":")
                return True
        return False
    
    def getFormalCharge(self):
        #atoms with incomplete octets can lower their formal charge in the future
        #so don't penalize
        if not self.hasOctet():
            return 0
        valence=self.valenceElectrons
        for domain in self.electronDomains:
            valence-=(2 if domain==":" else domain.electrons/2)
        return valence
    
    def rearrange(self):
        """rearrages atom's bonds and lp's for hashing"""
        bonds = [b for b in self.electronDomains if not isinstance(b, str)]
        bonds.sort(key=lambda bond: ((bond.electrons,
                bond.getOther(self).atomicNumber)),reverse=True)
        lone_pairs = [lp for lp in self.electronDomains if isinstance(lp, str)]
        self.electronDomains = bonds + lone_pairs
        return True
    
    def addUntilOctet(self):
        """adds lone pairs to an atom until it has octet"""
       #these atoms dont accept lps
        if self.symbol in {"Be","B","H","Al"}:
            return False
        
        for _ in range(4-self.countDomains("-")):
            self.addLonePair()

        return True
        
        
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
          if self.type=="-":
               self.electrons=2
          elif self.type=="=":
               self.electrons=4
          elif self.type=="#":
               self.electrons=6
          self.atomOne=atomOne
          self.atomTwo=atomTwo

     
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
         if self.type=="-":
             self.type="="
             self.electrons=4
             return True
         if self.type=="=":
             self.type="#"
             self.electrons=6
             return True
         return False
     
     def removePi(self):
         if self.type=="#":
             self.type="="
             self.electrons=4
             return True
         elif self.type=="=":
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
         return (self.atomOne.symbol+
                 self.type+self.atomTwo.symbol)
     
     def __eq__(self,other):
         return other is self
     
     def __hash__(self):
        return hash(id(self))





     

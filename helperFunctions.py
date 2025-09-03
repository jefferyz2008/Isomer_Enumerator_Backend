def getDistance(x1,y1,x2,y2):
    return ((x1-x2)**2+(y1-y2)**2)**0.5
    
def getMidPoint(x1,y1,x2,y2):
      return ((x1+x2)/2,(y1+y2)/2)

def checkOverBonding(atom):
      """check if an atom has too many bonds for the element"""
      #moleculeList is there to see if theres and alternative
      if atom.currentElectrons>=8 and not atom.canExpandOctet:
            return True
      atomSymbol=atom.symbol
      if not (atomSymbol=="O" or atomSymbol=="N" or atomSymbol=="H" 
              or atom.valenceElectrons==7):
            return False
      if atomSymbol=="H" and atom.electronDomains: 
            return True
      atomSingles=atom.countDomains("-")
      if (atomSymbol=="O" and atomSingles>=2):
            return True
      if (atomSymbol=="N" and atomSingles>=3):
            return True
      if ((atom.valenceElectrons==7 and atomSingles>=1) 
            and not atom.canExpandOctet):
            return True
      return False


def badCarbon(molecule):
        """returns true if carbon only has 1 bond"""
        #used for eliminating skeletal structures
        for atom in molecule.atoms:
            if atom.symbol=="C" and atom.countDomains("bond")<2:
                  return True
        return False

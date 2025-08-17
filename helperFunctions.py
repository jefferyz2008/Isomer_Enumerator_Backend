def getDistance(x1,y1,x2,y2):
    return ((x1-x2)**2+(y1-y2)**2)**0.5
def getMidPoint(x1,y1,x2,y2):
      return ((x1+x2)/2,(y1+y2)/2)
def hasOverOxyNitro(molecule):
      """returns true if oxygen or nitrogen have over their preffered bonds"""
      for atom in molecule.atoms:
            atomSymbol=atom.symbol
            numBonds=atom.countDomains("bond")
            if atomSymbol=="O":
                  if numBonds>2:
                        return True
            if atomSymbol=="N":
                  if numBonds>3:
                        return True
      return False

def badCarbon(molecule):
        """returns true if carbon only has 1 bond"""
        for atom in molecule.atoms:
            if atom.symbol=="C" and atom.countDomains("bond")<2:
                  return True
        return False

def badHalogen(molecule):
      """returns true if halogen has more than 2 bonds"""
      for atom in molecule.atoms:
            if atom.valenceElectrons==7 and atom.countDomains("bond")>1:
                  return True
      return False


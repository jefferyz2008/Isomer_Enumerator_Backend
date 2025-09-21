def getDistance(x1,y1,x2,y2):
    return ((x1-x2)**2+(y1-y2)**2)**0.5
    
def getMidPoint(x1,y1,x2,y2):
      return ((x1+x2)/2,(y1+y2)/2)
def pointsAreClose(p1, p2, tol=0.1):
    return abs(p1[0]-p2[0]) < tol and abs(p1[1]-p2[1]) < tol

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

def lineIntersection(p1, q1, p2, q2):
    """
    Returns the intersection point of two line segments p1q1 and p2q2
    or None if they don't intersect.
    """
    x1, y1 = p1
    x2, y2 = q1
    x3, y3 = p2
    x4, y4 = q2

    # Compute denominators
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denom == 0:
        return None  # parallel or collinear

    # Intersection point of the infinite lines
    px = ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2) * (x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2) * (y3 - y4) - (y1 - y2) * (x3*y4 - y3*x4)) / denom
    intersection = (px, py)

    # Check if intersection point lies on both segments
    if (min(x1, x2) <= px <= max(x1, x2) and
        min(y1, y2) <= py <= max(y1, y2) and
        min(x3, x4) <= px <= max(x3, x4) and
        min(y3, y4) <= py <= max(y3, y4)):
        return intersection

    return None

def pointOnSegment(p1, p2, q):
    """
    Returns True if point q lies on the line segment p1p2, otherwise False.
    p1, p2, q are tuples (x, y).
    """
    x1, y1 = p1
    x2, y2 = p2
    xq, yq = q

    # Check collinearity using cross product
    cross = (y2 - y1) * (xq - x1) - (x2 - x1) * (yq - y1)
    if cross != 0:
        return False  # not collinear

    # Check if q is within the bounding rectangle of p1p2
    if (min(x1, x2) <= xq <= max(x1, x2) and
        min(y1, y2) <= yq <= max(y1, y2)):
        return True

    return False



    




def prettyPrintList(moleculeList):
    for molecule in moleculeList:
        prettyPrint(molecule)
        for _ in range(5):
           print(" ")

def prettyPrint(molecule):
    for atom in molecule.atoms:
        print(atom.symbol,":")
        domainList=[]
        for domain in atom.electronDomains:
            if domain==":":
                domainList.append(":")
            else:
                domainList.append(domain.type+domain.getOther(atom).symbol)
        print(domainList)

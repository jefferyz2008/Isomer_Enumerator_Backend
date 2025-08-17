from chemistry import*
from moleculeClass import*
from tokenizer import*
from backtracker import*
from breadthFirst import*
def getAllStructures(molecule):
    "takes all skeletal structures and finds solutions for them"
    moleculeList=getSkeletalStructuresDFS(molecule)
    #debugging purposes
    print(str(len(moleculeList))+"this is how many skeletal structures there are")
    #prettyPrintList(moleculeList)
   # moleculeList=[moleculeList[12]]
    
    validStructs=[]
    minScore=[float("inf")]
    for molecule in moleculeList:
        #the function takes in single-bonded molecules that have an octet on 
        #every atom
       # molecule.addUntilOctet()
       # for atom in molecule.atoms:
           # assert(atom.hasOctet())
        result=backtrackSkeletalStructure(molecule,[],minScore)
        if result:
           minScore=[result[1]]
           validStructs+=result[0]
   # prettyPrintList(validStructs)
    return validStructs

def getBestStructures(molecule):
    """sorts molecules based on score and returns a list of the lowest score
    molecules"""
    validStructures=getAllStructures(molecule)
    validStructures.sort(key=lambda molecule: molecule.getScore())
    allBestStructs=[molecule for molecule in validStructures if
                     molecule.getScore()==validStructures[0].getScore()]
    return allBestStructs



def solveExpandedOctet(molecule):
    centralAtom=molecule.getCentralAtom()
    for atom in molecule.atoms:
        if not atom is centralAtom:
            atom.bondTo(centralAtom,"-")
            atom.addUntilOctet()

#start_time = time.perf_counter()

molecule=tokenizeMolecule("C3H4")

prettyPrintList((getBestStructures(molecule)))
#molecule.rearrange()
#print(str(molecule.getScore())+"score")


#end_time = time.perf_counter()

# Calculate the elapsed time
#elapsed_time = end_time - start_time
#print(elapsed_time)
    


           
               
               
     
   

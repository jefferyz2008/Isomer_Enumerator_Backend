#all DFS/backtracking algorithms
from chemistry import*
from moleculeClass import*
from tokenizer import*
from backtracker import*
from breadthFirst import*
import time
def getAllStructures(molecule):
    "takes all skeletal structures and finds solutions for them"
    moleculeList=getSkeletalStructuresDFS(molecule)
    print(len(moleculeList))
    validStructs=[]
    #default min score
    minScore=[float("inf")]
    for molecule in moleculeList:
        result=backtrackSkeletalStructure(molecule,[],[],[],minScore)
        if result:
           minScore=[result[1]]
           validStructs+=result[0]
    return validStructs

def testBacktrack(molecule):
    moleculeList=getSkeletalStructuresDFS(molecule)
  #  prettyPrintList(moleculeList)
    print(len(moleculeList))
    completeList=[]
    minScore=[float('inf')]
    start=time.perf_counter()
    for mol in moleculeList:
        newMolecule=backtrackPiBondsLonePairs(mol,[],None,minScore,set())
        minScore=newMolecule[1]
        completeList+=newMolecule[0]
    end=time.perf_counter()
  #  print(start-end)
    return completeList
        
def getBestStructures(molecule):
    """sorts molecules based on score and returns a list of the lowest score
    molecules"""
    print("getting structure")
    func=getAllStructures if molecule.expandedOctet else testBacktrack
    #func=testBacktrack
    #validStructures=getAllStructures(molecule)
    validStructures2=func(molecule)
    validStructures2.sort(key=lambda molecule: molecule.getScore())
    allBestStructs=[molecule for molecule in validStructures2 if
                     molecule.getScore()==validStructures2[0].getScore()]
    return allBestStructs
   # return validStructures2
   
    return validStructures2
    return validStructures






           
               
               
     
   

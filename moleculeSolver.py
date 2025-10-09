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

def getAllStructures2(molecule):
    moleculeList=getSkeletalStructuresDFS(molecule)
    completeList=[]
    minScore=[float('inf')]
    for mol in moleculeList:
        newMolecule=backtrackPiBondsLonePairs(mol,[],None,minScore,set())
        minScore=newMolecule[1]
        completeList+=newMolecule[0]
    return completeList
        
def getBestStructures(moleculeName):
    """sorts molecules based on score and returns a list of the lowest score
    molecules"""
    print("getting structure")
    molecule=parseMolecule(moleculeName)
    func=getAllStructures if molecule.expandedOctet else getAllStructures2
    validStructures=func(molecule)
    validStructures.sort(key=lambda molecule: molecule.getScore())
    allBestStructs=[molecule for molecule in validStructures if
                     (molecule.getScore()==validStructures[0].getScore() and 
                      molecule.getScore()<=40)]
    return allBestStructs

           
               
               
     
   

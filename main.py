from fastapi import FastAPI
from pydantic import BaseModel
from moleculeClass import*
from helperFunctions import*
from fastapi.middleware.cors import CORSMiddleware
from backtracker import*
from breadthFirst import*
from moleculeSolver import*
from tokenizer import*

def getInfo(molecule):
    pass
def getMolecule(molecule):
    """gets all the positions of the atoms, bonds, and  lp's in the molecule"""
    molecule.assignPositions(750,450)
    atoms = {f"{atom.centerX},{atom.centerY}": 
             atom.symbol for atom in molecule.atoms}
    
    bonds={}
    for bond in molecule.getBondList():
        atomOne=bond.atomOne
        atomTwo=bond.atomTwo

        #type of the bond and - if it is horizontal or | if it is vertical
        typeAndOrientation=bond.type+("|" if 
                            atomOne.centerX==atomTwo.centerX else "-")
        
        bonds[str((atomOne.centerX+atomTwo.centerX)/2)+","
              +str((atomOne.centerY+atomTwo.centerY)/2)]=typeAndOrientation
    lonePairs=molecule.getLonePairs()
    return {"atoms": atoms, "bonds": bonds,"lonePairs":lonePairs}

api=FastAPI()
# allow requests from Electron (localhost)
api.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_methods=["*"],
    allow_headers=["*"],
)


@api.get("/molecule")
def getAllMolecules(moleculeName):
    result={}
    moleculeList=getBestStructures(tokenizeMolecule(moleculeName))
    for index in range(len(moleculeList)):
        result[str(index)]=getMolecule(moleculeList[index])
    print(result)
    return result



if __name__ == "__main__":
    import uvicorn
    uvicorn.run(api, host="127.0.0.1", port=8000, reload=False)

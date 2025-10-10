from fastapi import FastAPI
from pydantic import BaseModel
from moleculeClass import*
from helperFunctions import*
from fastapi.middleware.cors import CORSMiddleware
from backtracker import*
from breadthFirst import*
from moleculeSolver import*
from tokenizer import*
import asyncio
import json
import concurrent.futures
executor = concurrent.futures.ProcessPoolExecutor()



def getMoleculeInfo(molecule):
    """returns info about the molecule like polarity and molar mass"""
    result={} 
    result["molar mass"] = str(round(molecule.getMolarMass(), 2)) + "g"
    result["polarity"]=molecule.getPolarity()
    result["σ bonds"]=str(molecule.countSigma())
    result["π bonds"]=str(molecule.countPi())
    return result

def getAtomsInfo(molecule):
    """gets info of ATOMS like VSEPR ,bond angles, and Hybirdization"""
    result={}
    for atom in molecule.atoms:
        info={}
        info["Hybirdization"]= atom.getHybirdization()
        vsepr=atom.getVSEPR()
        info["VSEPR"]=vsepr[0]
        info["bond angles"]=vsepr[1]
        fc=(int(atom.getFormalCharge()))
        info["formal charge"]=str(fc) if fc<=0 else "+"+str(fc)
        result[f"{atom.centerX},{atom.centerY}"]=info
    return result



def getMolecule(molecule):
    """gets all the positions of the atoms, bonds, and  lp's in the molecule"""
    molecule.assignPositions(1000,600)
    atoms = {f"{atom.centerX},{atom.centerY}": 
             atom.symbol for atom in molecule.atoms}
    bonds=molecule.assignBonds()
    lonePairs=molecule.getLonePairs()
    atomsInfo=getAtomsInfo(molecule)
    molInfo=getMoleculeInfo(molecule)
    return {"atoms": atoms, "bonds": bonds,"lonePairs":lonePairs,
            "atomsInfo":atomsInfo,"molInfo":molInfo}

api=FastAPI()
api.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_methods=["*"],
    allow_headers=["*"],
)


@api.get("/molecule")
async def getAllMolecules(moleculeName:str):
    try:
        testMol=parseMolecule(moleculeName)
    except:
        return {}
    #moleculeStr is a unique string of the molecule
    moleculeStr=testMol.uniqueString()
    with open("cache.json", "r") as f:
        data = json.load(f)
    if moleculeStr in data:
        return data[moleculeStr]
    result={}
    try:
       moleculeList=await asyncio.to_thread(getBestStructures,moleculeName)
       print("molecule received")
    except:
        return {}
    for index in range(len(moleculeList)):
        result[str(index)]=getMolecule(moleculeList[index])
    data[moleculeStr]=result
    with open("cache.json", "w") as f:
        json.dump(data, f, indent=4)
    return result







if __name__ == "__main__":
    import uvicorn
    uvicorn.run(api, host="127.0.0.1", port=8000, reload=False)








if __name__ == "__main__":
    import uvicorn
    uvicorn.run(api, host="127.0.0.1", port=8000, reload=False)

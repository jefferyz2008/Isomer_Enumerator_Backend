from chemistry import*
from moleculeClass import Molecule
def parseAtoms(name):
    """parses all atoms in molecule"""
    try:
        assert(isinstance(name,str))
        assert(name.count("(")==name.count(")"))
        atomList=[]
        index=0
        while index<len(name):
            if name[index].islower():
                index+=1
                continue
            if name[index]=="^":
                break
            if name[index].isupper():
                if index==len(name)-1 or (not name[index+1].islower()):
                   atomList.append(Atom(name[index]))
                else:
                   atomList.append(Atom(name[index:index+2]))
            if name[index].isnumeric():
                try:
                    if name[index-1].isnumeric() and index-1>=0:
                        index+=1
                        continue
                except:
                    timesToAdd=0
                try:
                    if name[index-1]==")":
                        index+=1
                        continue
                except:
                    index+=1
                    continue
                try:
                    if name[index+1].isnumeric():
                        timesToAdd=int(name[index:index+2])-1
                    else:
                        timesToAdd=int(name[index])-1
                except:
                    timesToAdd=int(name[index])-1
                prevAtom=atomList[-1].symbol
                for _ in range(timesToAdd):
                    atomList.append(Atom(prevAtom))
            if name[index]=="(":
                for index2 in range(index+1,len(name)):
                    if name[index2]==")":
                        timesAdded=0
                        try:
                           if name[index2+1].isnumeric():
                               try:
                                   if name[index2+2].isnumeric():
                                       timesAdded=int(name[index2+1:index2+3])
                                   else:
                                       timesAdded=int(name[index2+1:index2+2])
                               except:
                                   timesAdded=int(name[index2+1:index2+2])
                           else:
                             timesAdded=1
                        except:
                            timesAdded=1
                        assert(timesAdded)
                        for _ in range(timesAdded):
                            atomList+=parseAtoms(name[index+1:index2])
                        break
                index=index2+1
                continue
            index+=1
        return atomList

    except:
        return []
    

def parseMolecule(name):
        commonNames={
    "hydrochloric acid": "HCl",
    "ammonia": "NH3",
    "water": "H2O",
    "carbon dioxide": "CO2",
    "carbon monoxide": "CO",
    "methane": "CH4",
    "ethane": "C2H6",
    "propane": "C3H8",
    "butane": "C4H10",
    "ethene": "C2H4",
    "propene": "C3H6",
    "ethyne": "C2H2",
    "formaldehyde": "CH2O",
    "methanol": "CH3OH",
    "ethanol": "C2H5OH",
    "glucose": "C6H12O6",
    "ozone": "O3",
    "hydrogen": "H2",
    "nitrogen": "N2",
    "oxygen": "O2",
    "hydrogen sulfide": "H2S",
    "nitric oxide": "NO",
    "nitrogen dioxide": "NO2",
    "sulfur dioxide": "SO2",
    "hydrogen cyanide": "HCN",
    "formic acid": "HCOOH",
    "acetic acid": "CH3COOH",
    "malic acid" :"C4H6O5",
    "maleic acid": "C4H4O4",
    "nitric acid":"HNO3",
    "hydrogen peroxide":"H2O2",
    "tartaric acid":"C4H6O6"

}
        if name in commonNames:
            atomList=parseAtoms(commonNames[name])
        else:
            atomList=parseAtoms(name)
        ###gets the charge
        if "^" in name:
            if "+"==name[-1]:
                charge=(1 if name[-2]=="^" else int(name[-2]))
            if "-"==name[-1]:
                charge=(-1 if name[-2]=="^" else -1*int(name[-2]))
            if charge>3:
                raise ValueError("charge too high")
            if charge<-4:
                raise ValueError("charge too low")
            
        else:
            charge=0
        newMolecule=Molecule(atomList,charge,name)
        if newMolecule.numElectrons%2==1:
           raise ValueError("Molecule has odd number of electrons")
        return newMolecule
    





        

        

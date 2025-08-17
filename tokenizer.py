from chemistry import Atom
from moleculeClass import Molecule
def tokenizeMolecule(name):
    #takes the name of the moleceule like NO3-
    #and turns it into a molecule object
    try:
        assert(isinstance(name,str))
        atomList=[]
        for index in range(len(name)):
            if name[index].islower():
                continue
            if name[index]=="^":
                break
            if name[index].isupper():
                if index==len(name)-1 or not name[index+1].islower():
                   atomList.append(Atom(name[index]))
                else:
                   atomList.append(Atom(name[index:index+2]))
            if name[index].isnumeric():
                prevAtom=atomList[-1].symbol
                number=int(name[index])-1
                for _ in range(number):
                    atomList.append(Atom(prevAtom))
    
        ###gets the charge
        if "^" in name:
            if "+"==name[-1]:
                charge=(1 if name[-2]=="^" else int(name[-2]))
            if "-"==name[-1]:
                charge=(-1 if name[-2]=="^" else -1*int(name[-2]))
        else:
            charge=0
        return Molecule(atomList,charge)

    except:
        print("This isnt a valid molecule")
        return False

def prettyPrintList(moleculeList):
    for molecule in moleculeList:
        prettyPrint(molecule)
        for _ in range(5):
           print(" ")

def prettyPrint(molecule):

    for atom in molecule.atoms:
        print(atom.symbol,":")
        print(atom.electronDomains)
        

        
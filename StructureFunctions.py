import FileFunctions as FF, conf
import os
from collections import defaultdict
import numpy as np
# Base Pairs from dot bracket Secondary structure

def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    for i in range(len(Struct)):  # sequence length
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))
    return lista

# Parse an RNAsubopt file to extract Base pairs
def GetBasePairsFromStructFile(faPath):   #return dic={structure:[liste de pairs de base ],....}
	#print faPath
	DicStruct={}
    	lines = FF.Parsefile(faPath)
	#print lines
    	SeqLen=len(lines[1])-1
	#print SeqLen,"seq length"
    	for j in range(len(lines)):
		DicStruct[j]=ListBasePairsFromStruct(lines[j].strip().split(' ')[0])
    	return len(lines),DicStruct
from collections import Counter
def OccurenceBasePairs(ListPairs, Threshold):
    return [(elem[0],elem[1],Counter(ListPairs)[elem]) for elem in Counter(ListPairs)  if Counter(ListPairs)[elem]>=Threshold]
def fromPairsToStruct(rna, Pairs):
    structure=["." for i in range(len(rna)-1)]
    for (i,j) in Pairs:
        structure[i]='('
        structure[j]=')'
    return "".join(structure)

#Calculate distance between each 2 structures and generate a matrix 2 by 2 , stored in the file SVMLFile
def DistanceTwoBPlist(Struct1,Struct2):
	return len(set(Struct1).symmetric_difference(set(Struct2)) )
	#return len(set(Struct1).intersection(set(Struct2)) )

def dd():
    return 0
def aa():
    return defaultdict(dd)

def DistanceStruct(StructFile, SVMlFile, numberofsruct, constrainte):
    Redondantestructure = defaultdict(aa)
    MatDist = defaultdict(aa)
    Redondantestructure1 = []
    Dicnumberofsruct = {}

    for i in range(len(constrainte)):
        Dicnumberofsruct[constrainte[i]] = numberofsruct


    nb, DicStruct = GetBasePairsFromStructFile(StructFile)

    for i in range(0, nb):
        for j in range(i + 1, nb):
            MatDist[i][j] = DistanceTwoBPlist(DicStruct[i], DicStruct[
                j])
            MatDist[j][i] = MatDist[i][j]
            ####### Check for redundancy
            if MatDist[i][j] == 0:
                jconstraint=int(j / numberofsruct)
                if j not in Redondantestructure1 and int(i / numberofsruct)==jconstraint:# To be sure that the redundant  structure belongs to the same probing condition
                    Dicnumberofsruct[constrainte[jconstraint]] -= 1
                    Redondantestructure1.append((j,jconstraint))

    for (elem,jconstraint) in Redondantestructure1:
        StructureNumber = elem - jconstraint * numberofsruct
        Redondantestructure[constrainte[jconstraint]][StructureNumber] = 1 # we mark redundant structures by value 1

    # store the distance matrix in the  SVMLFile
    if os.path.isfile(os.path.join("output",SVMlFile)):
        os.remove(os.path.join("output",SVMlFile)) # To clean the previous version
    o = open(os.path.join("output",SVMlFile), "w")
    for i in range(len(MatDist)):
        o.write("%i\t" % (i + 1))
        for j in range(len(MatDist)):
            if (i != j):
                o.write("%i:%.4f\t" % (j + 1, MatDist[i][j]))
        o.write("\n")
    o.close()

    if Redondantestructure!=0:
        print "Warning! redundant structures"
    FF.PickleVariable(MatDist, os.path.join(conf.PickledData,"dissmatrix.pkl"))
    FF.PickleVariable(Redondantestructure1, os.path.join(conf.PickledData,"Redondantestructures.pkl"))
    FF.PickleVariable(Redondantestructure, os.path.join(conf.PickledData, "Redondantestructures_Id.pkl"))
    FF.PickleVariable(Dicnumberofsruct,os.path.join(conf.PickledData,"Dicnumberofsruct.pkl"))
    return 0


def FromStructFiletoRNAEvalInput(StructFile, InputRNAeval, rna):
    lines = FF.Parsefile(StructFile)
    o = open(InputRNAeval,
             "w")  # geneate intermediate file with sequence+strcuture , seq+strcture .... as the input format  to use RNAeval
    # print "sdfqspkojr",len(lines)
    for i in range(1, len(lines)):
        o.write("%s%s\t" % (rna, lines[i]))
    o.close()

#StructFile contains the RNA sequence in the first line and list of correponding structures by line
def ENERGY_VALUES_STRUCTURES(StructFile,rna):
    	Energy=[]
    	#generate the rnaeval input file
    	FromStructFiletoRNAEvalInput(StructFile,"InputRNAeval",rna)
   	 # launch the RNaeval command
        os.system('RNAeval <' + "InputRNAeval" + '>' + "energyvalues")
    	# Parse the RNAevaloutput to extract energy values
    	lines=FF.Parsefile("energyvalues")
    	for i in xrange(1,len(lines),2):
        	# i is the stucture number and 'lines[i].split(" ")[1][1:-2]' is  the  corresponding  energy value
		#print 'holla',(lines[i].split(" ")[1][1:-2])
		Energy.append(lines[i].split(" ",1)[1][1:-2]) # TODO ,1 is to get the first occrence of the space !!!
    	return Energy

##Boltzmman energy   according to the formula B=exp^\frac{-e}{RT}
def BoltzmannFactor(Energy):
    T=37+273.15
    R=0.0019872370936902486
    #print np.exp(-float(Energy)/float(100.*R*T))
    return np.exp(-float(Energy)/float(100.*R*T))

def Boltzmann_Calc(constraintes, StructfileRepository, numberofsruct,  rna, Redondantestructure):
    Energy = defaultdict(aa)
    Boltzman = defaultdict(aa)
    ConditionalBoltzman = defaultdict(aa)
    ZBolzman = defaultdict(aa)
    # Calculate estructure energies in each condition sample
    for Condition in constraintes:
        FileStructure = os.path.join(StructfileRepository, Condition)
        Energy[Condition] = ENERGY_VALUES_STRUCTURES(FileStructure,rna)  # list of energy values for the structures present in the Condition

    for Condition in constraintes:
        ListwithoutRedundnacy = []
        for i in range(numberofsruct):
            Boltzman[Condition][i] = BoltzmannFactor(Energy[Condition][i])
            if Redondantestructure[Condition][i]==0:  # if the structure is not redundant
                ListwithoutRedundnacy.append(Boltzman[Condition][i])

        # Calculate the normalization term as the sum over all Boltzmann probabilities for one copy of each structure
        ZBolzman[Condition] = sum(ListwithoutRedundnacy)  # Partition function

    #FF.PickleVariable(Boltzman, os.path.join(conf.PickledData, "Boltzman.pkl"))
    listall = []
    for Condition in constraintes:  # to not count MFES
        lista = []
        for i in range(numberofsruct):
            if Redondantestructure[Condition][i]==0: # a npn redundnat structure
                lista.append(Boltzman[Condition][i] / ZBolzman[Condition])
            else:
                lista.append(0)  # Redundant structures have a conditional Boltzmann value NULL
        listall += lista
        ConditionalBoltzman[Condition] = lista
        #print "Condition \t  ConditionalBoltzman", Condition, ConditionalBoltzman[Condition]
    FF.PickleVariable(Boltzman, os.path.join(conf.PickledData, "Boltzman.pkl"))
    FF.PickleVariable(ConditionalBoltzman, os.path.join(conf.PickledData, "ConditionalBoltzman.pkl"))
    FF.PickleVariable(ZBolzman, os.path.join(conf.PickledData, "ZBolzman.pkl"))

    return ConditionalBoltzman

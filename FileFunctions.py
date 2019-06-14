import os,pickle
from itertools import islice
from os import listdir
from os.path import isfile, join
#create folder if it doesn't exist
def CreateFold(dir):
 	try:
    		os.stat(dir)
	except:
         	os.mkdir(dir)


    # get all files with a specific extension from a specific path
def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]

#Parse a file by returning lines it contains
def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = fileIn.readlines()
    fileIn.close()
    return lines


def GetlinefromFile(Path, Linenumber):
    return Parsefile(Path)[Linenumber]

def MergeFiles(Path,output,fileslist ,Startline):
	print "Sample construction from the Probing conditions", fileslist
	with open(output, 'w') as outfile:
    		for fname in[Path+'/'+i for i in fileslist]:
        		with open(fname) as infile:
            			for line in islice(infile, Startline, None):
                			outfile.write(line)
	return output


def MergespecificsFiles(Path,output,fileslist ,Startline):
	with open(output, 'w') as outfile:
    		for fname in[Path+'/'+i for i in fileslist]:
        		with open(fname) as infile:
            			for line in islice(infile, Startline, None):
                			outfile.write(line)
	return output

def parseReactivityfile(fileinput):
	Reactvities=[]
	lines=Parsefile(fileinput)
	for it in range(len(lines)):
                Reactvities.append(lines[it].split("\t")[1])
	return Reactvities

def PickleVariable(variable,file):
    fileOut = open(os.path.join(os.getcwd(), file), "wb")# 'wb' instead 'w' for binary file
    pickle.dump(variable, fileOut,-1) # -1 specifies highest binary protocol
    fileOut.close()

def UnpickleVariable(file):
    fileIn = open(os.path.join(os.getcwd(),file), "rb")
    unpickled = pickle.load(fileIn)
    fileIn.close()
    return unpickled
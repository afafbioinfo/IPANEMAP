import os, pickle
from itertools import islice
from os import listdir
from os.path import isfile, join
from Progress import progress
from conf import loadConfig

# create folder if it doesn't exist
def CreateFold(dir):
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)


# get all files with a specific extension from a specific path
def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if
            isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]


# Parse a file by returning lines it contains
def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = [l.strip() for l in fileIn.readlines()]
    fileIn.close()
    return lines


def GetlinefromFile(Path, Linenumber):
    return Parsefile(Path)[Linenumber]


def MergeFiles(Path, output, fileslist, Startline):
    progress.Print("Merging samples for conditions %s"%(fileslist))
    with open(output, 'w') as outfile:
        for fname in [Path + '/' + i for i in fileslist]:
            with open(fname) as infile:
                for line in islice(infile, Startline, None):
                    outfile.write(line)
    return output


def MergespecificsFiles(Path, output, fileslist, Startline):
    with open(output, 'w') as outfile:
        for fname in [Path + '/' + i for i in fileslist]:
            with open(fname) as infile:
                for line in islice(infile, Startline, None):
                    outfile.write(line)
    return output


def parseReactivityfile(fileinput):
    Reactvities = []
    lines = Parsefile(fileinput)
    for it in range(len(lines)):
        data = lines[it].split()
        if len(data)>1:
            Reactvities.append(data[1])
    return Reactvities


def PickleVariable(variable, file):
    conf = loadConfig()
    fileOut = open(os.path.join(conf.OutputFolder,"tmp", conf.PickledData, file), "wb")  # 'wb' instead 'w' for binary file
    pickle.dump(variable, fileOut, -1)  # -1 specifies highest binary protocol
    fileOut.close()


def UnpickleVariable(file):
    conf = loadConfig()
    fileIn = open(os.path.join(conf.OutputFolder,"tmp", conf.PickledData,  file), "rb")
    unpickled = pickle.load(fileIn)
    fileIn.close()
    return unpickled


def LocateSoftConstraintFile(rna, cond):
    conf = loadConfig()
    path = os.path.join(conf.PathConstraintsFileShape, rna + cond + ".txt")
    if not os.path.isfile(path):
        path = ""
    return path

class IPANEMAPError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg
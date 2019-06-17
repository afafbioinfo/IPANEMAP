import FileFunctions as FF
import os, subprocess
import conf as conf
from Progress import progress

NUM_HEADER_LINES = 2
THERMO_CONDITION = "Thermo"

def StructSampling(Pathconstraints, Conditions, numberStructures, T, m, b):
    dir = os.path.join(conf.OutputFolder,"tmp",'OutputSamples' + str(numberStructures))
    FF.CreateFold(dir)
    for filename in Conditions:
        lines = []
        header = []
        progress.StartTask("Processing %s"%(filename))
        while len(lines) - NUM_HEADER_LINES < numberStructures:
            Path = ""
            for p in Pathconstraints:
                Input = os.path.join(p, filename + '.' + conf.FASTAExtension)
                if os.path.isfile(Input):
                    Path = p
                    break
            output = os.path.join(dir, filename)
            Command = 'RNAsubopt  -p ' + str(numberStructures) + ' -s -T ' + str(T)
            if os.path.isfile(Input):  # to take into account the case with no constraints (free-constraint)
                if filename.endswith(THERMO_CONDITION):
                   pass
                elif Path == conf.PathConstraintsFile:
                    Command += ' -C --enforceConstraint'
                elif Path == conf.PathConstraintsFileShape:
                    ShapeFile = conf.PathConstraintsFileShape + os.sep + filename + 'Probing.txt'
                    if not os.path.isfile(ShapeFile):
                        raise Exception("Error: Missing reactivity file for condition %s" % (filename))
                    Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm' + str(m) + 'b' + str(b) + '"'
                else:
                    raise Exception("Error: Cannot process FASTA file %s (file outside of hard/soft directories, and not thermo)" % (filename))
            else:
                raise Exception("Error: Missing input FASTA file for condition %s"%(filename))
            #progress.Print(Command)
            subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(output, 'wb'),
                            stderr=open(os.devnull, 'w'), shell=True)
            with open(output, 'r') as f:
                nlines = f.readlines()
                header = nlines[:NUM_HEADER_LINES]
                lines += nlines[NUM_HEADER_LINES:]
        with open(output, 'w') as f:
            f.writelines(header+lines[:numberStructures])
        progress.EndTask()
    return dir

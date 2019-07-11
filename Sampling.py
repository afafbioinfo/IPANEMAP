import FileFunctions as FF
import os, subprocess
from conf import loadConfig
from Progress import progress
import IPANEMAP

NUM_HEADER_LINES = 2

def StructSampling(Pathconstraints, Conditions, numberStructures, T, m, b, defaultFasta):
    conf = loadConfig()
    dir = os.path.join(conf.OutputFolder, "tmp", 'OutputSamples' + str(numberStructures))
    FF.CreateFold(dir)
    thermoMsgShown = False
    for filename in Conditions:
        lines = []
        header = []
        progress.StartTask("Processing %s"%(filename))
        while len(lines) - NUM_HEADER_LINES < numberStructures:

            # If alternative sequence file found in constraints folder, use it rather than default
            Input = defaultFasta
            for p in Pathconstraints:
                tmpInput = os.path.join(p, filename + '.' + IPANEMAP.FASTA_EXTENSION)
                if os.path.isfile(tmpInput):
                    Input = tmpInput

            output = os.path.join(dir, filename)
            Command = 'RNAsubopt  -p ' + str(numberStructures) + ' -s -T ' + str(T)

            (hasHardConstraints, hasSoftConstraints) = (False, False)

            hardConstraintFile = os.path.join(conf.PathConstraintsFile, filename + '.txt')
            if os.path.isfile(hardConstraintFile):
                Command += ' -C --enforceConstraint '
                hasHardConstraints = True
                Input = hardConstraintFile

            ShapeFile = os.path.join(conf.PathConstraintsFileShape, filename + '.txt')
            if os.path.isfile(ShapeFile):
                Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm' + str(m) + 'b' + str(b) + '"'
                hasSoftConstraints = True

            if not (hasHardConstraints or hasSoftConstraints or thermoMsgShown):
                progress.Print("Warning: Did not find suitable constraint file for this condition, using purely thermodynamic sampling")
                thermoMsgShown = True
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

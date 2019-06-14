import FileFunctions as FF
import os, subprocess
import conf as conf


def StructSampling(Pathconstraints, numberofsruct, Tmpr,Extension,m,b):
    dir='OutputSamples' + str(numberofsruct)
    FF.CreateFold(dir)
    for Path in Pathconstraints:
        for filename in FF.GetListFile(Path, Extension):
	    print "Processing the sampling for ",filename
            Input = os.path.join(Path, filename+'.'+Extension)
            output = os.path.join(dir , filename)
            Command = 'RNAsubopt  -p ' + str(numberofsruct) + ' -s -T ' + str(Tmpr)
	    if os.path.isfile( Path): # to take into account the case whith no constraints (free-constraint)
            	if Path ==conf.PathConstrainteFile:
                	Command += ' -C  --enforceConstraint'
            	if Path == conf.PathConstrainteFileShape:
                	ShapeFile = conf.PathConstrainteFileShape + "/" + filename + 'Probing.txt'
               		Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm'+str(m)+'b'+str(b)+'"'
            subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(output, 'wb'), shell=True)

            #because the version 2.3 of rnaeval does not consider the rna, second line should be removed
            with open(output, 'r') as f:
                lines = f.readlines()
            with open(output, 'w') as f:
                f.writelines(lines[:1] + lines[2:])
            print "Sampling done  with success for ", filename
    return dir


            
           



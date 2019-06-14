from ConfigParser import SafeConfigParser 
import sys

#Connect to  the config file
config= SafeConfigParser()
config.read("IPANEMAP.Config")
#Get parameters
PathSequenceFasta=config.get("Paths","PathSequenceFasta")
PathConstrainteFile=config.get("Paths","PathConstrainteFile")
PathConstrainteFileShape=config.get("Paths","PathConstraintsProbing")
OutputLogfile=config.get("Paths","OutputLogfile")
FastaExtension= config.get("Paths","FastaExtension")
PickledData=config.get("Paths","PickledData")
#sampling
sampling=config.get("Probing","sampling")
numberofsruct=config.get("Probing","numberofsruct")
#MFEs=(config.get("Probing","MFEs")).split(',')

Temperature=config.get("Probing","Temperature")
m=config.get("Probing","m")
b=config.get("Probing","b")
SHAPEVis=config.get("Probing","SHAPEVis")
percent=int(config.get("Pareto","percent"))
CutoffZcondition=float(config.get("Pareto","CutoffZcondition"))
cutoff=config.get("Probing","cutoffBasePairs")
Psdotpath=config.get("Paths","Psdotpath")
Matrixproba=config.get("Paths","Matrixproba")
Fastaextenstion=config.get("Paths","Extension")
#constraintes=config.get("Conditions","Constraintes")
constraintes=(config.get("Probing","constraints")).split(',')
MiniBatchKMean=config.get("Clustering","MiniBatchKMean")
Diana=config.get("Clustering","DIANA")
maxDiameterThreshold = float(config.get("Clustering","maxDiameterThreshold"))
maxAverageDiameterThreshold = float(config.get("Clustering","maxAverageDiameterThreshold"))
#function to cretae a log file
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

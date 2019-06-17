from ConfigParser import SafeConfigParser 
import sys


class Logger(object):
    """ Class to manage basic logging """
    def __init__(self, filename="IPANEMAP.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)


# Connect to the config file
config = SafeConfigParser()
config.read("IPANEMAP.Config")


# General Input
PathConstraintsFile = config.get("Input", "HardConstraintsDir")
PathConstraintsFileShape = config.get("Input", "SoftConstraintsDir")
PathRNAFASTA = config.get("Input", "RNADir")
Conditions = [s.strip() for s in (config.get("Input", "Conditions")).split(',')]
PickledData = "PickledData"

# Path section
OutputFolder = config.get("Paths", "WorkingDir")
FASTAExtension = config.get("Paths", "FASTAExtension")
OutputLogfile = config.get("Paths", "LogFile")

# Sampling section
Sampling = config.get("Sampling", "DoSampling")
SampleSize = config.get("Sampling", "NumStructures")
Temperature = config.get("Sampling", "Temperature")
m = config.get("Sampling", "m")
b = config.get("Sampling", "b")
# VisualizedCondition = config.get("Sampling", "VisualizedCondition")
# cutoff = config.get("Sampling", "cutoffBasePairs")


# Pareto Front section
percent = int(config.get("Pareto", "Percent"))
CutoffZcondition = float(config.get("Pareto", "CutoffZCondition"))


maxDiameterThreshold = float(config.get("Clustering", "maxDiameterThreshold"))
maxAverageDiameterThreshold = float(config.get("Clustering", "maxAverageDiameterThreshold"))


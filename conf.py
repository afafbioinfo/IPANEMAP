from ConfigParser import SafeConfigParser 
import sys
import argparse
from argparse import Namespace


class Logger(object):
    """ Class to manage basic logging """
    def __init__(self, filename="IPANEMAP.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)


conf = None

def loadConfig():
    def override(init, new):
        for x in vars(new):
            vars(init)[x] = vars(new)[x]
    global conf
    if conf is None:
            conf = Namespace()
            # Connect to the config file
            config = SafeConfigParser()
            config.read("IPANEMAP.Config")

            # General Input
            conf.PathConstraintsFile = config.get("Input", "HardConstraintsDir")
            conf.PathConstraintsFileShape = config.get("Input", "SoftConstraintsDir")
            conf.PathRNAFASTA = config.get("Input", "RNADir")
            conf.Conditions = [s.strip() for s in (config.get("Input", "Conditions")).split(',')]
            conf.PickledData = "PickledData"

            # Path section
            conf.OutputFolder = config.get("Paths", "WorkingDir")
            conf.FASTAExtension = config.get("Paths", "FASTAExtension")
            conf.OutputLogfile = config.get("Paths", "LogFile")

            # Sampling section
            conf.Sampling = config.get("Sampling", "DoSampling")
            conf.SampleSize = config.get("Sampling", "NumStructures")
            conf.Temperature = config.get("Sampling", "Temperature")
            conf.m = config.get("Sampling", "m")
            conf.b = config.get("Sampling", "b")
            # VisualizedCondition = config.get("Sampling", "VisualizedCondition")
            # cutoff = config.get("Sampling", "cutoffBasePairs")


            # Pareto Front section
            conf.percent = int(config.get("Pareto", "Percent"))
            conf.CutoffZcondition = float(config.get("Pareto", "CutoffZCondition"))

            conf.maxDiameterThreshold = float(config.get("Clustering", "maxDiameterThreshold"))
            conf.maxAverageDiameterThreshold = float(config.get("Clustering", "maxAverageDiameterThreshold"))


            # Load additional command line options
            parser = argparse.ArgumentParser(
                description='IPANEMAP: Integrative Probing Analysis Informed by Multiple Accessibility Profiles')
            parser.add_argument('--conditions', metavar='Cond', nargs='+', dest='Conditions',
                                help='One or several conditions')
            args = parser.parse_args()
            override(conf, args)

    return conf

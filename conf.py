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
            if vars(new)[x] is not None:
                vars(init)[x] = vars(new)[x]
    global conf
    if conf is None:
            conf = Namespace()
            # Connect to the config file
            config = SafeConfigParser()
            config.read("IPANEMAP.cfg")

            # General Input
            conf.PathConstraintsFile = config.get("Input", "HardConstraintsDir")
            conf.PathConstraintsFileShape = config.get("Input", "SoftConstraintsDir")
            conf.RNA = config.get("Input", "RNA")
            conf.Conditions = [s.strip() for s in (config.get("Input", "Conditions")).split(',')]
            conf.PickledData = "PickledData"

            # Path section
            conf.OutputFolder = config.get("Paths", "WorkingDir")
            conf.OutputLogfile = config.get("Paths", "LogFile")

            # Sampling section
            conf.Sampling = config.get("Sampling", "DoSampling")
            conf.SampleSize = config.get("Sampling", "NumStructures")
            conf.Temperature = config.get("Sampling", "Temperature")
            conf.m = config.get("Sampling", "m")
            conf.b = config.get("Sampling", "b")

            # Pareto Front section
            conf.percent = int(config.get("Pareto", "Percent"))
            conf.CutoffZcondition = float(config.get("Pareto", "ZCutoff"))

            # Clustering section
            conf.maxDiameterThreshold = float(config.get("Clustering", "MaxDiameterThreshold"))
            conf.maxAverageDiameterThreshold = float(config.get("Clustering", "MaxAverageDiameterThreshold"))

            # Visualization section
            conf.DrawModels = config.get("Visualization", "DrawModels")
            conf.DrawCentroids = config.get("Visualization", "DrawCentroids")
            conf.ShowProbing = config.get("Visualization", "ShowProbing")

            # Load additional command line options
            parser = argparse.ArgumentParser(
                description='Integrative Probing Analysis Informed by Multiple Accessibility Profiles',
                epilog="Further options can be set by editing IPANEMAP.cfg")
            parser.add_argument('--RNA', metavar='r', dest='RNA',
                                help='path to the RNA considered for the analysis (FASTA file)')
            parser.add_argument('--cond', metavar='c', nargs='+', dest='Conditions',
                                help='space-separated list of probing conditions')
            args = parser.parse_args()
            override(conf, args)

    return conf

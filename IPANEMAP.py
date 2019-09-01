import time
import os
import sys
from collections import defaultdict
from conf import loadConfig, Logger
import FileFunctions as FF
import Sampling as SP
import StructureFunctions as SF
import StructureFunctions as SF
import VisualizationTools as VT
import ClustersTrait as CT
import Optimize_clustering as OC
from Progress import progress

FASTA_EXTENSION = "fa"

if __name__ == "__main__":
    conf = loadConfig()

    # Create folders
    FF.CreateFold(conf.OutputFolder)
    FF.CreateFold(os.path.join(conf.OutputFolder, "tmp"))
    FF.CreateFold(os.path.join(conf.OutputFolder, "tmp", conf.PickledData))

    # Redirects all the print to the output Log file
    sys.stdout = Logger(os.path.join(conf.OutputFolder,conf.OutputLogfile))

    # ******************************** Generate sample

    try:
        rna = os.path.split(conf.RNA)[-1]
        RNAName = rna[:-(len(FASTA_EXTENSION)+1)]
        progress.StartTask("Processing RNA %s" % (RNAName))
        if not os.path.isfile(conf.RNA):
            raise FF.IPANEMAPError("Input file '%s' not found"%(conf.RNA))

        # Get the rna sequence
        RNASequence = FF.Parsefile(conf.RNA)[1].strip()

        # Get probing conditions for the treated RNA
        ProbingConditions = [RNAName + state for state in conf.Conditions]

        # Specify  whether to generate new sample or use a previously  generated one
        OutputSamples = os.path.join(conf.OutputFolder, "tmp",'OutputSamples' )+ conf.SampleSize
        if str.lower(conf.Sampling) == "true" or not os.path.isdir(OutputSamples):
            progress.StartTask("Sampling %s structures for each condition" % (conf.SampleSize))
            OutputSamples = SP.StructSampling([conf.PathConstraintsFile, conf.PathConstraintsFileShape],
                                              ProbingConditions,
                                              int(conf.SampleSize),
                                              conf.Temperature, conf.m, conf.b, conf.RNA)
            progress.EndTask()
        else:
            progress.Print("Using existing sample")

        progress.Print("Probing conditions: %s" % (ProbingConditions))
        # Create a global file that contains structures sampled from the list of Probing conditions
        FF.MergeFiles(OutputSamples, os.path.join(OutputSamples, 'Samples.txt'), ProbingConditions, SP.NUM_HEADER_LINES)

        # Create a distance matrix file
        progress.StartTask("Computing dissimilarity matrix")
        SVMlFile = "DissimilarityMatrix" + conf.SampleSize
        # Calculate distance and identify redundant structures within the same condition
        SF.DistanceStruct(os.path.join(OutputSamples, 'Samples.txt'), SVMlFile, int(conf.SampleSize), ProbingConditions)
        progress.EndTask()

        ################################# Calculate Conditional Boltzmann probabilities
        # for each condition, calculate Z over all non redundant structures and return a conditional Boltzmann probability for all structures with null value for redundant ones.
        progress.StartTask("Computing Boltzmann probabilities")
        BoltzmannFactor = defaultdict(lambda: defaultdict())
        ConditionalBoltzmannProbability = defaultdict(lambda: defaultdict())
        Zprobabilities = defaultdict(lambda: defaultdict())
        Redondantestructure = FF.UnpickleVariable("Redondantestructures_Id.pkl")
        ConditionalBoltzmannProbability = SF.Boltzmann_Calc(ProbingConditions, OutputSamples, int(conf.SampleSize),
                                                            RNASequence,
                                                            Redondantestructure)
        progress.EndTask()

        ################################# Clustering of structures based on their base pair distance
        progress.StartTask("Iterative clustering")
        # Load the pickled dissimilarity matrix
        DM = FF.UnpickleVariable("dissmatrix.pkl")
        # Get the list of redundant structures
        Redundant = []
        Redundant = FF.UnpickleVariable("Redondantestructures.pkl")
        BoltzmannFactor = FF.UnpickleVariable("Boltzman.pkl")
        method = "MiniBatchKMean"
        Clusters, CentroidStructure = OC.DefineNumberCluster(os.path.join(conf.OutputFolder, "tmp", SVMlFile),
                                                             Redundant, method,
                                                             DM, BoltzmannFactor, ProbingConditions, RNASequence)
        # Get Clusters from Pickled data
        # Clusters = FF.UnpickleVariable("Clusters" + method + ".pkl")
        progress.EndTask()

        progress.StartTask("Analyzing clusters/Drawing centroids")
 

        centroidPath = os.path.join(conf.OutputFolder, "Centroids."+FASTA_EXTENSION)
        SF.StructsToRNAEvalInput(CentroidStructure, centroidPath, RNASequence)
        Centroids_Energies = SF.RunEval(centroidPath)

        CT.ClustersDistributions(Clusters, os.path.join(conf.OutputFolder, "Clusters_details"), ProbingConditions,
                                 int(conf.SampleSize))

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!Election of the best structures starting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Calculate cumulated  Boltzmaan Energy for each cluster

        CumulBE = {}
        CardinalConditions = {}

        CumulBE = CT.CumulatedConditionalBoltzmannbyCluster(Clusters, ConditionalBoltzmannProbability,
                                                            int(conf.SampleSize), ProbingConditions)
        # epsilon= 1/(n+1)
        Epsilon = 1. / float(len(Clusters) + 1)
        CardinalConditions = CT.GetCardinalConditions(Clusters, ConditionalBoltzmannProbability, ProbingConditions,
                                                      int(conf.SampleSize), Epsilon)

        ImageFolder = os.path.join(conf.OutputFolder, "img")
        if not os.path.isdir(ImageFolder):
            os.mkdir(ImageFolder)

        ProbingPath = ""
        if conf.ShowProbing:
            cond = conf.Conditions[0]
            ProbingPath = FF.LocateSoftConstraintFile(RNAName, cond)

        for index in (CentroidStructure.keys()):
            progress.Print("%s\t%s\t%s\t%s\t%s" % (index+1, CentroidStructure[index], Centroids_Energies[index],
                                                   CardinalConditions[index], CumulBE[index]))
            if conf.DrawCentroids:
                VT.drawStructure(RNASequence, CentroidStructure[index], ProbingPath, os.path.join(ImageFolder, "Centroid-%s.svg"%(index+1)))

        progress.Print('Pareto optimal structure(s):')
        
        Dict = {}  # a dictionary that contains the three variables characterizing clusters (Criteria of optimization)
        for ClusterNumber in Clusters:
            Dict[ClusterNumber] = [CardinalConditions[ClusterNumber],
                                   CumulBE[ClusterNumber]]  # only Two criterions are considered

        ListOptimalClusters = CT.Pareto(Dict)
        progress.Print("Structure\tdG\t#SupportingConditions\tBoltzmannProbability", output=True)
        i = 0
        for index in ListOptimalClusters:
            progress.Print("%s\t%s\t%s\t%s" % (CentroidStructure[index], Centroids_Energies[index],
                                               CardinalConditions[index], CumulBE[index]), output=True)
            if conf.DrawModels:
                VT.drawStructure(RNASequence, CentroidStructure[index], ProbingPath, os.path.join(ImageFolder, "Optimal-%s.svg"%(i+1)))
            i += 1
        progress.EndTask()
        
    except FF.IPANEMAPError as e:
        progress.Print("Error: %s" % (e))
        progress.Flush()

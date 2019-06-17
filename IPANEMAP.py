import time
import os
import sys
from collections import defaultdict
import conf
import FileFunctions as FF
import Sampling as SP
import StructureFunctions as SF
import VisualizationTools as VT
import ClustersTrait as CT
import Optimize_clustering as OC
from Progress import progress

# Create folders
FF.CreateFold(conf.OutputFolder)
FF.CreateFold(os.path.join(conf.OutputFolder, "tmp"))
FF.CreateFold(os.path.join(conf.OutputFolder, "tmp", conf.PickledData))

# Redirects all the print to the output Log file
sys.stdout = conf.Logger(os.path.join(conf.OutputFolder,conf.OutputLogfile))

# ******************************** Generate sample of structures for all fasta files in the fastaShape/fastaConstraint folders

# ******************************** Treat all RNA sequences in the fasta folder
for RNAName in FF.GetListFile(conf.PathRNAFASTA, conf.FASTAExtension):
    progress.StartTask("Processing RNA %s" % (RNAName))
    # Get the rna sequence
    RNASequence = FF.Parsefile(os.path.join(conf.PathRNAFASTA, RNAName + '.' + conf.FASTAExtension))[1].strip()

    # Get probing conditions for the treated RNA
    ProbingConditions = [RNAName + state for state in conf.Conditions]

    # Specify  whether to generate new sample or use a previously  generated one
    if str.lower(conf.Sampling) == "true":
        progress.StartTask("Sampling %s structures for each condition" % (conf.SampleSize))
        OutputSamples = SP.StructSampling([conf.PathConstraintsFile, conf.PathConstraintsFileShape],
                                          ProbingConditions,
                                          int(conf.SampleSize),
                                          conf.Temperature, conf.m, conf.b)
        progress.EndTask()
    else:
        progress.StartTask("Using existing sample of %s structures" % (conf.SampleSize))
        OutputSamples = 'OutputSamples' + conf.SampleSize
        progress.EndTask()

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
    BoltzmanFactor = defaultdict(lambda: defaultdict())
    ConditionalBoltzmannProbability = defaultdict(lambda: defaultdict())
    Zprobabilities = defaultdict(lambda: defaultdict())
    Redondantestructure = FF.UnpickleVariable("Redondantestructures_Id.pkl")
    ConditionalBoltzmannProbability = SF.Boltzmann_Calc(ProbingConditions, OutputSamples, int(conf.SampleSize), RNASequence,
                                                        Redondantestructure)
    progress.EndTask()

    ################################# Clustering of structures based on their base pair distance
    progress.StartTask("Iterative clustering")
    # Load the pickled dissimilarity matrix
    DM = FF.UnpickleVariable("dissmatrix.pkl")
    # Get the list of redundant structures
    Redundant = []
    Redundant = FF.UnpickleVariable("Redondantestructures.pkl")
    BoltzmanFactor = FF.UnpickleVariable("Boltzman.pkl")
    method = "MiniBatchKMean"
    Clusters, CentroidStructure = OC.DefineNumberCluster(os.path.join(conf.OutputFolder,"tmp", SVMlFile), Redundant, method,
                                                         DM, BoltzmanFactor, ProbingConditions, RNASequence)
    # Get Clusters from Pickled data
    # Clusters = FF.UnpickleVariable("Clusters" + method + ".pkl")
    progress.EndTask()

    progress.StartTask("Analyzing clusters")
    # We should create an intermediate  file to be sure that  RNAeval  works!!


    centroidPath = os.path.join(conf.OutputFolder,"Centroids.fa")
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

    for index in (CentroidStructure.keys()):
        progress.Print("%s\t%s\t%s\t%s\t%s" % (index, CentroidStructure[index], Centroids_Energies[index],
                                       CardinalConditions[index], CumulBE[index]))


    progress.Print('Pareto optimal structure(s):')

    # print 'mean_distance_in_a_given_cluster', "\t", E.values()
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Cluster analysis to elect the best cluster
    Dict = {}  # a dictionary that contains the three variables characterizing clusters (Criteria of optimization)
    for ClusterNumber in Clusters:
        Dict[ClusterNumber] = [CardinalConditions[ClusterNumber],
                               CumulBE[ClusterNumber]]  # only Two criterions are considered

    ListOptimalClusters = CT.Pareto(Dict)
    progress.Print("Structure\tdG\t#SupportingConditions\tBoltzmannProbability",output=True)
    for index in ListOptimalClusters:
        progress.Print("%s\t%s\t%s\t%s" % (CentroidStructure[index], Centroids_Energies[index],
                                       CardinalConditions[index], CumulBE[index]),output=True)

    progress.EndTask()

    #DominatedClusters = [clusteri for clusteri in Clusters if clusteri not in ListOptimalClusters]
    #Filealloptimalcentroides = os.path.join("Multiprobing/", str(ProbingConditions) + ".optimals")
    # VT.Drawvarna(Filealloptimalcentroides, ListOptimalClusters, CentroidStructure, conf.numberofsruct,
    #	         rna, Centroids_Energies, conf.SHAPEVis)
    progress.EndTask()

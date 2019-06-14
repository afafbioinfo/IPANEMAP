
import conf , FileFunctions as FF, Sampling as SP, StructureFunctions as SF, VisualizationTools as VT, ClustersTrait as CT, Optimize_clustering as OC
import time,os,sys
from collections import defaultdict


#Redirect all the print to an output Log file
sys.stdout =conf.Logger(conf.OutputLogfile)


#******************************** Generate sample of structures for all fasta files in the fastaShape/fastaConstraint folders
# Specify  whether to generate new sample or use a previously  generated one

if str.lower(conf.sampling)=="true":
    print("Sampling Process for % s Structures" % (conf.numberofsruct))
    OutputSamples = SP.StructSampling([conf.PathConstrainteFile, conf.PathConstrainteFileShape],conf.numberofsruct, conf.Temperature, conf.Fastaextenstion,conf.m,conf.b)
else:
    print ("Use of an existing sample for % s Structures" % (conf.numberofsruct))
    OutputSamples = 'OutputSamples' + conf.numberofsruct

#******************************** Treat all RNA sequences in the fasta folder
for file in FF.GetListFile(conf.PathSequenceFasta, conf.FastaExtension):
	print "RNAId ", "\t",file
	startimebig = time.time()
    # Get the rna sequence
	rna = FF.Parsefile(os.path.join(conf.PathSequenceFasta, file + '.' + conf.FastaExtension))[1]

    # Get probing conditions for the treated RNA
	Probingconditions = [file + state for state in conf.constraintes]
	print "The considered probing conditions are ", Probingconditions
	# Create a global file that contains structures sampled from the list of Probing conditions
	FF.MergeFiles(OutputSamples, os.path.join(OutputSamples, 'Samples.txt'), Probingconditions, 1)

	# Create a distance matrix file
	SVMlFile = "DissimilarityMatrix" + conf.numberofsruct
	startime = time.time()

	print("Distance Matrix generation for % d Structures" % (int(conf.numberofsruct) * len(Probingconditions)))
	# Calculate distance and identify redundant structures within the same condition
	SF.DistanceStruct(os.path.join(OutputSamples, 'Samples.txt'), SVMlFile, int(conf.numberofsruct), Probingconditions)
	endtime = time.time()
	print("End of distance calculation between the structures in the sample  %53f\t" % (endtime - startime))

	################################# Calculate Conditional Boltzmann probabilities
	# for each condition, calculate Z over all non redundant structures and return a conditional Boltzmann probability for all structures with null value for redundant ones.
	startime = time.time()
	print ("Start Conditional  Boltzmann Probabilities calculation")
	startime = time.time()
	BoltzmanFactor = defaultdict(lambda: defaultdict())
	ConditionalBoltzmaanProbability = defaultdict(lambda: defaultdict())
	Zprobabilities = defaultdict(lambda: defaultdict())
	Redondantestructure = FF.UnpickleVariable(os.path.join(conf.PickledData, "Redondantestructures_Id.pkl"))
	ConditionalBoltzmaanProbability = SF.Boltzmann_Calc(Probingconditions, OutputSamples, int(conf.numberofsruct), rna,Redondantestructure)
	endtime = time.time()
	print ("Boltzmann probabilities Calculated with success in :  %53f\t" % (endtime - startime))

	################################# Clustering of structures based on their base pair distance
	# Load the pickled dissimilarity matrix
	DM = FF.UnpickleVariable(os.path.join(conf.PickledData, "dissmatrix.pkl"))
	# Get the list of redundant structures
	Redundant = []
	Redundant=FF.UnpickleVariable(os.path.join(conf.PickledData,"Redondantestructures.pkl"))
	BoltzmanFactor=FF.UnpickleVariable(os.path.join(conf.PickledData, "Boltzman.pkl"))
	print "Start of the clustering process"
	method= "MiniBatchKMean"
	Clusters,CentroidStructure=OC.define_number_cluster(os.path.join("output", SVMlFile), Redundant, method, BoltzmanFactor, Probingconditions, rna)
	#Get Clusters from Pickled data
	#Clusters = FF.UnpickleVariable(os.path.join(conf.PickledData, "Clusters" + method + ".pkl"))

	# We should create an intermediate  file to be sure that  RNAeval  works!!
	with open("filecentroide", "w")as Ctrdfile:
		Ctrdfile.write(">RNA")
		print CentroidStructure
		for elem in (CentroidStructure.keys()):
			Ctrdfile.write("\n%s" % (CentroidStructure[elem]))
			print CentroidStructure[elem]
	Centroids_Energies = SF.ENERGY_VALUES_STRUCTURES("filecentroide", rna)
	print "Centroids_Energies", "\t", Centroids_Energies
	CT.ClustersDistributions(Clusters, os.path.join("output", "Clusters_details"),Probingconditions, int(conf.numberofsruct))


	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!Election of the best structures strating!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# Calculate cumulated  Boltzmaan Energy for each cluster

	CumulBE = {}
	CardinalConditions = {}


	CumulBE = CT.CumulatedConditionalBoltzmannbyCluster(Clusters,  ConditionalBoltzmaanProbability, int(conf.numberofsruct), Probingconditions)
	# epsilon= 1/(n+1)
	Epsilon = 1 / float(len(Clusters)+1)
	CardinalConditions = CT.GetCardinalConditions(Clusters,  ConditionalBoltzmaanProbability, Probingconditions,
		                                                             int(conf.numberofsruct), Epsilon)
	print 'verification pareto values'
	print 'Cardinal condition values', CardinalConditions.values()
	print 'Cumulated Boltzmann', "\t", CumulBE.values()
	#print 'mean_distance_in_a_given_cluster', "\t", E.values()
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Cluster analysis to elect the best cluster
	Dict = {}  # a dictionary that contains the three variables characterizing clusters (Criteria of optimization)
	for ClusterNumber in Clusters:
		Dict[ClusterNumber] = [CardinalConditions[ClusterNumber],
		                       CumulBE[ClusterNumber]]  #  only Two criterions are considered



	ListOptimalClusters = CT.Pareto(Dict)
	print "The elected clusters figuring in the Pareto front", ListOptimalClusters

	Dominatedclusters = [clusteri for clusteri in Clusters if clusteri not in ListOptimalClusters]
	Filealloptimalcentroides = os.path.join("Multiprobing/", str(Probingconditions)+ ".optimals")
	#VT.Drawvarna(Filealloptimalcentroides, ListOptimalClusters, CentroidStructure, conf.numberofsruct,
	#	         rna, Centroids_Energies, conf.SHAPEVis)

	print " IPANEMAP has been run successfully"
	print ("End of IPANEMAP analysis for a sampling of %s rna   % s structure for %s  conditions: %53f\t" % (file, int(conf.numberofsruct), len(Probingconditions), time.time() - startimebig))

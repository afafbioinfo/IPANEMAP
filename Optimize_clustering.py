import numpy as np
from collections import defaultdict
from sklearn.datasets import load_svmlight_file
from sklearn import cluster
import FileFunctions as FF, Clustering as CL, StructureFunctions as SF
import sys, math, re, os
from conf import loadConfig
import ClustersTrait as CT
from Progress import progress

def CumulatedBoltzmannsbyCluster(clusters, BoltzmanCluster, numberofsruct, constraintes):
    cBE = {}
    for ClusterNumber in clusters:
        l = 0.
        for structure in clusters[ClusterNumber]:
            ConditionNumber = int((structure - 1) / numberofsruct)
            StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
            l += BoltzmanCluster[constraintes[ConditionNumber]][StructureNumber]
        cBE[ClusterNumber] = l
    return cBE


def DefineNumberCluster(SVMLMatrix, Redundant, method, DM, BoltzmanFactor, Probingconditions, rna):
    conf = loadConfig()

    epsilon = 1  # Cetroid base pair distance threshold

    Cluster = defaultdict(lambda: defaultdict(CL.a))
    Clust = defaultdict(lambda: defaultdict(CL.a))
    CumulBE = defaultdict(lambda: defaultdict(CL.a))
    Centroids = defaultdict(lambda: defaultdict(CL.a))

    progress.StartTask("Initialization step")
    # Initialization step
    Cluster[1] = CL.MiniBatchKMeans(SVMLMatrix, 1)
    Centroids[1], BZ, X, Y, Z, IntradistanceStop = CT.CentroidBycluster(Cluster[1],
                                                                        os.path.join(conf.OutputFolder,"tmp",'OutputSamples' + str(conf.SampleSize), 'Samples.txt'),
                                                                        BoltzmanFactor,
                                                                        int(conf.SampleSize),
                                                                        Probingconditions, rna)
    CumulBE[1] = CumulatedBoltzmannsbyCluster(Cluster[1], BZ, int(conf.SampleSize), Probingconditions)
    #print  "***************************************verification bz", "Cluster  \t Centroids  \t CumulBE \t ", Centroids[1], CumulBE[1]
    progress.EndTask()
    for nb in range(2, 21):
        progress.StartTask("Clustering with %s clusters"%nb)
        progress.StartTask("Run MBKM")
        Clust[nb] = CL.MiniBatchKMeans(SVMLMatrix, nb)
        progress.EndTask()
        Centroids[nb], BZ, X, Y, Z, Intradistance = CT.CentroidBycluster(Clust[nb],
                                                                        os.path.join(conf.OutputFolder,"tmp",'OutputSamples' + str(conf.SampleSize), 'Samples.txt'),
                                                                         BoltzmanFactor,
                                                                         int(conf.SampleSize),
                                                                         Probingconditions, rna)
        CumulBE[nb] = CumulatedBoltzmannsbyCluster(Clust[nb], BZ, int(conf.SampleSize),
                                                   Probingconditions)

        lista = []

        '''
        ####***************************************************First crierion:
        if len([ elem for elem in IntradistanceStop if elem <= epsilon_intradist ] )!=0:
            print "************************************* Clustering done with ", nb ," as the optimal number of clusters using the first criterion  intradistance*********************************************************************"
            break
        # ************************************* second criterion
        '''
        for elem1 in Centroids[nb - 1].keys():
            rep = []
            '''
            print "distance to all elements"
            print "Ref \t i \t i+1 \t BPdist \t CumulatedBz i \t CumulatedBz i+1 \t  CumulatedBzdist"
            '''
            for elem2 in Centroids[nb].keys():
                rep.append((elem2, SF.DistanceTwoStructs(SF.BasePairsFromStruct(Centroids[nb - 1][elem1]),
                                                        SF.BasePairsFromStruct(Centroids[nb][elem2]))))

            minima = np.min([item[1] for item in rep])
            pos = [elem[0] for elem in rep if elem[1] == minima][0]

            l1 = CumulBE[nb - 1][elem1]
            l2 = CumulBE[nb][pos]
            # print "what s wrong!", l1,l2
            Dist = l1 - l2

            lista.append((minima, (l1, l2, Dist)))
        ########## The new criterion i about the existence of probable cluster
        Bzmepsilon = 0.3 * CumulBE[1][0]

        BP_All_probable_centroids = [BPdist for BPdist, Bzmandist in lista if Bzmandist[0] >= Bzmepsilon]
        progress.EndTask()

        if (len([elem for elem in Intradistance if elem <= epsilon]) != 0 or len(
                [distance for distance in BP_All_probable_centroids if distance <= epsilon]) == len(
                BP_All_probable_centroids)):
            FF.PickleVariable(Cluster[nb],  "Clusters" + method + ".pkl")
            progress.Print("Choosing %s as the optimal number of clusters"%nb)
            break

        # for the entire clusters while keeping redundancy
    return Clust[nb], Centroids[nb]

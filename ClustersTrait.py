from collections import defaultdict
import StructureFunctions as SF, FileFunctions as FF, VisualizationTools as VT
import numpy as np
import scipy
import time
import sys
from Progress import progress


def PointStructure(structure, numberofsruct):
    ConditionNumber = int((structure - 1) / numberofsruct)
    StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
    return ConditionNumber, StructureNumber


def CumulatedConditionalBoltzmannbyCluster(clusters, ConditionalBoltzman, numberofsruct, constraintes):
    cBE = {}
    for ClusterNumber in clusters:
        l = 0.
        for structure in clusters[ClusterNumber]:
            ConditionNumber, StructureNumber = PointStructure(structure, numberofsruct)
            l += ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
        cBE[ClusterNumber] = l
    return cBE


# *************************Get cardinal for conditions verifying Epsilon test
def GetCardinalConditions(clusters, ConditionalBoltzman, constraintes, numberofsruct, Epsilon):
    CardinalConditions = {}

    for ClusterNumber in clusters:
        # for each cluster, a sum overall boltzmann proba for a given condition is calculated
        # a condition counts for the cardinal of existing constrainte if the sum overall its strcutures is greater than epsilon[condition]
        l = defaultdict(lambda: 0)

        for structure in clusters[ClusterNumber]:
            ConditionNumber = int((structure - 1) / numberofsruct)
            StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
            # print "structure, ConditionNumber, StructureNumber ",structure, ConditionNumber, StructureNumber
            temp = ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
            l[ConditionNumber] += temp  # sum of Boltzm probabilites for a given condition
        for ConditionNumber in range(len(constraintes)):
            # print ClusterNumber ," : ", constraintes[ConditionNumber],":",l[ConditionNumber]
            if l[
                ConditionNumber] < Epsilon:  # if cumulated boltzmann is less than 1/(n+1), we do not consider the condition
                del l[ConditionNumber]
        CardinalConditions[ClusterNumber] = len(l)
    return CardinalConditions


# Cluster Diameters calculation
def ClustersDiameter(clusters, BPStructs):
    progress.StartTask("Computing cluster diameters")
    eliminated_clusters = []
    lista = []

    for ClusterNumber in clusters:
        if len(clusters[ClusterNumber]) > 1:  # not unique structure
            d = max([SF.DistanceTwoStructs(BPStructs[ClusterNumber][structure1],
                                           BPStructs[ClusterNumber][structure2])
                     for structure1 in clusters[ClusterNumber] for structure2 in clusters[ClusterNumber]])
        else:
            d = 0
            eliminated_clusters.append(ClusterNumber)
        lista.append(d)
    progress.EndTask()
    return lista, eliminated_clusters


def BasePairsbyCluster(clusters, Structures, Boltzmann, numberofsruct, constrainte):
    ClusterBoltzmann = defaultdict(lambda: defaultdict())
    ListBPbyCluster = defaultdict()
    ListBPbyStruct = defaultdict(lambda: defaultdict())
    Zcluster = defaultdict(lambda: defaultdict())
    BoltzmannOverPairsbyCluster = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
    for ClusterNumber in clusters:
        liste = []
        # Calculate Z
        acc = 0.
        for structure in clusters[ClusterNumber]:
            ConditionNumber, StructureNumber = PointStructure(structure, numberofsruct)
            acc += Boltzmann[constrainte[ConditionNumber]][StructureNumber]
        Zcluster[ClusterNumber] = acc

        list1 = []
        for structure in clusters[ClusterNumber]:
            ConditionNumber, StructureNumber = PointStructure(structure, numberofsruct)
            ClusterBoltzmann[constrainte[ConditionNumber]][StructureNumber] = Boltzmann[constrainte[ConditionNumber]][
                                                                                  StructureNumber] / float(Zcluster[ClusterNumber])
            for (i, j) in Structures[structure-1]:
                BoltzmannOverPairsbyCluster[ClusterNumber][i][j] += ClusterBoltzmann[constrainte[ConditionNumber]][
                    StructureNumber]

            ListBPbyStruct[ClusterNumber][structure] = Structures[structure-1]

        ListBPbyCluster[ClusterNumber] = list1
    return ListBPbyStruct, ListBPbyCluster, BoltzmannOverPairsbyCluster, ClusterBoltzmann


def ClustersDistances(clusters, Boltzmanprobabilities, ListBPbystrcut, numberofsruct, condition):
    Emd = defaultdict()
    for ClusterNumber in clusters:
        liste = []
        acc = 0.
        for structure1 in clusters[ClusterNumber]:
            for structure2 in clusters[ClusterNumber]:
                # TODO Change the conversion to  be considered as fucntion!!!
                ConditionNumber1 = int((structure1 - 1) / numberofsruct)
                StructureNumber1 = (structure1 - 1) - ConditionNumber1 * numberofsruct
                ConditionNumber2 = int((structure2 - 1) / numberofsruct)
                StructureNumber2 = (structure2 - 1) - ConditionNumber2 * numberofsruct
                val = (Boltzmanprobabilities[condition[ConditionNumber1]][StructureNumber1] *
                             Boltzmanprobabilities[condition[ConditionNumber2]][StructureNumber2] *
                             SF.DistanceTwoStructs(ListBPbystrcut[ClusterNumber][structure1], ListBPbystrcut[ClusterNumber][structure2]))
                #liste.append(val)
            acc += val
        # print "verify", ClusterNumber,2 * sum(liste) / float(len(clusters[ClusterNumber]) * (len(clusters[ClusterNumber]) - 1)),liste
        if len(clusters[ClusterNumber]) != 0:
            Emd[ClusterNumber] = 2. * acc / float(len(clusters[ClusterNumber]) * (len(clusters[ClusterNumber]) - 1))
        else:
            Emd[ClusterNumber] = 0.
    # print "Mean distance in the clusters:", Emd
    return Emd


def backtracking(W, BPp, rna, P_ss, i, j, pair):
    if i < j:
        # print i, j ,"what is the problem???",P_ss[j],W[i][j-1]
        if W[i][j] == P_ss[i] + W[i + 1][j]:
            # print i,j,W[i][j],P_ss[i],W[i+1][j]
            backtracking(W, BPp, rna, P_ss, i + 1, j, pair)
            return pair
        elif W[i][j] == P_ss[j] + W[i][j - 1]:
            # print i,j,W[i][j],P_ss[j],W[i][j-1]
            backtracking(W, BPp, rna, P_ss, i, j - 1, pair)
            return pair
        elif W[i][j] == (2 * BPp[i][j] + W[i + 1][j - 1]):
            pair.append((i, j))
            backtracking(W, BPp, rna, P_ss, i + 1, j - 1, pair)
            return pair
        else:
            for k in range(i, j):
                if W[i][j] == W[i][k] + W[k + 1][j]:
                    backtracking(W, BPp, rna, P_ss, i, k, pair)
                    backtracking(W, BPp, rna, P_ss, k + 1, j, pair)
                    return pair
    '''
    elif abs(i-j) <= 1:
        return pair
    raise Exception("Missing case in backtrack")
    '''

def MEA(BPp, rna):
    pair = []
    P_ss = defaultdict(lambda: 0)
    W = defaultdict(lambda: defaultdict(lambda: 0))
    for i in range(len(rna)):
        P_ss[i] = 1 - sum([BPp[min(i, j)][max(i, j)] for j in range(len(rna))])
    W = MEA_algo(P_ss, BPp, rna)
    endtime = time.time()
    # print backtracking(W,BPp,rna,P_ss,0,len(rna)-1)
    # print fromPairsToStruct(rna, backtracking(W,BPp,rna,P_ss,0,len(rna)-1,pair))
    MEAStruct = backtracking(W, BPp, rna, P_ss, 0, len(rna) - 1, pair)
    return SF.fromPairsToStruct(rna, MEAStruct), set(MEAStruct)


def MEA_algo(P_ss, BPp, rna):
    nb = len(rna)
    W = [[-sys.maxint for j in range(nb)] for i in range(nb)]
    for i in range(0, nb):
        W[i][i] = P_ss[i]
    for d in range(1, nb):
        for j in range(d, nb):
            i = j - d
            W[i][j] = max(W[i][j], P_ss[i] + W[i + 1][j])
            W[i][j] = max(W[i][j], P_ss[j] + W[i][j - 1])
            W[i][j] = max(W[i][j], 2 * BPp[i][j] + W[i + 1][j - 1])
            for k in range(i, j):
                W[i][j] = max(W[i][j], W[i][k] + W[k + 1][j])

        # W[i][j]=np.max([P_ss[i]+W[i+1][j],P_ss[j]+W[i][j-1],2*BPp[i][j]+W[i+1][j-1],np.max([W[i][k]+W[k+1][j] for k in range(i,j)])])
    return W


def CentroidBycluster(clusters, StructFile, Boltzmann, numberofsruct, constrainte, rna):
    progress.StartTask("Computing centroids")
    dim_clustering = len(clusters)
    E = defaultdict()
    mycentroid = defaultdict()
    Intradistance = []
    centroids = defaultdict(lambda: defaultdict())
    Myproba = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))

    ListStructures = [SF.BasePairsFromStruct(Struct) for Struct in FF.Parsefile(StructFile)]

    progress.StartTask("Gathering base pairs")
    ListBPbystructure, ListBP, Myproba, Boltzmancluster = BasePairsbyCluster(clusters, ListStructures, Boltzmann,
                                                                             numberofsruct,
                                                                             constrainte)
    # Eliminate cluster reporting one structure
    ListDiameters, Listeliminated_clusers = ClustersDiameter(clusters, ListBPbystructure)
    for elem in Listeliminated_clusers:
        del clusters[elem]
    progress.EndTask()

    progress.StartTask("Computing cluster distance distribution")
    E = ClustersDistances(clusters, Boltzmann, ListBPbystructure, numberofsruct, constrainte)
    progress.EndTask()
    progress.StartTask("Computing MEA centroids")
    for ClusterNumber in clusters:
        mycentroid[ClusterNumber], centroids[ClusterNumber] = MEA(Myproba[ClusterNumber], rna)
    progress.EndTask()

    MatriceDistanceCentroids = scipy.zeros([dim_clustering, dim_clustering])
    MatriceDistanceClustersEucld = scipy.zeros([dim_clustering, dim_clustering])
    for ClusterNumber in clusters.keys():
        for ClusterNumber2 in clusters.keys():
            if ClusterNumber2 > ClusterNumber:
                l = SF.DistanceTwoStructs(centroids[ClusterNumber], centroids[ClusterNumber2])
                #print "BP_centoid_distances", "\t", ClusterNumber, "\t", ClusterNumber2, "\t", l
                Intradistance.append(l)
                # print "distance between clusters comparing the centroide's distances",l
                MatriceDistanceCentroids[ClusterNumber][ClusterNumber2] = l
                MatriceDistanceCentroids[ClusterNumber2][ClusterNumber] = l
                # print "distance between clusters comparing the means distances", ClusterNumber, ClusterNumber2, np.abs(E[ClusterNumber]-E[ClusterNumber2]),np.sqrt(abs(pow(E[ClusterNumber],2)-pow(E[ClusterNumber2],2)))
                # print E
                l = np.sqrt(abs(pow(E[ClusterNumber], 2) - pow(E[ClusterNumber2], 2)))
                MatriceDistanceClustersEucld[ClusterNumber][ClusterNumber2] = l
                MatriceDistanceClustersEucld[ClusterNumber2][ClusterNumber] = l
            # print "distance between clusters compring the centroide's distances", ClusterNumber, ClusterNumber2, DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][listCentroidStructure[ClusterNumber][0]],ListBPbystrcut[ClusterNumber2][listCentroidStructure[ClusterNumber2][0]])
    # VT.plotDistanceClusters(MatriceDistanceCentroids, clusters, "blue", " Base pair distance between centroids")
    # VT.plotDistanceClusters(MatriceDistanceClustersEucld, clusters, "red", "Eucledian distance between structures")
    #print "BZ_distance_btw_clusters", "\t", E
    progress.EndTask()
    return mycentroid, Boltzmancluster, E, MatriceDistanceCentroids, ListDiameters, Intradistance


# count the occurence of present conditions in a given cluster
def SamplePresentInClusters(origine, occ, clusters, numberofsruct):
    for ClusterNumber in clusters:
        for StructureNumber in clusters[ClusterNumber]:
            origine[ClusterNumber][StructureNumber] = GetOrigineofStructure(StructureNumber, numberofsruct)
        # calculate occurence of conditions present whithin a cluster
    for ClusterNumber in clusters:
        for ConditionInCluster in origine[ClusterNumber].values():
            occ[ClusterNumber][ConditionInCluster] = origine[ClusterNumber].values().count(ConditionInCluster)

        return occ


# for a given structure characterized by a 'StructureNumber' return the condition represented by this structure
def GetOrigineofStructure(StructureNumber, numberofsruct):
    if (StructureNumber % numberofsruct != 0):
        return int(StructureNumber / numberofsruct) + 1
    else:
        return int(StructureNumber / numberofsruct)


def ClustersDistributions(clusters, Filename, constraintes, numberofsruct):
    origine = defaultdict(lambda: defaultdict(lambda: 0))
    occ = defaultdict(lambda: defaultdict(lambda: 0))
    it = defaultdict(lambda: 0)
    numberssamples = len(constraintes)
    o = open(Filename, "w")
    o.write("Cluster \t structures \t")
    for j in range(0, numberssamples):
        o.write("constraint %i = %s \t" % (j + 1, constraintes[j]))
    o.write("Number of structures \t Number of groups\n")
    occ = SamplePresentInClusters(origine, occ, clusters, numberofsruct)

    for elem in clusters:
        o.write("%i \t %s \t" % (elem + 1, clusters[elem]))
        for j in range(1, numberssamples + 1):
            if (occ[elem][j] != 0):
                it[elem] += 1
            o.write("%i\t" % (occ[elem][j]))
        o.write("%i\t%i\t" % (len(clusters[elem]), it[elem]))
        o.write("\n")
    o.write("Cluster(s) with  high number of  present conditions is(are) : %s" % (
        [v + 1 for v in it.keys() if it[v] == max(it.values())]))


# *************************!!!!!!!!!!!!!Pareto front*************************!!!!!!!!!!!!!!!!!!
def dominates(row, rowCandidate):
    return all(r >= rc for r, rc in zip(row, rowCandidate))


def Pareto(Dico):
    cleareCluster = []
    remaining = Dico.keys()

    while remaining:
        candidateKey = remaining.pop()
        candidate = Dico[candidateKey]

        new_remaining = []
        for other in remaining:
            if not dominates(candidate, Dico[other]):
                new_remaining.append(other)

        if not any(dominates(Dico[other], candidate) for other in new_remaining):
            cleareCluster.append(candidateKey)
        remaining = new_remaining
    # print len(remaining)
    # return cleareCluster,cleared
    return cleareCluster

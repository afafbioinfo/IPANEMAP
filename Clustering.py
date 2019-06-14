import numpy as np
from collections import defaultdict
from sklearn.datasets import load_svmlight_file
from sklearn import cluster
import FileFunctions as FF
import sys, math, re,os
import conf


def FilterCluster(cluster,redundant):# function to eliminate redundant elements
	return [elem for elem in cluster if elem not in redundant]
def a():
    return []
def MiniBatchKMeans(SVMLMatrix,n_clusters):
    clusters = defaultdict(a)
    #clusters_clean= defaultdict(a)
    X, y = np.array(load_svmlight_file(SVMLMatrix))
    #print X,y
    algorithm = cluster.MiniBatchKMeans(n_clusters=n_clusters)  # The cas eof 2 clusters is calculated
    algorithm.fit(X)
    if hasattr(algorithm, 'labels_'):
        y_pred = algorithm.labels_.astype(np.int)
    else:
        y_pred = algorithm.predict(X)
    for i in range(len(y_pred)):
        clusters[y_pred[i]].append(i + 1)
    #print clusters
    '''
    for elem in clusters:
        #print "clusters_properties", "\t", elem, "\t", len(clusters[elem])
        clusters_clean[elem] = FilterCluster(clusters[elem], Redundant)
    '''
    return clusters


#The DIANA implementation is extracted from Vienna-Package code and adapted to our treatment
class DIANA:
    @staticmethod
    def doClustering(dissmatrix, structs, maxDiameter, maxAverageDiameter):
        """
        Computes the DIANA clustering

        Args:
            maxDiameter = threshold for clustering
            maxAverageDiameter = "
        """

        structs = list(structs)
        # start DIANA
        # 1. initialization
        # dissMatrix = computeBasePairDistanceMatrix(structs)

        c_root = Cluster()
        c_root.structures.extend([x for x in range(0, len(structs))])

        # start the real DIANA algorithm.
        createClusterTree(c_root, dissmatrix, maxDiameter, maxAverageDiameter)
        # end DIANA

        clusters = convertTreeToListOfClusters(c_root, structs)

        return clusters



class Cluster:
    """
    diana internal tree-like cluster structure.
    """

    def __init__(self):
        self.structures = []
        self.representative = ""
        self.childNodes = []


def maxDiameterInLeafNodes(clusterTree, dissMatrix):
    maxDiameter = -1
    maxCluster = None
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            dia, cluster = maxDiameterInLeafNodes(cn, dissMatrix)
            if dia > maxDiameter:
                maxDiameter = dia
                maxCluster = cluster
    else:
        maxDiameter = diameter(clusterTree.structures, dissMatrix)
        maxCluster = clusterTree

    return maxDiameter, maxCluster


def averageDiameterInLeafNodes(clusterTree, dissMatrix):
    """
    also known as Divisive Coefficient (DC).
    """
    sumDiameters = 0
    numberOfClusters = 0
    stackChildNodes = []
    stackChildNodes.append(clusterTree)
    while (len(stackChildNodes) > 0):
        cn = stackChildNodes.pop()
        if len(cn.childNodes) == 0:
            # is leaf node/cluster --> sum
            sumDiameters += diameter(cn.structures, dissMatrix)
            numberOfClusters += 1
        else:
            for c in cn.childNodes:
                stackChildNodes.append(c)

    averageDiameter = sumDiameters / float(numberOfClusters)
    return averageDiameter


def createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold):
    """
    The core of the diana algorithm (recursive function).
    c_root = the rootnode of the clusterTree. It contains the main cluster as childnode.
    """
    # 1. select cluster with the largest diameter from all leafnodes.
    maxDiameter, c_m = maxDiameterInLeafNodes(c_root, dissMatrix)
    dc = averageDiameterInLeafNodes(c_root, dissMatrix)
    if dc < maxAverageDiameterThreshold:
        return
    if maxDiameterThreshold >= 0:
        if maxDiameter <= maxDiameterThreshold:
            return

    if c_m == None:
        return
    if len(c_m.structures) > 1:
        # 2. object with highest dissimilarity to all others defines the new cluster (c_newA).
        o = objectIndexWithHighestDissimilarity(c_m.structures, dissMatrix)
        c_newA = Cluster()
        c_newA.structures.append(o)
        # B contains all other elements
        c_newB = Cluster()
        c_newB.structures.extend([x for x in c_m.structures if x != o])

        # 3. select best elements for each cluster (move closest elements from B to A).
        # for each object outside the splinter group:
        positiveDiExists = True
        while (positiveDiExists):
            positiveDiExists = False
            largestD_i = 0.0
            bestElement = None
            for j in c_newB.structures:
                d_i = averageDissimilarity(j, c_newB.structures, dissMatrix) - averageDissimilarity(j,
                                                                                                    c_newA.structures,
                                                                                                    dissMatrix)
                if (d_i > largestD_i):
                    largestD_i = d_i
                    bestElement = j
            if (bestElement != None) & (len(c_newB.structures) > 1):
                positiveDiExists = True
                c_newA.structures.append(bestElement)
                c_newB.structures.remove(bestElement)

        c_m.childNodes.append(c_newB)
        c_m.childNodes.append(c_newA)

        # build the subtrees for each child
        createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold)


def objectIndexWithHighestDissimilarity(cluster, dissMatrix):
    maxAvgDiss = 0.0
    maxIndex = -1
    for i1 in cluster:
        avgDiss = averageDissimilarity(i1, cluster, dissMatrix)
        if (avgDiss > maxAvgDiss):
            maxAvgDiss = avgDiss
            maxIndex = i1
    return maxIndex


def averageDissimilarity(index, cluster, dissMatrix):
    """
    computes the average dissimilarity
    index is a structure index, which is conform the the index in the dissMatrix.
    cluster is a list of structure indices.
    """
    avgDiss = 0.0
    subtract = 0.0
    if index in cluster:
        subtract = 1.0
    for i in cluster:
        avgDiss += float(dissMatrix[index][i])
    NumberOfElementsWithoutIndexElement = float(len(cluster) - subtract)
    if NumberOfElementsWithoutIndexElement > 0.0:
        avgDiss = avgDiss / NumberOfElementsWithoutIndexElement
    else:
        avgDiss = 0.0
    return avgDiss


def objectIndexWithHighestDissimilarity(cluster, dissMatrix):
    maxAvgDiss = 0.0
    maxIndex = -1
    for i1 in cluster:
        avgDiss = averageDissimilarity(i1, cluster, dissMatrix)
        if (avgDiss > maxAvgDiss):
            maxAvgDiss = avgDiss
            maxIndex = i1
    return maxIndex


def convertTreeToListOfClusters(clusterTree, structs, clusterList=[]):
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            convertTreeToListOfClusters(cn, structs, clusterList)
        return clusterList
    else:
        cluster = [structs[c] for c in clusterTree.structures]
        clusterList.append(cluster)


def diameter(cluster, dissMatrix):
    """
    computes the largest basepairdistance between all pairs of structures in the cluster.
    """
    maxDist = 0
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):
            dist = dissMatrix[cluster[i]][cluster[j]]
            if dist > maxDist:
                maxDist = dist
    return maxDist




    @staticmethod
    def printClusters(clusters):
        cid = 0
        for c in clusters:
            cid += 1
            print "ClusterID:", cid
            for s in c:
                print s

import numpy as np
from collections import defaultdict
from sklearn.datasets import load_svmlight_file
from sklearn import cluster
import FileFunctions as FF, Clustering as CL, StructureFunctions as SF
import sys, math, re,os
import conf
import ClustersTrait as CT


def CumulatedBoltzmannsbyCluster(clusters, BoltzmanCluster, numberofsruct, constraintes):
    cBE = {}
    #print "test", clusters #,BoltzmanCluster
    for ClusterNumber in clusters:

	    l = 0.
	    for structure in clusters[ClusterNumber]:
		    ConditionNumber = int((structure - 1) / numberofsruct)
		    StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
		    l += BoltzmanCluster[constraintes[ConditionNumber]][StructureNumber]
            #print structure,  "llll",StructureNumber ,BoltzmanCluster[constraintes[ConditionNumber]][StructureNumber]
	    cBE[ClusterNumber] = l
    return cBE


def define_number_cluster(SVMLMatrix, Redundant,method, BoltzmaanFactor,Probingconditions,rna):

    epsilon = 1 # Cetroid base pair distance threshold

    Cluster= defaultdict(lambda: defaultdict(CL.a))
    Clust = defaultdict(lambda: defaultdict(CL.a))
    CumulBE= defaultdict(lambda: defaultdict(CL.a))
    Centroids= defaultdict(lambda: defaultdict(CL.a))

    # Initialization step
    Cluster[1]=CL.MiniBatchKMeans(SVMLMatrix, 1)
    Centroids[1],BZ,X,Y,Z,IntradistanceStop = CT.CentroidBycluster(Cluster[1],os.path.join('OutputSamples' + str(conf.numberofsruct), 'Samples.txt'),
                                         BoltzmaanFactor,
                                         int(conf.numberofsruct),
                                        Probingconditions, rna)
    CumulBE[1] = CumulatedBoltzmannsbyCluster(Cluster[1], BZ, int(conf.numberofsruct),Probingconditions)
    print  "***************************************verification bz","Cluster  \t Centroids  \t CumulBE \t ", Centroids[1], CumulBE[1]

    for nb in range(2,21):
        print "************ Start clustering with", nb, " clusters*****************"
        print "nbr_clusters \t",nb
        Clust[nb] = CL.MiniBatchKMeans(SVMLMatrix, nb)
        Centroids[nb],BZ,X,Y,Z,Intradistance=CT.CentroidBycluster(Clust[nb], os.path.join('OutputSamples' + str(conf.numberofsruct),'Samples.txt'),
                                                                                                       BoltzmaanFactor,
                                                                                                       int(conf.numberofsruct),
                                                                                              Probingconditions,rna)
        CumulBE[nb] = CumulatedBoltzmannsbyCluster( Clust[nb], BZ, int(conf.numberofsruct),
                                                        Probingconditions)
        
        lista=[]
        
        '''
        ####***************************************************First crierion:
        if len([ elem for elem in IntradistanceStop if elem <= epsilon_intradist ] )!=0:
            print "************************************* Clustering done with ", nb ," as the optimal number of clusters using the first criterion  intradistance*********************************************************************"
            break
        # ************************************* second criterion
        '''
        for elem1 in Centroids[nb-1].keys():
            #print "YOUPIIIII",k1
            rep = []
            '''
            print "distance to all elements"
            print "Ref \t i \t i+1 \t BPdist \t CumulatedBz i \t CumulatedBz i+1 \t  CumulatedBzdist"
            '''
            for elem2 in Centroids[nb].keys():
                rep.append((elem2, SF.DistanceTwoBPlist(SF.ListBasePairsFromStruct(Centroids[nb-1][elem1]),SF.ListBasePairsFromStruct(Centroids[nb][elem2]))))
            

            minima = np.min([item[1] for item in rep])
            pos = [elem[0] for elem in rep if elem[1] == minima][0]
            

            l1 = CumulBE[nb-1][elem1]
            l2 = CumulBE[nb][pos]
        #print "what s wrong!", l1,l2
            Dist = l1- l2

            lista.append((minima, (l1,l2,Dist)))
        # to emove to count for the condition
        # Base pair distance is less then epsilon detected in all clusters from i + Boltzmann distance less then Bzmepsilon for threshhold number of clusters from step i to step i+1.4
        # distance intraclusters and inter-clusters cheeck
        '''
        if len([BPdist for BPdist,Bzmandist in lista if BPdist<= epsilon]) == nb-1 and len([Bzmandist for BPdist,Bzmandist in lista if
                                                                                        abs(Bzmandist) <= Bzmepsilon ]) == nb-1: # and nb-1 != 1:
        '''
        '''
        second criterion
        len([Bzmandist for BPdist,Bzmandist in lista if
                                                                                        abs(Bzmandist) <= Bzmepsilon ]) == nb-1) or (nb >2 and len([distance for distance in BP_All_probable_centroids if distance <= epsilon]) ==len(BP_All_probable_centroids)):
        '''
        ########## The new criterion i about the existance of probable cluster
        Bzmepsilon=0.3*CumulBE[1][0]

        BP_All_probable_centroids=[BPdist for BPdist, Bzmandist in lista if  Bzmandist[0]>=Bzmepsilon]
        
        if ( len([elem for elem in Intradistance if elem <=epsilon])!=0  or  len([distance for distance in BP_All_probable_centroids if distance <= epsilon]) ==len(BP_All_probable_centroids)):

            FF.PickleVariable(Cluster[nb], os.path.join(conf.PickledData, "Clusters" + method + ".pkl"))
            print "************************************* Clustering done with ", nb ," as the optimal number of clusters*********************************************************************"
            break

        # for the entire clusters while keeping redundancy
    return Clust[nb],Centroids[nb]

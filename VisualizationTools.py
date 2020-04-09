import conf as CF, FileFunctions as FF, StructureFunctions as SF
import scipy, numpy as np, os
import subprocess
from itertools import cycle


def ThreeD_MatDistance_Boltzmann(MatDist, Klust, Boltzmannprobabilty, numberofsruct, constrainte, MFEs):
    
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    fig_handle = plt.figure()
    ax1 = fig_handle.add_subplot(111, projection='3d')
    # plot 3 d
    elem = [i for i in range(len(MatDist))]
    D = scipy.zeros([len(MatDist), len(MatDist)])
    for i in range(len(MatDist)):
        for j in range(len(MatDist)):
            D[i][j] = MatDist[i][j]
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    plt.subplots_adjust(bottom=0.1)
    for k, col in zip(range(len(Klust)), colors):
        pos = -1
        for i in Klust[k]:
            pos += 1
            if i - 1 < numberofsruct * (len(constrainte) - 1):
                ConditionNumber = int((i - 1) / numberofsruct)
                StructureNumber = i - 1 - ConditionNumber * numberofsruct
                ax1.scatter(coords[i - 1, 0], coords[i - 1, 1],
                            Boltzmannprobabilty[constrainte[ConditionNumber]][StructureNumber], c=col, marker='.')
            else:
                ConditionNumber = len(constrainte) - 1
                StructureNumber = i - 1 - ConditionNumber * numberofsruct
                ax1.text(coords[i - 1, 0], coords[i - 1, 1],
                         Boltzmannprobabilty[constrainte[ConditionNumber]][StructureNumber],
                         '*%s' % (MFEs[StructureNumber]), color=col)
    ax1.set_xlabel('MDS axis1')
    ax1.set_ylabel('MDS axis2')
    ax1.set_zlabel('Boltzmann probability')
    ax1.set_title(' Secondary Structures Multidimensional Scaling ')
    fig_handle.savefig('output_plots/MDS_Structures_MFES.svg')

    return 0

# To eminitae the MFEs structures
def FilterClusters(Klust, lenconst):
    for C in range(len(Klust)):
        # to eliminate MFEs
        Klust[C] = [v for v in Klust[C] if v < lenconst]
        # if the cluster becomes empty , delete it
        if Klust[C] == ' ':
            del (Klust[C])
    return Klust

def threedcentoids(MatDist, Centroids_Energies, ListDiameters):
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    fig_handle = plt.figure()
    ax = fig_handle.add_subplot(111, projection='3d')

    D = scipy.zeros([len(MatDist), len(MatDist)])
    for i in range(len(MatDist)):
        for j in range(len(MatDist)):
            D[i][j] = MatDist[i][j]
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    plt.subplots_adjust(bottom=0.1)
    for k, col in zip(range(len(MatDist)), colors):
        ax.scatter(coords[k, 0], coords[k, 1], Centroids_Energies[k], c=col, marker='.')
        # ax.add_patch(mpatches.Circle((coords[k, 0], coords[k, 1]),ListDiameters/2,color=col,edgecolor="black"))
        ax.text(coords[k, 0], coords[k, 1], Centroids_Energies[k], 'C%s' % (k + 1), color=col)

    ax.set_xlabel('MDS axis1 (Base pairs distance)')
    ax.set_ylabel('MDS axis2')
    ax.set_zlabel('Boltzmann energy Centroids')
    ax.set_title('Clusters distances with centoid s Boltzmann energies')
    fig_handle.savefig('centroids_distribution.svg')
    
def plotDistanceClusters(D, clusters, coloro, title):
    Dic = {}
    for elem in clusters:
        Dic[elem] = elem
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    # plot results
    fig = plt.figure()
    plt.subplots_adjust(bottom=0.1)
    plt.scatter(coords[:, 0], coords[:, 1], marker='o')
    for label, x, y in zip(Dic.values(), coords[:, 0], coords[:, 1]):
        plt.annotate(
            label,
            xy=(x, y), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc=coloro, alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    fig.savefig('Distance_clusters_' + title + '.png')

def plotPareto(paretoPoints, dominatedPoints):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    dp = np.array(list(dominatedPoints))
    pp = np.array(list(paretoPoints))
    print(pp.shape, dp.shape)
    ax.scatter(dp[:, 0], dp[:, 1], dp[:, 2])
    ax.scatter(pp[:, 0], pp[:, 1], pp[:, 2], color='red')

    import matplotlib.tri as mtri
    triang = mtri.Triangulation(pp[:, 0], pp[:, 1])
    ax.plot_trisurf(triang, pp[:, 2], color='mediumvioletred')
    plt.show()

def plotPairs(lista, n):
    fig = PLT.figure()
    x = [elem[0] for elem in lista]
    y = [elem[1] for elem in lista]
    z = [elem[2] / float(n) for elem in lista]
    gridsize = 60
    PLT.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
    PLT.axis([min(x) - 1, max(x) + 1, min(y) - 1, max(y) + 1])
    cb = PLT.colorbar()
    cb.set_label('Probability value  in all  optimal centroids')

def plotClustercBECard(clusternumber, cBE, cardinal, xlabelo, ylabelo, output):
    Labels = clusternumber
    X = cBE
    Y = cardinal
    fig = plt.figure()
    fig.suptitle('Pareto front for clusters', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.set_title('Cluster distribution')
    ax.set_xlabel(xlabelo)
    ax.set_ylabel(ylabelo)
    ax.set_ylim(bottom=0)
    plt.axis([0, max(X) + 1, 0, max(Y) + np.mean(Y)])
    ax.grid(True)
    for i in Labels:
        ax.text(X[i] + 0.2, Y[i], i + 1, fontsize=10, horizontalalignment='center', color='b')
    plt.plot(X, Y, 'r*')
    fig.savefig(output)

def plotsecodnarystructures(rnaString, lista, lista2, n, reactivities):
    fig = plt.figure()
    fg, ax = plt.subplots(1, 1)
    import matplotlib as mp

    min_val = 0
    max_val = 1
    my_cmap = cm.get_cmap('Greys')
    norm = matplotlib.colors.Normalize(min_val, max_val)  #
    x = [elem[0] for elem in lista]
    y = [elem[1] for elem in lista]
    z = [elem[2] / float(n) for elem in lista]

    x2 = [elem[0] for elem in lista2]
    y2 = [elem[1] for elem in lista2]
    z2 = [elem[2] / float(n) for elem in lista2]

    rna = list(rnaString)
    cc = ["black" for i in range(len(rna))]
    for elem in range(len(rna)):
        if float(reactivities[elem]) < 0.2:
            cc[elem] = "black"  # "#00509d"#blue
        if float(reactivities[elem]) >= 0.2 and float(reactivities[elem]) < 0.4:
            cc[elem] = "gray"  # "#00c200"#green
        if float(reactivities[elem]) >= 0.4 and float(reactivities[elem]) < 0.7:
            cc[elem] = "mediumvioletred"  # "#f28f00"#yellow
        if float(reactivities[elem]) >= 0.7:
            cc[elem] = "plum"  # "#f20000"#red

    p0 = Rectangle((0, 0), 1, 1, fc="black")
    p1 = Rectangle((0, 0), 1, 1, fc="gray")
    p2 = Rectangle((0, 0), 1, 1, fc="mediumvioletred")
    p3 = Rectangle((0, 0), 1, 1, fc="plum")
    ax.legend([p0, p1, p2, p3], ["Reactivity <0.2", "0.2< <0.4", "0.4<  <0.7", ">0.7"])

    pac = [mpatches.Arc([x[i] + 0.5 + (y[i] - x[i] - 1) / float(2), 0], y[i] - x[i] - 1, (y[i] - x[i] - 1) / float(2),
                        angle=0, theta1=0, theta2=180, color=my_cmap(norm(z[i])), linewidth=1) for i in
           range(len(x))]  # linestyle='dotted', linestyle='dashed
    pac2 = [mpatches.Arc([x2[i] + 0.5 + (y2[i] - x2[i] - 1) / float(2), 0], y2[i] - x2[i] - 1,
                         (y2[i] - x2[i] - 1) / float(2), angle=0, theta1=180, theta2=360, color=my_cmap(norm(z[i])),
                         linewidth=1) for i in range(len(x2))]
    for arc, arc2 in zip(pac, pac2):
        ax.add_patch(arc)
        ax.add_patch(arc2)
    cmmapable = cm.ScalarMappable(norm, my_cmap)
    cmmapable.set_array(range(min_val, max_val))
    colorbar(cmmapable, fraction=0.046, pad=0.04, ticks=[0, 0.5, 1])

    fontProp = mp.font_manager.FontProperties(family="monospace", style="normal", weight="bold", size="8")
    ax.axis([0, max(x) + 20, -max(y) / 3, max(y) / 3])
    for i in range(len(rna)):
        nuc = rna[i]
        ax.add_patch(
            mpatches.Circle((i + 0.5, 0), 0.5, color=cc[i], edgecolor="black"))  # circle at center (x,y), radius 0.5
        ax.annotate(nuc, (i + 0.5, 0), color='white', weight='bold', fontsize=6, ha='center', va='center')
        ax.annotate(i + 1, (i + 0.5, -1), color='black', weight='bold', fontsize=6, ha='center', va='center')
    ax.set_aspect("equal")
    ax.get_yaxis().set_visible(False)
    fg.canvas.draw()

    ax.set_title('Combined optimal centroids ')

    # plt.show()
    plt.savefig("res.eps", format='eps', dpi=1000)
    fig.savefig('Arcs_structures.svg')

def plotClusteringDistribution(lenconstraint, Folder_name, Lenrna):
    D = scipy.zeros([lenconstraint, lenconstraint])
    Dic, B = SF.Load_Probabilities(Folder_name)

    # calculate the Eucledian distance between different matrix
    D = SF.Eucledian_distance(B, Lenrna)
    # D = SF.Absolute_distance(B, Lenrna)
    # tril for lower triangular matrix
    print Dic.values()
    print D
    # Clustering process with th plot
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    # plot results

    fig = plt.figure()
    plt.subplots_adjust(bottom=0.1)
    plt.scatter(coords[:, 0], coords[:, 1], marker='o')
    for label, x, y in zip(Dic.values(), coords[:, 0], coords[:, 1]):
        plt.annotate(
            label,
            xy=(x, y), xytext=(0, 20),
            textcoords='offset points', ha='left', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='green', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    # plt.show()
    fig.savefig('Euclidian_distance_dot_plot_Matrix.png')

COLOR_MAP = ' -colorMapStyle "$-0.5:#A0A0A0;$-0.499:#FFFFFF;$0:#FFFFFF;3:#FF0000"'

def drawStructure(Sequence, Structure, Shapefile, OutFile):
    conf = CF.loadConfig()
    cmopt = ""
    #print "shape",Shapefile
    if os.path.isfile(Shapefile):
        vals = FF.parseReactivityfile(Shapefile)
        cmopt = ' -colorMap "' + ";".join(["%.3f" % float(v) for v in vals]) + '"' + COLOR_MAP
    dummyout = os.path.join(conf.OutputFolder, "tmp", "varnamsg.txt")
    cmd = 'java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -bpStyle simple -sequenceDBN "%s" -structureDBN "%s" '%(Sequence, Structure) + cmopt + ' -algorithm line -o ' + OutFile
    #print cmd
    subprocess.call(cmd, stdin=None, stdout=open(dummyout, 'wb'),
                    stderr=open(dummyout, 'w'), shell=True)
    
def Convert2DDict_npArray(dic):
    return np.array([[dic[i][j] for j in sorted(dic[i])] for i in sorted(dic)])


def HeatMapplot(Distance, labels, ConvertDist):
    # Plot it out
    fig, ax = plt.subplots()
    if ConvertDist == 'True':
        nba_sort = Convert2DDict_npArray(Distance)

    else:
        nba_sort = Distance
    print 'conversion done'
    heatmap = ax.pcolor(nba_sort, cmap=plt.cm.Blues, alpha=1)  # alpha float (0.0 transparent through 1.0 opaque)

    # Format
    fig = plt.gcf()
    fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False


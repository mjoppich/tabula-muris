import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

import os,sys

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))
sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/nameConvert/nameconvert/")))

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE
from itertools import chain, combinations


from stores import UniprotStore, GeneIdentity
from collections import defaultdict
import re
import io
from natsort import natsorted

def perpareExpNames():

    """
    DROPLET CELLS
    """

    dropletAortaEpi = []
    dropletLiverEpi = ["Liver.duct epithelial cell"]
    dropletLungEpi = []
    dropletKidneyEpi = []
    dropletBladderEpi = []

    dropletAortaEndo = ["Heart_and_Aorta.endothelial cell"]
    dropletLiverEndo = ["Liver.endothelial cell of hepatic sinusoid"]
    dropletLungEndo = ["Lung.lung endothelial cell"]
    dropletKidneyEndo = ["Kidney.kidney capillary endothelial cell"]
    dropletBladderEndo = ["Bladder.endothelial cell"]


    dropletAortaCir = []
    dropletLiverCir = ["Liver.leukocyte"]
    dropletLungCir = ["Lung.myeloid cell", "Lung.B cell", "Lung.natural killer cell", "Lung.T cell", "Lung.non-classical monocyte", "Lung.leukocyte", "Lung.classical monocyte"]
    dropletKidneyCir = ["Kidney.leukocyte"]
    dropletBladderCir = ["Bladder.leukocyte"]

    """
    FACS CELLS
    """

    facsAortaEndo = ["Aorta.endothelial cell"]
    facsAortaEpi = []
    facsAortaCir = ["Aorta.erythrocyte"]

    facsLungEndo = ["Lung.lung endothelial cell"]
    facsLungEpi = ["Lung.epithelial cell of lung"]
    facsLungCir = ["Lung.myeloid cell", "Lung.monocyte", "Lung.classical monocyte", "Lung.B cell", "Lung.T cell", "Lung.natural killer cell", "Lung.leukocyte"]

    facsLintEndo = []
    facsLintEpi = ["Large_Intestine.epithelial cell of large intestine"]
    facsLintCir = []

    facsKidneyEndo = ["Kidney.endothelial cell"]
    facsKidneyEpi = []
    facsKidneyCir = ["Kidney.leukocyte"]

    facsBladderEndo = []
    facsBladderEpi = []
    facsBladderCir = []

    facsLiverEndo = ["Liver.endothelial cell of hepatic sinusoid"]
    facsLiverEpi = []
    facsLiverCir = ["Liver.natural killer cell", "Liver.B cell"]

    facsBrainEndo = ["Brain_Non-Myeloid.endothelial cell"]
    facsBrainEpi = ["Brain_Non-Myeloid.brain pericyte"]
    facsBrainCir = []

    """
    WELLS CELLS
    """
    wellsAortaEndo = ["Vascular endothelial cell(Neonatal-Heart)",
    "Endothelial cell_Igfbp5 high(Neonatal-Heart)",
    "Endothelial cell_Eln high(Neonatal-Heart)",
    "Endothelial cell_Enpp2 high(Neonatal-Heart)"]
    wellsLungEndo = ["Endothelial cell_Kdr high(Lung)",
    "Endothelial cell_Tmem100 high(Lung)",
    "Endothelial cells_Vwf high(Lung)"]
    wellsLintEndo = []
    wellsKidneyEndo = ["Endothelial cell(Kidney)",
    "Fenestrated endothelial cell_Plvap high(Kidney)",
    "Fenestrated endothelial cell_Tm4sf1 high(Kidney)"]
    wellsBladderEndo = ["Vascular endothelial cell(Bladder)",
    "Endothelial cell_Ly6c1 high(Bladder)"]
    wellsLiverEndo = ["Endothelial cell(Liver)"]

    wellsAortaEpi = ["Epithelial cell(Neonatal-Heart)"]
    wellsLungEpi = []
    wellsLintEpi = ["Columnar epithelium(Small-Intestine)",
    "Epithelial cell_Kcne3 high(Small-Intestine)",
    "Epithelial cell_Sh2d6 high(Small-Intestine)"]
    wellsKidneyEpi = ["Epithelial cell_Cryab high(Kidney)"]
    wellsBladderEpi = ["Basal epithelial cell(Bladder)",
    "Epithelial cell_Upk3a high(Bladder)",
    "Epithelial cell_Gm23935 high(Bladder)"]
    wellsLiverEpi = ["Epithelial cell(Liver)",
    "Epithelia cell_Spp1 high(Liver)"]

    wellsAortaCir = ["Neutrophil_Ngp high(Neonatal-Heart)",
    "Neutrophil_Retnlg high(Neonatal-Heart)"]
    wellsLungCir = ["B Cell(Lung)",
    "NK Cell(Lung)",
    "Dividing T cells(Lung)",
    "Eosinophil granulocyte(Lung)",
    "Neutrophil granulocyte(Lung)",
    "T Cell_Cd8b1 high(Lung)",
    "Basophil(Lung)"]

    wellsLintCir = ["T cell(Fetal_Intestine)",
    "B cell_Igkv12-46 high(Small-Intestine)",
    "B cell_Ms4a1 high(Small-Intestine)",
    "T cell_Icos high(Small-Intestine)",
    "B cell_Jchain high(Small-Intestine)",
    "T cell_Cd7 high(Small-Intestine)",
    "T cell_Ccl5 high(Small-Intestine)",
    "T cell_Ms4a4b high(Small-Intestine)",
    "B cell_Ighd high(Small-Intestine)",
    "Erythroblast(Small-Intestine)"]
    wellsKidneyCir = ["B cell(Kidney)",
    "T cell(Kidney)",
    "Neutrophil progenitor_S100a8 high(Kidney)"]
    wellsBladderCir = ["NK cell(Bladder)"]
    wellsLiverCir = ["B cell_Jchain high(Liver)",
    "Granulocyte(Liver)",
    "T cell_Trbc2 high(Liver)",
    "Erythroblast_Hbb-bt high(Liver)",
    "T cell_Gzma high(Liver)",
    "Erythroblast_Hbb-bs high(Liver)",
    "B cell_Fcmr high(Liver)",
    "Neutrophil_Ngp high(Liver)"]



    method2lists = defaultdict(lambda: defaultdict(lambda: dict()))

    """
    DROPLET
    """

    method2lists["droplet"]["aorta"]["epi"] = dropletAortaEpi
    method2lists["droplet"]["liver"]["epi"] = dropletLiverEpi
    method2lists["droplet"]["lung"]["epi"] = dropletLungEpi
    method2lists["droplet"]["kidney"]["epi"] = dropletKidneyEpi
    method2lists["droplet"]["bladder"]["epi"] = dropletBladderEpi

    method2lists["droplet"]["aorta"]["endo"] = dropletAortaEndo
    method2lists["droplet"]["liver"]["endo"] = dropletLiverEndo
    method2lists["droplet"]["lung"]["endo"] = dropletLungEndo
    method2lists["droplet"]["kidney"]["endo"] = dropletKidneyEndo
    method2lists["droplet"]["bladder"]["endo"] = dropletBladderEndo


    method2lists["droplet"]["aorta"]["cir"] = dropletAortaCir
    method2lists["droplet"]["liver"]["cir"] = dropletLiverCir
    method2lists["droplet"]["lung"]["cir"] = dropletLungCir
    method2lists["droplet"]["kidney"]["cir"] = dropletKidneyCir
    method2lists["droplet"]["bladder"]["cir"] = dropletBladderCir

    """
    FACS
    """
    method2lists["facs"]["aorta"]["epi"] = facsAortaEpi
    method2lists["facs"]["liver"]["epi"] = facsLiverEpi
    method2lists["facs"]["lung"]["epi"] = facsLungEpi
    method2lists["facs"]["kidney"]["epi"] = facsKidneyEpi
    method2lists["facs"]["bladder"]["epi"] = facsBladderEpi
    method2lists["facs"]["brain"]["epi"] = facsBrainEpi

    method2lists["facs"]["aorta"]["endo"] = facsAortaEndo
    method2lists["facs"]["liver"]["endo"] = facsLiverEndo
    method2lists["facs"]["lung"]["endo"] = facsLungEndo
    method2lists["facs"]["kidney"]["endo"] = facsKidneyEndo
    method2lists["facs"]["bladder"]["endo"] = facsBladderEndo
    method2lists["facs"]["brain"]["endo"] = facsBrainEndo


    method2lists["facs"]["aorta"]["cir"] = facsAortaCir
    method2lists["facs"]["liver"]["cir"] = facsLiverCir
    method2lists["facs"]["lung"]["cir"] = facsLungCir
    method2lists["facs"]["kidney"]["cir"] = facsKidneyCir
    method2lists["facs"]["bladder"]["cir"] = facsBladderCir
    method2lists["facs"]["brain"]["cir"] = facsBrainCir

    """
    WELLS
    """

    method2lists["wells"]["aorta"]["epi"] = wellsAortaEpi
    method2lists["wells"]["liver"]["epi"] = wellsLiverEpi
    method2lists["wells"]["lung"]["epi"] = wellsLungEpi
    method2lists["wells"]["kidney"]["epi"] = wellsKidneyEpi
    method2lists["wells"]["bladder"]["epi"] = wellsBladderEpi
    method2lists["wells"]["lint"]["epi"] = wellsLintEpi

    method2lists["wells"]["aorta"]["endo"] = wellsAortaEndo
    method2lists["wells"]["liver"]["endo"] = wellsLiverEndo
    method2lists["wells"]["lung"]["endo"] = wellsLungEndo
    method2lists["wells"]["kidney"]["endo"] = wellsKidneyEndo
    method2lists["wells"]["bladder"]["endo"] = wellsBladderEndo
    method2lists["wells"]["lint"]["endo"] = wellsLintEndo


    method2lists["wells"]["aorta"]["cir"] = wellsAortaCir
    method2lists["wells"]["liver"]["cir"] = wellsLiverCir
    method2lists["wells"]["lung"]["cir"] = wellsLungCir
    method2lists["wells"]["kidney"]["cir"] = wellsKidneyCir
    method2lists["wells"]["bladder"]["cir"] = wellsBladderCir
    method2lists["wells"]["lint"]["cir"] = wellsLintCir

    for method in method2lists:
        for tissue in method2lists[method]:
            for ltype in method2lists[method][tissue]:
                new_names = [x.replace(" ", "_") for x in method2lists[method][tissue][ltype]]
                method2lists[method][tissue][ltype] = new_names


    tname2cats = defaultdict(lambda: dict())

    for method in method2lists:
        for tissue in method2lists[method]:
            for ltype in method2lists[method][tissue]:
                for tname in method2lists[method][tissue][ltype]:
                    tname2cats[(tname, method)] = (ltype, tissue, method)

    return method2lists, tname2cats

def makeSubcellularLocations(subcellStr):

    if subcellStr == None:
        return set()

    sublocstr = subcellStr.replace("SUBCELLULAR LOCATION: ", "")

    asubloc = [x.strip() for x in sublocstr.split(";")]

    asubloc = [re.sub(r'\{[^}]*\}', '', x) for x in asubloc]
    asubloc = [re.sub(r'Note\=.*', '', x) for x in asubloc]

    asublocs = set()

    for elem in asubloc:

        for sl in elem.split("."):

            cleanSL = sl.strip()
            
            if "Note=" in cleanSL:
                continue

            if len(cleanSL) == 0:
                continue

            asublocs.add(sl.strip())

    return asublocs


cacheGene2Subloc = defaultdict(set)
cacheGene2Pfam = defaultdict(set)
cacheGene2Interpro = defaultdict(set)


def fetchUniprotInfo(geneList):

    gene2subloc = defaultdict(set)
    gene2pfam = defaultdict(set)
    gene2interpro = defaultdict(set)


    if len(geneList) == 0:
        return gene2subloc, gene2pfam, gene2interpro


    newGeneList = [x for x in geneList if not x.upper() in cacheGene2Subloc]

    print("Calling uniprot for genes: ", len(newGeneList))

    if len(newGeneList) > 0:

        try:
            up = UniprotStore()
            allConvertedIDs = up.fetch(GeneIdentity.GENE_NAME, newGeneList, toEntities=[GeneIdentity.TAXID, GeneIdentity.ENSEMBL, GeneIdentity.SUBCELLULAR_LOCATION, GeneIdentity.INTERPRO, GeneIdentity.PFAM, GeneIdentity.GENE_SYMBOL], customParams={"taxon": "Mouse [10090]"})

            #print(allConvertedIDs)
            header = allConvertedIDs.getIdxHeader()

            for row in allConvertedIDs:
                geneID = row[header[GeneIdentity.GENE_SYMBOL]].upper()
                sublocs = row[header[GeneIdentity.SUBCELLULAR_LOCATION]]
                sublocs = makeSubcellularLocations(sublocs)


                pfams = row[header[GeneIdentity.PFAM]]
                interpros = row[header[GeneIdentity.INTERPRO]]

                if pfams != None:
                    pfams = [x.strip() for x in pfams.split(";")]
                    pfams = [x for x in pfams if len(x) > 0]

                    cacheGene2Pfam[geneID] = pfams

                if interpros != None:
                    interpros = [x.strip() for x in interpros.split(";")]
                    interpros = [x for x in interpros if len(x) > 0]

                    cacheGene2Interpro[geneID] = interpros

                cacheGene2Subloc[geneID] = sublocs

                gene2subloc[geneID] = sublocs
                gene2pfam[geneID] = pfams
                gene2interpro[geneID] = interpros
        except:
            print("Uniprot call failed")
            pass
        
        print("Uniprot call ended")

    oldGenes = [x for x in geneList if not x in newGeneList]

    for x in oldGenes:
        geneID=x.upper()
        gene2subloc[geneID] = cacheGene2Subloc.get(geneID, None)
        gene2pfam[geneID] = cacheGene2Pfam.get(geneID, None)
        gene2interpro[geneID] = cacheGene2Interpro.get(geneID, None)

    return gene2subloc, gene2pfam, gene2interpro


def plotVennDiag(datadict, title, outfile):
    allSets = []
    allDescr= []

    for x in datadict:
        allSets.append(set(datadict[x]))
        allDescr.append(x)

        print("Venn Diag", x, len(datadict[x]))

    for x in allSets:
        print(type(x), len(x))


    while len(allSets) < 3:
        allSets.append(set())
        allDescr.append("")

    plt.figure()
    plt.axis('off')

    venn3(allSets, allDescr)
    ax = plt.gca()
    ax.set_axis_off()
    plt.title(title)

    if outfile != None:
        
        outname = outfile.name.replace(".html", ".png")
        plt.savefig(outname)
        return outname

    else:
        plt.show()
        return None

def loadTMResults(infile):

    gene2locev = defaultdict(lambda: defaultdict(set))

    if infile==None:
        return gene2locev

    for line in infile:

        aline = line.strip().split("\t")

        location = aline[0]
        gene = aline[3]
        pmid = aline[6]

        gene2locev[gene][location].add(pmid)

    for gene in gene2locev:
        for location in gene2locev[gene]:
            gene2locev[gene][location] = natsorted(gene2locev[gene][location], reverse=True)

    return gene2locev


def makeShortTMList(gene2locev, geneID):

    geneID = geneID.upper()

    allLocEvTuples = []
    
    curRoundIdx = 0
    roundAdded = 1
    while len(allLocEvTuples) < 10 and roundAdded > 0:

        roundAdded = 0
        
        for location in gene2locev[geneID]:

            if len(gene2locev[geneID][location]) > curRoundIdx:

                pmid = gene2locev[geneID][location][curRoundIdx]
                allLocEvTuples.append( (location, pmid) )
                roundAdded += 1

        curRoundIdx += 1

    return sorted(allLocEvTuples, key=lambda x: x[0])

def makeShortedList(inputlist):

    allelems = natsorted(inputlist)

    if len(allelems) < 10:
        return allelems

    return allelems[-10:]


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

import upsetplot
import pandas as pd
import numpy as np

def makeUpsetPlot(data, title, outfile):

    def make_data(indata):

        n_sets = len(indata)

        unionElems = set()
        for x in indata:
            unionElems = unionElems.union(set(indata[x]))

        unionElems = tuple(unionElems)
        n_samples = len(unionElems)

        setLabels = [x for x in indata]

        df = pd.DataFrame({'value': np.zeros(n_samples)})
        for i in range(len(setLabels)):

            r = []
            for x in unionElems:
                if x in indata[setLabels[i]]:
                    r.append(True)
                else:
                    r.append(False)

            df[setLabels[i]] = r
            df['value'] += r

        df.set_index([setLabels[i] for i in range(n_sets)], inplace=True)

        return df.value.groupby(level=list(range(n_sets))).count()



    ndata = make_data(data)
    upsetplot.plot( ndata , show_counts='%d')
    plt.title(title)

    if outfile != None:
        plt.savefig(outfile)
    else:
        plt.show()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-dp', '--pairwise', nargs='+', type=argparse.FileType("r"))
    parser.add_argument('-ds', '--simple', nargs='+', type=argparse.FileType("r"))
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=None)
    parser.add_argument('-p', '--maxpval', type=float, required=False, default=0.05)
    parser.add_argument('-n', '--name', type=str, required=False, default="DE comparison")
    parser.add_argument('-t', '--tm', type=argparse.FileType('r'), default=None)

    args = parser.parse_args()

    method2lists, tname2cats = perpareExpNames()
    gene2locev = loadTMResults(args.tm)


    """
    Here we prepare DE data

    """

    
    compType2DE = defaultdict(set)
    compType2line = defaultdict(lambda: dict())
    compType2idx = {}

    compType2Class = {}
    seenCompClasses = set()
    compClass2Tissue = defaultdict(set)


    for infile in args.pairwise:

        infile_category = None

        if 'facs' in infile.name:
            infile_category = 'facs'
        elif 'droplet' in infile.name:
            infile_category = 'droplet'
        elif 'wells' in infile.name:
            infile_category = 'wells'

        print(infile.name, infile_category)


        baseName = os.path.splitext(os.path.basename(infile.name))[0]
        aBaseName = baseName.split("___")

        targetTissue = aBaseName[1]
        sourceTissue = aBaseName[2]
        infile_target_tissue = aBaseName[0].split("_")[1]

        catsTarget = tname2cats[(targetTissue, infile_category)]

        defaultSource = [x for x in catsTarget]
        defaultSource[0] = "rest"
        defaultSource = tuple(defaultSource)

        catsSource = tname2cats.get((sourceTissue, infile_category), defaultSource)

        #print(catsTarget, "vs.", catsSource)

        seenCompClasses.add( (catsTarget, catsSource) )
        compType2Class[ (targetTissue, sourceTissue) ] = (catsTarget, catsSource)
        compClass2Tissue[ (catsTarget, catsSource) ].add((targetTissue, sourceTissue))

        print(targetTissue, sourceTissue, infile_target_tissue, catsTarget, catsSource)

        col2idx = {}
        for lidx, line in enumerate(infile):

            aline = [x.strip() for x in line.split("\t")] 

            
            if lidx == 0:

                for cidx, elem in enumerate(aline):
                    col2idx[elem] = cidx+1

                compType2idx[(targetTissue, sourceTissue)] = col2idx
                continue


            geneID = aline[0]
            genePVal = float(aline[col2idx["p_val_adj"]])

            compType2line[(targetTissue, sourceTissue)][geneID] = aline

            if genePVal < args.maxpval:
                compType2DE[(targetTissue, sourceTissue)].add(geneID)


    simpleCompType2DE = defaultdict(set)
    simpleCompType2line = defaultdict(lambda: dict())
    simpleCompType2idx = {}

    simpleCompType2nonDE = defaultdict(set)

    for sfname in args.simple:


        filename = sfname.name
        fname = os.path.splitext(os.path.basename(filename))[0]
        afilename = fname.split("_")

        print(afilename)

        # should become something like (ltype, tissue, method)
        tissue = afilename[1]
        ltype = afilename[0]

        if 'facs' in filename:
            infile_category = 'facs'
        elif 'droplet' in filename:
            infile_category = 'droplet'
        elif 'wells' in filename:
            infile_category = 'wells'

        targetTissue = (ltype, tissue, infile_category)
        sourceTissue = (ltype, 'all', infile_category)

        col2idx = {}
        for lidx, line in enumerate(sfname):

            aline = [x.strip() for x in line.split("\t")] 
            
            if lidx == 0:

                for cidx, elem in enumerate(aline):
                    col2idx[elem] = cidx+1

                simpleCompType2idx[fname] = col2idx
                continue


            geneID = aline[0]
            genePVal = float(aline[col2idx["p_val_adj"]])

            simpleCompType2line[fname][geneID] = aline

            if genePVal < args.maxpval:
                simpleCompType2DE[fname].add(geneID)
            else:
                simpleCompType2nonDE[fname].add(geneID)



    outDF = DataFrame()
    (head, _) = outDF._makeHTMLString("init")

    outDF = None

    args.output.write("<html><head>"+head + "</head><body>")

    allTargets = set([x[0] for x in seenCompClasses])
    allSources = set([x[1] for x in seenCompClasses])

    dfcount = 0

    unionDF = DataFrame()
    unionDF.addColumns(["GeneID"])


    targetPSet =  sorted([x for x in powerset(allTargets)], key=lambda x: len(x), reverse=True)
    sourcePSet =  sorted([x for x in powerset(allSources)], key=lambda x: len(x), reverse=True)


    targetPSet = [tuple(sorted(x, key=lambda x: (x[1], x[0]))) for x in targetPSet]
    sourcePSet = [tuple(sorted(x, key=lambda x: (x[1], x[0]))) for x in sourcePSet]


    print("All Targets", allTargets)
    print("All Sources", allSources)


    plannedComparisons = []
    performedComparisons = {}

    allTests = []
    for tlist in targetPSet:
        for slist in sourcePSet:

            if len(tlist) == 0 or len(slist) == 0:
                continue

            print("Added comparison", tlist, slist)
            plannedComparisons.append((tlist, slist))

    savedComparisonOutput=[]
    upsetData = {}
    upsetDataUpDown = defaultdict(set)


    for tlist in targetPSet:
        for slist in sourcePSet:

            if len(tlist) == 0 or len(slist) == 0:
                continue

            print("Performing Comparison", tlist, slist)

            tliststr = ["-".join(x) for x in tlist]
            sliststr = ["-".join(x) for x in slist]

            comparisonOutput = io.StringIO("")

            comparisonOutput.write("<h2>Comparison: {cl1} vs. {cl2} and {simples}</h2>".format(cl1="&".join(tliststr), cl2="&".join(sliststr), simples=",".join([x for x in simpleCompType2DE])))


            compTissuesDE = []
            
            for telem in tlist:
                for selem in slist:
                    compTissuesDE += compClass2Tissue.get((telem, selem), [])

            #print(compTissuesDE)

            simpleNonDE = set()
            for simpleCase in simpleCompType2nonDE:
                simpleNonDE = simpleNonDE.union(simpleCompType2nonDE[simpleCase])

            simpleUnion = set()
            for simpleCase in simpleCompType2DE:
                simpleUnion = simpleUnion.union(simpleCompType2DE[simpleCase])

            totalUnion = set()
            for method in compTissuesDE:
                totalUnion = totalUnion.union(compType2DE[method])
                #print(method, "simple non de, de", len(simpleNonDE.intersection(compType2DE[method])), "simple de, non de", len(simpleUnion.difference(compType2DE[method])), "simple non de, de genes", simpleNonDE.intersection(compType2DE[method]))

            totalUnion = totalUnion.union(simpleUnion)

            totalIntersect = set([x for x in totalUnion])

            for method in compTissuesDE:
                totalIntersect = totalIntersect.intersection(compType2DE[method])

            if len(slist) == 1:
                upsetData[(tuple(tliststr), tuple(sliststr))] = [x for x in totalIntersect]


            simpleIntersect = set([x for x in totalUnion])
            for simpleCase in simpleCompType2DE:
                simpleIntersect = simpleIntersect.intersection(simpleCompType2DE[simpleCase])

            #print("Sparc in simple", "Sparc" in simpleIntersect)

            # is there any simple gene which is 
            #print("Simple Intersect Analysis")
            #print("Number of intersect genes", len(totalIntersect))
            #print("Number of simple intersect genes", len(simpleIntersect))

            #print("intersect genes simple-intersect", len(simpleIntersect.difference(totalIntersect)))
            #print("intersect genes intersect-simple", len(totalIntersect.difference(simpleIntersect)), totalIntersect.difference(simpleIntersect))
            #print("non de intersect total", simpleNonDE.intersection(totalIntersect))


            """

            ABOVE: PREPARE DATA

            BELOW: UPDATE UNIONDF

            """
            print(len(totalUnion))


            columns = ['Direction']

            method2pref = {}
            for method in compTissuesDE:
                methodPrefix = str(method[0]) +" vs. "+ str(method[1])
                method2pref[method] = methodPrefix
                columns.append(methodPrefix + "Adjusted P-Value")
                columns.append(methodPrefix + "Average logFC")

            for simpleCase in simpleCompType2DE:
                columns.append("".join(simpleCase) + "Adjusted P-Value")
                columns.append("".join(simpleCase) + "Average logFC")

            unionDF.addColumns(columns, ignoreDuplicates=True)

            updateDRs = []
            for geneID in totalUnion:

                dfDict = {
                    "GeneID": geneID,
                }

                for method in compTissuesDE:

                    geneMethodLine = compType2line[method].get(geneID, None)
                    methodIdx = compType2idx[method]

                    if geneMethodLine == None:
                        dfDict[method2pref[method] + "Adjusted P-Value"] = 1.0
                        dfDict[method2pref[method] + "Average logFC"] = 0
                    else:
                        dfDict[method2pref[method] + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                        dfDict[method2pref[method] + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]


                for simpleCase in simpleCompType2DE:
                    geneMethodLine = simpleCompType2line[simpleCase].get(geneID, None)
                    methodIdx = simpleCompType2idx[simpleCase]

                    if geneMethodLine == None:
                        dfDict["".join(simpleCase) + "Adjusted P-Value"] = 1.0
                        dfDict["".join(simpleCase) + "Average logFC"] = 0
                    else:
                        dfDict["".join(simpleCase) + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                        dfDict["".join(simpleCase) + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]


                sameDirection = set()

                for x in dfDict:
                    if "Average logFC" in x:

                        foldChange = float(dfDict[x])
                        if foldChange < 0.0:
                            sameDirection.add("down")
                        elif foldChange > 0.0:
                            sameDirection.add("up")

                dfDict["Direction"] = "#".join(sameDirection)

                if len(sameDirection) > 1:
                    upsetDataUpDown[(tuple(tliststr), tuple(sliststr))].add(geneID)

                dr = DataRow.fromDict(dfDict)
                updateDRs.append(dr)
            
            unionDF.updateRowIndexed("GeneID", updateDRs, addIfNotFound=True)


            """

            BELOW: ADD COMPARISON IF ENOUGH CELLS

            """
            (head, body) = ("", "<h4>Not enough cell-types for comparison</h4>")
            outDF = None

            if len(slist)+1 >= len(allSources) or len(slist) == 1:

                print(len(totalIntersect))
                vennDiagPath = ""

                columns = ["GeneID", "Subcellular Locations", "TM Locations", "PFAM", "Interpro", "Direction"]

                method2pref = {}
                for method in compTissuesDE:
                    methodPrefix = str(method[0]) +" vs. "+ str(method[1])+" "
                    method2pref[method] = methodPrefix
                    columns.append(methodPrefix + "Adjusted P-Value")
                    columns.append(methodPrefix + "Average logFC")


                for simpleCase in simpleCompType2DE:
                    columns.append("".join(simpleCase) + "Adjusted P-Value")
                    columns.append("".join(simpleCase) + "Average logFC")

                tiLocation, tiPfam, tiInterpro = fetchUniprotInfo(totalIntersect)


                outDF = DataFrame()
                outDFnoHTML = DataFrame()

                outDF.addColumns(columns)
                outDFnoHTML.addColumns(columns)

                (head, body) = ("", "<h4>No Intersection</h4>")
                genesNotDiffInSimple = set()

                if len(totalIntersect) > 0:

                    performedComparisons[(tlist, slist)] = len(totalIntersect)

                    for geneID in totalIntersect:

                            dfDict = {
                                "GeneID": "<a target=\"_blank\" href=\"https://www.uniprot.org/uniprot/?query="+geneID+"%20and%20organism:10090\">"+geneID+"</a>"
                            }

                            dfDict["Subcellular Locations"] = "; ".join(tiLocation.get(geneID.upper(), set()))
                            dfDict["TM Locations"] = ""

                            dfDict["PFAM"] = ""
                            dfDict["Interpro"] = ""

                            dfNoHTMLDict = {
                                "GeneID": geneID,
                                "Subcellular Locations": "; ".join(tiLocation.get(geneID.upper(), set())),
                                "TM Locations": "",
                                "PFAM": "",
                                "Interpro": ""
                            }
                            
                            tmPfams = tiPfam.get(geneID.upper(), None)
                            tmInterpros = tiInterpro.get(geneID.upper(), None)
                            tmLocation = gene2locev.get(geneID.upper(), None)

                            if tmLocation != None:

                                geneLocEvs = makeShortTMList(gene2locev, geneID.upper())
                                allLinkStr = []

                                for location, pmid in geneLocEvs:
                                    linkstr = "<a target=\"_blank\" href=\"{link}\">{desc}</a>".format(link="https://www.ncbi.nlm.nih.gov/pubmed/"+str(pmid), desc=location + " ("+str(pmid)+")")
                                    allLinkStr.append(linkstr)

                                dfDict["TM Locations"] = "</br>".join(allLinkStr)
                                dfNoHTMLDict["TM Locations"] = ";".join([str(x[0]) + ", " + str(x[1]) for x in geneLocEvs])


                            if tmPfams != None:
                                allLinkStr = []

                                for pfamID in tmPfams:
                                    linkstr = "<a target=\"_blank\" href=\"{link}\">{desc}</a>".format(link="https://pfam.xfam.org/family/"+str(pfamID), desc=pfamID)
                                    allLinkStr.append(linkstr)

                                dfDict["PFAM"] = "</br>".join(allLinkStr)
                                dfNoHTMLDict["PFAM"] = ";".join([pfamID for pfamID in tmPfams])

                            if tmInterpros != None:
                                allLinkStr = []

                                for interproID in tmInterpros:
                                    linkstr = "<a target=\"_blank\" href=\"{link}\">{desc}</a>".format(link="https://www.ebi.ac.uk/interpro/entry/"+str(interproID), desc=interproID)
                                    allLinkStr.append(linkstr)

                                dfDict["Interpro"] = "</br>".join(allLinkStr)
                                dfNoHTMLDict["Interpro"] = ";".join([interproID for interproID in tmInterpros])

                            for method in compTissuesDE:

                                geneMethodLine = compType2line[method][geneID]
                                methodIdx = compType2idx[method]

                                dfDict[method2pref[method] + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                                dfDict[method2pref[method] + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]

                                dfNoHTMLDict[method2pref[method] + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                                dfNoHTMLDict[method2pref[method] + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]

                            for simpleCase in simpleCompType2DE:
                                geneMethodLine = simpleCompType2line[simpleCase].get(geneID, None)
                                methodIdx = simpleCompType2idx[simpleCase]

                                if geneMethodLine == None:
                                    dfDict["".join(simpleCase) + "Adjusted P-Value"] = 1.0
                                    dfDict["".join(simpleCase) + "Average logFC"] = 0

                                    dfNoHTMLDict["".join(simpleCase) + "Adjusted P-Value"] = 1.0
                                    dfNoHTMLDict["".join(simpleCase) + "Average logFC"] = 0
                                else:
                                    dfDict["".join(simpleCase) + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                                    dfDict["".join(simpleCase) + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]

                                    dfNoHTMLDict["".join(simpleCase) + "Adjusted P-Value"] = geneMethodLine[methodIdx['p_val_adj']]
                                    dfNoHTMLDict["".join(simpleCase) + "Average logFC"] = geneMethodLine[methodIdx['avg_logFC']]

                                    if float(geneMethodLine[methodIdx['p_val_adj']]) > args.maxpval:
                                        genesNotDiffInSimple.add(geneID)


                            sameDirection = set()

                            for x in dfDict:
                                if "Average logFC" in x:
                                    foldChange = float(dfDict[x])
                                    if foldChange < 0.0:
                                        sameDirection.add("down")
                                    elif foldChange > 0.0:
                                        sameDirection.add("up")

                            dfDict["Direction"] = "#".join(sameDirection)
                            dfNoHTMLDict["Direction"] = "#".join(sameDirection)

                            dr = DataRow.fromDict(dfDict)
                            outDF.addRow(dr)
                
                            drNoHTML = DataRow.fromDict(dfNoHTMLDict)
                            outDFnoHTML.addRow(drNoHTML)



            print("Total union", len(totalUnion))

            dfcount += 1

            if len(totalIntersect) > 0 and outDF != None:
                (head, body) = outDF._makeHTMLString("diffdf" + str(dfcount))
                
                testname = os.path.splitext(args.output.name)[0] + "___" + "_".join(tliststr) + "_vs_" + "_".join(sliststr) + ".intersect"

                print("Exporting tsv") 
                outDF.export(testname + ".tsv", ExportTYPE.TSV)
                outDFnoHTML.export(testname + ".nohtml.tsv", ExportTYPE.TSV)

                if len(slist) == len(allSources) and len(tlist) == len(allTargets):
                    outbase = os.path.splitext(args.output.name)[0] + ".full.tsv"
                    outDF.export(outbase, ExportTYPE.TSV)

                    outDFnoHTML.export(os.path.splitext(args.output.name)[0] + ".full.nohtml.tsv", ExportTYPE.TSV)



                print("Exporting xlsx") 
                outDF.export(testname + ".xlsx", ExportTYPE.XLSX)
                outDFnoHTML.export(testname + ".nohtml.xlsx", ExportTYPE.XLSX)
                
                comparisonOutput.write("<h4>Additional Files</h4><ol>")
                comparisonOutput.write("<li><a href=\"{loc}\" target=\"_blank\">TSV File</a></li>".format(loc=os.path.basename(testname) + ".tsv"))
                comparisonOutput.write("<li><a href=\"{loc}\" target=\"_blank\">XLSX File</a></li>".format(loc=os.path.basename(testname) + ".xlsx"))
                comparisonOutput.write("<li><a href=\"{loc}\" target=\"_blank\">TSV File (no HTML)</a></li>".format(loc=os.path.basename(testname) + ".nohtml.tsv"))
                comparisonOutput.write("<li><a href=\"{loc}\" target=\"_blank\">XLSX File (noHTML)</a></li>".format(loc=os.path.basename(testname) + ".nohtml.xlsx"))
                comparisonOutput.write("</ol>")
                if len(slist) > 1:
                    comparisonOutput.write("<h5>Genes differential < {pval} in intersection but not in simple</h5>".format(pval=args.maxpval))
                comparisonOutput.write("<p>{genes}</p>".format(genes=", ".join(genesNotDiffInSimple)))
                
            comparisonOutput.write("<h3>Intersection Overview</h3></br><h4>Intersection Genes: {geneCount}</h4>".format(geneCount=len(totalIntersect)))
            comparisonOutput.write(body)

            savedComparisonOutput.append(comparisonOutput)

            
            #dfcount += 1
            #(headCombined, bodyCombined) = combinedDF._makeHTMLString("combdf" + str(dfcount))

            #args.output.write("<h3>Combined Overview</h3>")
            #args.output.write(bodyCombined)


    dfcount += 1

    print("Found a union DF", len(unionDF))

    (headUnion, bodyUnion) = unionDF._makeHTMLString("uniondf" + str(dfcount))

    unionOutBase = os.path.splitext(args.output.name)[0] + "_union"
    unionDF.export(unionOutBase + ".xlsx", ExportTYPE.XLSX)
    unionDF.export(unionOutBase + ".tsv", ExportTYPE.TSV)


    args.output.write("<h1>All Comparisons</h1>")
    args.output.write("<ul>")

    allComparisonsLI = []
    for tlist, slist in plannedComparisons:
            tliststr = ["-".join(x) for x in tlist]
            sliststr = ["-".join(x) for x in slist]

            allComparisonsLI.append( "<li>" +  " & ".join(tliststr) + " vs. " +  " & ".join(sliststr) + ": {count}</li>".format(count=performedComparisons.get((tlist, slist), 0)))

    args.output.write( "\n".join(allComparisonsLI))
    args.output.write("</ul>")

    args.output.write("<h1>Differential Gene Overlap</h1>")

    for compOut in savedComparisonOutput:
        args.output.write(str(compOut.getvalue()))
    
    args.output.write("<h3>Gene Union Overview</h3>")



    args.output.write("Overlap of tissues")


    niceUpsetData = {}

    for x in upsetData:
        print(x)
        name="&".join(x[0]) + " vs. " + "&".join(x[1])

        niceUpsetData[name] = upsetData[x]

    makeUpsetPlot(niceUpsetData, "Tissue Overlap", unionOutBase + ".issue.upset.png")
    args.output.write("<h4>Tissue Overlap</h4>")
    args.output.write("<img src=\"{loc}\"/>".format(loc=os.path.basename(unionOutBase) + ".issue.upset.png"))

    niceUpsetData = {}

    for x in upsetData:
        name="&".join(x[0]) + " vs. " + "&".join(x[1])

        genes = upsetData[x]
        genes = genes.difference(upsetDataUpDown[x])

        niceUpsetData[name] = genes


    makeUpsetPlot(niceUpsetData, "Tissue Overlap (no up/down)", unionOutBase + ".issue.upset.noupdown.png")
    args.output.write("<h4>Tissue Overlap</h4>")
    args.output.write("<img src=\"{loc}\"/>".format(loc=os.path.basename(unionOutBase) + ".issue.upset.noupdown.png"))

    args.output.write(bodyUnion)

    args.output.write("<a href=\"{loc}\" target=\"_blank\">TSV File</a><br/>".format(loc=os.path.basename(unionOutBase) + ".tsv"))
    args.output.write("<a href=\"{loc}\" target=\"_blank\">XLSX File</a><br/>".format(loc=os.path.basename(unionOutBase) + ".xlsx"))

    args.output.write("</body></html>")


    #file:///D:/dev/git/tabula-muris/00_data_ingest/sanne_analysis/droplet_analysis_excl_overview/lung.html
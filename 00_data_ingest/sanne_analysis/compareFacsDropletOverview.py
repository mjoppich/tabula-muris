


import glob
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

import os,sys

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


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
        
        outname = outfile.replace(".html", ".png")
        plt.savefig(outname)
        return outname

    else:
        plt.show()
        return None

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-f', '--facs', type=argparse.FileType("r"))
    parser.add_argument('-d', '--droplet', type=argparse.FileType("r"))
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=None)

    args = parser.parse_args()



    #for facsFile in glob.glob("facs_analysis_excl_overview/*.intersect.tsv"):

    #    dropletFile = facsFile.replace("facs", "droplet")

    #    facsFileExists = os.path.exists(facsFile)
    #    dropletFileExists = os.path.exists(dropletFile)

    #    if not facsFileExists or not dropletFileExists:
    #        continue
    
    facsDF = DataFrame.parseFromFile(args.facs)
    dropletDF = DataFrame.parseFromFile(args.droplet)

    #metaName = facsFile.replace("facs_analysis_excl_overview/", "droplet_facs_meta/")
    #metaName = metaName.replace("-facs", "")

    method2genes = defaultdict(set)

    method2row = defaultdict(lambda: dict())

    for drow in facsDF:
        geneID = drow["GeneID"]
        method2genes["FACS"].add(geneID)

        method2row["FACS"][geneID] = drow


    dropletRows = []
    for drow in dropletDF:
        geneID = drow["GeneID"]
        method2genes["DROPLET"].add(geneID)

        method2row["DROPLET"][geneID] = drow



    facsOnly = method2genes["FACS"].difference(method2genes["DROPLET"])
    dropletOnly = method2genes["DROPLET"].difference(method2genes["FACS"])
    bothFacsDroplet = method2genes["FACS"].intersection(method2genes["DROPLET"])


    vennDiagPath = plotVennDiag(method2genes, "", args.output.name)


    combinedDF = DataFrame()

    newColumns = ["GeneID"]
    for colName in facsDF.getHeader():
        if not colName in newColumns:
            newColumns.append(colName)

    for colName in dropletDF.getHeader():
        if not colName in newColumns:
            newColumns.append(colName)

    combinedDF.addColumns(newColumns)

    combinedDF.updateRowIndexed("GeneID", [method2row["FACS"][x] for x in bothFacsDroplet], addIfNotFound=True)
    combinedDF.updateRowIndexed("GeneID", [method2row["DROPLET"][x] for x in bothFacsDroplet], addIfNotFound=True)


    (head, body ) = combinedDF._makeHTMLString("combdf")

    args.output.write("<html><head>"+head + "</head><body><h1>Differential Gene Overlap</h1><img src=\"./"+os.path.basename(vennDiagPath)+"\"/>")

    args.output.write("<p>")
    args.output.write("Original Filename FACS: ")
    args.output.write(args.facs.name)
    args.output.write("</p>")

    args.output.write("<p>")
    args.output.write("Original Filename Droplet: ")
    args.output.write(args.droplet.name)
    args.output.write("</p>")

    args.output.write("<h1>FACS only genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(facsOnly))))

    args.output.write("<h1>DROPLET only genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(dropletOnly))))

    args.output.write("<h1>FACS/DROPLET genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(bothFacsDroplet))))

    args.output.write("<h1>Combined Overview</h1>"+body)
    args.output.write("</body></html>")




import glob
import argparse
from collections import defaultdict

import matplotlib as mpl
mpl.use('Agg')

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
    parser.add_argument('-f', '--facs', type=argparse.FileType("r"), default=None)
    parser.add_argument('-d', '--droplet', type=argparse.FileType("r"), default=None)
    parser.add_argument('-w', '--wells', type=argparse.FileType("r"), default=None)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()



    #for facsFile in glob.glob("facs_analysis_excl_overview/*.intersect.tsv"):

    #    dropletFile = facsFile.replace("facs", "droplet")

    #    facsFileExists = os.path.exists(facsFile)
    #    dropletFileExists = os.path.exists(dropletFile)

    #    if not facsFileExists or not dropletFileExists:
    #        continue
    
    if args.facs != None:
        facsDF = DataFrame.parseFromFile(args.facs)
    else:
        facsDF = None

    dropletDF = DataFrame.parseFromFile(args.droplet)
    wellsDF = DataFrame.parseFromFile(args.wells)

    #metaName = facsFile.replace("facs_analysis_excl_overview/", "droplet_facs_meta/")
    #metaName = metaName.replace("-facs", "")

    method2genes = defaultdict(set)

    method2row = defaultdict(lambda: dict())


    def oldrow2new(drow, pref):
        drowdict = drow.to_dict()
        ndict = {}

        ncolumns = set()
        for k in drowdict:
            if 'P-Value' in k or 'logFC' in k:
                ndict[pref + " " + k] = drowdict[k]
            else:
                ndict[k] = drowdict[k]

        ncolumns = [x for x in ndict]
        return DataRow.fromDict(ndict), ncolumns

    newColumns = set()

    if facsDF!=None:
        for drow in facsDF:
            geneID = drow["GeneID"]
            method2genes["FACS"].add(geneID)

            method2row["FACS"][geneID], ncolumns = oldrow2new(drow, "FACS")

            newColumns = newColumns.union(ncolumns)


    dropletRows = []
    for drow in dropletDF:
        geneID = drow["GeneID"]
        method2genes["DROPLET"].add(geneID)

        method2row["DROPLET"][geneID], ncolumns = oldrow2new(drow, "DROPLET")
        newColumns = newColumns.union(ncolumns)

    wellsRows = []
    for drow in wellsDF:
        geneID = drow["GeneID"]
        method2genes["WELLS"].add(geneID)

        method2row["WELLS"][geneID], ncolumns = oldrow2new(drow, "WELLS")
        newColumns = newColumns.union(ncolumns)

    facsOnly = method2genes["FACS"].difference(method2genes["DROPLET"]).difference(method2genes["WELLS"])
    dropletOnly = method2genes["DROPLET"].difference(method2genes["FACS"]).difference(method2genes["WELLS"])
    wellsOnly = method2genes["WELLS"].difference(method2genes["FACS"]).difference(method2genes["DROPLET"])

    bothFacsDroplet = method2genes["FACS"].intersection(method2genes["DROPLET"])
    bothDropletWells = method2genes["WELLS"].intersection(method2genes["DROPLET"])

    facsDropletDropletWells = bothFacsDroplet.union(bothDropletWells)


    allFacsDropletWells = method2genes["FACS"].intersection(method2genes["DROPLET"]).intersection(method2genes["WELLS"])

    bothFacsDropletWithWells = method2genes["FACS"].intersection(method2genes["DROPLET"]).union(allFacsDropletWells)

    allGenes = method2genes["FACS"].union(method2genes["DROPLET"]).union(method2genes["WELLS"])

    vennDiagPath = plotVennDiag(method2genes, "", args.output.name)


    combinedDF = DataFrame()

    sortedColumns = ["GeneID"]
    for colName in dropletDF.getHeader():
        if colName in newColumns:
            sortedColumns.append(colName)
    
    for x in sorted(newColumns):
        if not x in sortedColumns:
            sortedColumns.append(x)


    combinedDF.addColumns(sortedColumns)

    displayGenes = None
    overviewName = ""

    if args.wells != None and args.facs != None and args.droplet != None:

        if len(allFacsDropletWells) > 0:
            displayGenes = allFacsDropletWells
            overviewName = "FACS&DROPLET&WELLS"

        elif len(facsDropletDropletWells) > 0:
            displayGenes = facsDropletDropletWells
            overviewName = "FACS&DROPLET + DROPLET&WELLS"

        else:
            displayGenes = bothFacsDroplet
            overviewName = "FACS&DROPLET"

        if displayGenes == None:
            print("NO GENES")

    elif args.wells == None and args.facs != None and args.droplet != None:
        displayGenes = bothFacsDropletWithWells
        overviewName = "FACS&DROPLET+WELLS"
    elif args.wells != None and args.facs == None and args.droplet != None:
        print("PRINT DROPLET WELLS")
        displayGenes = bothDropletWells
        overviewName = "WELLS&DROPLET"
    else:
        print("ALL GENES DISPLAYED!")
        displayGenes = allGenes
        overviewName = "(ALL GENES)"

    if len(displayGenes) == 0:
        print("NO DISPLAY GENES")
        displayGenes = allGenes
        overviewName = "(ALL GENES)"


    combinedDF.updateRowIndexed("GeneID", [method2row["FACS"][x] for x in displayGenes if x in method2row["FACS"]], ignoreMissingCols=True, addIfNotFound=True)
    combinedDF.updateRowIndexed("GeneID", [method2row["DROPLET"][x] for x in displayGenes if x in method2row["DROPLET"]], ignoreMissingCols=True, addIfNotFound=True)
    combinedDF.updateRowIndexed("GeneID", [method2row["WELLS"][x] for x in displayGenes if x in method2row["WELLS"]], ignoreMissingCols=True, addIfNotFound=True)



    def rowFunc(rowdata):

        sameDirection = set()
        for x in combinedDF.column2idx:
            if "Average logFC" in x:
                idx = combinedDF.column2idx[x]

                if rowdata[idx] == None:
                    continue

                foldChange = float(rowdata[idx])
                if foldChange < 0.0:
                    sameDirection.add("down")
                elif foldChange > 0.0:
                    sameDirection.add("up")

        return "#".join(sameDirection)


    #now we must somehow update the Direction
    combinedDF.applyByRow("Direction", rowFunc)

    (head, body ) = combinedDF._makeHTMLString("combdf")

    args.output.write("<html><head>"+head + "</head><body><h1>Differential Gene Overlap</h1><img src=\"./"+os.path.basename(vennDiagPath)+"\"/>")

    args.output.write("<p>")
    args.output.write("Original Filename FACS: ")
    args.output.write(args.facs.name if args.facs != None else "")
    args.output.write("</p>")

    args.output.write("<p>")
    args.output.write("Original Filename Droplet: ")
    args.output.write(args.droplet.name)
    args.output.write("</p>")


    args.output.write("<p>")
    args.output.write("Original Filename Wells: ")
    args.output.write(args.wells.name)
    args.output.write("</p>")

    """

    facsOnly
    dropletOnly
    wellsOnly

    bothFacsDroplet
    allFacsDropletWells
    bothFacsDropletWithWells

    """


    args.output.write("<h1>FACS only genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(facsOnly))))

    args.output.write("<h1>DROPLET only genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(dropletOnly))))

    args.output.write("<h1>WELLS only genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(wellsOnly))))

    args.output.write("<h1>FACS&DROPLET genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(bothFacsDroplet))))
    args.output.write("<h1>WELLS&DROPLET genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(bothDropletWells))))

    args.output.write("<h1>FACS&DROPLET&WELLS genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(allFacsDropletWells))))

    args.output.write("<h1>FACS&DROPLET + FACS&DROPLET&WELLS genes</h1>")
    args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(bothFacsDropletWithWells))))

    if displayGenes != None:
        args.output.write("<h1>"+overviewName+"</h1>")
        args.output.write("<p>{fog}</p>".format(fog=", ".join(sorted(displayGenes))))


    args.output.write("<h1>Combined Overview "+overviewName+"</h1>"+body)
    args.output.write("</body></html>")

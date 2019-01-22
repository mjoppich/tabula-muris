import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

import os,sys

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))

from porestat.utils.DataFrame import DataFrame, DataRow


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

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-de', '--diffreg', nargs='+', type=argparse.FileType("r"))
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=None)
    parser.add_argument('-p', '--maxpval', type=float, required=False, default=0.05)
    parser.add_argument('-n', '--name', type=str, required=False, default="DE comparison")

    args = parser.parse_args()

    print([x.name for x in args.diffreg])

    method2de = defaultdict(set)
    method2line = defaultdict(lambda: dict())
    method2idx = {}
    foundMethods = set()

    for infile in args.diffreg:

        infile_category = None

        if 'facs' in infile.name:
            infile_category = 'facs'
        elif 'droplet' in infile.name:
            infile_category = 'droplet'
        elif 'wells' in infile.name:
            infile_category = 'wells'

        foundMethods.add(infile_category)

        col2idx = {}
        for lidx, line in enumerate(infile):

            aline = [x.strip() for x in line.split("\t")] 

            
            if lidx == 0:

                for cidx, elem in enumerate(aline):
                    col2idx[elem] = cidx+1

                #print(col2idx)
                method2idx[infile_category] = col2idx
                continue


            geneID = aline[0]
            genePVal = float(aline[col2idx["p_val_adj"]])

            method2line[infile_category][geneID] = aline

            if genePVal < args.maxpval:
                method2de[infile_category].add(geneID)


            
    totalUnion = set()
    for method in method2de:
        totalUnion = totalUnion.union(method2de[method])

    totalIntersect = totalUnion

    for method in method2de:
        totalIntersect = totalIntersect.intersection(method2de[method])

    #print(totalIntersect)

    vennDiagPath = plotVennDiag(method2de, args.name, args.output)

    outDF = DataFrame()
    outDF.addColumns(["GeneID", "Method", "Adjusted P-Value", "Average logFC"])

    for geneID in totalIntersect:
        for method in foundMethods:

            geneMethodLine = method2line[method][geneID]
            methodIdx = method2idx[method]

            dfDict = {
                "GeneID": geneID,
                "Method": method,
                "Adjusted P-Value": geneMethodLine[methodIdx['p_val_adj']],
                "Average logFC": geneMethodLine[methodIdx['avg_logFC']]
            }

            dr = DataRow.fromDict(dfDict)
            outDF.addRow(dr)

            #print(geneID, method, geneMethodLine[methodIdx['p_val_adj']], geneMethodLine[methodIdx['avg_logFC']])
        

    combinedDF = DataFrame()
    combinedDF.addColumns(["GeneID", "Method", "Adjusted P-Value", "Average logFC"])


    for geneID in totalUnion:
        for method in foundMethods:

            geneMethodLine = method2line[method].get(geneID, None)
            methodIdx = method2idx[method]

            if geneMethodLine == None:
                continue

            dfDict = {
                "GeneID": geneID,
                "Method": method,
                "Adjusted P-Value": geneMethodLine[methodIdx['p_val_adj']],
                "Average logFC": geneMethodLine[methodIdx['avg_logFC']]
            }

            dr = DataRow.fromDict(dfDict)
            combinedDF.addRow(dr)

    (head, body) = outDF._makeHTMLString("diffdf")
    (headCombined, bodyCombined) = combinedDF._makeHTMLString("combdf")


    args.output.write("<html><head>"+head + "</head><body><h1>Differential Gene Overlap</h1><img src=\"./"+os.path.basename(vennDiagPath)+"\"/><h1>Tabular Overview</h1>"+body + "<h1>Combined Overview</h1>"+bodyCombined + "</body></html>")
    
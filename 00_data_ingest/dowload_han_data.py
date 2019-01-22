import urllib.request
import os

with open("han_data/han_metadata.csv") as fin:

    for line in fin:

        aline = line.strip().split(",")

        gsm = aline[0]

        if gsm == None or not gsm.startswith("GSM"):
            continue

        print(gsm)

        agsm = gsm.split("_")[0]

        urlBase = "https://www.ncbi.nlm.nih.gov/geo/download/?acc="
        
        
        if os.path.isfile('han_data/' + gsm + ".gz"):
            continue


        gsmUrl = urlBase + agsm + "&format=file&file=" + gsm + ".gz"

        try:
            urllib.request.urlretrieve(gsmUrl, 'han_data/' + gsm + ".gz")
        except:
            print("ERROR", "file not found", gsm, gsmUrl)

        
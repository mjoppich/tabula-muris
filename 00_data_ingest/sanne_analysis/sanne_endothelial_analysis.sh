set -x
SCRIPT=combineSimplePairwise.py

if true; then
    for TECH in droplet wells
    do

        for TISSUE in kidney liver bladder lung
        do
            CMD="python3 $SCRIPT --expr ./${TECH}_expression/ --exprAll ./${TECH}_expression_all/ --tm gene_location.combined.pmid --simple ${TECH}_analysis/*endo_${TISSUE}_all.tsv --pairwise `ls ${TECH}_analysis_excl/endo_excl_${TISSUE}_* | grep -v Aorta | grep -v Heart` --output ${TECH}_endo_analysis/$TISSUE.html"
            echo $CMD
            $CMD || exit
        done

    done
fi

if true; then

    for TECH in facs
    do
        for TISSUE in kidney liver lung aorta brain
        do
            CMD="python3 $SCRIPT --expr ./${TECH}_expression/ --exprAll ./${TECH}_expression_all/ --tm gene_location.combined.pmid --simple ${TECH}_analysis/*endo_${TISSUE}_*.tsv --pairwise `ls ${TECH}_analysis_excl/endo_excl_${TISSUE}_*` --output ${TECH}_endo_analysis/$TISSUE.html"
            echo $CMD

            $CMD || exit
        done
    done

fi


for TISSUE in kidney liver lung
do
python3 compareTechnologies.py --facs facs_endo_analysis/$TISSUE.full.tsv --droplet droplet_endo_analysis/$TISSUE.full.tsv --wells wells_endo_analysis/$TISSUE.full.tsv --output droplet_facs_meta/$TISSUE.html
done

TISSUE=bladder
python3 compareTechnologies.py --droplet droplet_endo_analysis/$TISSUE.full.tsv --wells wells_endo_analysis/$TISSUE.full.tsv --output droplet_facs_meta/$TISSUE.html
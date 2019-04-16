
for METHOD in facs droplet wells;
do

rm -rf ./$METHOD\_expression/
rm -rf ./$METHOD\_expression_all/

mkdir -p ./$METHOD\_expression_all
mkdir -p ./$METHOD\_expression


scp markus@weismann:/usr/local/ssd/sanne/tabula-muris/00_data_ingest/sanne_analysis/$METHOD\_expression/* ./$METHOD\_expression
scp markus@weismann:/usr/local/ssd/sanne/tabula-muris/00_data_ingest/sanne_analysis/$METHOD\_expression_all/* ./$METHOD\_expression_all

done
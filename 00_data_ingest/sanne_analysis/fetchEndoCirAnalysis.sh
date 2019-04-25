
for METHOD in facs droplet wells;
do

rm -rf ./${METHOD}_endo_vs_cir/
mkdir -p ./${METHOD}_endo_vs_cir/

scp markus@weismann:/usr/local/ssd/sanne/tabula-muris/00_data_ingest/sanne_analysis/${METHOD}_endo_vs_cir/* ${METHOD}_endo_vs_cir/

done
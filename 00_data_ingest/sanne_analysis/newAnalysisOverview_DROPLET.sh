
mkdir -p droplet_analysis_excl_overview


python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/cirepi_excl_tissue_kidney.html --name "Kidney CirEpiEndo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cirepi_excl_kidney* droplet_analysis_excl/endo_excl_kidney*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/cirepi_excl_tissue_lung.html --name "Lung CirEpiEndo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cirepi_excl_lung* droplet_analysis_excl/endo_excl_lung*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/cirepi_excl_tissue_aorta.html --name "Aorta CirEpiEndo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cirepi_excl_aorta* droplet_analysis_excl/endo_excl_aorta*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/cirepi_excl_tissue_bladder.html --name "Bladder CirEpiEndo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cirepi_excl_bladder* droplet_analysis_excl/endo_excl_bladder*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/cirepi_excl_tissue_liver.html --name "Liver CirEpiEndo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cirepi_excl_liver* droplet_analysis_excl/endo_excl_liver*

exit()

python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_tissue_kidney.html --name "Kidney" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/excl_kidney*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_tissue_lung.html --name "Lung" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/excl_lung*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_tissue_aorta.html --name "Aorta" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/excl_aorta*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_tissue_bladder.html --name "Bladder" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/excl_bladder*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_tissue_liver.html --name "Liver" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/excl_liver*



python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_kidney.html --name "Kidney Endos" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_excl_kidney_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_lung.html --name "Lung Endos" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_excl_lung_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_aorta.html --name "Aorta Endos" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_excl_aorta_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_bladder.html --name "Bladder Endos" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_excl_bladder_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_liver.html --name "Liver Endos" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_excl_liver_*


python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_epi_kidney.html --name "Kidney Endo-Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_epi_excl_kidney_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_epi_lung.html --name "Lung Endo-Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_epi_excl_lung_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_epi_aorta.html --name "Aorta Endo-Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_epi_excl_aorta_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_epi_bladder.html --name "Bladder Endo-Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_epi_excl_bladder_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_epi_liver.html --name "Liver Endo-Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_epi_excl_liver_*


python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_kidney.html --name "Kidney Endo-Cir" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_kidney_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_lung.html --name "Lung Endo-Cir" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_lung_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_aorta.html --name "Aorta Endo-Cir" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_aorta_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_bladder.html --name "Bladder Endo-Cir" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_bladder_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_liver.html --name "Liver Endo-Cir" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_liver_*


python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_kidney.html --name "Kidney Endo-Cir/Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_kidney_* droplet_analysis_excl/endo_epi_excl_kidney_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_lung.html --name "Lung Endo-Cir/Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_lung_* droplet_analysis_excl/endo_epi_excl_lung_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_aorta.html --name "Aorta Endo-Cir/Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_aorta_* droplet_analysis_excl/endo_epi_excl_aorta_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_bladder.html --name "Bladder Endo-Cir/Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_bladder_* droplet_analysis_excl/endo_epi_excl_bladder_*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_liver.html --name "Liver Endo-Cir/Epi" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_liver_* droplet_analysis_excl/endo_epi_excl_liver_*



x_endo against all cir, all epi and y_endo
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_kidney.html --name "Kidney Endo-Cir/Epi/Endo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_kidney_* droplet_analysis_excl/endo_epi_excl_kidney_* droplet_analysis_excl/excl_kidney*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_lung.html --name "Lung Endo-Cir/Epi/Endo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_lung_* droplet_analysis_excl/endo_epi_excl_lung_* droplet_analysis_excl/excl_lung*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_aorta.html --name "Aorta Endo-Cir/Epi/Endo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_aorta_* droplet_analysis_excl/endo_epi_excl_aorta_* droplet_analysis_excl/excl_aorta*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_bladder.html --name "Bladder Endo-Cir/Epi/Endo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_bladder_* droplet_analysis_excl/endo_epi_excl_bladder_* droplet_analysis_excl/excl_bladder*
python3 prepareExclusiveAnalysis.py --output droplet_analysis_excl_overview/excl_endo_cir_epi_liver.html --name "Liver Endo-Cir/Epi/Endo" --tm gene_location.combined.pmid --diffreg droplet_analysis_excl/endo_cir_excl_liver_* droplet_analysis_excl/endo_epi_excl_liver_* droplet_analysis_excl/excl_liver*
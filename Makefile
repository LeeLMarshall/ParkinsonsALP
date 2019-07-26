objects=\
	IndexPage/index.html \
	Test/index.html \
	Mice_CecalPatch_Padlock/index.html \
	Mice_DSS_Padlock/index.html \
	Appendix_PDvsControls_Padlock/index.html \
	Appendix_PDvsControls_RNAseq/index.html \
	Appendix_AgeAcceleration_Padlock/index.html \
	Appendix_AgeAcceleration_RNAseq/index.html \
	Brain_PFC_Padlock_CGonly/index.html \
	Brain_OFB_Padlock_CGonly/index.html \
	Brain_PFCRep_Padlock_withGLU/index.html \
	Brain_PFCRep_Padlock_withGLU_Braak/index.html \
	Brain_AgeAcceleration_Padlock \
	Discover_OverlapAppendixMice/index.html \
	Discover_OverlapAppendixBrain/index.html \
	Discover_GenomicEnrichment/index.html \
	Discover_Aging/index.html \
	Discover_Appendix_Meth2Expression/index.html \
	Results_Appendix/index.html \
	Results_Brain/index.html \
	Results_Aging/index.html \
	Results_Mice/index.html \
	Results_Proteomics/index.html \
	Results_etc/index.html \
	Figure1/index.html \
	Figure2/index.html \
	Figure3/index.html \
	Figure4/index.html \
	Figure5/index.html
# 	Brain_PFC_hMeDIPseq/index.html

all: navbar $(objects)


navbar: include/before_body.html

# Prepare navigation bar
include/before_body.html: code/generateNavigationBar.R index.json
	Rscript code/generateNavigationBar.R

# Rendering a Rmd file
# has to be redone if navbar has changed
%.html: %.Rmd include/before_body.html
	Rscript -e 'rmarkdown::render("$<", knit_root_dir="./")'
	touch $(dir $<)restart.txt

# deployment to development server
devel:
	rsync -rav \
	--include '*/www/*' \
	--exclude '*.RDS' \
	--exclude '*.RData' \
	--exclude '*.bed' \
	--exclude '*.bdg' \
	--exclude '.*' \
	--exclude 'input' \
	--exclude '*.txt.gz' \
	. minge.ibt.lt:ShinyApps/VAI/PD2018
	ssh minge.ibt.lt chmod -R ugo-w ShinyApps/VAI/PD2018/*

deploy:
	rsync -rav \
	--include '*/www/*' \
	--exclude '*.RDS' \
	--exclude '*.RData' \
	--exclude '*.bed' \
	--exclude '*.bdg' \
	--exclude '.*' \
	--exclude 'input' \
	--exclude '*.txt.gz' \
	. shiny@minge.ibt.lt:VAI/
	ssh shiny@minge.ibt.lt chmod -R ugo-w VAI/*

# $(objects): code rsrc FORCE
# 	make -C $@



# code: FORCE
# 	rsync -rav ./code minge.ibt.lt:ShinyApps/VAI/PD2018/

# rsrc: FORCE
# 	rsync -rav ./rsrc minge.ibt.lt:ShinyApps/VAI/PD2018/

# FORCE:

# TODO deployment script
# Make sure restart txt is deployed too

# dashboard: FORCE
# 	rsync -rav ./dashboard shiny@minge.ibt.lt:VAI/	
# 	ssh shiny@minge.ibt.lt touch VAI/dashboard/restart.txt	




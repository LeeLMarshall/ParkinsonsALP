objects=index.Rmd \
	../code/common.R

all: index.html


index.html: $(objects)
	R -e 'rmarkdown::render("index.Rmd")'
	touch restart.txt

devel:
	rsync -rav ../code minge.ibt.lt:ShinyApps/VAI/PD2018/
	rsync -rav \
	--include 'www/*' \
	--include 'index.Rmd' \
	--include 'index.html' \
	--include 'index_data' \
	--exclude '*.*' \
	. minge.ibt.lt:ShinyApps/VAI/PD2018/Test
	ssh minge.ibt.lt touch ShinyApps/VAI/PD2018/Test/restart.txt

deploy:
	rsync -rav --exclude '*.RDS' . shiny@minge.ibt.lt:/home/shiny/VAI/Test


clean:
	rm -rf index_*
	rm index.html

rmrds:
	rm -rf *.RDS
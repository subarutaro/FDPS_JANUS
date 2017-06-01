FILE=manual
TARGET=$(FILE).pdf
TEX=platex

all	:$(TARGET)

$(TARGET)	:$(FILE).dvi 
	dvipdfmx $(FILE).dvi
	open $(TARGET)

$(FILE).dvi	:$(FILE).tex
	platex   $(FILE).tex
	bibtex 	 $(FILE)
	platex   $(FILE).tex

back	:
	mkdir bkup
	cp -r *.tex Makefile figure bkup/.
	tar zcvf bkup$(shell date +"%Y%m%d").tar.gz bkup
	scp bkup$(shell date +"%Y%m%d").tar.gz tote:master_thesis/.
	rm -r bkup

clean	:
	rm $(FILE).dvi $(FILE).pdf $(FILE).aux $(FILE).bbl

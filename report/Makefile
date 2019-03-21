# Makefile 

# Master file:
MASTER = MSc_Report

# File lists:
RNW = $(wildcard chapter*.Rnw)
R = $(patsubst %.Rnw,%.R,$(RNW))
TEX = $(patsubst %.Rnw,%.tex,$(RNW))

# Some technical details
.SUFFIXES: .Rnw .R .tex .pdf
.SILENT: *.pdf *.tex
.PHONY: FORCE_LATEXMK all short clean



# Targets:
all: $(MASTER).pdf

$(MASTER).pdf: $(MASTER).tex FORCE_LATEXMK
	latexmk -bibtex -pdf -use-make -synctex=1 \
		 -pdflatex="pdflatex -interaction=nonstopmode" $(MASTER).tex 

$(MASTER).tex: $(RNW) $(MASTER).Rnw title.tex 
	 Rscript -e "library(knitr); knitr::knit('$(MASTER).Rnw')"

%.tex: %.Rnw
	Rscript -e "library(knitr); knitr::knit('$<')"
	
%.R: %.Rnw
	Rscript -e "library(knitr); knitr::purl('$<',documentation=0L)"

%.pdf: %.tex
	latexmk -bibtex -pdf -use-make -synctex=1 \
		 -pdflatex="pdflatex -interaction=nonstopmode" $< 


clean:
	rm -f *~   *.out Rplots.pdf comment.cut \
              *.idx *.ilg *.brf *.blg *.spl  $(MASTER).dvi \
              *.backup *.toc *.fls  *fdb_latexmk *.synctex.gz  *-concordance.tex

cleanreallyall: clean
	rm -f  *.out Rplots.pdf comment.cut
	rm -f  *.idx *.ilg *.brf *.blg *.spl  $(MASTER).dvi
	rm -f  *.backup *.toc *.fls  *fdb_latexmk *.synctex.gz  *.concordance.tex
	rm -f *.aux *.log *.ind *.ist
	rm -f figure/ch??_fig*.pdf
	rm -f chapter??.tex $(MASTER).bbl  $(MASTER).tex $(MASTER).pdf
	rm -f cache/*.* cache/__packages
	echo "\nRemoved all files to precompilation status\n"


# 2018 Reinhard Furrer, Florian Gerber, Josef Stocker, Roman Flury

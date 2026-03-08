EPS = $(wildcard *.eps)

%.pdf: %.tex $(EPS)
	( \
	\pdflatex $<; \
	while \grep -q "Rerun to get cross-references right." $(<:.tex=.log); \
	do \
		\pdflatex $<; \
	done \
	)

TEX = $(wildcard *.tex)
PDF = $(patsubst %.tex,%.pdf, $(TEX))
LOG = $(patsubst %.tex,%.log, $(TEX))
AUX = $(patsubst %.tex,%.aux, $(TEX))
NAV = $(patsubst %.tex,%.nav, $(TEX))
SNM = $(patsubst %.tex,%.snm, $(TEX))
OUT = $(patsubst %.tex,%.out, $(TEX))
TOC = $(patsubst %.tex,%.toc, $(TEX))
SPL = $(patsubst %.tex,%.spl, $(TEX))

all: $(PDF)

clean:
	rm -f $(PS) $(PDF) $(LOG) $(AUX) $(NAV) $(OUT) $(TOC) $(SNM) $(SPL)

#! /bin/bash

fname=report

shopt -s extglob
rm -v phd_report!(.tex) 2>/dev/null
shopt -u extglob

pdflatex $fname
biber $fname
#bibtex $fname
pdflatex $fname
makeglossaries $fname
makeindex ${fname}.nlo -s nomencl.ist -o ${fname}.nls
makeindex ${fname}.idx
pdflatex $fname
pdflatex $fname

#!/bin/bash

cat header.tex > results.tex
sh printSum.sh >> results.tex
echo "\\hline" >> results.tex 
echo "\\end{tabular} " >> results.tex 
echo "\\end{sidewaystable}" >> results.tex 
echo "\\end{document}" >> results.tex 
pdflatex results.tex 
evince results.pdf &


#!/bin/bash

gnuplot lbImprSprint.gp
gnuplot lbImprLong.gp
gnuplot lbImprMed.gp

epstopdf lbImprSprint.eps
epstopdf lbImprLong.eps
epstopdf lbImprMed.eps

cp lbImpr*.pdf ~/artigos/patat2012/img/

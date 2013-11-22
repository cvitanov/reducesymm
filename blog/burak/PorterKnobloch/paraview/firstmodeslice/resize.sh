#!/bin/bash

#source: http://www.ehow.com/how_6823473_reduce-pdf-file-size-linux.html

gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -sOutputFile=gorbitsandslicelofi.pdf gorbitsandslice.pdf

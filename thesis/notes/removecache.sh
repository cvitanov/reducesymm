# Downloaded from Michele Vallisneri's webpage, www.vallis.org
#!/bin/sh
gs -q -dBATCH -dNOPAUSE -dNOCACHE -sDEVICE=epswrite -sOutputFile=$2 $1 -c quit

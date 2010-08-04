#!/bin/bash
# Overwrites all .php files in current directory with ``translated'' php files from src directory
# in local server(my laptop). That is, the downloaded files will not be dynamic. Use with caution
# to avoid overwritting edited files. For the moment images, etc are not updated by this script.

ls *.php | while read i
do
	        curl -O  http://127.0.0.1/siminos/www/src/$i 
done

# needed to work on cns server
cp index.php index.html
# cp auxiliary files
cp ../src/css/*.css css/
cp ../src/css/*.js css/

ES 					May-16-2010

ES homepage php sources. Require local server with php support to generate html content,
or can be run directly in a php enabled server. Edit the sources here, then 

cd  ../public_html 
./update

Manually copy images, etc. 

contains:

.	: All content pages in top level, to make use of php includes easier (not optimal)
include : banners, headers, etc. Write/edit once, include in all content pages.
css	: cns css style, local additions, css(+js) style for collapsible text.
bib	: bibtex files, to be read in publication list. Future: convert to xml, read by php.
images	: all images

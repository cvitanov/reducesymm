reducesymm/bibtex/00ReadMe.txt
===========================================================

following files home is here
    reducesymm/bibtex/QFT.bib

    reducesymm/bibtex/DasGroup.bib  % include, harmonize with:
                      [x] dasgroup/thepaper/groups-bibliography.bib
                      [x] 2016-12-01 Alcock-Zeilinger GroupTheory.bib
                      [x] 2016-12-01 extracted reducesymm/bibtex/QFT.bib parts
                      [.] dasgroup/book/chapter/refs.tex
                      [ ] dasgroup/book/toInclude/refsInsert.tex

    reducesymm/bibtex/diffuse.bib   (Zhang entries)
    reducesymm/bibtex/linresp.bib   (McInroe  entries)

following files are copied to here, and NOT TO be edited here
    always edit siminos/bibtex/siminos.bib
then copy to here
    always edit dasbuch/bibtex/ChaosBook.bib
then copy to here

bibtex/bibreduce
	perl script that examines article.aux, creates a list of
	the papers that are cited in article.tex and extracts them
	from large.bib to small.bib.

NOTES
-----

reducesymm/bibtex/QFT.bib
2013-11-23 Predrag % Partially extracted from
    siminos/bibtex/siminos.bib

see also svn
    dasgroup/thepaper/groups-bibliography.bib
which already includes
    dasgroup/WWW/PUPress/WestburyVogel.bib

----------

Zoteromania: do not rename any existing citation, either
in siminos.bib or when you are adding a bibTeX entry to
a citation copied from ChaosBook.org.

When submitting to arXiv, journals,
please DO NOT generate yet another *.bib:
	submit article_name.bbl file

why:
	1) purging main BibTeX file is waste of time
	2) it's tedious to keep all versions updated

Fix these:
----------
                    Predrag Mar 10 2007

-- do not disable hyperlinking in bibtex/siminos.bib
   improve it, not remove it


==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= was siminos/bibtex/00ReadMe.txt until copied to git 2013-07-07

= DONE:
check http://arxiv.org/hypertex/bibstyles/ for hyperref compatible
BiBTeX styles

= DONE:

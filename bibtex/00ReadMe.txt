reducesymm/bibtex/00ReadMe.txt
===========================================================
Burak Budanur and Evangelos Siminos

    remember to diff and merge with siminos/bibtex/siminos.bib
    every so often...

              ***
Maintain only ONE bib file for all dynamical systems publications
New entries at THE TOP of this file
              ***

	bibtex/siminos.bib
		all references should go here
    bibtex/05bibTools.txt
        notes about BibTeX tools, web resources, etc
    bibtex/zotero.txt
        if we return to zotero...

bibtex/bibreduce
	perl script that examines article.aux, creates a list of
	the papers that are cited in article.tex and extracts them
	from large.bib to small.bib.

NOTES
-----

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

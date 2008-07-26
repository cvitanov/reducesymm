/bibtex/00ReadMe.txt
$Author$ $Date$
===========================================================

              ***
Maintain only ONE bib file for all dynamical systems publications
              ***

Evangelos Siminos 
	bibtex/siminos.bib

bibtex/05bibTools.txt
	notes about BibTeX tools, web resources, etc

bibtex/bibreduce: 
	perl script that examines article.aux, creates a list of 
	the papers that are cited in article.tex and extracts them
	from large.bib to small.bib.


NOTES
-----

When submitting rpo.tex, remember to
    1) Submit rpo.bib file. 
       rpo_ks/rpo.bib is always a bibreduce subset of siminos.bib

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

= DONE: 
check http://arxiv.org/hypertex/bibstyles/ for hyperref compatible
BiBTeX styles

= DONE:  


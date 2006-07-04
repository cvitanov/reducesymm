/home/predrag/articles/siminos/thesis/notes/00ReadMe.txt
	Evangelos Siminos			jul  3 2006

===========================================================

----------
Processing of eps figures with mathematica [and other] fonts

PROBLEM:
> dvips rpo
dvips: Font Mathematica1Mono used in file figs/rpo22-55-4.eps is not in the mapping file.
dvips: Font Mathematica2Mono used in file figs/rpo22-55-4.eps is not in the mapping file.

SOLUTION:

./removecache.sh rpo22-55-4.eps rpo22-55-4-clean.eps

Mathematica uses its own fonts for any mathematical symbol, even + and -. 
Script removecache.sh fixes the problem by converting the fonts to curves.
The resulting figure is much smaller

> ls -s
168 rpo22-55-4.eps
 40 rpo22-55-4-clean.eps

but have to inspect it - fonts do look different


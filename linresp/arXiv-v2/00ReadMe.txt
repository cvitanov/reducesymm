lippolis/hyperb/arXiv-v2/00ReadMe.txt
$Author: predrag $ $Date: 2009-11-16 17:02:17 -0500 (Mon, 16 Nov 2009) $
===========================================================
                    arXiv resubmission: Domenico Nov ?? 2009
                    finished rewriting: Predrag  Jan  5 2009
                    finished rewriting: Domenico Jan  5 2009
                    started writing:  Predrag  Aug 24 2008

How well can one resolve the \statesp\ of a chaotic map?
Domenico Lippolis and Predrag Cvitanovi\'c

-----------------------------------------------------------

arXiv re-submission notes
----------------------

in the directory above:

[NO] this was not done by Domenico:
cp hyperb.tex pruned_prl.tex
	prune all extraneous text, comments, etc,
	insert this as preface to pruned_prl.tex:
			 %
% hyperb.tex
% $Author: predrag $ $Date: 2009-11-16 17:02:17 -0500 (Mon, 16 Nov 2009) $
%
        \newif\ifdraft \draftfalse  % draft version is default
        \newif\ifpreparepdf \preparepdftrue % hyperlinked pdf default
%% uncomment for final version: %%%%%%%%%%%%%%%%%%%%%%%
% \draftfalse
% \preparepdffalse % for B&W, print version
%
%% ------------------ for arXiv resubmission ----------------------------
%
% Title:    How well can one resolve the state space of a chaotic flow?
% Authors:  Domenico Lippolis and Predrag Cvitanovic
% Comments: 4 pages, 3 postscript figures, uses revtex4
% Files:    hyperb.tex hyperbDefs.tex hyperb.bbl
%           hyperbMarkov3lbld.eps  repOverlap.eps  synth.eps
%
%% ------------------ cut here ----------------------------------------
    [NO] mark here when done [Domenico Nov ?? 2009]

cp pruned_prl.tex arXiv-v2/hyperb.tex
	then in arXiv-v2/hyperb.tex cut everything above line marked
%% ------------------ cut here ----------------------------------------
    [ ] mark here when done [Domenico Nov ?? 2009]

remove extraneous macros from hyperbDefs.tex before submission
    [ ] mark here when done [Domenico Nov ?? 2009]

resubmit arXiv, version 2
    [X] mark here when done [Domenico Nov 16 2009]

domenico@gatech.edu
0902.4269
user/password combination for this article is

 User-ID: 0902.4269
 Password: xmpqn

save arXiv confirmation email as arXiv-v2/arXivAckn.txt
    [ ] mark here when done [Predrag Nov 12 2009]

To submit this article to an overlay journal based on arXiv you
may need to supply the following identifier:

 arXiv:tracking/59b1339b13da4505

-----------------------------------------------------------------

PRLett submission notes
---------------------

see lippolis/hyperb/PRlett-v2/00ReadMe.txt

===============================================================

NOTES
-----

Preparation of figures
----------------------

Fix these:
----------

==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= DONE:

svn rm  [all unused figures and all *.pdf]
    [ ] mark here when done [Predrag Nov 12 2009]

svn rm [all unused *.tex and other files]
    [ ] mark here when done [Domenico Nov ?? 2009]

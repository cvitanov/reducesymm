siminos/CLE/00ReadMe.txt
$Author$ $Date$
===========================================================
Continuous symmetry reduction
    and return maps for higher-dimensional flows
    Evangelos Siminos and Predrag Cvitanovi\'c

% Predrag, Evangelos, John last push        2010-01-30
% Predrag reorganized siminos/CLE/          2009-10-09
% Vaggelis created siminos/CLE/CLE.tex      2009-09-21
% Vaggelis created siminos/CLE              2009-02-05
            submission:       Vaggelis      2010-01-30


					Predrag                         2009-12-02
To get a move on this article, I propose that we respond to
the invitation, submit it to

    05PhysicaD.txt

and finalize it before Dec 21, 2009.
After that date Predrag teaches and things slow down.

-----------
            submission DEADLINE                     2010-01-31
for submission, read and enter all details in
    siminos/CLE/PhysD-v1/00ReadMe.txt

PhysicaD manuscript #???

arXiv.org/abs/1209.???
arXiv:????.???? (submitted v1):                     2010-02-??
    siminos/CLE/arxiv-v1/   low resolution figs
    siminos/CLE/PhysD-v1/   contains the submission
    siminos/CLE/PhysD-v1/reviews/   referee edits

Processing
----------
    make sure that you are in ../siminos/CLE, then
> ./update

Directories
-----------
until submission
    ../figs
    defCLE.tex   redefinitions specific to CLE.tex.
    ../bibtex/siminos   submit only rpo.bib file

please do not create a `current' subdirectory; work here,
save versions in subdirectories...

Working titles:
-----------
- Continuous symmetry reduction
    and return maps for high dimensional flows
- Continuous symmetry reduced return maps
    for high dimensional flows
- Returm maps of high dimensional flows with
    continuous symmetry: I. Symmetry reduction
- Towards continuous symmetry reduction
    in high dimensional flows



Submission
----------

=== submit arXiv - see arxiv-v1/01ReadMe.txt
    [X] check here when done, enter date:           Dec ?? 2009


=== submit APS  - see APS-v1/00ReadMe.txt
    [ ] check here when done, enter date:           Dec ?? 200?

    APS Journal of Blah
    www.???.org/journals/
    epubs.???.org/APS


NOTES
-----
                    Predrag                         Dec 29 2009
=== acknowledge U Chicago support
    Thomas Witten:
    Gene will tell you about the funding acknowledgment.
    Mention the james franck institute.

CLE.tex now formated for Physica D          20 dec 2009
==================================
                    Vaggelis                        Mar 10 200?
reformatted using APS macros
    www.???.org/journals/

check http://arxiv.org/hypertex/bibstyles/ for hyperref compatible
BiBTeX styles


When submitting, remember to
    1) Submit rpo.bib file.
    2) All PostScript figures should be sent in separate files.
    3) PostScript figures should be generated with sufficient
       line thickness.
							

------------------------------------------------------------------



Fix these:
----------

- bluh

==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= DONE:
  At submission, merge with defs.tex, prune

= DONE: - Use next command to get References header appear.
  The problem is a conflict between elsarticle and amsref
  over \bibsection command.
    \newcommand\bibsection{%
    \section*{\bibname\markright{\MakeUppercase{\bibname}}}}

= DONE: \input ../inputs/def        % TEMPORARY
        \input inputs/defsCLE.tex
 % eventually merge used commanDs into inputs/defsCLE.tex

siminos/atlas/Chaos-v1/00ReadMe.txt
$Author$
$Date$
------------------------------------------------------------
submit before     April 16, 2012

Slice & chart pedagogical introduction
====================
	Reduction of continuous symmetries of chaotic flows
	by the method of slices

Predrag Cvitanovic et al

-------------------------------------------------------------
Predrag ver. 0.2, 						Feb 22 2012
Predrag ver. 0.1, 						Dec 10 2011

Focus issue for Chaos: "50 Years of Chaos"
chaos.peerx-press.org website is currently open to accept papers for this
focus issue, and you may submit your papers at any point through the
manuscript DUE DATE,

                April 16, 2012.

-------------------------------------------------------------

- ../figs/
	all figs
- ../figsSrc/
	all fig source programs
- ../bibtex/siminos.bib

----- NOTES        ------------------------------------------

[ ] www.aip.org/pubservs/compuscript.html#prepare
[ ] style files work with REVTeX 4.1 and are included in REVTeX 4.1 package
    tex.stackexchange.com/tags/revtex/info
    [ ] udpated 12/27/2011, install this version or later
[ ] The Manuscript, including the abstract, references, and captions,
    should be set up for 21.6 x 28 cm (8-1/2 x 11 in. or A4) paper with
    ample margins. It is essential that the motivations, central results,
    and conclusion be stated in a nontechnical manner that is
    intelligible to a broad audience.
[ ] Be sure that your manuscript contains page numbers.
[ ] The title page should contain the title of the article, the names of
    the authors, a suitable byline, and a short abstract.
[ ] Parts of the manuscript should be arranged in the following order:
    title, author, affiliation, abstract, text, acknowledgments,
    appendices, and references.
[ ] The title of a paper should be as concise as possible but informative
    enough to facilitate information retrieval.
[x] The Abstract should be self-contained (contain no footnotes). It
    should be adequate as an index (giving all subjects, major and minor,
    about which new information is given) and as a summary (giving the
    conclusions and all results of general interest in the article). It
    should be about 5% of the length of the article, but less than 500
    words. The abstract should be written as
    [ ] one paragraph
    [ ] should not contain displayed mathematical equations
[x] Authors' names should preferably be written in a standard form for
    all publications to facilitate indexing and avoid ambiguities.
[ ] Authors with Chinese names may choose to have their names published
    in their own language alongside the English versions of their names
    in the author list of their publications, in either Simplified or
    Traditional characters: www.aip.org/pubservs/cjk_instructions.html

[x] The first paragraph of the article should be a Lead Paragraph and
    will be highlighted in the journal in boldface type. This paragraph,
    which essentially advertises the main points of the article, must
    describe in terms accessible to the nonspecialist reader the context
    and significance of the research problem studied and the importance
    of the results. The Editors will pay special attention to the clarity
    and accessibility of this paragraph, and in many cases may rewrite it

[ ] References and footnotes should be in the form shown in recent issues
    of this journal. They should be given in a double-spaced list at the
    end of the text.
    [ ] The names, including initials, of all authors in
    each reference should be given (in the text, the use of et al. is
    permissible).
    [ ] By number, in the order of first appearance, giving the names of
    the authors, the journal name, volume, year, and first page number
    only, as in:
    V. Bargmann, Proc. Natl. Acad. Sci. USA 38, 961 (1952).
[ ] For footnotes to title and bylines use a), b), c), etc.
    Avoid lengthy footnotes by inserting them in the text, except for the
    references.

[ ] adhere to these guidelines when preparing your illustrations for
    submission:
    chaos.aip.org/authors/information_for_contributors/how_to_prepare_illustrations
    [ ] figures must be in the final published size, not oversized.

[ ] Multimedia files can be included in the online version of published papers.
    chaos.aip.org/authors/information_for_contributors/guidelines_for_multimedia


-------------------------------------------------------------

Chaos resubmission  (here for future use)
----------------
                                        2011-04-28
Chaos referee responses, edits
----------------------------
Predrag:						      Apr 1 2009
    generate colored comments for referees,
        PCedit\{...\}, etc in defslice.tex,
    then go back to B&W for the final resubmission

Chaos submission
----------------

Submit  via http://chaos.peerx-press.org.
[ ] Focus issue for Chaos: "50 Years of Chaos"
[ ] Cover letter
        siminos/atlas/letters/submLett.doc
    title, authors, contact information, and mention this Focus issue
    Journal, corresponding author's e-mail address, any special requests.
[ ] a double-spaced manuscript/article file in LaTeX/REVTeX.
    [ ] a TeX file must be one file - no support for call-ins, BibTeX.
[ ] create your own .bbl
    [ ] include it as the bibliography section in the main .tex file.
[ ] Commands to include figures may be used. Ensure that the figure
    filename cited in the command matches that of the actual file upload;
    use only the simple filename, not a complete directory path;
    [ ] invitation letter: individual figures in high-quality PDF format
    [x] website: Color will appear in the online journal free
        of charge. A usable color graphics file: [...] or *.pdf
    [ ] the printed version: black-and-white
        [ ] a descriptive term other than a color is needed in the
        caption to support the data of discussion: instead of `the red
        and blue symbols' write `the red circles and blue squares,' etc.
    [?] website: only include figures in .eps format (not .tif or .pdf).
        If an EPS file contains a thumbnail preview, it may not render in
        the article PDF generated via the TeX (Predrag - ignore this)
[ ] Do not use the "\footnote" command; use the standard "\cite" command
    for footnotes as well as references.
[ ] http://chaos.peerx-press.org, fill out Transfer of Copyright Agreement.

To redo the TeX processing, e.g., after uploading a new or
    revised figure file, do a Replace on the .tex Article File. Doing a
    Rebuild of the Merged PDF file does NOT initiate another TeX
    processing.

For any regular published article that exceeds 12 journal pages, a
mandatory $150 page fee will be added for each page in excess of 12
pages. Read chaos.aip.org/chaos/authors/publication_charges



arXiv submission
----------------
    read siminos/atlas/arxiv-v1

===========================================================
Processing
----------

make sure that you are in repository siminos/atlas/ ,
	then, if in linux:
./update
    or
pdflatex atlas; bibtex atlas; pdflatex atlas; pdflatex atlas

Things to fix
-------------
    [when fixed, move the line to the FIXED section at the end]

= see 00edits.txt


Style files
-----------

Preparation of figures
----------------------

==============================================================
FORMERLY OUTSTANDING ITEMS, NOW DISPOSED OFF:

= Working titles
-- Frank Zappa and The Mothers of Invention
= Submission of custom .sty files is not supported.

%% Paket zum erstellen von Pr�entationspostern
%% Thomas Bartsch 07.04.2004:
%%    Titellayout und Platzierung des Logos veraendert
%%    Laengenparameter logowidth eingefuehrt

%% Ausgabe auf RUS-Posterdrucker:
%%
%% latex <datei>.tex
%% dvips -Ppdf <datei>.dvi
%% ps2ps -sPAPERSIZE=a4 -r2400 <datei>.ps <datei>-2.ps
%%
%% <datei>-2.ps dann mittels sendweb (s. u.) an den RUS-Drucker
%% HP DesignJet 2800 CP/PosterJet
%% schicken mit einer Skalierung von 400%
%%
%% sendweb: 
%%   http://www.uni-stuttgart.de/rus/drucken
%%   auf `sendweb' clicken
%%   Seitenformate gr�er A3 anw�len (Gruppe 3?)
%%   Papier ausw�len
%%   Drucker ausw�len (s. o.)

\ProvidesClass{theo1poster}[2003/04/25 0.4a - LaTeX2e presentation-poster class]

\LoadClass[10pt]{scrartcl}[1999/12/29]
%% todo: Weitergabe von Optionen an Klassen und Pakete
%%       Vorgabe des Abstandes zum Boxenrand durch Benutzer
%%          und automatisches Berechnen der Minipage-Breite

%% Layout-Befehle:
%%
%% \posterbgcolor{}  : Setzen der Hintergrundfarbe des Posters
%%                     Default: Normale Hintergrundfarbe
%% \headingfont{}    : Setzen des Fonts fr die Abschnittsberschriften
%%                     Default: \sffamily\bfseries
%% \headingfontsize{}: Gr�e der Sheet-�erschriften
%% \headingcolor{}   : Farbe der Sheet-�erschriften
%% \sheetframecolor{}: Farbe der Sheet-Rahmen
%% \sheetframewidth{}: Breite der Sheet-Rahmen
%% \sheetbgcolor{}   : Farbe der Sheet-Hintergrnde
%% \titlecolor{}     : Farbe fr den Titel
%% \logo{}           : Angabe eines Logos, wird links neben Autor gesetzt
%% \logowidth        : Laengenparameter fuer Breite des Logos
%% \address{}        : Angabe einer Adresse
%% \(no)colorsheets  : Farbboxen fr die Sheets oder nicht
%% \roundsheets      : Abgerundete Ecken, ohne Farbe
%% \roundcolorsheets : Abgerundete Ecken, farbig

%% Weitere Befehle
%%
%% \subheading{}     : Unterberschrift
%% \highlight[]{}    : Farbige Hinterlegung
%% \highlightcolor{} : Festlegung der Default-Hervorhebungsfarbe

%% Umgebungen:
%%
%% poster            : 1 optionales Argument: Breite (default: Textbreite)
%%                     1 obligatorisches Argument: Anzahl an Spalten
%% sheet             : 1 optionales Argument: Breite (default: Linienbreite)
%%                     1 obligatorisches Argument: �erschrift der Seite
%%
%% cequation[*]      : Farbig hinterlegte Formel
%%                     1 optionales Argument: Hintergrundfarbe


%% Laden von ben�igten Paketen

\RequirePackage{vmargin}
\RequirePackage{multicol}
%\RequirePackage{pstricks}
%\RequirePackage{color}
\RequirePackage{pstcol}
%\RequirePackage{graphicx}

\setmarginsrb{6mm}{6mm}{6mm}{6.5mm}{0pt}{0pt}{0pt}{0pt}


%% Einstellen der Schriftgr�en auf 6pt

\fontsize{6}{7}\selectfont

\renewcommand{\normalsize}{\fontsize{7}{8pt}\selectfont
   \abovedisplayskip 6\p@ \@plus1\p@ \@minus3\p@
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 3.5\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini %% Added 22 Dec 87
   \parsep 2pt plus 0.5pt minus 0.5pt 
         \itemsep 0.5\parsep
   \topsep 5.5pt plus 1pt minus 1.5pt}}

\renewcommand{\small}{\fontsize{5}{7pt}\selectfont
   \abovedisplayskip 5\p@ \@plus1\p@ \@minus2.5\p@
   \belowdisplayskip \abovedisplayskip
   \abovedisplayshortskip \z@ \@plus1.5\p@
   \belowdisplayshortskip 3\p@ \@plus1.5\p@ \@minus1.5\p@
   \def\@listi{\leftmargin\leftmargini %% Added 22 Dec 87
   \parsep 1pt plus 0.5pt minus 0.5pt 
         \itemsep \parsep
   \topsep 3.5pt plus 1pt minus 1.5pt}}

\renewcommand{\footnotesize}{\fontsize{5}{6pt}\selectfont}
%\renewcommand{\scriptsize}{\fontsize{3}{5pt}\selectfont}
\renewcommand{\scriptsize}{\fontsize{3}{4pt}\selectfont}
%\renewcommand{\tiny}{\fontsize{3}{4pt}\selectfont}
\renewcommand{\tiny}{\fontsize{3}{3.6pt}\selectfont}
%\renewcommand{\large}{\fontsize{10}{14pt}\selectfont}
\renewcommand{\large}{\fontsize{10}{12pt}\selectfont}
%\renewcommand{\Large}{\fontsize{14.4}{17pt}\selectfont}
\renewcommand{\Large}{\fontsize{14.4}{17.3pt}\selectfont}
%\renewcommand{\huge}{\fontsize{14.4}{20pt}\selectfont}
\renewcommand{\huge}{\fontsize{14.4}{18pt}\selectfont}
%\renewcommand{\Huge}{\fontsize{20.74}{28pt}\selectfont}
\renewcommand{\Huge}{\fontsize{20.74}{25pt}\selectfont}

%%% �derung von Gliederungsbefehlen

\renewcommand\section{\@startsection{section}{1}{\z@}%
  {-3.5ex \@plus -1ex \@minus -.2ex}%
  {0.6ex \@plus.2ex}%
  {\centering\@headingcolor\sectfont\@headingfontsize}}
      %\raggedsection\normalfont\sectfont\nobreak\size@section\nobreak}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {0.6ex \@plus .2ex}%
  {\centering\@headingcolor\@headingfontsize\centering\itshape}}

%\newcommand*\subheading[1]{\vspace{.4\baselineskip}%
%  {\@headingcolor\@headingfontsize%
%  \centering\itshape #1\vspace{.2\baselineskip}\par}}


%\renewcommand\section{\@startsection{section}{1}{\z@}%
%  {-3.5ex \@plus -1ex \@minus -.2ex}%
%%  {2.3ex \@plus.2ex}%
%  {0.2ex}%
%  {\mathversion{bold}\centering\normalfont\@headingfontsize\sectfont}}
%%    \mathversion{bold}\raggedsection\normalfont\@headingfontsize\sectfont}}

%%
%% Kein Section fr das Literaturverzeichnis
%%
\renewenvironment{thebibliography}[1]{%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \setlength{\topsep}{0pt}%
            \setlength{\partopsep}{0pt}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand*\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}


%% Logo, Adresse

\renewcommand*\@title{}
\renewcommand*\@author{}

\newcommand*\@logo{}
\newcommand*\logo[1]{%
  \renewcommand*\@logo{\includegraphics[width=\logowidth]{#1}}}

\newcommand*\@address{}
\newcommand*\address[1]{%
  \renewcommand*\@address{#1}}


%% Farbe fr den Hintergrund

\newcommand*\@posterbgcolor{}
\newcommand*\posterbgcolor[1]{%
  \renewcommand*\@posterbgcolor{\pagecolor{#1}}}


%% Farbe und Font fr die Sheet-�erschriften

%\newcommand*\@headingfont{\sffamily\bfseries}
\renewcommand*\sectfont{\sffamily\bfseries\mathversion{bold}}
\newcommand*\@headingfontsize{\large}
\newcommand*\@headingcolor{}
\newcommand*\headingfont[1]{%
  \renewcommand*\sectfont{#1}}
\newcommand*\headingfontsize[1]{%
  \renewcommand*\@headingfontsize{#1}}
\newcommand*\headingcolor[1]{%
  \renewcommand*\@headingcolor{\color{#1}}}


%% Farbe und Font fr den Titel

\newcommand*\@titlefontsize{\Huge}
\newcommand*\@titlecolor{}
\newcommand*\titlecolor[1]{%
  \renewcommand*\@titlecolor{\color{#1}}}


%% Farbbox oder nicht fr die Sheets?

\newcommand\@makesheet[3]{%
  \psframebox[linewidth=\fboxrule,%
  framesep=\@sheetframesep]{#3}}      %% Default: Schwarzwei� Sheets

\newcommand*\colorsheets{%            %% Auf Wunsch Farbig
  \renewcommand\@makesheet[3]{%
    \psframebox[linewidth=\fboxrule,framesep=\@sheetframesep,%
    linecolor=##1,fillstyle=solid,fillcolor=##2]{##3}}}

\newcommand*\roundcolorsheets{%
  \renewcommand\@makesheet[3]{%
    \psframebox[linearc=2.5mm,cornersize=absolute,linewidth=\fboxrule,%
    linecolor=##1,fillstyle=solid,fillcolor=##2,%
    framesep=\@sheetframesep]{##3}}}

\newcommand*\roundsheets{%
  \renewcommand\@makesheet[3]{%
    \psframebox[linearc=2.5mm,cornersize=absolute,linewidth=\fboxrule,%
    framesep=\@sheetframesep]{##3}}}

\newcommand*\nocolorsheets{           %% Oder wieder schwarzwei�und eckig
  \renewcommand\@makesheet[3]{%
    \psframebox[linewidth=\fboxrule,framesep=\@sheetframesep]{##3}}}
%% todo: Klassenoption nocolor, die alle Farben abstellt


%% Titelei

\renewcommand\maketitle{%
  \settowidth{\@authorwidth}{\large\@address}%
  {\centering\sectfont\@titlefontsize\@titlecolor\@title\par}%
  \vspace{2.5ex}%
  \hfill\parbox[c]{0cm}{\makebox[0cm][c]{\@logo}}
  \hfill%
  \parbox[c]{\@authorwidth}{%
     {\centering\large\@author\par\vspace{1ex}\@address\par}%
  }%
  \hfill\hfill%
  \vspace*{-2ex}%
}


%% Farbe fr Sheet-Rahmen und -Hintergrund

\newcommand*\@sheetframecolor{black}
\newcommand*\@sheetbgcolor{white}
\newcommand*\sheetframecolor[1]{%
  \renewcommand*\@sheetframecolor{#1}%
  \colorsheets}
\newcommand*\sheetbgcolor[1]{%
  \renewcommand*\@sheetbgcolor{#1}
  \colorsheets}


%% Dicke fr die Sheet-Rahmen und Abstand zum Rahmen

\setlength{\fboxrule}{1mm}
\newcommand*\sheetframewidth[1]{%
  \setlength{\fboxrule}{#1}%
  \setlength{\@autosheetinnerwidth}{\linewidth}%
  \addtolength{\@autosheetinnerwidth}{-2\fboxrule}%
  \addtolength{\@autosheetinnerwidth}{-2\@sheetframesep}}
\newcommand*\sheetframesep[1]{%
  \setlength{\@sheetframesep}{#1}%
  \setlength{\@autosheetinnerwidth}{\linewidth}%
  \addtolength{\@autosheetinnerwidth}{-2\fboxrule}%
  \addtolength{\@autosheetinnerwidth}{-2\@sheetframesep}}


%% Ben�igte L�gen und Boxen

\newlength{\@authorwidth}
\newlength{\logowidth}
\setlength{\logowidth}{1cm}
\newlength{\@parindentlength}

\newlength{\@sheetframesep}
\setlength{\@sheetframesep}{3pt}
\newlength{\@autosheetinnerwidth}  % Berechnete L�ge fr den Innenraum

\newsavebox{\@sheetsavebox}


%% Umgebung sheet

\newenvironment{sheet}[2][\@autosheetinnerwidth]{%
  %% todo: automatische Berechnung der Breite der minipage
  \begin{lrbox}{\@sheetsavebox}
    \begin{minipage}[t]{#1}
      \section*{#2}
%      {\centering\@headingcolor\sectfont\@headingfontsize #2\par}%
%      \section{#2}
%      \vspace{.3\baselineskip}
      \setlength{\parindent}{.5\@parindentlength}%
}{%
    \end{minipage}
  \end{lrbox}
  \@makesheet{\@sheetframecolor}{\@sheetbgcolor}{%
    \usebox{\@sheetsavebox}}\par
  \vspace{1ex plus 0.75ex minus 0.5ex}
  %% todo: Wahlm�lichkeit zwischen fester und Gummil�ge
}


%% Umgebung poster

\newenvironment{poster}[2][\textwidth]{%
  \setlength{\@parindentlength}{\parindent}%
  \setlength{\parindent}{0pt}%
  \newpage\thispagestyle{empty}\@posterbgcolor%
  %% todo: Bessere Bestimmung der Breite
  \begin{minipage}[c]{#1}%
    \maketitle%
    \begin{multicols}{#2}%
      \sheetframewidth{\fboxrule}  % Innenraum der Sheets berechnen
}{%
  \end{multicols}\end{minipage}%
}


%%
%% Zus�zliche Befehle, um das Setzen von Text auf Postern zu vereinfachen
%%

%% Unterberschriften

%\newcommand*\subheading[1]{\vspace{.4\baselineskip}%
%  {\@headingcolor\@headingfontsize%
%  \centering\itshape #1\vspace{.2\baselineskip}\par}}
\newcommand*\subheading[1]{\subsection*{#1}}


%% Befehl \highlight

\newcommand\highlight[2][\@defaulthighlightcolor]{%
  \colorbox{#1}{#2}}

%% Befehle fr colorenv und \highlight: Standard-Hervorhebung

\newcommand*\@defaulthighlightcolor{yellow}
\newcommand*\highlightcolor[1]{%
  \renewcommand*\@defaulthighlightcolor{#1}}

%%
%% Zus�zliche Umgebungen, um das Setzen von Text auf Postern zu vereinfachen
%%

%% Umgebung colorenv

\newsavebox{\@colorenvsavebox}
\newcommand*{\@colorenvsavecolor}{}

\newenvironment{colorenv}[1][\@defaulthighlightcolor]{%
  \renewcommand*{\@colorenvsavecolor}{#1}%
  \begin{lrbox}{\@colorenvsavebox}%
}{%
  \end{lrbox}%
  \colorbox{\@colorenvsavecolor}{\usebox{\@colorenvsavebox}}}

%% Umgebung cequation[*]

%\newsavebox{\@cequationsavebox}
%\newcommand*{\@cequationsavecolor}{}

%\newenvironment{cequation}[1]{%
%  \renewcommand{\@cequationsavecolor}{#1}
%  \begin{lrbox}{\@cequationsavebox}
%    \ensuremath\displaystyle}
%{%
%  \end{lrbox}
%  \begin{equation}
%    \colorbox{\@cequationsavecolor}{\usebox{\@cequationsavebox}}
%  \end{equation}
%}

\newenvironment{cequation}[1][\@defaulthighlightcolor]{%
  \begin{equation}%
    \begin{colorenv}[#1]%
      \begin{math}\displaystyle}%
{%
  \end{math}\end{colorenv}\end{equation}\@ignoretrue}

%\newenvironment{cequation*}[1]{%
%  \renewcommand{\@cequationsavecolor}{#1}
%  \begin{lrbox}{\@cequationsavebox}
%    \begin{math}\displaystyle}
%{%
%  \end{math}\end{lrbox}
%  \begin{equation*}
%    \colorbox{\@cequationsavecolor}{\usebox{\@cequationsavebox}}
%  \end{equation*}
%}

\newenvironment{cequation*}[1][\@defaulthighlightcolor]{%
  \begin{equation*}
    \begin{colorenv}[#1]
      \begin{math}\displaystyle}
{\end{math}\end{colorenv}\end{equation*}\@ignoretrue}


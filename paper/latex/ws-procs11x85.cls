%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ws-procs11x85.cls   :   19 Oct 2017
%% Class file for Proceedings Trim Size [11in x 8.5in] written in Latex2E. 
%% The content, structure, format and layout of this style file is the
%% property of World Scientific Publishing Co. Pte. Ltd.
%%
%% Trim Size: 11in x 8.5in
%% Text Area: 9.25in (include runningheads) x 6.6in
%% Main Text: 12/15pt
%%
%% For Technical support please email to Rajesh <rajesh@wspc.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{amsfonts,amssymb,amsmath}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ws-procs11x85}[2000/12/15 v1.2u Standard LaTeX document class]
\@twosidetrue\@mparswitchtrue%\@twocolumnfalse
\newif\if@openright
\@openrighttrue
\newif\if@restonecol
%
\newif\if@addchapnum
\newif\if@wsdraft
\newif\if@wssquare
%
\newdimen\wsdraftrule
\wsdraftrule=0pt
%
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{splissue}{\@splissuetrue}
\DeclareOption{wsdraft}{\@wsdrafttrue\setlength\overfullrule{5pt}\wsdraftrule0.01pt}%
%
\hoffset=-2.6cm
\voffset=-2cm
%
\DeclareOption{square}{\@wssquaretrue}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
%
\ExecuteOptions{a4paper,10pt,twoside,onecolumn,final}
\ProcessOptions
%
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\newcommand\bibname{References}%
%
\newcommand\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
    \space\number\day, \number\year}
\newcount\minute
\newcount\hour
\def\currenttime{%
    \minute\time
    \hour\minute
    \divide\hour60
    \the\hour:\multiply\hour60\advance\minute-\hour\the\minute}
%
\newdimen\trimheight
\newdimen\trimwidth
\newdimen\typeheight
\newdimen\typewidth
\newdimen\tempdimen
\newdimen\tablewidth
\newdimen\normaltextheight
\newbox\tempbox
\newdimen\tablewd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\foliofont{\fontsize{11}{13}\selectfont}
\def\bibfont{\fontsize{11}{13}\selectfont}
\def\rhfont{\footnotesize\itshape{}}
\def\catchlinefont{\footnotesize}
\def\titlefont{\fontsize{12}{15}\bfseries\boldmath\selectfont\centering{}}
\def\authorfont{\footnotesize}
\def\abstractfont{\footnotesize}
\def\addressfont{\small\itshape{}}
\let\affiliationfont\rhfont
\def\datefont{\footnotesize}
\def\sectionfont{\bfseries\boldmath\rightskip2pc plus1fill{}}
\def\subsectionfont{\bfseries\boldmath\itshape\rightskip2pc plus1fill{}}
\def\subsubsectionfont{\itshape\rightskip2pc plus1fill{}}
\def\paragraphfont{\bfseries\boldmath}
\def\subparagraphfont{\itshape}
\def\figurecaptionfont{\footnotesize}
\def\tablefont{\footnotesize}
\def\tablecaptionfont{\footnotesize}
%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xiipt{15}%
   \abovedisplayskip6\p@\@plus4\p@\@minus2\p@
   \abovedisplayshortskip\abovedisplayskip
   \belowdisplayshortskip\abovedisplayshortskip
   \belowdisplayskip\abovedisplayskip
   \let\@listi\@listI}
\normalsize
%
\newcommand\small{%
   \@setfontsize\small\@xipt{13}%
   \abovedisplayskip 8\p@ \@plus3\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus2\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}
%
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@xipt{13}%
   \abovedisplayskip 6\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 3\p@ \@plus\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 3\p@ \@plus\p@ \@minus\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip}
%
\newcommand\scriptsize{\@setfontsize\scriptsize\@viipt\@viiipt}
\newcommand\tiny{\@setfontsize\tiny\@vpt\@vipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}%
\newcommand\Large{\@setfontsize\Large\@xivpt{16}}%%
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
%
\DeclareOldFontCommand{\rm}{\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\scshape}{\@nomath\sc}
\DeclareOldFontCommand{\bi}{\bfseries\itshape}{\bfseries\itshape}
\renewcommand{\rmdefault}{cmr}%
%
\DeclareMathSizes{28}{28}{25}{20}
\DeclareMathSizes{24}{24}{21}{18}
\DeclareMathSizes{22}{22}{19}{17}
\DeclareMathSizes{20}{20}{17}{15}
\DeclareMathSizes{18}{18}{15}{10}
\DeclareMathSizes{12}{11}{8}{5}
\DeclareMathSizes{11}{11}{8}{5}
\DeclareMathSizes{10}{10}{7}{5}
\DeclareMathSizes{9.5}{9}{7}{5}
\DeclareMathSizes{9}{9}{7}{5}
\DeclareMathSizes{8}{8}{6}{5}
\DeclareMathSizes{7}{7}{5}{5}
\DeclareMathSizes{6}{6}{5}{5}
\DeclareMathSizes{5}{5}{5}{5}
\def\bm#1{\mathchoice
    {\mbox{\boldmath$#1$}}%
          {\mbox{\boldmath$#1$}}%
          {\mbox{\boldmath$\scriptstyle#1$}}%
          {\mbox{\boldmath$\scriptscriptstyle#1$}}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Initializations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\trimheight11.0truein
\trimwidth8.5truein
\typeheight9.2in
%
\setlength\parindent{0.25in}
\setlength\headheight{8\p@}%8
\setlength\headsep   {0.3in}%25
\setlength\topskip   {.7\p@}%{10\p@}
\setlength\footskip{0.17in}%
\setlength\maxdepth{4\p@}
\setlength\@maxdepth\maxdepth
\setlength\textwidth{6.6in}
\setlength\textheight\typeheight
\advance\textheight-\headheight
\advance\textheight-\headsep
\setlength\normaltextheight{\textheight}
\setlength\oddsidemargin{2cm}
\setlength\evensidemargin{2cm}
%\addtolength\evensidemargin{-\textwidth}
%\addtolength\evensidemargin{-\oddsidemargin}
\setlength\marginparwidth{.25in}
\setlength\marginparsep{7\p@}
\setlength\marginparpush{5\p@}
\setlength\topmargin{.75in}%
\setlength\footnotesep{6\p@}%6pt
\setlength{\skip\footins}{18\p@ \@plus 5\p@ \@minus 2\p@}
\setlength\floatsep    {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{18\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\setlength\partopsep{0\p@}
%
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}%
\renewcommand\baselinestretch{}%{}
\setlength\parskip{0\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
%
\newenvironment{contents}{%
    \cleardoublepage
    \vspace*{0.92in}
    {\centering{
    {\bf {\uppercase{Contents}}\par}\par}}
    \vskip 0.05in             %VSPACE BEFORE CONTENT TEXT
    \vskip -3pt
    \normalsize
}{\par\vfill\pagebreak\par}%
%
\newcommand\mastertoc{%
    \cleardoublepage
    \title{CONTENTS}
    \normalsize
    \@starttoc{toc}%
    \par\vfill\pagebreak\par
    \pagenumbering{arabic}
    \setcounter{page}{1}}
%
\newenvironment{titlepages}{%
    \thispagestyle{plain}
    \renewcommand\thepage{\roman{page}}
    \par\ \par
    \vfill
     \begin{flushleft}
     \large
     {\bfseries Proceedings Trim Size: 11in x 8.5in \\
     Text Area: 9.25in (include runningheads) x 6.6in \\
     Main Text is 10/13pt} \\[29pt]
     {\large For Half-Title Page (prepared by publisher)}
    \end{flushleft}
    \vfill
    \pagebreak\par
    \thispagestyle{plain}
    \ \par
    \vfill
    {\begin{center}\large Publishers' page --- (Blank page)\end{center} }\par
    \vfill
    \pagebreak\par
    \thispagestyle{plain}
    \ \par
    \vfill
    {\begin{center}\large For Full Title Page (prepared by publisher)
    \end{center} }\ \par
    \vfill
    \pagebreak
    \ \par    
    \vfill
    {\begin{center}\large For Copyright Page (prepared by publisher)
    \end{center} }\ \par
    \vfill
    \pagebreak\par}{\par}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\CNfont{\large\bfseries}
\def\CTfont{\Large\bfseries\boldmath{}}
\newcounter{part}
\renewcommand\thepart{\Alph{part}}
\newcommand\part{\cleardoublepage
                 \renewcommand\thepage{\arabic{page}}
                 \thispagestyle{pblind}%
                 \global\@topnum\z@
                 \@afterindentfalse
                 \secdef\@part\@spart}
%
\def\@part[#1]#2#3{%
 {\vspace*{3pc}
  \vbox{\parindent\z@ \reset@font
         \refstepcounter{part}
%    \addcontentsline{toc}{part}{#2}
    \addcontentsline{toc}{part}{Part~\thepart\quad #1}
         \vskip 1.4in
         {\centering{\CNfont PART \thepart}
         \par
         \vskip 12\p@
    \vbox{\CTfont #1}\par}
    \vskip 106\p@
    \vbox{\noindent #3}\par
    \nobreak}}
  \@endpart}
%
\def\@spart#1{%
    {\PTfont #1\par}%
    \@endpart}
%
\def\@endpart{\vfil\newpage
              \if@twoside
                \hbox{}%
                \thispagestyle{pblind}%
                \newpage
              \fi}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@clinebuf{}
%
\def\title#1{
    \setcounter{section}{0}
    \setcounter{subsection}{0}
    \setcounter{subsubsection}{0}
    \setcounter{footnote}{0}
    \setcounter{table}{0}
    \setcounter{figure}{0}
    \setcounter{equation}{0}
    \setcounter{appendixc}{0}
    \@addtoreset{theorem}{section}
    \renewcommand\theequation{\arabic{equation}}
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \long\def\@makefntext##1{\noindent
    \@textsuperscript{\normalfont\@thefnmark}{{##1}}}%
   %\cleardoublepage %to start a new paper in odd page
    \thispagestyle{plain}\@thanks
    \vspace*{-36pt}
    \@clinebuf
    \vspace*{30pt}
    {\centering{\titlefont #1\par}}%
    \vspace*{28pt}
}
\def\subtitle#1{\titlefont :~{{#1}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
\newcounter{ptitl}
\def\@empty{}

\def\wstoc#1#2{\refstepcounter{ptitl}
    \addtocontents{toc}
    {\protect\contentsline{chapter}
        \protect{\noindent#1}{\thepage}
    }
    \ifx#2\@empty\else
    \addcontentsline{toc}{authors}{\noindent#2}\fi
    }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Definitions %%%%%%%%%%%%%%%%%%%%%%
%
\def\cdash{$^{\raisebox{-0.5pt}{\hbox{--}}}$}
\def\mtb{\fboxsep=0pt\fbox{\phantom{e}}}
\def\email#1{\it #1\par}
\def\http#1{\it http://#1\par}

\def\breakhere{\par\vskip-\lastskip
            \pagebreak}
\def\forcebreak{\vadjust{\pagebreak}}
%

\newcommand{\cal}{\protect\pcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
%
\def\author#1{\par
    {\centering{\authorfont#1}\par\vspace*{0.05in}}
}

\def\address#1{\par
    {\centering{\affiliationfont#1\par}}\par\vspace*{13pt}
}
%
\newenvironment{abstract}{\par%
    \vspace*{6pt}
    \abstractfont
    \noindent\leftskip.3in\rightskip.3in
}{%
  \par}
%
\def\keywords#1{\par
    \vspace*{8pt}
    {\authorfont{\leftskip.3in\rightskip\leftskip
    \noindent{\it Keywords}\/:\ #1\par}}\vskip12pt}
%
%
\def\category#1#2{\par
    \vspace*{-8pt}
    {\authorfont{\leftskip.3in\rightskip\leftskip
    \noindent{\it Categories}\/:\ #1:\ #2\par}}\vspace*{18pt}\par}
%
\newif\if@cghi
\def\cite{\@cghitrue\@ifnextchar [{\@tempswatrue
    \@citex}{\@tempswafalse\@citex[]}}
\def\citelow{\@cghifalse\@ifnextchar [{\@tempswatrue
    \@citex}{\@tempswafalse\@citex[]}}
\def\@cite#1#2{{$^{#1}$\if@tempswa\typeout
    {WSPC warning: optional citation argument
    ignored: `#2'} \fi}}
\newcommand{\citeup}{\cite}
%
\newcount\@tempcntc
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi
 \@tempcnta\z@\@tempcntb\m@ne\def\@citea{}\@cite{\@for\@citeb:=#2\do
  {\@ifundefined
   {b@\@citeb}{\@citeo\@tempcntb\m@ne\@citea\def\@citea{,}%
     {\mbox{\bfseries ?}}\@warning
   {Citation `\@citeb' on page \thepage \space undefined}}%
  {\setbox\z@\hbox{\global\@tempcntc0\csname b@\@citeb\endcsname\relax}%
   \ifnum\@tempcntc=\z@ \@citeo\@tempcntb\m@ne
    \@citea\def\@citea{;\ }\hbox{\csname b@\@citeb\endcsname}%
   \else
    \advance\@tempcntb\@ne
    \ifnum\@tempcntb=\@tempcntc
    \else\advance\@tempcntb\m@ne\@citeo
    \@tempcnta\@tempcntc\@tempcntb\@tempcntc\fi\fi}}\@citeo}{#1}}
\def\@citeo{\ifnum\@tempcnta>\@tempcntb\else\@citea\def\@citea{,}%
 \ifnum\@tempcnta=\@tempcntb\the\@tempcnta\else
  {\advance\@tempcnta\@ne\ifnum\@tempcnta=\@tempcntb \else \def\@citea{\mbox{--\/}}\fi
   \advance\@tempcnta\m@ne\the\@tempcnta\@citea\the\@tempcntb}\fi\fi}
\def\@refcitex[#1]#2{\if@filesw\immediate\write\@auxout
    {\string\citation{#2}}\fi
\def\@citea{}\@refcite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{, }\@ifundefined
    {b@\@citeb}{{\bf ?}\@warning
    {Citation `\@citeb' on page \thepage \space undefined}}
    \hbox{\csname b@\@citeb\endcsname}}}{#1}}
 \def\@refcite#1#2{{#1\if@tempswa\typeout
        {WSPC warning: optional citation argument
    ignored: `#2'} \fi}}
 \def\refcite{\@ifnextchar[{\@tempswatrue
    \@refcitex}{\@tempswafalse\@refcitex[]}}
%
%%%%%%%%%%%%%%%%%%%%%%%% Table of Contents %%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@listI{\leftmargin\leftmargini
            \listparindent\parindent
            \parsep \z@\labelsep.5em
            \topsep 6.5\p@ \@plus3\p@% \@minus2\p@
            \itemsep0\p@}
%
\let\@listi\@listI
\@listi
%
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    3.5\p@ \@plus2\p@ \@minus\p@
              \parsep    \z@
              \itemsep\z@}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep 3.5\p@ \@plus\p@\@minus\p@
              \parsep \z@
             % \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep0\p@}%
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}
%
\setlength\leftmargini{3pc}
\setlength\leftmarginii{2.2em}
\setlength\leftmarginiii{1.87em}
\setlength\leftmarginiv{1.7em}
\setlength\leftmarginv{1em}
\setlength\leftmarginvi{1em}
\setlength\leftmargin{\leftmargini}
\setlength\listparindent{\parindent}
\setlength\labelsep{.5em}
\setlength\labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\alph{enumii}}
\renewcommand\theenumiii{\roman{enumiii}}
\renewcommand\theenumiv{\Alph{enumiv}}
\newcommand\labelenumi{(\theenumi)}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}
%
\def\enummax#1{\setbox\tempbox=\hbox{#1\hskip\labelsep}%
  \expandafter\global\csname leftmargin\romannumeral\the\@enumdepth\endcsname\wd\tempbox}
%
\def\enumerate{\@ifnextchar[{\@enumerate}{\@enumerate[\csname label\@enumctr\endcsname]}}
%
\def\@enumerate[#1]{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne\edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \enummax{#1}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\def\makelabel##1{\hss\llap{##1}}}\fi}
%
\let\Item\item
\newenvironment{enumeroman}{%
  \def\theenumi{\roman{enumi}}\def\theenumii{\alph{enumii}}%
  \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
        \let\item\Item
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{alphlist}{%
  \def\theenumi{\alph{enumi}}\def\theenumii{\alph{enumii}}%
  \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
        \let\item\Item
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{Alphlist}{%
  \def\theenumi{\Alph{enumi}}\def\theenumii{\Alph{enumii}}%
  \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
        \let\item\Item
  \begin{enumerate}%
}{%
  \end{enumerate}}
%
\newenvironment{arabiclist}{%
    \def\theenumi{\arabic{enumi}}\def\theenumii{\arabic{enumii}}
   \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}% %(2)
%    \def\labelenumi{\theenumi.}\def\labelenumii{(\theenumii)}%  %2. 
    \let\item\Item
    \begin{enumerate}
    }{%
    \end{enumerate}}
%
\newenvironment{romanlist}{%
    \def\theenumi{\roman{enumi}}\def\theenumii{\roman{enumii}}
    \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
    \let\item\Item
    \begin{enumerate}
    }{%
    \end{enumerate}}
%
\newenvironment{Romanlist}{%
    \def\theenumi{\Roman{enumi}}\def\theenumii{\Roman{enumii}}
    \def\labelenumi{(\theenumi)}\def\labelenumii{(\theenumii)}%
    \let\item\Item
    \begin{enumerate}
    }{%
    \end{enumerate}}
%
\newenvironment{itemlist}{%
    \def\labelenumi{\labelitemi}
    \let\item\Item
    \begin{enumerate}
    }{%
    \end{enumerate}}
%
\newenvironment{unnumlist}{%
    \let\item\Item
    \leftmargini2pc
    \ifnum \@enumdepth >3 \@toodeep\else
        \advance\@enumdepth \@ne
        \list{}{\itemindent-2pc\topsep6pt
              \def\makelabel##1{\hss\llap{##1}}}%
   \fi
}{%
  \endlist}
%
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{secnumdepth}{3}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thesection      {\arabic{section}}
\renewcommand\thesubsection   {\thesection.\arabic{subsection}}
\renewcommand\thesubsubsection{\thesubsection.\arabic{subsubsection}}
\renewcommand\theparagraph    {\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph {\theparagraph.\arabic{subparagraph}}
%
\newcommand\section{\@startsection{section}{1}{\z@}{-18pt \@plus-3pt minus-6pt}{5pt}{\sectionfont}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}{-18pt \@plus-3pt minus-6pt}{5pt}{\subsectionfont}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{-18pt \@plus-3pt minus-3pt}{5pt}{\subsubsectionfont}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}{12pt \@plus -3pt}{-.5em}{\paragraphfont}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}{3.25ex \@plus1ex}{-1em}{\subparagraphfont}}
%
\def\@seccntformat#1{{\upshape{\csname the#1\endcsname.}}\hskip .5em}
%
\def\body{
\setcounter{footnote}{0}
\def\thefootnote{\alph{footnote}}
\def\@makefnmark{{$^{\rm \@thefnmark}$}}
}
\def\bodymatter{\body}
%
\def\copyrightinfo#1{\vspace*{-12pt}\footnotetext[0]{\\[-6pt]#1}\body}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Floats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setcounter{topnumber}{3}
\def\topfraction{.98}
\setcounter{bottomnumber}{3}
\def\bottomfraction{.98}
\setcounter{totalnumber}{6}%4
\def\textfraction{.01}
\def\floatpagefraction{.98100}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.981}
\def\dblfloatpagefraction{.99}
%
%
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{6\p@}% plus3pt minus3pt}
\setlength\belowcaptionskip{-1\p@}
%
\def\FigName{figure}%
%
\newbox\captionbox
%
\long\def\@makecaption#1#2{%
  \ifx\FigName\@captype
        \vskip\abovecaptionskip
    \setbox\captionbox\hbox{%
    {\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2\par}}
    \ifdim\wd\captionbox>\hsize
    {\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2\par}\par
    \else
    \centerline{\figurecaptionfont#1.\unskip\ignorespaces\hskip1em#2}\par
    \fi\else\global\setbox\captionbox\hbox{%
    {\tablecaptionfont#1.\ignorespaces\hskip1em#2\par}}
    \fi
}
%
\def\figurename{Fig.}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{htbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}  %Figure 1
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
%
\def\ArtDir{\jobname/}
%
\newcounter{table}
\def\tablename{\upshape Table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{htbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}  %Table 1.
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
%
\long\def\tbl#1#2{%
    \parindent\z@\ignorespaces\noindent\tablecaptionfont
    \caption{#1}%
    \par\setbox\tempbox\hbox{\tablefont #2}%
    \tablewd\hsize\advance\tablewd-\wd\tempbox\global\divide\tablewd\tw@
    \ifdim\wd\captionbox<\wd\tempbox\centerline{\unhbox\captionbox}
    \else\leftskip\tablewd\rightskip\leftskip{\unhbox\captionbox}\par
    \fi\vskip5pt\centerline{\box\tempbox}
}%
\long\def\ttbl#1#2#3{%
    \parindent\z@
    \noindent
    \tablecaptionfont
    \caption{#2}%
    \tablewd\hsize
    \advance\tablewd-#1
    \global\divide\tablewd\tw@
    \ifdim \wd\captionbox<#1
    \centerline{\unhbox\captionbox}
    \else
    {\leftskip\tablewd% plus1fill{}
    \rightskip\tablewd% plus1fill{}
    \unhbox\captionbox\par}
    \fi
    \vskip-8pt
    \centerline{#3}
}%
%
\newenvironment{InTextTable}{%
    \par\addvspace{16pt plus6pt minus6pt}%
    \begin{center}
    \tablecaptionfont\noindent\ignorespaces
}{%
    \end{center}
    \par\addvspace{17pt plus6pt minus6pt}}
%
\newenvironment{tabnote}{%
    \par\vskip5pt\footnotesize%
    \@ifnextchar[{\@tabnote}{\@tabnote}}{%
    \par\vskip-5pt}
\def\@tabnote{%
    \leftskip\tablewd\rightskip\leftskip}
%
\def\Hline{%
  \noalign{\ifnum0=`}\fi\hrule \@height 2\arrayrulewidth \futurelet
   \@tempa\@xhline}

\def\x{@{\extracolsep{\fill}}}
\def\toprule{\\[-6pt]\hline\\[-5.5pt]}
\def\colrule{\\[-7.5pt]\hline\\[-5.5pt]}
\def\botrule{\\[-7pt]\hline\\[-8.5pt]}
\def\crule#1{\\[-6pt]\cline{#1}\\[-6.5pt]}
%
\newenvironment{tablehere}
  {\parskip0pt\topskip 13pt plus 2pt minus 1pt\def\@captype{table}\begin{small}\begin{center}}
  {\end{center}\end{small}\vskip 13pt plus 2pt minus 1pt}
  
\newenvironment{figurehere}
  {\parskip0pt\topskip 13pt plus 2pt minus 1pt\def\@captype{figure}}
  {\vskip 13pt plus 2pt minus 1pt}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Maths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\setlength\arraycolsep{1.5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}

\renewcommand\theequation{\arabic{equation}}

%
\setlength\columnsep{0.3in}
\setlength\columnseprule{0\p@}
%
\def\bstrut{\vrule width0pt depth5pt}
\def\tstrut{\vrule width0pt height10pt}
\jot=3pt
%
\renewenvironment{equation}{\@equation}{\end@equation}
%
\def\refstepeqcnt{\refstepcounter{equation}}
\def\eqnnum{\displaywidth\linewidth\hbox{\@eqnnum}}
%
\def\@equation{\let\topsep\abovedisplayskip\@beginparpenalty\predisplaypenalty\@endparpenalty\postdisplaypenalty
     \refstepeqcnt
     \trivlist \item[]\leavevmode
       \hbox to\linewidth\bgroup\hfill $\m@th% $
       \displaystyle}%
\def\end@equation{$\hfill\eqnnum\egroup\endtrivlist}
%
\renewenvironment{eqnarray}{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip \tw@\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr}
        {%
      \@@eqncr
      \egroup
      \global\advance\c@equation\m@ne
   $$\global\@ignoretrue
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Footnote %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
     \@makefnmark#1}
%
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
    \mathchar "278\or \mathchar "27B\or \|\or **\or \dagger\dagger
    \or \ddagger\ddagger \else\@ctrerr\fi}}
%
\renewcommand\thefootnote{\@fnsymbol\c@footnote}
\def\bodymatter{\renewcommand\thefootnote{\alph{footnote}}
\def\@makefnmark{{$^{\rm \@thefnmark}$}}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Page Style %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\mark{{}{}}
%
\def\ps@plain{%
     \let\@mkboth\@gobbletwo
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
    \def\@evenhead{\rhfont{\upshape{\thepage}}\hfill}%
    \def\@oddhead{\hfill\rhfont{\upshape{\thepage}}}}%
%
\def\RightMark{}
%
\def\ps@headings{%
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
    \def\@evenhead{\rhfont{\upshape{\thepage}}\hskip1em\leftmark\hfill}%
    \def\@oddhead{\hfill\rhfont\rightmark\hskip1em{\upshape{\thepage}}}%
    \let\@mkboth\markboth
    \def\authormark##1{\markboth{##1}{\RightMark}}%
    \def\chaptermark##1{\gdef\RightMark{##1}\markright{##1}}}%
    %\def\sectionmark##1{\gdef\RightMark{##1}\markright{##1}}}%
%
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\rlap{\foliofont \thepage}\hfill{\rhfont\leftmark}\hfill}%
    \def\@oddhead{\hfill{\rhfont\rightmark}\hfill{\foliofont\thepage}}%
    \let\@mkboth\@gobbletwo
    \let\chaptermark\@gobble
    \let\sectionmark\@gobble}
%
\def\ps@pblind{%
    \let\@mkboth\@gobbletwo
    \def\@evenhead{}%
    \def\@oddhead{}
    \def\@oddfoot{\hfill\mbox{\lower8pc\hbox{}}\hfill}
   \let\@evenfoot\@oddfoot}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Index Shortcuts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\aindx#1{\index{author}{#1}}    %For Author index
\def\sindx#1{\index{subject}{#1}}   %For Subject index
\def\gindx#1{\index{general}{#1}}   %For General index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CITATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand{\eref}[1]{Eq.~(\ref{#1})}
\newcommand{\sref}[1]{Sec.~\ref{#1}}
\newcommand{\cref}[1]{Chap.~\ref{#1}}
\newcommand{\fref}[1]{Fig.~\ref{#1}}
\newcommand{\tref}[1]{Table~\ref{#1}}
\newcommand{\Eref}[1]{Equation (\ref{#1})}
\newcommand{\Sref}[1]{Section~\ref{#1}}
\newcommand{\Cref}[1]{Chapter~\ref{#1}}
\newcommand{\Fref}[1]{Figure~\ref{#1}}
\newcommand{\Tref}[1]{Table~\ref{#1}}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bibliography %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newdimen\bibindent
\bibindent=2em
\def\@biblabel#1{#1.}
\newenvironment{thebibliography}[1]
     {\section*{\bibname}\bibfont%
      \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
        \labelsep6pt%.33em
            \advance\leftmargin\labelsep
        \itemsep0pt\parsep0pt
            \usecounter{enumiv}%
            \renewcommand\theenumiv{\arabic{enumiv}}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcounter{appendixc}
\newcounter{subappendixc}[appendixc]
\newcounter{subsubappendixc}[subappendixc]
\renewcommand{\theappendixc}{\Alph{appendixc}}
\renewcommand{\thesubappendixc}{\Alph{appendixc}.\arabic{subappendixc}}
\renewcommand{\thesubsubappendixc}
   {\Alph{appendixc}.\arabic{subappendixc}.\arabic{subsubappendixc}}

\newcommand{\appendix}[1]{
       \vspace{14pt}
       \refstepcounter{appendixc}
       \setcounter{figure}{0}
       \setcounter{table}{0}
       \setcounter{lemma}{0}
       \setcounter{theorem}{0}
       \setcounter{corollary}{0}
       \setcounter{definition}{0}
       \setcounter{equation}{0}
       \renewcommand{\thefigure}{\Alph{appendixc}\arabic{figure}}
       \renewcommand{\thetable}{\Alph{appendixc}\arabic{table}}
       \renewcommand{\theappendixc}{\Alph{appendixc}}
       \renewcommand{\thelemma}{\Alph{appendixc}.\arabic{lemma}}
       \renewcommand{\thetheorem}{\Alph{appendixc}.\arabic{theorem}}
       \renewcommand{\thedefinition}{\Alph{appendixc}.\arabic{definition}}
       \renewcommand{\thecorollary}{\Alph{appendixc}.\arabic{corollary}}
       \renewcommand{\theequation}{\Alph{appendixc}.\arabic{equation}}
       \noindent{\sectionfont\appendixname~\theappendixc.{~\boldmath~#1}}
       \par\vspace{5pt}\@afterheading}

\newcommand{\subappendix}[1]{\vspace{14pt}
       \refstepcounter{subappendixc}
       \noindent{\sectionfont\appendixname~\thesubappendixc.~{~\bfseries\boldmath\itshape #1}}
       \par\vspace{8pt}\@afterheading}

\newcommand{\subsubappendix}[1]{\vspace{14pt}
       \refstepcounter{subsubappendixc}
       \noindent{\appendixname~\thesubsubappendixc.~{~\subsubsectionfont #1}}
       \par\vspace{8pt}\@afterheading}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\bib{B\kern-.05em{I}\kern-.025em{B}\kern-.08em}
\def\btex{B\kern-.05em{I}\kern-.025em{B}\kern-.08em\TeX}
%
\def\enddocument{%
   \label{eof}
   \@enddocumenthook
   \@checkend{document}%
   \clearpage
   \begingroup
     \if@filesw
       \immediate\closeout\@mainaux
       \let\@setckpt\@gobbletwo
       \let\@newl@bel\@testdef
       \@tempswafalse
       \makeatletter \input\jobname.aux
     \fi
     \@dofilelist
     \ifdim \font@submax >\fontsubfuzz\relax
       \@font@warning{Size substitutions with differences\MessageBreak
                  up to \font@submax\space have occured.\@gobbletwo}%
     \fi
     \@defaultsubs
     \@refundefined
     \if@filesw
       \ifx \@multiplelabels \relax
         \if@tempswa
           \@latex@warning@no@line{Label(s) may have changed.
               Rerun to get cross-references right}%
         \fi
       \else
         \@multiplelabels
       \fi
     \fi
   \endgroup
   \deadcycles\z@\@@end}
%
%%%%%%%%%%%%%%%%%%%%%%%%%% INDEX FOR ROOT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newif\if@restonecol
\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
    \columnseprule \z@
    \columnsep 35pt
    \thispagestyle{plain}
    \twocolumn[\vspace*{0.15in}
          \title{Index}
          \maketitle
          ]
    \bibfont
    \@mkboth{Index}{Index}
    \parindent\z@
    \parskip\z@
    \let\item\@idxitem}
\def\@idxitem{\par\hangindent 12pt}
\def\subitem{\par\hangindent 24pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par \vskip 15pt plus 5pt minus 3pt\relax}
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Trim Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\wsdraftnote{\hbox to \trimwidth{\scriptsize\sf\ \today
\quad\currenttime\hfil WSPC - Proceedings Trim Size: 11in x 8.5in
\quad\jobname\hfil}}%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\cropmarks{\nointerlineskip%
     \vbox to 0pt{\vskip-8.6pt
     \parindent0pt\hspace*{6pt}\wsdraftnote
     \par\advance\trimheight by 2pc
     \vbox to \trimheight{\offinterlineskip\parindent0pt\vskip-2pc
     \hbox to \trimwidth{\raise2.6pc\hbox{\vrule height2pc width0.3pt}\raisebox{2pc} %lefttop vert
   {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill\raisebox{2pc} %lefttop hori
   {\vrule width2pc height0.3pt}\hspace*{-2.6pc}\raise2.6pc\hbox{\vrule height2pc width0.3pt}}\vfill %righttop hori & vert
     \hbox to \trimwidth{\lower1.2pc\hbox{\vrule height2pc width0.3pt}\raisebox{1.4pc} %leftbot vert
        {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill  %leftbot hor
    \lower1.2pc\hbox{\vrule height2pc width0.3pt}\rlap{\hspace*{0.6pc}\raisebox{1.4pc}{\vrule width2pc height0.3pt}}}}} %rightbot vert & hori
     \vbox to 0pt{\offinterlineskip\parindent0pt\vskip\topmargin%\vskip1.3pt
     \tempdimen\normaltextheight
     \advance\tempdimen27pt
     \moveright\@themargin
     \vbox to \tempdimen{\hrule height\wsdraftrule 
     \vbox to 0pt{\vskip\headheight\vskip\headsep%\vskip3pt
                \vrule height\wsdraftrule width\textwidth}% %2nd hrule
      \hbox to \textwidth{\vrule width\wsdraftrule \hfill\vrule height9.2in width\wsdraftrule}%
           \hrule height\wsdraftrule}}}%%
\def\xcropmarks{\nointerlineskip%
     \vbox to 0pt{\vskip-8.6pt
    \parindent0pt\hspace*{6pt}\wsdraftnote
     \par\advance\trimheight by 2pc
     \vbox to \trimheight{\offinterlineskip\parindent0pt\vskip-2pc
     \hbox to \trimwidth{\raise2.6pc\hbox{\vrule height2pc width0.3pt}\raisebox{2pc}
   {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill\raisebox{2pc}
   {\vrule width2pc height0.3pt}\hspace*{-2.6pc}\raise2.6pc\hbox{\vrule height2pc width0.3pt}}\vfill
     \hbox to \trimwidth{\lower1.2pc\hbox{\vrule height2pc width0.3pt}\raisebox{1.4pc}
        {\hspace*{-2.6pc}\vrule width2pc height0.3pt}\hfill 
    \lower1.2pc\hbox{\vrule height2pc width0.3pt}\rlap{\hspace*{0.6pc}\raisebox{1.4pc}{\vrule width2pc height0.3pt}}}}}
     \vbox to 0pt{\offinterlineskip\parindent0pt\vskip\topmargin%\vskip1.3pt
     \tempdimen\normaltextheight
     \advance\tempdimen27pt
     \moveright\@themargin
     \vbox to \tempdimen{\hrule height\wsdraftrule 
     \vbox to 0pt{\vskip\headheight\vskip\headsep%\vskip3pt
                \vrule height\wsdraftrule width\textwidth}%
      \hbox to \textwidth{\vrule width\wsdraftrule \hfill\vrule height\tempdimen width\wsdraftrule}%
           \hrule height\wsdraftrule}}}%%
%
\let\trimmarks\cropmarks
\def\trimmarks{\relax}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output Routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\@outputpage{%
\begingroup
  \let \protect \noexpand
  \@resetactivechars
  \let\-\@dischyph
  \let\'\@acci\let\`\@accii\let\=\@acciii
  \let\\\@normalcr
  \let\par\@@par
  \shipout \vbox{%
    \set@typeset@protect
    \aftergroup \endgroup
    \aftergroup \set@typeset@protect
                                % correct? or just restore by ending
                                % the group?
  \if@specialpage
    \global\@specialpagefalse\@nameuse{ps@\@specialstyle}%
  \fi
  \if@twoside
    \ifodd\count\z@ \let\@thehead\@oddhead \let\@thefoot\@oddfoot
         \let\@themargin\oddsidemargin
    \else \let\@thehead\@evenhead
       \let\@thefoot\@evenfoot \let\@themargin\evensidemargin
    \fi
  \fi
  \reset@font
  \normalsize
  \baselineskip\z@skip \lineskip\z@skip \lineskiplimit\z@
    \@begindvi\trimmarks
    \vskip \topmargin
    \moveright\@themargin \vbox {%
      \setbox\@tempboxa \vbox to\headheight{%
        \vfil
        \color@hbox
          \normalcolor
          \hb@xt@\textwidth {%
            \let \label \@gobble
            \let \index \@gobble
            \let \glossary \@gobble
            \@thehead
            }%
        \color@endbox
        }%
      \dp\@tempboxa \z@
      \box\@tempboxa
      \vskip \headsep
      \box\@outputbox
      \baselineskip \footskip
      \color@hbox
      \normalcolor
        \hb@xt@\textwidth{%
          \let \label \@gobble
          \let \index \@gobble    
          \let \glossary \@gobble 
          \@thefoot
          }%
      \color@endbox
      }%
    }%
\global \@colht \textheight
\stepcounter{page}%
\let\firstmark\botmark
}
%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process Options %%%%%%%%%%%%%%%%%%%%%%%%
%
\pagestyle{headings}
\pagenumbering{arabic}
\onecolumn
\frenchspacing
\flushbottom
\def\sloppy{\tolerance=100000\hfuzz=\maxdimen \vfuzz=\maxdimen}
\sloppy
\vbadness=12000
\hbadness=12000
\RequirePackage[figuresright]{ws-rotating}
\RequirePackage{url}
\RequirePackage{graphicx}
\RequirePackage{chapterbib}
\if@wssquare
\RequirePackage[space]{cite}
\def\@biblabel#1{[#1]}
\else
\RequirePackage[super,nospace]{cite}
\def\@biblabel#1{#1.}
\fi
\let\refcite\citen
\def\thepage{}
\endinput
%
%% End of file
%% For Support: ykoh@wspc.com.sg
%%              rajesh@wspc.com.sg
%% ws-procs11x85.cls (C) World Scientific Publishing Co.

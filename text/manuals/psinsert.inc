% For use with ArborText DVILASER/PS.
% insertplot parameters:
%           #1: PostScript file name (Example: PGPLOT /VPS file)
%           #2: figure height (in)
%           #3: figure width (in)
%           #4: lower left corner x (inches from llc of paper)
%           #5: lower left corner y (inches from llc of paper)
%           #6: desired magnification
%           #7: 1 => original in landscape mode, 0 => portrait
% Parameters 2-5 refer to the dimensions of the original figure
% when printed directly.
%
\def\PSinsert#1#2#3#4#5#6#7{
\vskip 10pt\nobreak\hbox to \hsize{\hss\dimen0=#3in\hbox to #6\dimen0{%
\dimen0=#2in\vbox to #6\dimen0{\vss
\special{ps::[asis, begin]
  0 SPB
  /PSavenv save def
  /showpage {} def 
  /initgraphics {} def 
  Xpos Ypos translate
  #6 dup scale 
  #7 1 eq { -90 rotate -8.15 72 mul -0.25 72 mul translate} if
  #4 72 mul neg #5 72 mul neg translate
}
\special{ps: plotfile #1 asis}
\special{ps::[asis,end]
  PSavenv restore
  0 SPE
}
}\hss}\hss}\vskip 10pt}

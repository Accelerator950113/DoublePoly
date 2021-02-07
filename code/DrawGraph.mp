verbatimtex
%&latex
\documentclass{article}
\usepackage{amsfonts}
\usepackage{bm}
\usepackage{amsmath}
\begin{document}
etex

prologues:=3;
filenametemplate "%j-%c.eps";
u=1mm;

def dp(expr x, s, c) =
	fill fullcircle scaled s shifted x withcolor c;
	draw fullcircle scaled s shifted x withcolor black;
enddef;

vardef COL(expr r, g, b) = ((r/255)*red + (g/255)*green + (b/255)*blue) enddef;

beginfig(1);

color nodeColor, curveColor;
nodeColor := green;
curveColor := red;

numeric Len;

Len := scantokens(readfrom "DrawSettings.txt")*u;
closefrom "nodeInfo.txt";

numeric nr, l, maxc;
nr := 4u;
maxc := 0;

numeric n;
numeric p[][];
numeric m, pn, ms;
numeric prex, prey, nowx, nowy;

n := scantokens(readfrom "nodeInfo.txt");
for i = 1 upto n :
	p[i][1] := scantokens(readfrom "nodeInfo.txt");
	if p[i][1] > maxc : maxc := p[i][1]; fi;
	p[i][2] := scantokens(readfrom "nodeInfo.txt");
	if p[i][2] > maxc : maxc := p[i][2]; fi;
endfor;
closefrom "nodeInfo.txt";

numeric nc[];
for i = 1 upto n :
	nc[i] = scantokens(readfrom "colorInfo.txt");
endfor;
closefrom "colorInfo.txt";

string lab[];
for i = 1 upto n :
	lab[i] = (readfrom "labelInfo.txt");
endfor;
closefrom "labelInfo.txt";

l := Len / maxc;

pickup pencircle scaled 5;

ms := scantokens(readfrom "SublineInfo.txt");
for i = 1 upto ms :
	pn := scantokens(readfrom "SublineInfo.txt");
	prex := scantokens(readfrom "SublineInfo.txt");
	prey := scantokens(readfrom "SublineInfo.txt");
	for j := 2 upto pn :
		nowx := scantokens(readfrom "SublineInfo.txt");
		nowy := scantokens(readfrom "SublineInfo.txt");
		draw (prex*l, prey*l)--(nowx*l, nowy*l) withcolor (0.35blue + 0.65white);
		prex := nowx;
		prey := nowy;
	endfor;
endfor;
closefrom "SublineInfo.txt";

pickup pencircle scaled 1;

m := scantokens(readfrom "edgeInfo.txt");
for i = 1 upto m :
	pn := scantokens(readfrom "edgeInfo.txt");
	prex := scantokens(readfrom "edgeInfo.txt");
	prey := scantokens(readfrom "edgeInfo.txt");
	for j := 2 upto pn :
		nowx := scantokens(readfrom "edgeInfo.txt");
		nowy := scantokens(readfrom "edgeInfo.txt");
		draw (prex*l, prey*l)--(nowx*l, nowy*l);
		prex := nowx;
		prey := nowy;
	endfor;
endfor;
closefrom "edgeInfo.txt";

pickup pencircle scaled 1;

pair zp[];
for i = 1 upto n :
	zp[i] := (p[i][1]*l, p[i][2]*l);
endfor;

for i = 1 upto n :
	dp(zp[i], nr, nc[i]*nodeColor+(1-nc[i])*white);
	%dp(zp[i], nr, nodeColor);
	label(lab[i] infont defaultfont scaled sqrt(1/length(lab[i])), zp[i]);
endfor;

endfig;

end;

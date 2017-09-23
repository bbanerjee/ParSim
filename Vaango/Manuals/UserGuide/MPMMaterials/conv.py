#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
#
# Mike Meylan is hacking this to include more latex commands
# and Malte has added a few things that get rid of his personal latex commands.
#
# This code has been modified by Anthony Miller for handling of inline 
# mathematics and more sophisticated documents.
#
# Original idea from : 
#       Maxime Biais <maxime@biais.org>
#     but has been nearly all rewritten since...
# A good fraction of this code was written by
# Marc Poulhiès <marc.poulhies@epfl.ch>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# $Id: latex2twiki.py,v 1.2 2005/07/27 12:40:53 poulhies Exp $

import sys, re

bullet_level=0
enum_level = 0
#bdoc = None
bdoc = 1
end_line = 1
count = 0;
countfig = 0;
label_list = [];
labelfig_list = [];
lines = [];

verbatim_mode = 0
math_mode = 0
eqnarry_mode = 0
label_mode = 0

def dummy():
    pass

def toggle_label():
	global label_mode
	label_mode += 1

def inc_bullet():
	global bullet_level
	bullet_level += 1

def dec_bullet():
	global bullet_level
	bullet_level -= 1

def inc_enum():
	global enum_level
	enum_level += 1
	
def dec_enum():
	global enum_level
	enum_level -= 1
	
def start_doc():
	global bdoc;
	bdoc = 1

def do_not_el():
	global end_line
	end_line = None

def do_el():
	global end_line;
	end_line = 1

def decide_el():
	global end_line
	if bullet_level == 0:
		return "\n"
	else:
		return " "

def decide_math_replace():
	global math_mode
	if math_mode == 1:
		return r"\1"
	else:
		return " "

def decide_math():
	global math_mode
	if math_mode == 1:
		return "$$"
	else:
		return "$$"
		
def start_verbatim():
	global verbatim_mode
	verbatim_mode = 1

def end_verbatim():
	global verbatim_mode
	verbatim_mode = 0

def start_eqnarry():
	global eqnarry_mode
	eqarry_mode = 1

def end_eqnarry():
	global eqnarry_mode
	eqnarry_mode = 0

def toggle_math():
	global math_mode
	math_mode = 1 - math_mode

conv_table = { '>':'&gt;',
			   '<':'&lt;'}

def translate_to_html(char):
	global verbatim_mode
	global conv_table
	if verbatim_mode == 0:
		return conv_table[char]
	else:
		return char

NONE = "__@NONE@__"

tr_list2 = [
	(r"\\lstset", None, dummy),
	(r"\\begin\{lstlisting}", (lambda : r"{% highlight xml %}"), dummy),
	(r"\\end\{lstlisting}", (lambda : r"{% endhighlight %}"), dummy),
	(r"\\footnotesize", None, dummy),
	(r"\\begin\{abstract}", None, dummy),
	(r"\\begin\{article}", None, dummy),
	(r"\\end\{abstract}", None, dummy),
	(r"\\end\{article}", None, dummy),
	(r"\\end\{document}", None, dummy),
	(r"\\protect", None, dummy),
	(r"\\small", None, dummy),
	(r"\\func", None, dummy),
	(r"\\begin\{document}", None, start_doc),
	(r"\\cite{(.*?)}", (lambda :r"[[\1]]"), dummy),
#	(r"\\label{(.*?)}", (lambda :r" (\1)"), dummy),
#	(r"\\ref{(.*?)}", (lambda :r"(\1)"), dummy),
	(r"\\citep{(.*?)}", (lambda :r"[[\1]]"), dummy),
	(r"\\citet{(.*?)}", (lambda :r"[[\1]]"), dummy),
	(r"\\emph{(.*?)}", (lambda :r"''\1'' "), dummy),
	(r"\\textit{(.*?)}", (lambda :r"''\1'' "), dummy),
	(r"\\texttt{(.*?)}", (lambda : r"=\1= "), dummy),
#	(r"\\text{(.*?)}", (lambda : r"=\1= "), dummy),
	(r"\\textbf{(.*?)}", (lambda : r"'''\1''' "), dummy),
	(r"\\begin{verbatim}", (lambda : "<br><code>"), start_verbatim),
	(r"\\end{verbatim}", (lambda : "</code><br>"), end_verbatim),
	(r"\\begin{itemize}", (lambda : "\n"), inc_bullet),
	(r"\\end{itemize}", None, dec_bullet),
	(r"\\begin{enumerate}", (lambda : "\n"), inc_enum),
	(r"\\end{enumerate}", None, dec_enum),
	(r"\\item (.*?)", (lambda :  (r"#" * enum_level) +(r"\n*" * bullet_level) + r"\1"), dummy),
	(r"\\begin{equation[*]*}", (lambda :"<div>\n$$"), toggle_math),
	(r"\\end{equation[*]*}", (lambda :"$$\n</div>"), toggle_math),
	(r"\\\[", (lambda :"<div>\n$$"), toggle_math),
	(r"\\dfrac", (lambda :r"\\frac"), dummy),
	(r"\\\]", (lambda :"$$\n</div>"), toggle_math),
	(r"\\begin{eqnarray[*]?}", (lambda :r"<div>\n$$\\begin{matrix}"), toggle_math),
	(r"\\begin{array[*]?}", (lambda :r"\\begin{matrix}"), toggle_math),
	(r"\\end{eqnarray[*]?}", (lambda :r"\\end{matrix}$$\n</div>"), toggle_math),
	(r"\\end{array[*]?}", (lambda :r"\\end{matrix}"), toggle_math),
#	(r"(\\begin{.*?})", decide_math_replace, dummy),
#	(r"(\\end{.*?})",decide_math_replace, dummy),
#	(r"~\\ref{([^}]*)}",(lambda : r" ---\1---"),dummy),
	(r"``(.*?)''", (lambda :r'"\1"'), dummy),
	(r"\\paragraph{(.*?)}", (lambda : r"##### \1"), dummy),
	(r"\\subsubsection{(.*?)}", (lambda : r"#### \1"), dummy),
	(r"\\subsection{(.*?)}", (lambda : r"### \1"), dummy),
	(r"\\section{(.*?)}", (lambda : r"## \1"), dummy),
	(r"\\_", (lambda :"_"), dummy),
	(r"\\title{(.*)}", (lambda :r"# \1 "),dummy),
        (r"\\author{(.*)}", (lambda :r"\1"),dummy),
        (r"\\date{(.*)}", (lambda :r"\1"),dummy),
	(r"\\tableofcontents",None, dummy),
	(r"\\null",None, dummy),
	(r"\\newpage",None, dummy),
	(r"\\thispagestyle{.*?}", None, dummy),
	(r"\\maketitle", None, dummy),
	(r"\n$", decide_el, dummy),
#	(r"[^\\]?\{", None, dummy),
#	(r"[^\\]?\}", None, dummy),
#	(r"\$(.*?)\$",(lambda :r"$$\1$$"),dummy),
	(r"\$",decide_math,toggle_math),
#	(r"%.*$",None, dummy),
	(r"\\r{(.*?)}", (lambda : r"\\mathrm{\1}"), dummy),
	(r"\\d ", (lambda : r"\\,\mathrm{d} "), dummy),
	(r"\\i ", (lambda : r"\\mathrm{i} "), dummy),
	(r"\\i\\", (lambda : r"\\mathrm{i}\\"), dummy),
	(r"\\e\^", (lambda : r"\\mathrm{e}^"), dummy),
	(r"\\begin{align[*]?}", (lambda :r"<div>\n$$\\begin{align}"), toggle_math),
	(r"\\end{align[*]?}", (lambda :r"\\end{align}$$\n</div>"), toggle_math),
	(r"\\begin{aligned[*]?}", (lambda :r"\\begin{align}"), toggle_math),
	(r"\\end{aligned[*]?}", (lambda :r"\\end{align}"), toggle_math),
	(r"\\begin{subequations[*]?}", None, dummy),
	(r"\\end{subequations[*]?}", None, dummy),
	(r"\\nonumber", None, dummy),
	(r"\\displaystyle", None, dummy),
      (r"\\hdots", None, dummy),
#      (r"\{c\}", None, dummy),
      (r"\{cc\}", None, dummy),
      (r"\{ccc\}", None, dummy),
      (r"\{cccc\}", None, dummy),
      (r"\{ccccc\}", None, dummy),
      (r"\{cccccc\}", None, dummy),
      (r"\{ccccccc\}", None, dummy),
      (r"\{cccccccc\}", None, dummy),
	(r"\\dint", (lambda : r"\\int"), dummy),
	(r"\\hspace{.*?}", None, dummy),        
	(r"\\vspace{.*?}", None, dummy),        
#	(r"\\bf", None, dummy),
	(r"\\Large", None, dummy),
	(r"\\boxed", None, dummy),
	(r"\\em", None, dummy),
	(r"\\newpage", None, dummy),
	(r"\\ast", (lambda : r"*"), dummy),
	(r"\\tfrac", (lambda : r"\\frac"), dummy),
	(r"\\leqslant", (lambda : r"\\leq"), dummy),
	(r"\\geqslant", (lambda : r"\\geq"), dummy),
	(r"\\notag", None, dummy),

   (r"\\Blue", None, dummy),
   (r"\\begin{thebibliography}.*$", None, dummy),
   (r"\\end{thebibliography}", None, dummy),
   (r"\\bibitem.*$", None, dummy),
   (r"\\newblock {(.*?)}", (lambda : r"''\1''"), dummy),
   (r"\\newblock", None, dummy),
	(r"{\\it(.*?)}", (lambda : r"''\1''"), dummy),
	(r"{\\bf(.*?)}", (lambda : r"'''\1'''"), dummy),
	(r"\\begin{figure}", (lambda :"{|align=\"center\" |[[image:figure.png|thumb|400px|Figure title]]"), dummy),
	(r"\\end{figure}", (lambda :"|}"), dummy),
	(r"\\begin{center}", (lambda :"{|align=\"center\" |[[image:figure.png|thumb|400px|Figure title]]"), dummy),
	(r"\\end{center}", (lambda :"|}"), dummy),

   (r"\\Dela", (lambda : r"\\Delta a"), dummy),
   (r"\\Delb", (lambda : r"\\Delta b"), dummy),
   (r"\\DA", (lambda : r"\\text{dA}"), dummy),
   (r"\\DAvec", (lambda : r"\\text{d}\mathbf{A}"), dummy),
   (r"\\Da", (lambda : r"\\text{da}"), dummy),
   (r"\\Davec", (lambda : r"\\text{d}\mathbf{a}"), dummy),
   (r"\\DV", (lambda : r"\\text{dV}"), dummy),
   (r"\\DSq", (lambda : r"\\text{d}\\square"), dummy),
   (r"\\Dtau", (lambda : r"\\text{d}\\tau"), dummy),
   (r"\\domega", (lambda : r"\\text{d}\\omega"), dummy),
   (r"\\dOmega", (lambda : r"\\text{d}\\Omega"), dummy),
   (r"\\dGamma", (lambda : r"\\text{d}\\Gamma"), dummy),
   (r"\\dzeta", (lambda : r"\\text{d}\\zeta"), dummy),
   (r"\\Ds", (lambda : r"\\text{d}s"), dummy),
   (r"\\Dt", (lambda : r"\\text{d}t"), dummy),
   (r"\\Dx", (lambda : r"\\text{d}\\mathbf{x}"), dummy),
   (r"\\dr", (lambda : r"\\text{d}r"), dummy),
   (r"\\dx", (lambda : r"\\text{d}x"), dummy),
   (r"\\dy", (lambda : r"\\text{d}y"), dummy),
   (r"\\dz", (lambda : r"\\text{d}z"), dummy),
   (r"\\dk", (lambda : r"\\text{d}k"), dummy),
   (r"\\dBr", (lambda : r"\\text{d}\\mathbf{r}"), dummy),
   (r"\\dBx", (lambda : r"\\text{d}\\mathbf{x}"), dummy),
   (r"\\dBk", (lambda : r"\\text{d}\\mathbf{k}"), dummy),
   (r"\\BCe", (lambda : r"\\mathcal{E}"), dummy),
   (r"\\Conj{(.*?)}",    (lambda : r"\\overline{\1}"), dummy),
   (r"\\Gradbar{(.*?)}", (lambda : r"\\overline{\\boldsymbol{\\nabla}} \1"), dummy),
   (r"\\Grady{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}_y \1"), dummy),
   (r"\\Grads{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}_s \1"), dummy),
   (r"\\Gradp{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}' \1"), dummy),
   (r"\\GradX{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}_{\\circ} \1"), dummy),
   (r"\\Divbar{(.*?)}",  (lambda : r"\\overline{\\boldsymbol{\\nabla}} \cdot\1"), dummy),
   (r"\\Divy{(.*?)}",  (lambda : r"\\boldsymbol{\\nabla}_y \cdot\1"), dummy),
   (r"\\Divp{(.*?)}",  (lambda : r"\\boldsymbol{\\nabla}'\cdot\1"), dummy),
   (r"\\DivX{(.*?)}",  (lambda : r"\\boldsymbol{\\nabla}_{\\circ}\cdot\1"), dummy),
   (r"\\Curly{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}_y\\times\1"), dummy),
   (r"\\Curls{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}_s\\times\1"), dummy),
   (r"\\Curlp{(.*?)}", (lambda : r"\\boldsymbol{\\nabla}'\\times\1"), dummy),
   (r"\\Edot", (lambda : r"\\dot{e}"), dummy),
   (r"\\Gdot", (lambda : r"\\dot{g}"), dummy),
   (r"\\Jdot", (lambda : r"\\dot{J}"), dummy),
   (r"\\Qdot", (lambda : r"\\dot{q}"), dummy),
   (r"\\Tdot", (lambda : r"\\dot{T}"), dummy),
   (r"\\Etadot", (lambda : r"\\dot{\eta}"), dummy),
   (r"\\Bxdot", (lambda : r"\\dot{\\mathbf{x}}"), dummy),
   (r"\\BEdot", (lambda : r"\\dot{\\boldsymbol{E}}"), dummy),
   (r"\\BFdot", (lambda : r"\\dot{\\boldsymbol{F}}"), dummy),
   (r"\\BQdot", (lambda : r"\\dot{\\boldsymbol{Q}}"), dummy),
   (r"\\BSdot", (lambda : r"\\dot{\\boldsymbol{S}}"), dummy),
   (r"\\BUdot", (lambda : r"\\dot{\\boldsymbol{U}}"), dummy),
   (r"\\BVdot", (lambda : r"\\dot{\\boldsymbol{V}}"), dummy),
   (r"\\IntInfT", (lambda : r"\\int_{-\\infty}^t"), dummy),
   (r"\\IIIntInfInf", (lambda : r"\\int_{-\\infty}^{\\infty}\\int_{-\\infty}^{\\infty}\\int_{-\\infty}^{\\infty}"), dummy),
   (r"\\IIntInfInf", (lambda : r"\\int_{-\\infty}^{\\infty}\\int_{-\\infty}^{\\infty}"), dummy),
   (r"\\IntInfInf", (lambda : r"\\int_{-\\infty}^{\\infty}"), dummy),
   (r"\\IntInfZero", (lambda : r"\\int_{-\\infty}^{0}"), dummy),
   (r"\\IntZeroInf", (lambda : r"\\int_0^{\\infty}"), dummy),
   (r"\\IntSq", (lambda : r"\\int_{\\square}"), dummy),
   (r"\\IntDOmega", (lambda : r"\\int_{\\partial\\Omega}"), dummy),
   (r"\\IntDomegaint", (lambda : r"\\int_{\\partial \\Omega_{\\text{int}}}"), dummy),
   (r"\\IntDomega", (lambda : r"\\int_{\\partial \\Omega}"), dummy),
   (r"\\IntOmegac", (lambda : r"\\int_{\\Omega_c}"), dummy),
   (r"\\IntOmegapc", (lambda : r"\\int_{\\Omega_p\\cap\\Omega_c}"), dummy),
   (r"\\IntOmegap", (lambda : r"\\int_{\\Omega'}"), dummy),
   (r"\\IntOmegaq", (lambda : r"\\int_{\\Omega_q\\cap\\Omega}"), dummy),
   (r"\\IntGammat", (lambda : r"\\int_{\\Gamma_t}"), dummy),
   (r"\\IntGammau", (lambda : r"\\int_{\\Gamma_u}"), dummy),
   (r"\\IntGammaq", (lambda : r"\\int_{\\Gamma_q}"), dummy),
   (r"\\IntGammaT", (lambda : r"\\int_{\\Gamma_T}"), dummy),
   (r"\\IntGamma", (lambda : r"\\int_{\\Gamma}"), dummy),
   (r"\\IntOmegat", (lambda : r"\\int_{\\Omega(t)}"), dummy),
   (r"\\IntDomegat", (lambda : r"\\int_{\\partial\\Omega(t)}"), dummy),
   (r"\\IntOmegatDelt", (lambda : r"\\int_{\\Omega(t+\\Delta t)}"), dummy),
   (r"\\IntOmegar", (lambda : r"\\int_{\\Omega_0}"), dummy),
   (r"\\IntDomegar", (lambda : r"\\int_{\\partial \\Omega_0}"), dummy),
   (r"\\Balphahat", (lambda : r"\widehat{\\boldsymbol{\\alpha}}"), dummy),
   (r"\\Balpha", (lambda : r"\\boldsymbol{\\alpha}"), dummy),
   (r"\\Bbeta", (lambda : r"\\boldsymbol{\\beta}"), dummy),
   (r"\\Bdelta", (lambda : r"\\boldsymbol{\\delta}"), dummy),
   (r"\\Bpsi", (lambda : r"\\boldsymbol{\\psi}"), dummy),
   (r"\\LimDelt", (lambda : r"\\lim_{\\Delta t \\rightarrow 0}"), dummy),
   (r"\\Beq", (lambda : r"<div>\n $$"), toggle_math),
   (r"\\Eeq", (lambda : r"$$\n</div>"), toggle_math),
   (r"\\Bal", (lambda : r"\\begin{align}"), dummy),
   (r"\\Eal", (lambda : r"\\end{align}"), dummy),
   (r"\\Truesdell{(.*?)}", (lambda : r"\\overset{\\circ}{\1}"), dummy),
   (r"\\GreenNaghdi{(.*?)}", (lambda : r"\\overset{\\square}{\1}"), dummy),
   (r"\\Jaumann{(.*?)}", (lambda : r"\\overset{\\triangle}{\1}"), dummy),
   (r"\\Oldroyd{(.*?)}", (lambda : r"\\overset{\\triangledown}{\1}"), dummy),
   (r"\\Convective{(.*?)}", (lambda : r"\\overset{\\diamond}{\1}"), dummy),
   (r"\\Tand",  (lambda : r"\\text{and}"), dummy),
   (r"\\Tbar",  (lambda : r"\\text{bar}"), dummy),
   (r"\\Tball", (lambda : r"\\text{ball}"), dummy),
   (r"\\Tbody", (lambda : r"\\text{body}"), dummy),
   (r"\\Teff",  (lambda : r"\\text{eff}"), dummy),
   (r"\\Text",  (lambda : r"\\text{ext}"), dummy),
   (r"\\Tint",  (lambda : r"\\text{int}"), dummy),
   (r"\\Tkin",  (lambda : r"\\text{kin}"), dummy),
   (r"\\Tmin",  (lambda : r"\\text{min}"), dummy),
   (r"\\Tmax",  (lambda : r"\\text{max}"), dummy),
   (r"\\Tor",   (lambda : r"\\text{or}"), dummy),
   (r"\\Trial", (lambda : r"\\text{trial}"), dummy),
   (r"\\Tr",    (lambda : r"\\text{tr}"), dummy),
   (r"\\TAi",    (lambda : r"\\text{Ai}"), dummy),
   (r"\\TBi",    (lambda : r"\\text{Bi}"), dummy),
   (r"\\sgn",    (lambda : r"\\text{sgn}"), dummy),
   (r"\\Half",  (lambda : r"\\frac{1}{2}"), dummy),
   (r"\\SThr",  (lambda : r"\\sqrt{3}"), dummy),
   (r"\\STT",   (lambda : r"\\frac{\\sqrt{3}}{2}"), dummy),
   (r"\\Third", (lambda : r"\\frac{1}{3}"), dummy),
   (r"\\Gradu", (lambda : r"\\boldsymbol{\\nabla}\\mathbf{u}"), dummy),
   (r"\\Gradv", (lambda : r"\\boldsymbol{\\nabla}\\mathbf{v}"), dummy),
   (r"\\Divu",  (lambda : r"\\boldsymbol{\\nabla}\\cdot\\mathbf{u}"), dummy),
   (r"\\Divv",  (lambda : r"\\boldsymbol{\\nabla}\\cdot\\mathbf{v}"), dummy),
   (r"\\Curlu", (lambda : r"\\boldsymbol{\\nabla}\\times\\mathbf{u}"), dummy),
   (r"\\Curlv", (lambda : r"\\boldsymbol{\\nabla}\\times\\mathbf{v}"), dummy),
   (r"\\Dualn", (lambda : r"\\Bn\\otimes\\Bn}"), dummy),
   (r"\\Inner{(.*?)}{(.*?)}",  (lambda : r"\\left\\langle\1,~\2\\right\\rangle"), dummy),
   (r"\\Bcross{(.*?)}{(.*?)}", (lambda : r"\1\\times\2"), dummy),
   (r"\\Bdot{(.*?)}{(.*?)}",   (lambda : r"\1\\cdot\2"), dummy),
   (r"\\Dyad{(.*?)}{(.*?)}",   (lambda : r"\1\\otimes\2"), dummy),
   (r"\\GradO{(.*?)}",          (lambda : r"\\boldsymbol{\\nabla}_{\\circ} \1"), dummy),
   (r"\\Grad{(.*?)}",          (lambda : r"\\boldsymbol{\\nabla} \1"), dummy),
   (r"\\Lap{(.*?)}",           (lambda : r"\\nabla^2 \1"), dummy),
   (r"\\Biharm{(.*?)}",        (lambda : r"\\nabla^4 \1"), dummy),
   (r"\\Div{(.*?)}",           (lambda : r"\\boldsymbol{\\nabla} \\cdot \1"), dummy),
   (r"\\Curl{(.*?)}",          (lambda : r"\\boldsymbol{\\nabla} \\times \1"), dummy),
   (r"\\Over{(.*?)}",            (lambda : r"\\frac{1}{\1}"), dummy),
   (r"\\Diff{(.*?)}{(.*?)}",     (lambda : r"\\frac{d \1}{d \2}"), dummy),
   (r"\\Partial{(.*?)}{(.*?)}",  (lambda : r"\\frac{\\partial \1}{\\partial \2}"), dummy),
   (r"\\PPartial{(.*?)}{(.*?)}", (lambda : r"\\frac{\\partial^2 \1}{\\partial \2^2}"), dummy),
   (r"\\FPartial{(.*?)}{(.*?)}", (lambda : r"\\frac{\\partial^4 \1}{\\partial \2^4}"), dummy),
   (r"\\PPartialA{(.*?)}{(.*?)}{(.*?)}", (lambda : r"\\frac{\\partial^2 \1}{\\partial \2 \\partial \3}"), dummy),
   (r"\\FPartialA{(.*?)}{(.*?)}{(.*?)}", (lambda : r"\\frac{\\partial^4 \1}{\\partial \2^2 \\partial \3^2}"), dummy),
   (r"\\DotMbT",   (lambda : r"\\dot{\\mathbf{T}}"), dummy),
   (r"\\TildeMbT", (lambda : r"\\widetilde{\\mathbf{T}}"), dummy),
   (r"\\BarT",     (lambda : r"\\overline{T}"), dummy),
   (r"\\Barq",     (lambda : r"\\overline{q}"), dummy),
   (r"\\Jump{(.*?)}",          (lambda : r"\\llbracket\1\\rrbracket"), dummy),
   (r"\\Blimitx{(.*?)}",       (lambda : r"\\left[\1\\right]_{x_a}^{x_b}"), dummy),
   (r"\\Deriv{(.*?)}{(.*?)}",  (lambda : r"\\cfrac{d \1}{d \2}"), dummy),
   (r"\\MDeriv{(.*?)}{(.*?)}", (lambda : r"\\cfrac{D \1}{D \2}"), dummy),
   (r"\\DDDDeriv{(.*?)}{(.*?)}", (lambda : r"\\cfrac{d^4 \1}{d \2^4}"), dummy),
   (r"\\DDDeriv{(.*?)}{(.*?)}", (lambda : r"\\cfrac{d^3 \1}{d \2^3}"), dummy),
   (r"\\DDeriv{(.*?)}{(.*?)}", (lambda : r"\\cfrac{d^2 \1}{d \2^2}"), dummy),
   (r"\\Intx",      (lambda : r"\\int_{x_a}^{x_b}"), dummy),
   (r"\\IntX",      (lambda : r"\\int_{X_a}^{X_b}"), dummy),
   (r"\\Intiso",    (lambda : r"\\int_{-1}^{1}"), dummy),
   (r"\\IntOmegaA", (lambda : r"\\int_{\\Omega_0}"), dummy),
   (r"\\IntOmega",  (lambda : r"\\int_{\\Omega}"), dummy),
   (r"\\Domega",  (lambda : r"\\partial \\Omega "), dummy),
   (r"\\Norm{(.*?)}{(.*?)}", (lambda : r"\\lVert\1\\rVert_{\2}"), dummy),
   (r"\\Var{(.*?)}",         (lambda : r"\\delta \1"), dummy),
   (r"\\DelTwo", (lambda : r"\\Delta/2"), dummy),
   (r"\\DelT", (lambda : r"\\Delta t"), dummy),
   (r"\\BCalD", (lambda : r"\\boldsymbol{\\mathcal{D}}"), dummy),
   (r"\\CalD", (lambda : r"\\mathcal{D}"), dummy),
   (r"\\CalE", (lambda : r"\\mathcal{E}"), dummy),
   (r"\\CalF", (lambda : r"\\mathcal{F}"), dummy),
   (r"\\CalH", (lambda : r"\\mathcal{H}"), dummy),
   (r"\\CalJ", (lambda : r"\\mathcal{J}"), dummy),
   (r"\\CalL", (lambda : r"\\mathcal{L}"), dummy),
   (r"\\CalN", (lambda : r"\\mathcal{N}"), dummy),
   (r"\\BCalM", (lambda : r"\\boldsymbol{\\mathcal{M}}"), dummy),
   (r"\\CalM", (lambda : r"\\mathcal{M}"), dummy),
   (r"\\CalP", (lambda : r"\\mathcal{P}"), dummy),
   (r"\\CalS", (lambda : r"\\mathcal{S}"), dummy),
   (r"\\BCalS", (lambda : r"\\boldsymbol{\\mathcal{S}}"), dummy),
   (r"\\CalT", (lambda : r"\\mathcal{T}"), dummy),
   (r"\\CalU", (lambda : r"\\mathcal{U}"), dummy),
   (r"\\CalV", (lambda : r"\\mathcal{V}"), dummy),
   (r"\\CalW", (lambda : r"\\mathcal{W}"), dummy),
   (r"\\CalX", (lambda : r"\\mathcal{X}"), dummy),
   (r"\\AvPowerInf", (lambda : r"\\left\\langle \\boldsymbol{\\sigma}:\\dot{\\boldsymbol{\\varepsilon}}\\right\\rangle"), dummy),
   (r"\\AvPowerPF", (lambda : r"\\left\\langle \\boldsymbol{P}^T:\\dot{\\boldsymbol{F}\\right\\rangle"), dummy),
   (r"\\AvPower", (lambda : r"\\left\\langle\\boldsymbol{\\sigma}:\\boldsymbol{\\nabla}\\mathbf{v} \\right\\rangle"), dummy),
   (r"\\AvP", (lambda : r"\\left\\langle \\boldsymbol{P} \\right\\rangle"), dummy),
   (r"\\AvWorkInf", (lambda : r"\\left\\langle \\boldsymbol{\\sigma}:\\boldsymbol{\\varepsilon}\\right\\rangle"), dummy),
   (r"\\AvGradudot", (lambda : r"\\left\\langle\\boldsymbol{\\nabla} \\dot{\\mathbf{u}} \\right\\rangle"), dummy),
   (r"\\AvGradu", (lambda : r"\\left\\langle \\boldsymbol{\\nabla} \\mathbf{u} \\right\\rangle"), dummy),
   (r"\\AvGradv", (lambda : r"\\left\\langle\\boldsymbol{\\nabla}\\mathbf{v} \\right\\rangle"), dummy),
   (r"\\AvSigBar", (lambda : r"\\overline{\\boldsymbol{\\sigma}}"), dummy),
   (r"\\AvSig", (lambda : r"\\left\\langle \\boldsymbol{\\sigma}} \\right\\rangle"), dummy),
   (r"\\AvTauBar", (lambda : r"\\overline{\\boldsymbol{\\tau}}"), dummy),
   (r"\\AvTau", (lambda : r"\\left\\langle \\boldsymbol{\\tau} \\right\\rangle"), dummy),
   (r"\\AvEpsdot", (lambda : r" \\left\\langle \\dot{\\boldsymbol{\\varepsilon}}\\right\\rangle"), dummy),
   (r"\\AvEps", (lambda : r"\\left\\langle \\boldsymbol{\\varepsilon} \\right\\rangle"), dummy),
   (r"\\AvDisp", (lambda : r"\\left\\langle\\mathbf{u}\\right\\rangle"), dummy),
   (r"\\AvFdot", (lambda : r"\\left\\langle\\dot{\\boldsymbol{F}}\\right\\rangle"), dummy),
   (r"\\AvF", (lambda : r"\\left\\langle \\boldsymbol{F}\\right\\rangle"), dummy),
   (r"\\Avl", (lambda : r"\\overline{\\boldsymbol{l}}"), dummy),
   (r"\\AvOmega", (lambda : r"\\left\\langle \\boldsymbol{\\omega} \\right\\rangle"), dummy),
   (r"\\Av{(.*?)}", (lambda : r"\\left\\langle \1 \\right\\rangle"), dummy),
   (r"\\Comp{(.*?)}{(.*?)}",         (lambda : r"\1 \\circ \2"), dummy),
   (r"\\Map{(.*?)}{(.*?)}{(.*?)}",   (lambda : r"\1 : \2 \\rightarrow \3"), dummy),
   (r"\\MapTo{(.*?)}{(.*?)}{(.*?)}", (lambda : r"\1 : \2 \\mapsto \3"), dummy),
   (r"\\Real{(.*?)}",                (lambda : r"\\mathbb{R}^{\1}"), dummy),
   (r"\\Rea", (lambda : r"\\text{Re}"), dummy),
   (r"\\Img", (lambda : r"\\text{Im}"), dummy),
   (r"\\BGammahat", (lambda : r"\\boldsymbol{\\mathit{\\widehat{\\Gamma}}}"), dummy),
   (r"\\BGamma", (lambda : r"\\boldsymbol{\\mathit{\\Gamma}}"), dummy),
   (r"\\Bveps", (lambda : r"\\boldsymbol{\\varepsilon}"), dummy),
   (r"\\Veps", (lambda : r"\\varepsilon"), dummy),
   (r"\\BHat{(.*?)}", (lambda : r"\\hat{\\boldsymbol{\1}}"), dummy),
   (r"\\BTx", (lambda : r"\\tilde{\\mathbf{x}}"), dummy),
   (r"\\Beh", (lambda : r"\\hat{\\mathbf{e}}"), dummy),
   (r"\\BHex", (lambda : r"\\hat{\\mathbf{e}}_1"), dummy),
   (r"\\BHey", (lambda : r"\\hat{\\mathbf{e}}_2"), dummy),
   (r"\\BHez", (lambda : r"\\hat{\\mathbf{e}}_3"), dummy),
   (r"\\BHn{(.*?)}", (lambda : r"\\hat{\\mathbf{n}}_{\1}"), dummy),
   (r"\\BHe{(.*?)}", (lambda : r"\\hat{\\mathbf{e}}_{\1}"), dummy),
   (r"\\BHg{(.*?)}", (lambda : r"\\hat{\\mathbf{g}}_{\1}"), dummy),
   (r"\\BHG{(.*?)}", (lambda : r"\\hat{\\mathbf{G}}_{\1}"), dummy),
   (r"\\Hn", (lambda : r"\\hat{\\mathbf{n}}"), dummy),
   (r"\\Mbatilde", (lambda : r"\\tilde{\\mathbf{a}}"), dummy),
   (r"\\Mba", (lambda : r"\\mathbf{a}"), dummy),
   (r"\\Mbb", (lambda : r"\\mathbf{b}"), dummy),
   (r"\\Mbd", (lambda : r"\\mathbf{d}"), dummy),
   (r"\\Mbf", (lambda : r"\\mathbf{f}"), dummy),
   (r"\\Mbntilde", (lambda : r"\\tilde{\\mathbf{n}}"), dummy),
   (r"\\Mbn", (lambda : r"\\mathbf{n}"), dummy),
   (r"\\Mbr", (lambda : r"\\mathbf{r}"), dummy),
   (r"\\Mbu", (lambda : r"\\mathbf{u}"), dummy),
   (r"\\Mbv", (lambda : r"\\mathbf{v}"), dummy),
   (r"\\Mbx", (lambda : r"\\mathbf{x}"), dummy),
   (r"\\MbA", (lambda : r"\\mathbf{A}"), dummy),
   (r"\\MbB", (lambda : r"\\mathbf{B}"), dummy),
   (r"\\MbC", (lambda : r"\\mathbf{C}"), dummy),
   (r"\\MbD", (lambda : r"\\mathbf{D}"), dummy),
   (r"\\MbF", (lambda : r"\\mathbf{F}"), dummy),
   (r"\\MbHbar", (lambda : r"\\overline{\\mathbf{H}}"), dummy),
   (r"\\MbH", (lambda : r"\\mathbf{H}"), dummy),
   (r"\\MbI", (lambda : r"\\mathbf{I}"), dummy),
   (r"\\MbKbar", (lambda : r"\\overline{\\mathbf{K}}"), dummy),
   (r"\\MbKtilde", (lambda : r"\\tilde{\\mathbf{K}}"), dummy),
   (r"\\MbK", (lambda : r"\\mathbf{K}"), dummy),
   (r"\\MbM", (lambda : r"\\mathbf{M}"), dummy),
   (r"\\MbN", (lambda : r"\\mathbf{N}"), dummy),
   (r"\\MbPbar", (lambda : r"\\overline{\\mathbf{P}}"), dummy),
   (r"\\MbP", (lambda : r"\\mathbf{P}"), dummy),
   (r"\\MbR", (lambda : r"\\mathbf{R}"), dummy),
   (r"\\MbT", (lambda : r"\\mathbf{T}"), dummy),
   (r"\\MbU", (lambda : r"\\mathbf{U}"), dummy),
   (r"\\MbV", (lambda : r"\\mathbf{V}"), dummy),
   (r"\\MbX", (lambda : r"\\mathbf{X}"), dummy),
   (r"\\Mbone", (lambda : r"\\mathbf{1}"), dummy),
   (r"\\Mbzero", (lambda : r"\\mathbf{0}"), dummy),
   (r"\\MbSig", (lambda : r"\\boldsymbol{\\sigma}"), dummy),
   (r"\\Mb", (lambda : r"\\left[\\mathsf{b}\\right]"), dummy),
   (r"\\Mu", (lambda : r"\\left[\\mathsf{u}\\right]"), dummy),
   (r"\\Mv", (lambda : r"\\left[\\mathsf{v}\\right]"), dummy),
   (r"\\Mw", (lambda : r"\\left[\\mathsf{w}\\right]"), dummy),
   (r"\\Mx", (lambda : r"\\left[\\mathsf{x}\\right]"), dummy),
   (r"\\MA", (lambda : r"\\left[\\mathsf{A}\\right]"), dummy),
   (r"\\MC", (lambda : r"\\left[\\mathsf{C}\\right]"), dummy),
   (r"\\MD", (lambda : r"\\left[\\mathsf{D}\\right]"), dummy),
   (r"\\MH", (lambda : r"\\left[\\mathsf{H}\\right]"), dummy),
   (r"\\MI", (lambda : r"\\left[\\mathsf{I}\\right]"), dummy),
   (r"\\ML", (lambda : r"\\left[\\mathsf{L}\\right]"), dummy),
   (r"\\MM", (lambda : r"\\left[\\mathsf{M}\\right]"), dummy),
   (r"\\MN", (lambda : r"\\left[\\mathsf{N}\\right]"), dummy),
   (r"\\MP", (lambda : r"\\left[\\mathsf{P}\\right]"), dummy),
   (r"\\MR", (lambda : r"\\left[\\mathsf{R}\\right]"), dummy),
   (r"\\MT", (lambda : r"\\left[\\mathsf{T}\\right]"), dummy),
   (r"\\MV", (lambda : r"\\left[\\mathsf{V}\\right]"), dummy),
   (r"\\Mone", (lambda : r"\\left[\\mathsf{1}\\right]"), dummy),
   (r"\\Mzero", (lambda : r"\\left[\\mathsf{0}\\right]"), dummy),
   (r"\\SfZero", (lambda : r"\\boldsymbol{\\mathsf{0}}"), dummy),
   (r"\\SfOne", (lambda : r"\\boldsymbol{\\mathsf{1}}"), dummy),
   (r"\\Sfa", (lambda : r"\\boldsymbol{\\mathsf{a}}"), dummy),
   (r"\\Sfc", (lambda : r"\\boldsymbol{\\mathsf{c}}"), dummy),
   (r"\\SfA", (lambda : r"\\boldsymbol{\\mathsf{A}}"), dummy),
   (r"\\SfC", (lambda : r"\\boldsymbol{\\mathsf{C}}"), dummy),
   (r"\\SfD", (lambda : r"\\boldsymbol{\\mathsf{D}}"), dummy),
   (r"\\SfI", (lambda : r"\\boldsymbol{\\mathsf{I}}"), dummy),
   (r"\\SfK", (lambda : r"\\boldsymbol{\\mathsf{K}}"), dummy),
   (r"\\SfL", (lambda : r"\\boldsymbol{\\mathsf{L}}"), dummy),
   (r"\\SfS", (lambda : r"\\boldsymbol{\\mathsf{S}}"), dummy),
   (r"\\SfT", (lambda : r"\\boldsymbol{\\mathsf{T}}"), dummy),
   (r"\\Msig", (lambda : r"\\left[\\boldsymbol{\\sigma}\\right]"), dummy),
   (r"\\Meps", (lambda : r"\\left[\\boldsymbol{\\varepsilon}\\right]"), dummy),
   (r"\\Ex", (lambda : r"\\mathbf{e}_1"), dummy),
   (r"\\Ey", (lambda : r"\\mathbf{e}_2"), dummy),
   (r"\\Ez", (lambda : r"\\mathbf{e}_3"), dummy),
   (r"\\Exp", (lambda : r"\\mathbf{e}^{'}_1"), dummy),
   (r"\\Eyp", (lambda : r"\\mathbf{e}^{'}_2"), dummy),
   (r"\\Ezp", (lambda : r"\\mathbf{e}^{'}_3"), dummy),
   (r"\\Epsxx", (lambda : r"\\varepsilon_{11}"), dummy),
   (r"\\Epsyy", (lambda : r"\\varepsilon_{22}"), dummy),
   (r"\\Epszz", (lambda : r"\\varepsilon_{33}"), dummy),
   (r"\\Epsyz", (lambda : r"\\varepsilon_{23}"), dummy),
   (r"\\Epszx", (lambda : r"\\varepsilon_{31}"), dummy),
   (r"\\Epsxy", (lambda : r"\\varepsilon_{12}"), dummy),
   (r"\\Sigxx", (lambda : r"\\sigma_{11}"), dummy),
   (r"\\Sigyy", (lambda : r"\\sigma_{22}"), dummy),
   (r"\\Sigzz", (lambda : r"\\sigma_{33}"), dummy),
   (r"\\Sigyz", (lambda : r"\\sigma_{23}"), dummy),
   (r"\\Sigzx", (lambda : r"\\sigma_{31}"), dummy),
   (r"\\Sigxy", (lambda : r"\\sigma_{12}"), dummy),
   (r"\\Eps{(.*?)}", (lambda : r"\\varepsilon_{\1}"), dummy),
   (r"\\Sig{(.*?)}", (lambda : r"\\sigma_{\1}"), dummy),
   (r"\\X", (lambda : r"X_1"), dummy),
   (r"\\Y", (lambda : r"X_2"), dummy),
   (r"\\Z", (lambda : r"X_3"), dummy),
   (r"\\Dhat", (lambda : r"\\widehat{D}"), dummy),
   (r"\\Ehat", (lambda : r"\\widehat{E}"), dummy),
   (r"\\Fhat", (lambda : r"\\widehat{F}"), dummy),
   (r"\\Phat", (lambda : r"\\widehat{P}"), dummy),
   (r"\\Uhat", (lambda : r"\\widehat{U}"), dummy),
   (r"\\Vhat", (lambda : r"\\widehat{V}"), dummy),
   (r"\\fhat", (lambda : r"\\widehat{f}"), dummy),
   (r"\\ghat", (lambda : r"\\widehat{g}"), dummy),
   (r"\\phat", (lambda : r"\\widehat{p}"), dummy),
   (r"\\uhat", (lambda : r"\\widehat{u}"), dummy),
   (r"\\vhat", (lambda : r"\\widehat{v}"), dummy),
   (r"\\xhat", (lambda : r"\\widehat{x}"), dummy),
   (r"\\yhat", (lambda : r"\\widehat{y}"), dummy),
   (r"\\rhohat", (lambda : r"\\widehat{\\rho}"), dummy),
   (r"\\phihat", (lambda : r"\\widehat{\\varphi}"), dummy),
   (r"\\sighat", (lambda : r"\\widehat{\\sigma}"), dummy),
   (r"\\epshat", (lambda : r"\\widehat{\\epsilon}"), dummy),
   (r"\\kappahat", (lambda : r"\\widehat{\\kappa}"), dummy),
   (r"\\ktilde", (lambda : r"\\tilde{k}"), dummy),
   (r"\\Etilde", (lambda : r"\\tilde{E}"), dummy),
   (r"\\Htilde", (lambda : r"\\tilde{H}"), dummy),
   (r"\\Rtilde", (lambda : r"\\tilde{R}"), dummy),
   (r"\\Ttilde", (lambda : r"\\tilde{T}"), dummy),
   (r"\\epseff", (lambda : r"\\epsilon_{\\text{eff}}"), dummy),
   (r"\\mueff", (lambda : r"\\mu_{\\text{eff}}"), dummy),
   (r"\\BPi", (lambda : r"\\boldsymbol{\\Pi}"), dummy),
   (r"\\Btheta", (lambda : r"\\boldsymbol{\\theta}"), dummy),
   (r"\\Brho", (lambda : r"\\boldsymbol{\\rho}"), dummy),
   (r"\\Bpi", (lambda : r"\\boldsymbol{\\pi}"), dummy),
   (r"\\Bchi", (lambda : r"\\boldsymbol{\\chi}"), dummy),
   (r"\\BVeps", (lambda : r"\\boldsymbol{\\varepsilon}"), dummy),
   (r"\\Bepseff", (lambda : r"\\boldsymbol{\\epsilon}_{\\text{eff}}"), dummy),
   (r"\\Bepshat", (lambda : r"\\widehat{\\boldsymbol{\\epsilon}}"), dummy),
   (r"\\Beps", (lambda : r"\\boldsymbol{\\epsilon}"), dummy),
   (r"\\Bkappa", (lambda : r"\\boldsymbol{\\kappa}"), dummy),
   (r"\\Bbeps", (lambda : r"\\bar{\\boldsymbol{\\varepsilon}}"), dummy),
   (r"\\Bnabla", (lambda : r"\\boldsymbol{\\nabla}"), dummy),
   (r"\\BOmega", (lambda : r"\\boldsymbol{\\Omega}"), dummy),
   (r"\\Bomega", (lambda : r"\\boldsymbol{\\omega}"), dummy),
   (r"\\Bsighat", (lambda : r"\\widehat{\\boldsymbol{\\sigma}}"), dummy),
   (r"\\Bsig", (lambda : r"\\boldsymbol{\\sigma}"), dummy),
   (r"\\Btau", (lambda : r"\\boldsymbol{\\tau}"), dummy),
   (r"\\Bvarphi", (lambda : r"\\boldsymbol{\\varphi}"), dummy),
   (r"\\Blambda", (lambda : r"\\boldsymbol{\\lambda}"), dummy),
   (r"\\Bmueff", (lambda : r"\\boldsymbol{\\mu}_{\\text{eff}}"), dummy),
   (r"\\Bmu", (lambda : r"\\boldsymbol{\\mu}"), dummy),
   (r"\\Bxi", (lambda : r"\\boldsymbol{\\xi}"), dummy),
   (r"\\Bonev", (lambda : r"\\boldsymbol{1}"), dummy),
   (r"\\Bone", (lambda : r"\\boldsymbol{\\mathit{1}}"), dummy),
   (r"\\BzeroT", (lambda : r"\\boldsymbol{\\mathit{0}}"), dummy),
   (r"\\Bzero", (lambda : r"\\boldsymbol{0}"), dummy),
   (r"\\Bahat", (lambda : r"\\widehat{\\mathbf{a}}"), dummy),
   (r"\\Ba", (lambda : r"\\mathbf{a}"), dummy),
   (r"\\Bbhat", (lambda : r"\\widehat{\\mathbf{b}}"), dummy),
   (r"\\BbT", (lambda : r"\\boldsymbol{b}"), dummy),
   (r"\\Bb", (lambda : r"\\mathbf{b}"), dummy),
   (r"\\Bc", (lambda : r"\\mathbf{c}"), dummy),
   (r"\\BdT", (lambda : r"\\boldsymbol{d}"), dummy),
   (r"\\Bd", (lambda : r"\\mathbf{d}"), dummy),
   (r"\\BeT", (lambda : r"\\boldsymbol{e}"), dummy),
   (r"\\Be", (lambda : r"\\mathbf{e}"), dummy),
   (r"\\BfT", (lambda : r"\\boldsymbol{f}"), dummy),
   (r"\\Bf", (lambda : r"\\mathbf{f}"), dummy),
   (r"\\BgT", (lambda : r"\\boldsymbol{g}"), dummy),
   (r"\\Bg", (lambda : r"\\mathbf{g}"), dummy),
   (r"\\Bh", (lambda : r"\\mathbf{h}"), dummy),
   (r"\\Bk", (lambda : r"\\mathbf{k}"), dummy),
   (r"\\BlT", (lambda : r"\\boldsymbol{l}"), dummy),
   (r"\\Bl", (lambda : r"\\mathbf{l}"), dummy),
   (r"\\Bm", (lambda : r"\\mathbf{m}"), dummy),
   (r"\\Bn", (lambda : r"\\mathbf{n}"), dummy),
   (r"\\Bo", (lambda : r"\\mathbf{o}"), dummy),
   (r"\\Bp", (lambda : r"\\mathbf{p}"), dummy),
   (r"\\Bq", (lambda : r"\\mathbf{q}"), dummy),
   (r"\\Br", (lambda : r"\\mathbf{r}"), dummy),
   (r"\\Bs", (lambda : r"\\mathbf{s}"), dummy),
   (r"\\Bt", (lambda : r"\\mathbf{t}"), dummy),
   (r"\\Buhat", (lambda : r"\\widehat{\\mathbf{u}}"), dummy),
   (r"\\Bu", (lambda : r"\\mathbf{u}"), dummy),
   (r"\\Bv", (lambda : r"\\mathbf{v}"), dummy),
   (r"\\BwT", (lambda : r"\\boldsymbol{w}"), dummy),
   (r"\\Bw", (lambda : r"\\mathbf{w}"), dummy),
   (r"\\Bx", (lambda : r"\\mathbf{x}"), dummy),
   (r"\\By", (lambda : r"\\mathbf{y}"), dummy),
   (r"\\BAeff", (lambda : r"\\boldsymbol{A}_{\\text{eff}}"), dummy),
   (r"\\BAhat", (lambda : r"\\widehat{\\boldsymbol{A}}"), dummy),
   (r"\\BAv", (lambda : r"\\mathbf{A}"), dummy),
   (r"\\BA", (lambda : r"\\boldsymbol{A}"), dummy),
   (r"\\BBhatv", (lambda : r"\\widehat{\\mathbf{B}}"), dummy),
   (r"\\BBhat", (lambda : r"\\widehat{\\boldsymbol{B}}"), dummy),
   (r"\\BBv", (lambda : r"\\mathbf{B}"), dummy),
   (r"\\BB", (lambda : r"\\boldsymbol{B}"), dummy),
   (r"\\BC", (lambda : r"\\boldsymbol{C}"), dummy),
   (r"\\BDhatv", (lambda : r"\\widehat{\\mathbf{D}}"), dummy),
   (r"\\BDtildev", (lambda : r"\\mathbf{\\tilde{D}}"), dummy),
   (r"\\BDv", (lambda : r"\\mathbf{D}"), dummy),
   (r"\\BD", (lambda : r"\\boldsymbol{D}"), dummy),
   (r"\\BEhatv", (lambda : r"\\widehat{\\mathbf{E}}"), dummy),
   (r"\\BEv", (lambda : r"\\mathbf{E}"), dummy),
   (r"\\BE", (lambda : r"\\boldsymbol{E}"), dummy),
   (r"\\BFhatv", (lambda : r"\\widehat{\\mathbf{F}}"), dummy),
   (r"\\BFv", (lambda : r"\\mathbf{F}"), dummy),
   (r"\\BF", (lambda : r"\\boldsymbol{F}"), dummy),
   (r"\\BG", (lambda : r"\\boldsymbol{G}"), dummy),
   (r"\\BHhatv", (lambda : r"\\widehat{\\mathbf{H}}"), dummy),
   (r"\\BHhat", (lambda : r"\\widehat{\\boldsymbol{H}}"), dummy),
   (r"\\BHv", (lambda : r"\\mathbf{H}"), dummy),
   (r"\\BH", (lambda : r"\\boldsymbol{H}"), dummy),
   (r"\\BI", (lambda : r"\\boldsymbol{I}"), dummy),
   (r"\\BJv", (lambda : r"\\mathbf{J}"), dummy),
   (r"\\BJ", (lambda : r"\\boldsymbol{J}"), dummy),
   (r"\\BKbar", (lambda : r"\\boldsymbol{\\bar{K}}"), dummy),
   (r"\\BK", (lambda : r"\\boldsymbol{K}"), dummy),
   (r"\\BL", (lambda : r"\\boldsymbol{L}"), dummy),
   (r"\\BMv", (lambda : r"\\mathbf{M}"), dummy),
   (r"\\BM", (lambda : r"\\boldsymbol{M}"), dummy),
   (r"\\BNv", (lambda : r"\\mathbf{N}"), dummy),
   (r"\\BN", (lambda : r"\\boldsymbol{N}"), dummy),
   (r"\\BPhatv", (lambda : r"\\widehat{\\mathbf{P}}"), dummy),
   (r"\\BPv", (lambda : r"\\mathbf{P}"), dummy),
   (r"\\BP", (lambda : r"\\boldsymbol{P}"), dummy),
   (r"\\BQv", (lambda : r"\\mathbf{Q}"), dummy),
   (r"\\BQ", (lambda : r"\\boldsymbol{Q}"), dummy),
   (r"\\BRv", (lambda : r"\\mathbf{R}"), dummy),
   (r"\\BR", (lambda : r"\\boldsymbol{R}"), dummy),
   (r"\\BS", (lambda : r"\\boldsymbol{S}"), dummy),
   (r"\\BTv", (lambda : r"\\mathbf{T}"), dummy),
   (r"\\BT", (lambda : r"\\boldsymbol{T}"), dummy),
   (r"\\BU", (lambda : r"\\boldsymbol{U}"), dummy),
   (r"\\BVhatv", (lambda : r"\\widehat{\\mathbf{V}}"), dummy),
   (r"\\BVv", (lambda : r"\\mathbf{V}"), dummy),
   (r"\\BV", (lambda : r"\\boldsymbol{V}"), dummy),
   (r"\\BW", (lambda : r"\\boldsymbol{W}"), dummy),
   (r"\\BXT", (lambda : r"\\boldsymbol{X}"), dummy),
   (r"\\BX", (lambda : r"\\mathbf{X}"), dummy),
   (r"\\BY", (lambda : r"\\boldsymbol{Y}"), dummy),
   (r"\\pbar", (lambda : r"\\bar{p}"), dummy)
]
 
in_stream  = sys.stdin;
out_stream = sys.stdout

# precompile regular expressions
tr_list3 = map(lambda x: (re.compile(x[0]),x[1],x[2]),tr_list2)

label = re.compile(r"\\label{eq:(.*?)}")
labelfig = re.compile(r"\\label{fig:(.*?)}")
for line in in_stream.readlines():

	if label.search(line):
		mch = label.search(line)
#		print mch
#		print line
		count +=1
		label_list.append(mch.group(1))
		line = label.sub(r"\\text{("+'%d'%count+r")} \\qquad ",line)
	if labelfig.search(line):
		mch = labelfig.search(line)
#		print mch
#		print line
		countfig +=1
		labelfig_list.append(mch.group(1))
		line = labelfig.sub(r"\\text{("+'%d'%countfig+r")} ",line)
	for reg in tr_list3:
		if reg[0].search(line):
			reg[2]()
			if reg[1] != None:
                          if ('%' in reg[1]()):
#                            print reg[0].pattern, reg[1](), line
			    line = reg[0].sub(reg[1](),line)
#                            print line
                          else:
			    line = reg[0].sub(reg[1](),line)
#			    print line
			else:
				line = reg[0].sub("", line)
	if bdoc != None:
#  	        print line
		lines.append(line)


ref = re.compile(r"\\ref{eq:(.*?)}")
reffig = re.compile(r"\\ref{fig:(.*?)}")
for line in lines:
	if ref.search(line):
		list = ref.findall(line)
		for item in list:
			r=re.compile(r"\\ref{eq:" + item + r"}")
			if item in label_list:
				line = r.sub('%d'%(label_list.index(item)+1), line)
	if reffig.search(line):
		list = reffig.findall(line)
		for item in list:
			r=re.compile(r"\\ref{fig:" + item + r"}")
			if item in labelfig_list:
				line = r.sub('%d'%(labelfig_list.index(item)+1), line)
	print >> out_stream, line,


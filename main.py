import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp, arcsin, arccos, arctan2, sinh, cosh
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi


plt.rcParams["text.usetex"] = True
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
}

plt.rcParams.update(tex_fonts)
plt.rc('text', usetex=True)







def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# 1. teil: quecksilber

x11 = 1.99/1000 #
x21 = 1.67/1000 #
N1 = 201 # pm 10


x12 = 2.91/1000
x22 = 2.64/1000
N2 = 200


x13 = 4.865/1000
x23 = 4.52/1000
N3 = 200

#natrium
x11 = 4.485/1000
x21 = 4.19/1000
N1 = 200

x12 = 4.19/1000
x22 = 3.88
N2 = 206

#Schwebung:

s0 = 6.44
s1 = 4.96
s2 = 3.44
s3 = 1.97
s4 = 0.5


s12 = 7.88
s22 = 9.33
s32 = 10.79
s42 = 12.24
s52 = 13.63

# auswerteformeln: delta_lambda = 2*lambda_m^2/lambda_s

x1 = 4.23/1000
#platte einbringen, und neu suchen
x2 = 1.27/1000
dicke = 1.1/1000 #pm 0.05/1000

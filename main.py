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



matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)




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

x1 = 1.99/1000 #
sigma_x1 = 0.05/1000
x2 = 1.67/1000 #
sigma_x2 = 0.05/1000
N1 = 201 # pm 10
sigma_N1 = 40

ü1 = N1*546.1e-9/(x1 - x2)
sigma_ü1 = gauss("N1*546.1*10**-9/(x1 - x2)")



x1 = 2.91/1000
sigma_x1 = 0.05/1000
x2 = 2.64/1000
sigma_x2 = 0.05/1000
N1 = 200
sigma_N1 = 40

ü2 = N1*(546.1e-9)/(x1 - x2)
sigma_ü2 = gauss("N1*(546.1*10**-9)/(x1 - x2)")

x1 = 4.865/1000
sigma_x1 = 0.05/1000
x2 = 4.52/1000
sigma_x2 = 0.05/1000
N1 = 200
sigma_N1 = 40

ü3 = N1*(546.1e-9)/(x1 - x2)
sigma_ü3 = gauss("N1*(546.1*10**-9)/(x1 - x2)")

u_ges = mean(array([ü1, ü2, ü3]))
sigma_u_ges = max(array([sigma_ü1, sigma_ü2, sigma_ü3]))/sqrt(3)



#natrium
x1 = 4.485/1000
x2 = 4.19/1000
N1 = 200
sigma_N1 = 10

lamda_m1 = u_ges*(x1 - x2)/N1
sigma_lamda_m1 = gauss("u_ges*(x1 - x2)/N1")

x1 = 4.19/1000
x2 = 3.88/1000
N1 = 206

lamda_m2 = u_ges*(x1 - x2)/N1
sigma_lamda_m2 = gauss("u_ges*(x1 - x2)/N1")

lamda_m = mean(array([lamda_m1, lamda_m2]))
sigma_lamda_m = max(array([sigma_lamda_m1, sigma_lamda_m2]))/sqrt(2)

#Schwebung:

s4 = 6.44
s3 = 4.96
s2 = 3.44
s1 = 1.97
s0 = 0.5


s5 = 7.88
s6 = 9.33
s7 = 10.79
s8 = 12.24
s9 = 13.63

Y = array([s0, s1, s2, s3, s4, s5, s6, s7, s8, s9])
sigma_Y = 0.1*ones(10)
X = array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

def linear(x, a):
    return a*x + 0.5


errorbar(X, Y, sigma_Y, None,'x', label='Schwebung')
optimizedParameters1, s = opt.curve_fit(linear, X, Y)
plot(X, linear(X, *optimizedParameters1), label="fit1")
xlabel('Index des Kontrastminimums', fontsize=20)
ylabel('Abstand in mm', fontsize=20)
legend(fontsize=13, loc='center left')
grid()
plt.tight_layout()
savefig('SchwebNa')
show()

delta_x = optimizedParameters1[0]
sigma_delta_x = np.diag(s)[0]

lamda_s = u_ges*delta_x/1000
sigma_lamda_s = gauss("u_ges*delta_x/1000")

delta_lamda = 2*lamda_m**2/lamda_s
sigma_delta_lamda = gauss("2*lamda_m**2/lamda_s")


# auswerteformeln: delta_lambda = 2*lambda_m^2/lambda_s


#wolfram
x1 = 4.23/1000
#platte einbringen, und neu suchen
x2 = 1.27/1000
dicke = 1.5/1000 #pm 0.05/1000
sigma_dicke = 0.05/1000

n = (u_ges*(x1 - x2) + dicke)/dicke
sigma_n = gauss("(u_ges*(x1 - x2) + dicke)/dicke")


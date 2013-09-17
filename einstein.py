import sys
from math import *
#from numpy import coth
# Usage python einstein.py T [spec] [#] [K in eV/A^2]


class Spec: # Holds information regarding each species
    def __init__(self):
        self.name = []
        self.num = []
        self.K = []
        self.u = []
        self.eps = []
        self.U = []
        self.F = []
        self.Cv = []
# Constants
hbar = 6.5821192e-16 #eV*S
k = 8.6173324e-5 #eV/K
eVtoHa =   0.036749325

massdict = { 'Fe':55.845, 'Si':28.0855, 'O':15.9994, 'Mg':24.3050, 'S':32.065, 'C':12.0107, 'N':14.0067, 'H':1.00794, 'He':4.0026202 }
#massdict = { 'Fe':55.845, 'Si':28.0855, 'O':15.9994 }


# Read in arguments
args = sys.argv[1::]
T = float( args.pop(0))

spec = Spec()

i = 0
while i == 0 :
    try: 
        float(args[0])
        i = 1
    except:
        spec.name.append( args.pop(0))

while i == 1:
    try:
        int(args[0])
        spec.num.append( int(args.pop(0)))
    except:
        i = 2
while i == 2:
    try:
        float(args[0])
        spec.K.append( float(args.pop(0)))
    except:
        i = 3

#print spec.num
#print spec.name

nomatch = 'number of values do not much'
assert len(spec.name) == len(spec.num), nomatch
assert len(spec.name) == len(spec.K), nomatch

# Calculate values per species
for i in range(len(spec.name)):
    u = massdict[spec.name[i]]
    m = u*9.31494061e8 #eV/c^2
    omega = sqrt(float(spec.K[i])/m*8.98755179e36) # s^-1  [c^2/A^2] * 8.98755179e36 = [ Hz ]^2
    eps = hbar*omega #eV
    Cv = 3*k*pow(eps/(2*k*T),2) /(sinh( eps / (2*k*T)))
    U = 3 * .5 * eps * cosh( eps / (2*k*T)) / sinh( eps / (2*k*T)) 
    F = 3 * k * T * log(2*sinh(eps/(2*k*T)))

    spec.u.append(u)
    spec.eps.append(eps)
    spec.Cv.append(Cv)
    spec.U.append(U)
    spec.F.append(F)

# sum F for system
F = 0.0
for i in range(len(spec.name)):
    F = F + float(spec.num[i])*spec.F[i]
#    print spec.name[i], spec.num[i], spec.F[i] 

print F,' eV ', F*eVtoHa,' Ha ' 


from ngsolve import *

omega = 2*pi*50     # Frequenz 50 Hz
mu0 = 4*pi*1e-7     # Magnetische Feldkonstante

sigmaCu = 56e6      # S/m (Kupfer)
sigmaCore = 56e6    # S/m (Kupfer)

lamAir = 0.0262     # W/(m K)
lamCu = 400         # W/(m K)
lamCore = 400       # W/(m K)

# Reihenfolge: air, copper, core
sigma = CoefficientFunction([0 , sigmaCu , sigmaCore ])
mur = CoefficientFunction([1 ,1 ,1])
nu = CoefficientFunction(1/mu0 * 1/mur ) # Permeabilität des Materials
lam = CoefficientFunction([lamAir, lamCu, lamCore])
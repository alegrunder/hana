from ngsolve import *

omega = 2*pi*50       # Frequenz 50 Hz
mu0 = 4*pi*1e-7       # Magnetische Feldkonstante

sigmaCu = 56e6        # S/m (Kupfer)
sigmaCoreTopf = 56e6  # S/m (Core und Topf - Kupfer)
sigmaWater = 0.01     # S/m (Leitungswasser: 0.005 ... 0.05 S/m)

lamAir = 0.0262       # W/(m K)
lamCu = 400           # W/(m K)
lamCoreTopf = 400     # W/(m K) (Core und Topf)
lamWater = 0.6        # W/(m K)

# Reihenfolge: air, copper, core, water
sigma = CoefficientFunction([0 , sigmaCu , sigmaCoreTopf, sigmaWater])
mur = CoefficientFunction([1 ,1 ,1, 1])
nu = CoefficientFunction(1/mu0 * 1/mur ) # Permeabilität des Materials
lam = CoefficientFunction([lamAir, lamCu, lamCoreTopf, lamWater])
from ngsolve import *

### Definition FE-Raum für Wirbelstromproblem
V = H1(mesh,order=3, complex=True, dirichlet='rotsym')
u,v = V.TnT()

### Definition Wirbelstromproblem
# Wechsel in Zylinderkoordinaten (es wird nur ein Querschnitt dargestellt (phi = 0))
r = y
z = x
uz, ur = grad(u)
vz, vr = grad(v) 

# Schwache Gleichung als Bilinearform
a = BilinearForm(V)
# erster Bilinearterm (sigma=0 in Luft)
a += 1j*omega*sigma*u*v*dx(definedon=mesh.Materials('copper|core')) 
# zweiter Bilinearterm
a += (nu*(1/r*u+ur)*vr + nu*uz*vz)*dx 

# Richtung und Stärke der externen Stromquelle
Jimp = CoefficientFunction([0,ic/(r_wire**2*pi),0])

f = LinearForm(V)
f += Jimp*v*dx(definedon=mesh.Materials('copper'))

def solveWirbelstromproblem(gfu:GridFunction):
    global a, f
    a.Assemble()
    f.Assemble()
    gfu.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec
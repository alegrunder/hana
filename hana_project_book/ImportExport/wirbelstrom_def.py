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
a += 1j*omega*sigma*u*v*dx(definedon=mesh.Materials('copper|core|water')) 
# zweiter Bilinearterm
a += (nu*(1/r*u+ur)*vr + nu*uz*vz)*dx 

# Richtung und Stärke der externen Stromquelle in Jimp enthalten
f = LinearForm(V)
f += Jimp*v*dx(definedon=mesh.Materials('copper'))

# Assemling
a.Assemble()
f.Assemble()

# Lösen des Systems
gfu.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec
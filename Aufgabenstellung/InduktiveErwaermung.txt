# Projekt Code Induktive Erwärmung

Projekt Code als Text-File aus dem PDF **ohne funktionierenden** Zusammenhang.

```python
### Listing 1 : Definition FE-Raum für Wirbelstromproblem
V = H1(mesh,order=3, complex = True, dirichlet='rotsym')
u,v = V.TnT()
gfu = GridFunction(V)

### Listing 2 : Definition Wirbelstrom Problem
r = y
z = x

uz, ur = grad(u) # partielle Ableitung nach z, r
vz, vr = grad(v) 

a = BilinearForm(V)
a += (nu*(1/r*u+ur)*vr + nu*uz*vz)*dx
a += 1j*omega*sigma*u*v*
     dx(definedon=mesh.Materials('copper|core'))

Jimp = CoefficientFunction([0,ic/(r_wire**2*pi),0])

f = LinearForm(V)
f += Jimp*v*dx(definedon=mesh.Materials('copper'))

a.Assemble()
f.Assemble()

### Listing 3 : Berechnung der induktiven Wärmequelle
Ez = -1j*omega*gfu
Jz = sigma*Ez
Jtot = Jz + Jimp
Qe = 1/2*Norm(InnerProduct(Jtot,Conj(Ez)))

### Listing 4 : FE-Raum für die Temperaturverteilung
V2 = H1(mesh, order = 4, dirichlet = 'outer')
uT,vT = V2.TnT()
gfT = GridFunction(V2)

### Listing 5 : Bilinearform Temperaturfeld
aT = BilinearForm(V2)
aT += lam*grad(uT)*grad(vT)*dx

fT = LinearForm(V2)
fT += Qe*vT*dx

aT.Assemble()
fT.Assemble()

### Listing 6 : Geometrie Beispiel
from netgen.geom2d import SplineGeometry
l_air = 0.1
r_wire = 0.001
nx_wire = 10
ny_wire = 3
dxdy_wire = 0.5*r_wire

rA_core = 0.01
rI_core = 0.008
l_core = 0.08
epsZero = 1e-8

geo = SplineGeometry()
geo.AddRectangle(p1=(-l_air/2,epsZero),
                 p2=(l_air/2,l_air/2),
                 bcs=["rotsym","outer","outer","outer"],
                 leftdomain=1,
                 rightdomain=0)
pts = []
for i in range(nx_wire):
    for j in range(ny_wire):
        pts.append([
            i*(2*r_wire+dxdy_wire), 
            rA_core+dxdy_wire+r_wire+j*(2*r_wire+dxdy_wire)
        ])
        geo.AddCircle(
            c=(i*(2*r_wire+dxdy_wire), 
               rA_core+dxdy_wire+r_wire+j*(2*r_wire+dxdy_wire)),
            r=r_wire, bc="inner",
            leftdomain=2,
            rightdomain=1)
pts = np.array(pts)
geo.AddRectangle(p1=(-l_core/2,rI_core),
                 p2=(l_core/2,rA_core),
                 bc="inner",
                 leftdomain=3,
                 rightdomain=1)
geo.SetMaterial (1, "air")
geo.SetMaterial (2, "copper")
geo.SetMaterial (3, "core")
geo.SetDomainMaxH(3,(rA_core-rI_core)/2)

mesh = Mesh(geo.GenerateMesh(maxh=0.0025))

### Listing 7 : Definition der vom Gebiet abhängigen Materialparameter.
omega = 2*pi*50
mu0 = 4*pi*1e-7

sigmaCu = 56e6
sigmaCore = 56e6 #Kupfer
sigma = CoefficientFunction([0, sigmaCu, sigmaCore])
mur = CoefficientFunction([1,1,1])
nu = CoefficientFunction(1/mu0*1/mur)
```

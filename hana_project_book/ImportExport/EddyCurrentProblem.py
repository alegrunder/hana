from ngsolve import *
from ngsolve.webgui import Draw

class EddyCurrentProblem:
    def __init__(self, mymesh, Jimp_coeff, omega_val, sigma_val, nu_val):
        """
        Initializes the Eddy Current Problem.

        Args:
            mymesh: The NGSolve mesh object.
            Jimp_coeff: CoefficientFunction or GridFunction for the impressed current density.
            omega_val: Angular frequency (omega).
            sigma_val: Electrical conductivity (sigma), can be a CoefficientFunction.
            nu_val: Magnetic reluctivity (1/mu), can be a CoefficientFunction.
        """
        self.mesh = mymesh
        self.Jimp = Jimp_coeff
        self.omega = omega_val
        self.sigma = sigma_val
        self.nu = nu_val

        # These will be defined in setup_problem
        self.V = None
        self.u = None
        self.v = None
        self.a = None
        self.f = None
        self.gfu = None

        self._setup_fem_space()
        self._define_forms()

    def _setup_fem_space(self):
        """Defines the Finite Element space and GridFunction."""
        self.V = H1(self.mesh, order=3, complex=True, dirichlet='rotsym')
        self.u, self.v = self.V.TnT()
        self.gfu = GridFunction(self.V, name="A_phi") # Added a name for clarity

    def _define_forms(self):
        """Defines the bilinear and linear forms for the eddy current problem."""
        # NGSolve's symbolic coordinates for cylindrical systems
        # These are available in the context of integrators like dx
        r = y
        z = x

        # Gradients of trial and test functions
        uz, ur = grad(self.u) # grad(u) = (du/dx, du/dy) -> (du/dz, du/dr)
        vz, vr = grad(self.v)

        # Bilinear form
        self.a = BilinearForm(self.V)
        # Term 1: j*omega*sigma*u*v (only in conductive regions)
        self.a += (
            1j
            * self.omega
            * self.sigma
            * self.u
            * self.v
            * dx(definedon=self.mesh.Materials('copper|core'))
        )
        # Term 2
        self.a += (
            self.nu * (1 / r * self.u + ur) * vr + self.nu * uz * vz
        ) * dx

        # Linear form
        self.f = LinearForm(self.V)
        self.f += self.Jimp * self.v * dx(definedon=self.mesh.Materials('copper'))

    def solve(self):
        """Assembles the forms and solves the system."""
        if self.a is None or self.f is None or self.gfu is None:
            print("Problem not fully set up. Call setup methods first.")
            return None

        print("Assembling...")
        self.a.Assemble()
        self.f.Assemble()

        print("Solving system...")
        # Solve A * gfu.vec = f.vec for FreeDofs
        self.gfu.vec.data = self.a.mat.Inverse(
            self.V.FreeDofs()
        ) * self.f.vec
        print("System solved.")
        return self.gfu

    def draw_solution(self, draw_mesh=True, **kwargs):
        """Draws the solution (gfu)."""
        if self.gfu is None or self.gfu.vec.FV().Norm() == 0: # Check if solution exists
            print("No solution to draw or solution is zero.")
            return
        print("Drawing A_phi (Norm)")
        if draw_mesh:
            Draw(Norm(self.gfu), self.mesh, **kwargs)
        else:
            Draw(Norm(self.gfu), **kwargs)

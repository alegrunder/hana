from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.tri import Triangulation

class EddyCurrentProblem:
    def __init__(self, mymesh, Jimp_coeff, omega_val, sigma_val, nu_val):
        """
        Initializes the Eddy Current Problem.

        Args:
            mymesh: The NGSolve mesh object.
            Jimp_coeff: CoefficientFunction or GridFunction for the impressed current density (J_phi).
            omega_val: Angular frequency (omega).
            sigma_val: Electrical conductivity (sigma), can be a CoefficientFunction.
            nu_val: Magnetic reluctivity (1/mu), can be a CoefficientFunction.
        """
        self.mesh = mymesh
        self.Jimp = Jimp_coeff
        self.omega = omega_val
        self.sigma = sigma_val
        self.nu = nu_val

        self.V = None
        self.u = None
        self.v = None
        self.a = None
        self.f = None
        self.gfu = None # This will store A_phi

        # For B-field calculations and visualization
        self.B_r_coeff = None
        self.B_z_coeff = None
        self.B_r_gf = None
        self.B_z_gf = None
        self.time_param = Parameter(0.0) # Symbolic time for time-harmonic fields
        self.B_vector_cf_timed = None    # Time-dependent B-field CoefficientFunction

        self._setup_fem_space()
        self._define_forms()

    def _setup_fem_space(self):
        """Defines the Finite Element space and GridFunction for A_phi."""
        self.V = H1(self.mesh, order=3, complex=True, dirichlet='rotsym')
        self.u, self.v = self.V.TnT()
        self.gfu = GridFunction(self.V, name="A_phi")

    def _define_forms(self):
        """Defines the bilinear and linear forms for the eddy current problem."""
        r = y  # NGSolve's y-coordinate is radial (r)
        z = x  # NGSolve's x-coordinate is axial (z)

        uz, ur = grad(self.u) # du/dz, du/dr
        vz, vr = grad(self.v)

        self.a = BilinearForm(self.V)
        self.a += (
            1j * self.omega * self.sigma * self.u * self.v
            * dx(definedon=self.mesh.Materials('copper|core'))
        )
        self.a += (
            self.nu * (1 / r * self.u + ur) * vr + self.nu * uz * vz
        ) * dx

        self.f = LinearForm(self.V)
        self.f += self.Jimp * self.v * dx(definedon=self.mesh.Materials('copper'))

    def solve(self):
        """Assembles the forms and solves the system for A_phi."""
        if self.a is None or self.f is None or self.gfu is None:
            print("Error: Problem not fully set up.")
            return None

        print("Assembling BilinearForm 'a'...")
        self.a.Assemble()
        print("Assembling LinearForm 'f'...")
        self.f.Assemble()

        print("Solving system for A_phi...")
        self.gfu.vec.data = self.a.mat.Inverse(
            self.V.FreeDofs()
        ) * self.f.vec
        print("System solved.")
        return self.gfu

    def draw_solution(self, **kwargs):
        """Draws the solution A_phi."""
        if self.gfu is None or self.gfu.vec.FV().Norm() == 0:
            print("No A_phi solution to draw or solution is zero.")
            return
        print("Drawing |A_phi| (Magnitude of A_phi)...")
        Draw(Norm(self.gfu), self.mesh, name="A_phi (complex)", **kwargs)

    def calculate_b_field(self, epsZero=1e-8):
        """
        Calculates B_r and B_z from A_phi (self.gfu).
        B_r = -dA_phi/dz
        B_z = (1/r) * d(r*A_phi)/dr
        """
        if self.gfu is None or self.gfu.vec.FV().Norm() == 0:
            print("A_phi solution (gfu) is not available or is zero. Cannot calculate B-field.")
            return False

        print("Calculating B-field components...")
        r = y
        z = x # NGSolve's x-coordinate is axial (z)
        
        grad_gfu = grad(self.gfu)
        self.B_r_coeff = -grad_gfu[0]  # Negative derivative with respect to z (x-coordinate)
        rAphi = GridFunction(self.V)
        rAphi.Set(r * self.gfu)
        self.B_z_coeff = grad(rAphi)[1] / r
        self.B_z_coeff = IfPos(r-epsZero, self.B_z_coeff, 0) # Handle potential division by zero at r=0 axis

        # Create GridFunctions for visualization
        self.B_r_gf = GridFunction(self.V, name="B_r")
        self.B_z_gf = GridFunction(self.V, name="B_z")
        self.B_r_gf.Set(self.B_r_coeff)
        self.B_z_gf.Set(self.B_z_coeff)
        
        print("B-field components calculated.")
        return True

    def draw_b_field_components(self, **kwargs):
        """Draws B_r, B_z, and the norm of B."""
        if self.B_r_gf is None or self.B_z_gf is None:
            print("B-field components not calculated. Call calculate_b_field() first.")
            return

        print("Drawing B_z (complex)...")
        Draw(self.B_z_gf, self.mesh, "B_z (complex)", **kwargs)
        
        print("Drawing B_r (complex)...")
        Draw(self.B_r_gf, self.mesh, "B_r (complex)", **kwargs)

        B_cf = CoefficientFunction((self.B_z_gf, self.B_r_gf))
        print("Drawing |B| (Magnitude of B-field)...")
        Draw(Norm(B_cf), self.mesh, "|B| (Magnitude)", **kwargs)

    def draw_b_vector_field(self, time_val=0.0, **kwargs):
        """
        Draws the B-field (Re(B_z), Re(B_r)) as a vector field at a given time.
        Note: In NGSolve's 2D plots, the first component of the vector CF
              is along the x-axis (our z-axis), and the second is along
              the y-axis (our r-axis). So we plot (B_z, B_r).
        """
        if self.B_r_coeff is None or self.B_z_coeff is None:
            print("B-field components not calculated. Call calculate_b_field() first.")
            return

        self.time_param.Set(time_val) # Set the symbolic time parameter

        # Define the time-dependent complex B-field vector if not already defined
        if self.B_vector_cf_timed is None:
            self.B_vector_cf_timed = CoefficientFunction((
                self.B_z_gf * exp(1j * self.omega * self.time_param), # Component along z-axis (NGSolve x)
                self.B_r_gf * exp(1j * self.omega * self.time_param)  # Component along r-axis (NGSolve y)
            ))
        
        print(f"Drawing B-vector field (real part) at t = {time_val*1000} ms ...")
        Draw(self.B_vector_cf_timed.real, self.mesh,
             f"B-field (Re) at t={time_val:.2e}", vectors=True, **kwargs)

    def showAnimation(self, store_ani=False):
        # --- Mesh coordinates ---
        mesh_coords = [(v.point[0], v.point[1]) for v in self.mesh.vertices]
        x_vals = np.array([p[0] for p in mesh_coords])
        y_vals = np.array([p[1] for p in mesh_coords])

        # --- Triangulation ---
        triangles = []
        for el in self.mesh.Elements(VOL):
            if len(el.vertices) == 3:
                triangles.append([v.nr for v in el.vertices])
        triang = Triangulation(x_vals, y_vals, np.array(triangles))

        # --- Mask vertices in copper regions ---
        is_in_copper = np.zeros(len(self.mesh.vertices), dtype=bool)
        for el in self.mesh.Elements(VOL):
            if el.mat == "copper":
                for v in el.vertices:
                    is_in_copper[v.nr] = True
        # --- Additional mask: quiver hide box around wires ---
        # Compute bounding box of the copper wires
        r_wire=0.001
        x_min, y_min = [-0.002,   0.0095]  #pts.min(axis=0) - 2*r_wire
        x_max, y_max = [0.0245, 0.0185] #pts.max(axis=0) + 2*r_wire
        is_in_quiver_mask = np.zeros(len(self.mesh.vertices), dtype=bool)

        for i, (x, y) in enumerate(mesh_coords):
            # Already exclude copper region
            if is_in_copper[i]:
                is_in_quiver_mask[i] = True
            # Also exclude if inside bounding box of wire grid
            elif x_min <= x <= x_max and y_min <= y <= y_max:
                is_in_quiver_mask[i] = True

        # --- Prepare B-field storage ---
        Bx_all, By_all, Bmag_all = [], [], []
        T = 1 / (self.omega / (2 * pi))
        times = np.linspace(0, T, 100)

        # --- Grid function ---
        V = VectorH1(self.mesh, order=3)
        B_gf = GridFunction(V)

        # --- Evaluate B field over time ---
        for t_val in times:
            self.time_param.Set(t_val)
            B_gf.Set(self.B_vector_cf_timed.real)

            Bx_frame, By_frame, Bmag_frame = [], [], []
            for i, (x, y) in enumerate(mesh_coords):
                try:
                    val = B_gf(x, y)
                    Bx, By = float(val[0]), float(val[1])
                    if is_in_quiver_mask[i]:
                        Bx_frame.append(0)
                        By_frame.append(0)
                    else:
                        Bx_frame.append(Bx)
                        By_frame.append(By)
                    Bmag_frame.append(np.sqrt(Bx**2 + By**2))
                except:
                    Bx_frame.append(0)
                    By_frame.append(0)
                    Bmag_frame.append(0)

            Bx_all.append(Bx_frame)
            By_all.append(By_frame)
            Bmag_all.append(Bmag_frame)

        # --- Plotting ---
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.set_aspect('equal')
        ax.set_xlim(min(x_vals), max(x_vals))
        ax.set_ylim(min(y_vals), max(y_vals))

        # --- Draw material boundaries ---
        edges = {}
        for el in self.mesh.Elements(VOL):
            material = el.mat
            for i in range(len(el.vertices)):
                v1 = el.vertices[i]
                v2 = el.vertices[(i+1) % len(el.vertices)]
                edge = tuple(sorted([v1.nr, v2.nr]))
                if edge not in edges:
                    edges[edge] = []
                edges[edge].append(material)

        boundary_lines = []
        for edge, materials in edges.items():
            if len(set(materials)) > 1:
                v1, v2 = edge
                p1, p2 = mesh_coords[v1], mesh_coords[v2]
                boundary_lines.append([p1, p2])

        line_collection = LineCollection(boundary_lines, linewidths=1, colors='black')
        ax.add_collection(line_collection)

        # --- Background field magnitude ---
        cmap = plt.cm.turbo
        vmin = np.min(Bmag_all)
        vmax = np.max(Bmag_all)
        bg = ax.tripcolor(triang, np.array(Bmag_all[0]), shading='gouraud', cmap=cmap, vmin=vmin, vmax=vmax)
        cb = fig.colorbar(bg, ax=ax, label='|Re(B)| magnitude (T)')

        # --- Initial quiver plot (masked)
        Bx0 = np.ma.masked_array(Bx_all[0], mask=is_in_copper)
        By0 = np.ma.masked_array(By_all[0], mask=is_in_copper)
        quiver = ax.quiver(x_vals, y_vals, Bx0, By0, color='black', scale=0.4)

        # --- Animation update ---
        def update_quiver(frame):
            Bx = np.ma.masked_array(Bx_all[frame], mask=is_in_copper)
            By = np.ma.masked_array(By_all[frame], mask=is_in_copper)
            quiver.set_UVC(Bx, By)
            bg.set_array(np.asarray(Bmag_all[frame]))
            ax.set_title(f"B-field (Re(B) shown, no quiver in copper)\nTime: {times[frame]*1000:.1f} ms")
            return quiver, bg

        ani = animation.FuncAnimation(fig, update_quiver, frames=len(times), interval=50)
        plt.tight_layout()
        plt.show()
        if store_ani:
            writer_gif = PillowWriter(fps=20)
            ani.save("../ImportExport/animation.gif", writer=writer_gif)
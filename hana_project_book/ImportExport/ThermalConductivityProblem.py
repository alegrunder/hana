from ngsolve import *
from ngsolve.webgui import Draw
import numpy as np
import matplotlib.pyplot as plt

class ThermalConductivityProblem:
    def __init__(self, mesh, Qe_source, lambda_coeff, dirichlet_boundaries=''):
        """
        Initializes the Thermal Conductivity Problem.

        Args:
            mesh: The NGSolve mesh object.
            Qe_source: CoefficientFunction for the heat source term.
            lambda_coeff: CoefficientFunction for thermal conductivity.
            dirichlet_boundaries: String specifying Dirichlet boundaries,
                                  e.g., 'outer_right|outer_top|outer_left'.
                                  Can be empty if no Dirichlet BCs are set initially.
        """
        self.mesh = mesh
        self.Qe = Qe_source
        self.lam = lambda_coeff
        self.dirichlet_boundaries_str = dirichlet_boundaries
        self.inhomogDirichlet = False

        self.V_T = None      # H1 Finite Element space for Temperature
        self.u_T = None      # Trial function for Temperature
        self.v_T = None      # Test function for Temperature
        self.gf_T = None     # GridFunction for Temperature solution
        self.a_T = None      # BilinearForm for the stiffness matrix
        self.f_T = None      # LinearForm for the load vector

        self.x_min, self.x_max = None, None
        self.y_min, self.y_max = None, None
        self.z_min, self.z_max = None, None # For potential 3D meshes
        self.is_3d = False
        self._calculate_mesh_ranges()

        self._setup_fem_space()
        self._define_base_forms()

    def _calculate_mesh_ranges(self):
        """Calculates and stores the x, y, (and z) ranges of the mesh."""
        if self.mesh.nv > 0:
            v0 = self.mesh.vertices[0].point
            self.x_min, self.x_max = v0[0], v0[0]
            self.y_min, self.y_max = v0[1], v0[1]
            
            if len(v0) == 3:
                self.is_3d = True
                self.z_min, self.z_max = v0[2], v0[2]
            else:
                self.is_3d = False
                self.z_min, self.z_max = float('inf'), float('-inf')


            for v_idx in range(self.mesh.nv):
                p = self.mesh.vertices[v_idx].point
                self.x_min = min(self.x_min, p[0])
                self.x_max = max(self.x_max, p[0])
                self.y_min = min(self.y_min, p[1])
                self.y_max = max(self.y_max, p[1])
                if self.is_3d:
                    self.z_min = min(self.z_min, p[2])
                    self.z_max = max(self.z_max, p[2])
            
            #print(f"Mesh X range: [{self.x_min}, {self.x_max}]")
            #print(f"Mesh Y range: [{self.y_min}, {self.y_max}]")
            if self.is_3d:
                None
                #print(f"Mesh Z range: [{self.z_min}, {self.z_max}]")
        else:
            print("The mesh has no vertices.")
            # Set default ranges to avoid errors if methods are called on empty mesh
            self.x_min, self.x_max = 0,0
            self.y_min, self.y_max = 0,0
            self.z_min, self.z_max = 0,0


    def _setup_fem_space(self, order=4):
        """
        Defines the Finite Element space for Temperature.
        The dirichlet boundaries are taken from self.dirichlet_boundaries_str.
        """
        print(f"Setting up H1 space with Dirichlet boundaries: '{self.dirichlet_boundaries_str}'")
        self.V_T = H1(self.mesh, order=order, dirichlet=self.dirichlet_boundaries_str)
        self.u_T, self.v_T = self.V_T.TnT()
        self.gf_T = GridFunction(self.V_T, name="Temperature")

    def update_dirichlet_boundaries(self, new_dirichlet_boundaries_str, order=4):
        """
        Allows changing the Dirichlet boundaries and re-initializes the FE space
        and base forms.
        """
        self.dirichlet_boundaries_str = new_dirichlet_boundaries_str
        self.inhomogDirichlet = False
        self._setup_fem_space(order=order) # Re-setup space with new BCs
        self._define_base_forms()      # Redefine forms for the new space
        print(f"Dirichlet boundaries updated to: '{self.dirichlet_boundaries_str}'")


    def _define_base_forms(self):
        """
        Defines the base bilinear and linear forms.
        These can be extended by other methods for specific boundary conditions.
        """
        if self.V_T is None:
            print("FE space not set up. Call _setup_fem_space or update_dirichlet_boundaries first.")
            return

        # Bilinear form (stiffness matrix part)
        self.a_T = BilinearForm(self.V_T)
        self.a_T += self.lam * grad(self.u_T) * grad(self.v_T) * dx

        # Linear form (source term part)
        self.f_T = LinearForm(self.V_T)
        if self.Qe is not None:
            self.f_T += self.Qe * self.v_T * dx
        else:
            print("Warning: Qe_source is None. No volumetric heat source added to f_T.")

    def add_neumann_boundary_condition(self, boundary_name, dTdn, lamBorder):
        """
        Adds a Neumann boundary condition to the linear form f_T.
        flux_value_coeff: CoefficientFunction representing the heat flux.
                          Positive dTdn means heat flowing INTO the domain.
                          The term in linear form is integral(lam * dTdn * v_T * ds).
                          If heat flows OUTSIDE, dTdn should be negative.
        """
        if self.f_T is None:
            self._define_base_forms() # Ensure f_T is initialized
            if self.f_T is None: # Still None, something is wrong
                print("Error: LinearForm f_T could not be initialized.")
                return

        print(f"Adding Neumann BC on '{boundary_name}' with given dT/dn")
        # The term is integral(q_n * v * ds) where q_n is the prescribed flux.
        # If flux_value_coeff is heat flux INTO domain, it's + lam_air * dTdn * v_T * ds
        self.f_T += lamBorder * dTdn * self.v_T * ds(definedon=self.mesh.Boundaries(boundary_name))

    def add_robin_boundary_condition(self, boundary_name, alpha_coeff, lamBorder):
        """
        Adds a Robin (convective) boundary condition.
        dT/dn = - alpha * T
        Weak form contribution: -integral(-alpha_coeff*u_T*v_T*ds)
        So, add to a_T: integral(alpha_coeff * u_T * v_T * ds)

        Args:
            boundary_name: Name of the boundary.
            alpha_coeff: Heat transfer coefficient (CoefficientFunction).
        """
        if self.a_T is None or self.f_T is None:
            self._define_base_forms()
            if self.a_T is None or self.f_T is None:
                print("Error: BilinearForm a_T or LinearForm f_T could not be initialized.")
                return

        print(f"Adding Robin BC on '{boundary_name}' with given alpha")
        self.a_T += lamBorder * alpha_coeff * self.u_T * self.v_T * ds(definedon=self.mesh.Boundaries(boundary_name))

    def set_dirichlet_values(self, boundary_name, value_coeff):
        """
        Sets Dirichlet values on a specified boundary.
        Note: The boundary must have been included in 'dirichlet_boundaries_str'
              during FE space setup for this to take effect directly on gfu.Set.
        Args:
            boundary_name: Name of the boundary.
            value_coeff: CoefficientFunction or float for the temperature value.
        """
        if self.gf_T is None:
            print("GridFunction gf_T not initialized.")
            return
        print(f"Setting Dirichlet value on '{boundary_name}'.")
        self.gf_T.Set(0)
        self.gf_T.Set(value_coeff, definedon=self.mesh.Boundaries(boundary_name))
        self.inhomogDirichlet = True

    def solve(self):
        """Assembles the forms and solves the thermal system."""
        if self.a_T is None or self.f_T is None or self.gf_T is None:
            print("Error: Problem forms or GridFunction not fully set up.")
            return None

        print("Assembling BilinearForm a_T...")
        self.a_T.Assemble()
        print("Assembling LinearForm f_T...")
        self.f_T.Assemble()

        print("Solving thermal system...")
        if self.inhomogDirichlet:
            # Apply Dirichlet conditions (gf_T might have been set by set_dirichlet_values)
            # The solver handles the FreeDofs correctly based on the dirichlet flags in V_T
            self.gf_T.vec.data += self.a_T.mat.Inverse(
                self.V_T.FreeDofs()
            ) * (self.f_T.vec - self.a_T.mat * self.gf_T.vec) # Standard way to handle inhomogeneous Dirichlet
        else:
            self.gf_T.vec.data = self.a_T.mat.Inverse(self.V_T.FreeDofs()) * self.f_T.vec

        print("Thermal system solved.")
        return self.gf_T

    def draw_temperature_solution(self, **kwargs):
        """Draws the temperature solution gf_T."""
        if self.gf_T is None or self.gf_T.vec.FV().Norm() == 0:
            print("No temperature solution to draw or solution is zero.")
            return
        print("Drawing Temperature...")
        Draw(self.gf_T, self.mesh, name="Temperature", **kwargs)
        Redraw(blocking=True)

    def plot_boundary_conditions(self):
        """
        Plots Dirichlet and Neumann values on the boundaries of a 2D mesh.
        Assumes gf_T has been computed.
        """
        if self.gf_T is None:
            print("Temperature solution (gf_T) not available. Cannot plot boundary conditions.")
            return
        if self.is_3d:
            print("Boundary condition plotting is currently implemented for 2D meshes only.")
            return
        if self.x_min is None: # Mesh ranges not calculated
            print("Mesh ranges not calculated. Cannot plot boundary conditions.")
            return

        print("Plotting boundary conditions...")
        # Ensure gf_T has meaningful data, otherwise grad(gf_T) might be problematic
        if self.gf_T.vec.FV().Norm() < 1e-12: # Check if solution is non-trivial
            print("Warning: Temperature solution gf_T is close to zero. Gradient might be unreliable.")


        xp = np.linspace(self.x_min, self.x_max, 200)
        yp = np.linspace(self.y_min, self.y_max, 200)

        fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(10, 12), sharex=False)
        fig.subplots_adjust(hspace=0.45)

        # For grad(gf_T), it's better to project it to a suitable space if evaluating point-wise
        # However, for plotting along boundaries, direct evaluation might be acceptable
        # but can be noisy. Consider projecting grad(gf_T) to an L2 space for smoother results.
        grad_gfT_cf = grad(self.gf_T)

        # ========== First plot: left/right boundaries (constant x) ==========
        ax2 = ax1.twinx()
        try:
            dir_left_vals = np.array([self.gf_T(self.mesh(self.x_min, y_val)) for y_val in yp]).flatten()
            dir_right_vals = np.array([self.gf_T(self.mesh(self.x_max, y_val)) for y_val in yp]).flatten()
            
            # Normal vector on left boundary is (-1, 0), so -dT/dx = -grad(T)[0]
            # Normal vector on right boundary is (1, 0), so dT/dx = grad(T)[0]
            # Neumann value q_n = -lambda * dT/dn. We plot dT/dn here.
            # For simplicity, we plot components of grad(T) and interpret them.
            # dT/dn on left: -grad(T)[0]
            # dT/dn on right: +grad(T)[0]
            neu_left_vals = np.array([-grad_gfT_cf(self.mesh(self.x_min, y_val))[0] for y_val in yp]).flatten()
            neu_right_vals = np.array([grad_gfT_cf(self.mesh(self.x_max, y_val))[0] for y_val in yp]).flatten()

            ax1.plot(yp, dir_left_vals, '-', color='blue', label='Left Dirichlet (T)')
            ax1.plot(yp, dir_right_vals, '-', color='red', label='Right Dirichlet (T)')
            ax2.plot(yp, neu_left_vals, '--', color='cyan', label='Left dT/dn_x')
            ax2.plot(yp, neu_right_vals, '--', color='magenta', label='Right dT/dn_x')
        except Exception as e:
            print(f"Error evaluating on left/right boundaries: {e}")
            return

        dir_max = max(abs(dir_left_vals).max(), abs(dir_right_vals).max(), 1e-3) * 1.1
        neu_max = max(abs(neu_left_vals).max(), abs(neu_right_vals).max()) * 1.1
        ax1.set_ylim(-dir_max, dir_max)
        ax2.set_ylim(-neu_max, neu_max)
        ax1.set_ylabel("Temperature (Dirichlet)")
        ax2.set_ylabel("dT/dn (Normal Gradient)")
        ax1.set_xlabel("y-coordinate")
        ax1.set_title("Boundary Values on Left & Right Edges")
        ax1.grid(True)
        ax2.grid(False)
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

        # ========== Second plot: bottom/top boundaries (constant y) ==========
        ax4 = ax3.twinx()
        try:
            dir_bottom_vals = np.array([self.gf_T(self.mesh(x_val, self.y_min)) for x_val in xp]).flatten()
            dir_top_vals = np.array([self.gf_T(self.mesh(x_val, self.y_max)) for x_val in xp]).flatten()

            # Normal vector on bottom boundary is (0, -1), so -dT/dy = -grad(T)[1]
            # Normal vector on top boundary is (0, 1), so dT/dy = grad(T)[1]
            # dT/dn on bottom: -grad(T)[1]
            # dT/dn on top: +grad(T)[1]
            neu_bottom_vals = np.array([-grad_gfT_cf(self.mesh(x_val, self.y_min))[1] for x_val in xp]).flatten()
            neu_top_vals = np.array([grad_gfT_cf(self.mesh(x_val, self.y_max))[1] for x_val in xp]).flatten()

            ax3.plot(xp, dir_bottom_vals, '-', color='blue', label='Bottom Dirichlet (T)')
            ax3.plot(xp, dir_top_vals, '-', color='red', label='Top Dirichlet (T)')
            ax4.plot(xp, neu_bottom_vals, '--', color='cyan', label='Bottom dT/dn_y')
            ax4.plot(xp, neu_top_vals, '--', color='magenta', label='Top dT/dn_y')
        except Exception as e:
            print(f"Error evaluating on bottom/top boundaries: {e}")
            return
        dir_max = max(abs(dir_bottom_vals).max(), abs(dir_top_vals).max(), 1e-3) * 1.1
        neu_max = max(abs(neu_bottom_vals).max(), abs(neu_top_vals).max()) * 1.1
        ax3.set_ylim(-dir_max, dir_max)
        ax4.set_ylim(-neu_max, neu_max)
        ax3.set_ylabel("Temperature (Dirichlet)")
        ax4.set_ylabel("dT/dn (Normal Gradient)")
        ax3.set_xlabel("x-coordinate")
        ax3.set_title("Boundary Values on Bottom & Top Edges")
        ax3.grid(True)
        ax4.grid(False)
        lines1, labels1 = ax3.get_legend_handles_labels()
        lines2, labels2 = ax4.get_legend_handles_labels()
        ax3.legend(lines1 + lines2, labels1 + labels2, loc='best')

        plt.tight_layout()
        plt.show()
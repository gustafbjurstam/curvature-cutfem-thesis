#include "FESpace/FESpace.hpp"
#include "FESpace/finiteElement.hpp"
#include "FESpace/paraview.hpp"
#include "common/base_interface.hpp"
#include "common/cut_mesh.hpp"
#include "cutfem.hpp"
#include "num/util.hpp"
#include "problem/baseCutProblem.hpp"
#include "problem/itemVF.hpp"
#include "problem/testFunction.hpp"
#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>
#include <iomanip>

// Types for better readability
using mesh_t        = MeshQuad2; // Standard 2D Quad Mesh
using fun_test_t    = TestFunction<mesh_t>;
using fun_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using interface_t   = Interface<mesh_t>;
using space_t       = GFESpace<mesh_t>;
using cutspace_t    = CutFESpace<mesh_t>;
using cutfem_t      = CutFEM<mesh_t>;
using paraview_t    = Paraview<mesh_t>;
using Rd_t          = mesh_t::Rd;
using lagrange_t    = LagrangeQuad2; // Lagrange polynomials on quads

using namespace globalVariable;

// ==============================================================================
// Pure Double-Precision Geometry Namespaces (No Algoim)
// ==============================================================================

namespace peanut {
    const std::string example = "Peanut";

    const double offset_x = 0.0;
    const double offset_y = 0.0;

    const double X0 = -1.65 + offset_x;
    const double Y0 = -0.7 + offset_y;
    const double DX = 1.65 + 1.6;
    const double DY = 1.4;

    const double h0 = DX / 100;
    const int n_refinements = 12;
    const double tauF[] = {1e-3, 1e-2, 1e-1, 1e0, 1e1};

    const double a = 1.0;
    const double b = 1.02 * a;
    const double c = 0.5;

    inline double sq(double x) { return x * x; }

    double phi(double x, double y) {
        double factor1 = sq(x - a) + sq(y);
        double factor2 = sq(x + a) + sq(y);
        return factor1 * factor2 - (b*b*b*b) + c*x;
    }
    double phi_x(double x, double y) {
        double factor1 = sq(x - a) + sq(y);
        double factor2 = sq(x + a) + sq(y);
        return 2.0*(x - a)*factor2 + factor1*2.0*(x + a) + c;
    }
    double phi_y(double x, double y) {
        return 4.0 * y * (sq(x) + a*a + sq(y));
    }
    double phi_xx(double x, double y) {
        return 4.0 * (sq(x) + sq(a) + sq(y)) + 8.0 * (sq(x) - sq(a));
    }
    double phi_xy(double x, double y) {
        return 8.0 * x * y;
    }
    double phi_yy(double x, double y) {
        return 4.0 * (sq(x) + sq(a) + 3.0 * sq(y));
    }
}

namespace polar {
    inline double r(double x, double y) { return std::hypot(x, y); }
    inline double theta(double x, double y) { return std::atan2(y, x); }
    inline double rx(double x, double y) { return x / r(x,y); }
    inline double ry(double x, double y) { return y / r(x,y); }
    inline double thetax(double x, double y) { double r2 = x*x + y*y; return -y / r2; }
    inline double thetay(double x, double y) { double r2 = x*x + y*y; return x / r2; }
    inline double rxx(double x, double y) { double r_val = r(x,y); return (y*y) / (r_val*r_val*r_val); }
    inline double ryy(double x, double y) { double r_val = r(x,y); return (x*x) / (r_val*r_val*r_val); }
    inline double rxy(double x, double y) { double r_val = r(x,y); return -(x*y) / (r_val*r_val*r_val); }
    inline double thetaxx(double x, double y) { double r2 = x*x + y*y; return (2.0*x*y) / (r2*r2); }
    inline double thetayy(double x, double y) { double r2 = x*x + y*y; return -(2.0*x*y) / (r2*r2); }
    inline double thetaxy(double x, double y) { double r2 = x*x + y*y; return (y*y - x*x) / (r2*r2); }
}

namespace star {
    const std::string example = "Star";

    const double offset_x = 0.00333; // Needed to make sure interface doesn't go through node
    const double offset_y = 0.00333; // Needed to make sure interface doesn't go through node
    const double X0 = -1.4 + offset_x;
    const double Y0 = -1.5 + offset_y;
    const double DX = 1.4 * 2.0;
    const double DY = 3.0;

    const double h0 = DX / 80;
    const int n_refinements = 12;
    const double tauF[] = {1e-4, 1e-3, 1e-2, 1e-1, 1e0};

    const double R0 = 1.0; 
    const double A  = 0.3; 
    const double K  = 7.0; 

    inline double get_R(double t) { return R0 + A * std::cos(K*t); }
    inline double get_dR(double t) { return -A * K * std::sin(K*t); }
    inline double get_ddR(double t) { return -A * K * K * std::cos(K*t); }

    double phi(double x, double y) { return polar::r(x,y) - get_R(polar::theta(x,y)); }
    double phi_x(double x, double y) {
        double t = polar::theta(x,y);
        return polar::rx(x,y) - get_dR(t) * polar::thetax(x,y);
    }
    double phi_y(double x, double y) {
        double t = polar::theta(x,y);
        return polar::ry(x,y) - get_dR(t) * polar::thetay(x,y);
    }
    double phi_xx(double x, double y) {
        double t = polar::theta(x,y);
        return polar::rxx(x,y) - (get_ddR(t) * std::pow(polar::thetax(x,y), 2) + get_dR(t) * polar::thetaxx(x,y));
    }
    double phi_yy(double x, double y) {
        double t = polar::theta(x,y);
        return polar::ryy(x,y) - (get_ddR(t) * std::pow(polar::thetay(x,y), 2) + get_dR(t) * polar::thetayy(x,y));
    }
    double phi_xy(double x, double y) {
        double t = polar::theta(x,y);
        return polar::rxy(x,y) - (get_ddR(t) * polar::thetax(x,y) * polar::thetay(x,y) + get_dR(t) * polar::thetaxy(x,y));
    }
}

namespace amoeba {
    const std::string example = "Amoeba";

    const double offset_x = 0.00333; // Needed to make sure interface doesn't go through node
    const double offset_y = 0.00333; // Needed to make sure interface doesn't go through node
    const double X0 = -1.0 + offset_x;
    const double Y0 = -1.4 + offset_y;
    const double DX = 1.0 + 1.65;
    const double DY = 1.4 + 1.3;

    const double h0 = DX / 80;
    const int n_refinements = 12;
    const double tauF[] = {1e-3, 1e-2, 1e-1, 1e0, 1e1};

    const double R0 = 1.0;
    const double C1 = 0.3;
    const double S2 = 0.2;
    const double C3 = 0.1;

    inline double get_R(double t) { return R0 + C1*std::cos(t) + S2*std::sin(2.0*t) + C3*std::cos(3.0*t); }
    inline double get_dR(double t) { return -C1*std::sin(t) + 2.0*S2*std::cos(2.0*t) - 3.0*C3*std::sin(3.0*t); }
    inline double get_ddR(double t) { return -C1*std::cos(t) - 4.0*S2*std::sin(2.0*t) - 9.0*C3*std::cos(3.0*t); }

    double phi(double x, double y) { return polar::r(x,y) - get_R(polar::theta(x,y)); }
    double phi_x(double x, double y) { double t = polar::theta(x,y); return polar::rx(x,y) - get_dR(t) * polar::thetax(x,y); }
    double phi_y(double x, double y) { double t = polar::theta(x,y); return polar::ry(x,y) - get_dR(t) * polar::thetay(x,y); }
    double phi_xx(double x, double y) { double t = polar::theta(x,y); return polar::rxx(x,y) - (get_ddR(t) * std::pow(polar::thetax(x,y), 2) + get_dR(t) * polar::thetaxx(x,y)); }
    double phi_yy(double x, double y) { double t = polar::theta(x,y); return polar::ryy(x,y) - (get_ddR(t) * std::pow(polar::thetay(x,y), 2) + get_dR(t) * polar::thetayy(x,y)); }
    double phi_xy(double x, double y) { double t = polar::theta(x,y); return polar::rxy(x,y) - (get_ddR(t) * polar::thetax(x,y) * polar::thetay(x,y) + get_dR(t) * polar::thetaxy(x,y)); }
}

// Select Geometry here:
// #define PEANUT
// #define STAR
#define AMOEBA

#if defined(PEANUT)
    using namespace peanut;
#elif defined(STAR)
    using namespace star;
#elif defined(AMOEBA)
    using namespace amoeba;
#else
    #error "No example defined"
#endif

// ==============================================================================
// Analytical Curvature & Level Set Wrappers
// ==============================================================================
namespace common {
    double curvature(double x, double y) {
        double phix = phi_x(x,y);
        double phiy = phi_y(x,y);
        double phixx = phi_xx(x,y);
        double phixy = phi_xy(x,y);
        double phiyy = phi_yy(x,y);

        double kappa = phiy*phiy*phixx - 2.0*phix*phiy*phixy + phix*phix*phiyy;
        kappa /= std::pow(phix*phix + phiy*phiy, 1.5);

        return kappa;
    }

    R2 curvature_vector(double x, double y) {
        double phix = phi_x(x,y);
        double phiy = phi_y(x,y);
        double r    = std::hypot(phix, phiy);
        double kappa = curvature(x, y);
        return R2(-kappa * phix / r, -kappa * phiy / r);
    }

    double fun_levelset(double* P, int comp) {
        return phi(P[0], P[1]);
    }

    double fun_curvature(double* P, int direction) {
        R2 H = curvature_vector(P[0], P[1]);
        if (direction == 0) return H[0];
        if (direction == 1) return H[1];
        
        std::cout << "Error, no valid curvature direction specified";
        return std::numeric_limits<double>::max();
    }
}
using namespace common;

// ==============================================================================
// Solver Class
// ==============================================================================
class Curvature {
public:
    Curvature(const space_t &V, const interface_t &interface)
        : interface_(interface)
    {
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_ -> createSurfaceMesh(interface_);
        Vh_           = std::make_unique<cutspace_t>(*surface_mesh_, V);
        problem       = std::make_unique<cutfem_t>(*Vh_);
    }

    void assemble(const double face_stab, const double patch_stab){
        constexpr int D_mesh = mesh_t::Rd::d; // 2
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        Normal n;
        Rnm Id(D_mesh, D_mesh);
        Id = 0.0;
        for (int i = 0; i < D_mesh; ++i){
            Id(i,i) = 1.0;
        }

        const double h = surface_mesh_->Th.get_mesh_size();

        // B_h, unstabilized left hand side
        problem->addBilinear(+ innerProduct(H,v), interface_);

        // Face stabilization
        if (face_stab > 0.0){
            problem->addFaceStabilization(
                face_stab * innerProduct(jump(grad(H)*n), jump(grad(v)*n)),
                *surface_mesh_);
        }

        // Patch stabilization
        if (patch_stab > 0.0){
            problem->addPatchStabilization(
                patch_stab * pow(h, -3) * innerProduct(jump(H),jump(v)), 
                *surface_mesh_);
        }

        // RHS
        problem->addLinear( - contractProduct(Id, gradS(v)), interface_);
    }

    fun_t solve() {
        problem->solve();
        data = problem->rhs_;
        return fun_t(*Vh_, data);
    }

    const activemesh_t &mesh() const { return *surface_mesh_; }

private:
    const interface_t &interface_;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::vector<double> data;
    std::unique_ptr<cutfem_t> problem;
};


// ==============================================================================
// Main Loop
// ==============================================================================
int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    const std::string base_output_path = "../output_files/curvature_2d_standard/" + example + "/";
    std::string path_output_data(base_output_path + "data/");
    std::string path_output_figures(base_output_path + "paraview/");
    
    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
        std::cout << "Running example: " << example << " (P1 elements, No Algoim)\n";
    }

    double h = h0;
    int stab_cases = sizeof(tauF) / sizeof(tauF[0]);
    std::vector<std::vector<double>> error(n_refinements, std::vector<double>(stab_cases,0.0));
    std::vector<double> hs(n_refinements, 0.0);

    for (int i = 0; i < n_refinements; ++i) {
        hs[i] = h;
        const int nx = static_cast<int>(DX / h);
        const int ny = static_cast<int>(DY / h);
        mesh_t Th(nx, ny, X0, Y0, DX, DY);
        
        // FEM spaces:
        lagrange_t FEcurv(1);
        space_t CurvV(Th, FEcurv);
        space_t Lh(Th, DataFE<mesh_t>::P1); // Linear interpolation of the level set

        // Geometry
        fun_t levelset(Lh, fun_levelset);
        InterfaceLevelSet<mesh_t> interface(Th, levelset);

        for (int j = 0; j < stab_cases; ++j){
            // Solver setup
            Curvature curvature(CurvV, interface);
            
            // Assemble: pass tau for either patch or face stab depending on your preference
            // e.g. assemble(face_stab, patch_stab)
            curvature.assemble(0.0, tauF[j]); 
            
            fun_t H = curvature.solve(); // Mean curvature vector

            error[i][j] = L2normSurf(H, fun_curvature, interface, 0, 2);
            // Is this the right one?

            // Fetch exact curvature for Paraview
            cutspace_t CutCurvV(curvature.mesh(), CurvV);
            fun_t curvature_exact(CutCurvV, fun_curvature);

            // Reconstruct exact P0 normals for Paraview
            auto dphix = dx(levelset.expr(0));
            auto dphiy = dy(levelset.expr(0));
            auto mag = sqrt((dphix * dphix) + (dphiy * dphiy));
            auto n_x = dphix / mag;
            auto n_y = dphiy / mag;

            if (MPIcf::IamMaster()) {
                // Output curvature on the induced mesh
                paraview_t para_curv(curvature.mesh(), path_output_figures + "curvature_" + std::to_string(i+1)+ ".vtk");
                para_curv.add(H, "H", 0, 2);
                para_curv.add(curvature_exact, "H_exact", 0, 2);
                para_curv.add(levelset, "levelset", 0, 1);
                para_curv.add(n_x, "n_x");
                para_curv.add(n_y, "n_y");
            }
        }
        std::cout << "h = " << h << ", error[" << i << "] = " << error[i] <<"\n";
        h *= sqrt(0.5); 
    }

    // Output Summary
    if (MPIcf::IamMaster()) {
        std::cout << std::setprecision(16);

        std::cout << "error = [\n";
        for (size_t i = 0; i < n_refinements; ++i) {
            std::cout << "  ["; 
            for (size_t j = 0; j < stab_cases; ++j) {
                std::cout << error[i][j] << (j < stab_cases - 1 ? ", " : "");
            }
            std::cout << "]" << (i < n_refinements - 1 ? ";\n" : "\n");
        }
        std::cout << "];\n";
        
        std::cout << "---------------------\n";

        std::cout << "h = [";
        for (size_t i = 0; i < n_refinements; ++i) {
            std::cout << hs.at(i) << (i < n_refinements - 1 ? ", " : "");
        }
        std::cout << "];\n";
        std::cout << "---------------------\n";
    }

    return 0;
}
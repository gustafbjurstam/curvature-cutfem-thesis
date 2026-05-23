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



// Types for better readability

using mesh_t        = Mesh3;
using fun_test_t    = TestFunction<mesh_t>;
using fun_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using interface_t   = Interface<mesh_t>;
using space_t       = GFESpace<mesh_t>;
using cutspace_t    = CutFESpace<mesh_t>;
using cutfem_t      = CutFEM<mesh_t>;
using paraview_t    = Paraview<mesh_t>;
using Rd_t          = mesh_t::Rd;
using lagrange_t    = Lagrange3; // Lagrange polynomials

using namespace globalVariable;


// Begin with an ellipsoid

namespace Ellipsoid {
    const std::string example = "Ellipsoid";

    // axes x, y, z
    const double a = 1.0;
    const double b = 0.7;
    const double c = 0.92;

    //  Background domain
    const double X0 = -1.1 * a;
    const double Y0 = -1.1 * b;
    const double Z0 = -1.1 * c;
    const double DX =  2.2 * a;
    const double DY =  2.2 * b;
    const double DZ =  2.2 * c;

    // Initial mesh size (dependent on memory availability)
    const double h0 = DX / 20;
    const int n_refinements = 6;
    // const double tauF[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}; // patch
    const double tauF[] = {1e-2, 5e-2, 1e-1, 5e-1, 2e0}; // face

    double fun_levelset(double* P, int i) {
        // Normalized coordinates
        const double x = P[0] / a, y = P[1] / b, z = P[2] / c;
        return sqrt(x*x + y*y + z*z) - 1;
    }

    double fun_curvature(double* P, int direction) {
        const double x = P[0], y = P[1], z = P[2];

        // nabla rho
        const double gx = 2 * x / (a*a);
        const double gy = 2 * y / (b*b);
        const double gz = 2 * z / (c*c);

        const double g_norm = sqrt(gx*gx + gy*gy + gz*gz);

        const double lap = 2 / (a*a) + 2 / (b*b) + 2 / (c*c);
        const double gHg = 2 * (gx/a) * (gx/a) + 2 * (gy/b) * (gy/b) + 2 * (gz/c) * (gz/c);

        const double H_s = (gHg/(g_norm*g_norm) - lap) / g_norm;

        if (direction == 0) {
            return H_s * gx / g_norm;
        } else if (direction == 1) {
            return H_s * gy / g_norm;
        } else if (direction == 2) {
            return H_s * gz / g_norm;
        }

        std::cout << "Error, no valid curvature direction specified";
        return std::numeric_limits<double>::max();
    }
    
} // namespace Ellipsoid

namespace std_tor {
    const std::string example = "Torus";
    const double R = 1.0; // Standard torus
    const double a = 0.5; // Standard torus
    const double b = 0.5; // Standard torus
    const double tauF[] = {1.0e-2, 5.0e-2, 1.0e-1, 5.0e-1, 2.0e0}; // Standard torus
}
namespace ell_tor {
    const std::string example = "Elliptic_Torus";
    const double R = 1.0; // Flattened torus
    const double a = 0.4; // Flattened torus
    const double b = 0.3; // Flattened torus
    const double tauF[] = {1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2, 2.0e-1}; // Flattened torus
}

namespace Torus {
    // const std::string example = "Torus";
    // const std::string example = "Elliptic_Torus";
    using namespace std_tor;
    // using namespace  ell_tor;

    // Parameters
    // const double R = 1.0; // Standard torus
    // const double a = 0.5; // Standard torus
    // const double b = 0.5; // Standard torus

    const double R = 1.0; // Flattened torus
    const double a = 0.4; // Flattened torus
    const double b = 0.3; // Flattened torus

    // Domain (Matching Julia: Lx = (R+a)*1.3, Lz = b*2.5)
    const double Lx = (R + a) * 1.3;
    const double Lz = b * 2.5;

    const double X0 = -Lx;
    const double Y0 = -Lx;
    const double Z0 = -Lz;
    const double DX = 2.0 * Lx;
    const double DY = 2.0 * Lx;
    const double DZ = 2.0 * Lz;

    // Initial mesh size (dependent on memory availability)
    const double h0 = DX / 15; 
    const int n_refinements = 8; 
    // const double tauF[] = {1.0e-2, 5.0e-2, 1.0e-1, 5.0e-1, 2.0e0}; // Standard torus
    const double tauF[] = {1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2, 2.0e-1}; // Flattened torus

    double fun_levelset(double* P, int i) {
        const double x = P[0], y = P[1], z = P[2];
        const double r_xy = sqrt(x*x + y*y);
        
        const double term1 = (r_xy - R) / a;
        const double term2 = z / b;
        
        return term1 * term1 + term2 * term2 - 1.0;
    }

    double fun_curvature(double* P, int direction) {
        const double x = P[0], y = P[1], z = P[2];
        const double S = sqrt(x*x + y*y);

        // Safety for center singularity
        if (S < 1e-10) return 0.0;

        // Precompute constants
        const double A = 1.0 / (a*a);
        const double B = 1.0 / (b*b);

        // Gradient
        const double fac = 2.0 * A * (S - R) / S;
        const double gx = fac * x;
        const double gy = fac * y;
        const double gz = 2.0 * B * z;

        const double gn2 = gx*gx + gy*gy + gz*gz;
        const double g_norm = sqrt(gn2);

        // Laplacian
        const double lap = 2.0 * A * (2.0 - R / S) + 2.0 * B;

        // Hessian Quadratic Form: g^T H g
        // Terms derived analytically
        const double Hxx = 2.0 * A * (1.0 - R/S + x*x * R / (S*S*S));
        const double Hyy = 2.0 * A * (1.0 - R/S + y*y * R / (S*S*S));
        const double Hxy = 2.0 * A * (x*y * R / (S*S*S));
        const double Hzz = 2.0 * B;
        // Note: Hxz and Hyz are 0

        const double gHg = gx*(Hxx*gx + Hxy*gy) + gy*(Hxy*gx + Hyy*gy) + gz*(Hzz*gz);

        // Mean Curvature calculation (matching Ellipsoid sign convention)
        // H_s = (gHg / |g|^2 - lap) / |g|
        const double H_s = (gHg / gn2 - lap) / g_norm;

        if (direction == 0) {
            return H_s * gx / g_norm;
        } else if (direction == 1) {
            return H_s * gy / g_norm;
        } else if (direction == 2) {
            return H_s * gz / g_norm;
        }

        std::cout << "Error, no valid curvature direction specified";
        return std::numeric_limits<double>::max();
    }
} // namespace Torus

namespace Decocube {
    const std::string example = "Decocube";

    const double c = 0.86602540378; // sqrt(0.75)
    const double delta = 0.02;

    // Domain
    const double L = 1.4;
    const double X0 = -L;
    const double Y0 = -L;
    const double Z0 = -L;
    const double DX = 2.0 * L;
    const double DY = 2.0 * L;
    const double DZ = 2.0 * L;

    // Initial mesh size (dependent on memory availability)
    const double h0 = DX / 40;
    const int n_refinements = 5;
    // const double tauF[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}; // patch
    const double tauF[] = {2.0e-3, 5.0e-3, 1.0e-2, 5.0e-2, 2.0e-1}; // face


    // Helper to compute Term and its derivatives: T = ((u^2 + v^2 - c^2)^2 + (w^2 - 1)^2)
    struct TermProps {
        double val;
        double du, dv, dw;
        double duu, dvv, dww, duv; // mixed terms involving w are 0
    };

    TermProps compute_term(double u, double v, double w) {
        TermProps p;
        const double r2 = u*u + v*v;
        const double c2 = c*c;
        const double base1 = r2 - c2;
        const double base2 = w*w - 1.0;

        p.val = base1*base1 + base2*base2;

        // First derivs
        p.du = 4.0 * base1 * u;
        p.dv = 4.0 * base1 * v;
        p.dw = 4.0 * base2 * w;

        // Second derivs
        p.duu = 4.0 * base1 + 8.0 * u * u;
        p.dvv = 4.0 * base1 + 8.0 * v * v;
        p.dww = 4.0 * base2 + 8.0 * w * w;
        p.duv = 8.0 * u * v;

        return p;
    }

    double fun_levelset(double* P, int i) {
        const double x = P[0], y = P[1], z = P[2];
        
        // T1: u=x, v=y, w=z
        double t1 = (pow(x*x + y*y - c*c, 2) + pow(z*z - 1.0, 2));
        // T2: u=y, v=z, w=x
        double t2 = (pow(y*y + z*z - c*c, 2) + pow(x*x - 1.0, 2));
        // T3: u=z, v=x, w=y
        double t3 = (pow(z*z + x*x - c*c, 2) + pow(y*y - 1.0, 2));

        return t1 * t2 * t3 - delta;
    }

    double fun_curvature(double* P, int direction) {
        const double x = P[0], y = P[1], z = P[2];

        // 1. Compute properties of the three factor terms
        // Note careful mapping of (x,y,z) to (u,v,w) for each term
        TermProps T1 = compute_term(x, y, z); // T1(x, y, z)
        TermProps T2 = compute_term(y, z, x); // T2(y, z, x) -> u=y, v=z, w=x
        TermProps T3 = compute_term(z, x, y); // T3(z, x, y) -> u=z, v=x, w=y

        // 2. Gradient of L = T1*T2*T3
        // dL/dx = (dT1/dx)*T2*T3 + T1*(dT2/dx)*T3 + T1*T2*(dT3/dx)
        // Mapping derivatives:
        // T1(x,y,z): dx->du
        // T2(y,z,x): dx->dw
        // T3(z,x,y): dx->dv
        double gx = T1.du*T2.val*T3.val + T1.val*T2.dw*T3.val + T1.val*T2.val*T3.dv;
        double gy = T1.dv*T2.val*T3.val + T1.val*T2.du*T3.val + T1.val*T2.val*T3.dw;
        double gz = T1.dw*T2.val*T3.val + T1.val*T2.dv*T3.val + T1.val*T2.val*T3.du;

        const double g_norm = sqrt(gx*gx + gy*gy + gz*gz);
        if (g_norm < 1e-6) return 0.0;

        // 3. Hessian terms (Chain rule on product of 3 terms)
        // H(ABC) = H(A)BC + A H(B)C + AB H(C) + cross terms
        // Cross terms example: 2*( grad(A)*grad(B)^T * C + ... )

        // Helper to compute second derivative component d2L / d(var)2
        auto compute_H_diag = [&](double t1, double t2, double t3, 
                                  double d1, double d2, double d3, 
                                  double dd1, double dd2, double dd3) {
            return dd1*t2*t3 + t1*dd2*t3 + t1*t2*dd3 
                    + 2.0*(d1*d2*t3 + d1*t2*d3 + t1*d2*d3);
        };

        // Hxx
        double Hxx = compute_H_diag(T1.val, T2.val, T3.val, 
                                    T1.du,  T2.dw,  T3.dv, 
                                    T1.duu, T2.dww, T3.dvv);
        // Hyy
        double Hyy = compute_H_diag(T1.val, T2.val, T3.val, 
                                    T1.dv,  T2.du,  T3.dw, 
                                    T1.dvv, T2.duu, T3.dww);
        // Hzz
        double Hzz = compute_H_diag(T1.val, T2.val, T3.val, 
                                    T1.dw,  T2.dv,  T3.du, 
                                    T1.dww, T2.dvv, T3.duu);
        
        // Laplacian
        double lap = Hxx + Hyy + Hzz;

        // Mixed derivatives (d2L / dxdy)
        // Need gradients and mixed terms of individual Ts
        // T1(x,y): duv. T2(y,z): duw=0. T3(z,x): dvw=0.
        
        // This is getting verbose, so we compute gHg explicitly via sum of products
        // gHg = gx*gx*Hxx + gy*gy*Hyy + gz*gz*Hzz + 2*gx*gy*Hxy + 2*gx*gz*Hxz + 2*gy*gz*Hyz
        
        // We need Hxy, Hxz, Hyz.
        // d2L/dxdy = d/dy ( gx )
        // gx = T1_x T2 T3 + T1 T2_x T3 + T1 T2 T3_x
        // d/dy applied to each term...
        
        // Helper for cross derivatives
        auto compute_H_cross = [&](double val1, double val2, double val3,
                                   double d1_a, double d2_a, double d3_a, // d/da (e.g. x)
                                   double d1_b, double d2_b, double d3_b, // d/db (e.g. y)
                                   double dd1_ab, double dd2_ab, double dd3_ab) // d2/dadb
        {
             double term_ii = dd1_ab*val2*val3 + val1*dd2_ab*val3 + val1*val2*dd3_ab;
             double term_ij = d1_a*d2_b*val3 + d1_a*val2*d3_b +
                              d1_b*d2_a*val3 + val1*d2_a*d3_b +
                              d1_b*val2*d3_a + val1*d2_b*d3_a;
             return term_ii + term_ij;
        };

        // Hxy: x -> u,w,v | y -> v,u,w
        double Hxy = compute_H_cross(T1.val, T2.val, T3.val,
                                     T1.du, T2.dw, T3.dv,   // dx
                                     T1.dv, T2.du, T3.dw,   // dy
                                     T1.duv, 0.0,  0.0);    // d2/dxdy matches duv for T1 only

        // Hxz: x -> u,w,v | z -> w,v,u
        double Hxz = compute_H_cross(T1.val, T2.val, T3.val,
                                     T1.du, T2.dw, T3.dv,   // dx
                                     T1.dw, T2.dv, T3.du,   // dz
                                     0.0,   0.0,   T3.duv); // d2/dxdz matches duv for T3 only (v=x, u=z)

        // Hyz: y -> v,u,w | z -> w,v,u
        double Hyz = compute_H_cross(T1.val, T2.val, T3.val,
                                     T1.dv, T2.du, T3.dw,   // dy
                                     T1.dw, T2.dv, T3.du,   // dz
                                     0.0,   T2.duv, 0.0);   // d2/dydz matches duv for T2 only (u=y, v=z)


        const double gHg = gx*(gx*Hxx + gy*Hxy + gz*Hxz) +
                           gy*(gx*Hxy + gy*Hyy + gz*Hyz) +
                           gz*(gx*Hxz + gy*Hyz + gz*Hzz);

        const double H_s = (gHg/(g_norm*g_norm) - lap) / g_norm;

        if (direction == 0) return H_s * gx / g_norm;
        if (direction == 1) return H_s * gy / g_norm;
        if (direction == 2) return H_s * gz / g_norm;

        std::cout << "Error, no valid curvature direction specified";
        return std::numeric_limits<double>::max();
    }
} // namespace Decocube

// #define ELLIPSOID
// #define TORUS
#define DECOCUBE

#if defined(ELLIPSOID)
    using namespace Ellipsoid;
#elif defined(TORUS)
    using namespace Torus;
#elif defined(DECOCUBE)
    using namespace Decocube;
#else
    #error "No example defined"
#endif

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
    void assemble(const double face_stab, const double edge_stab, const double patch_stab){
        constexpr int D_mesh = mesh_t::Rd::d;
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

        // TODO: Edge stabilization
        // This might require quite a lot of additional work

        // TODO: Test patch stabilization on P1
        // Higher needs aditional term
        // Stab Higher order surface PDEs (2020) Sara, Mats
        // We assume that the scaling should be like the face stab
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

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    const std::string base_output_path = "../output_files/curvature/" + example + "/";
    std::string path_output_data(base_output_path + "data/");
    std::string path_output_figures(base_output_path + "paraview/");
    std::ofstream output_data(path_output_data + "output.dat", std::ofstream::out);

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
        std::cout << "Running example: " << example << "\n";
    }

    double h = h0;
    int stab_cases = sizeof(tauF) / sizeof(tauF[0]);
    std::vector<std::vector<double>> error(n_refinements, std::vector<double>(stab_cases,0.0));
    std::vector<double> hs(n_refinements, 0.0);

    for (int i = 0; i < n_refinements; ++i) {
        hs[i] = h;
        const int nx = static_cast<int>(DX / h);
        const int ny = static_cast<int>(DY / h);
        const int nz = static_cast<int>(DZ / h);
        mesh_t Th(nx, ny, nz, X0, Y0, Z0, DX, DY, DZ);
        // FEM spaces:
        lagrange_t FEcurv(1);
        space_t CurvV(Th,FEcurv);
        space_t Lh(Th, DataFE<mesh_t>::P1); // Linear interpolation of the level set

        // Geometry
        fun_t levelset(Lh, fun_levelset);
        InterfaceLevelSet<mesh_t> interface(Th, levelset);
        for (int j = 0; j < stab_cases; ++j){
            // Solver setup
            Curvature curvature(CurvV, interface);
            curvature.assemble(0.0, 0.0, tauF[j]);
            // curvature.assemble(tauF[j], 0.0, 0.0);
            fun_t H = curvature.solve(); // Mean curvature vector

            error[i][j] = L2normSurf(H, fun_curvature, interface, 0 ,2);

            cutspace_t CutCurvV(curvature.mesh(), CurvV);
            fun_t curvature_exact(CutCurvV, fun_curvature);
            auto dphix = dx(levelset.expr());
            auto dphiy = dy(levelset.expr());
            auto dphiz = dz(levelset.expr());
            auto mag = sqrt((dphix * dphix) + (dphiy * dphiy) + (dphiz * dphiz));
            auto nx = dphix / mag;
            auto ny = dphiy / mag;
            auto nz = dphiz / mag;


            if (MPIcf::IamMaster()) {
                // Output curvature on the induced mesh
                paraview_t para_curv(curvature.mesh(), path_output_figures + "curvature_" + std::to_string(i+1)+ ".vtk");
                para_curv.add(H, "H", 0, 3);
                para_curv.add(curvature_exact, "H_exact", 0, 3);
                para_curv.add(levelset, "levelset", 0, 1);
                para_curv.add(nx, "nx");
                para_curv.add(ny, "ny");
                para_curv.add(nz, "nz");
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
            std::cout << "  ["; // Start of row (refinement level i)
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
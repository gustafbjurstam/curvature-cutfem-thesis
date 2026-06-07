#include "FESpace/FESpace.hpp"
#include "FESpace/finiteElement.hpp"
#include "FESpace/paraview.hpp"
#include "common/cut_mesh.hpp"
#include "cutfem.hpp"
#include "num/util.hpp"
#include "problem/itemVF.hpp"
#include "problem/testFunction.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <malloc.h>
#include <memory>
#include <vector>

// function space
#define HEXA_Q1 // SCRIPT_MARKER3
// #define HEXA_Q2
// #define TET_P2

// Geometry. Exactly one should be defined.
#define ELL_TOR // SCRIPT_MARKER1
// #define DECO_CUBE
// #define SPHERE

// Stabilisation. Exactly one should be defined.
#define ELEM_STAB // SCRIPT_MARKER2
// #define ELEM_STAB
// #define PATCH_ONLY

#if defined(ELEM_STAB)
    const bool eleming = true;
#elif defined(PATCH_ONLY)
    const bool eleming = false;
#else
    #error "No stabilisation mode defined"
#endif

#if defined(HEXA_Q1)
    constexpr int fe_order = 1;
    const std::string fe_suffix = "Q1";
    using lagrange_t    = LagrangeHexa3;
    using mesh_t        = MeshHexa;
    const std::string mesh = "Hexa";

    const bool patching = true;
    const bool facing   = false;
#elif defined(HEXA_Q2)
    constexpr int fe_order = 2;
    const std::string fe_suffix = "Q2";
    using lagrange_t    = LagrangeHexa3;
    using mesh_t        = MeshHexa;
    const std::string mesh = "Hexa";

    const bool patching = true;
    const bool facing   = false;
#elif defined(TET_P2)
    constexpr int fe_order = 2;
    const std::string fe_suffix = "P2";
    using lagrange_t    = Lagrange3;
    using mesh_t        = Mesh3;
    const std::string mesh = "Tet";

    const bool patching = false;
    const bool facing   = true;
#else
    #error "No FE order defined"
#endif

using fun_test_t    = TestFunction<mesh_t>;
using fun_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using space_t       = GFESpace<mesh_t>;
using cutspace_t    = CutFESpace<mesh_t>;
using paraview_t    = Paraview<mesh_t>;
using Rd_t          = mesh_t::Rd;
using interface_t   = AlgoimInterface<mesh_t, fun_t>;
using cutfem_t      = AlgoimCutFEMUnified<mesh_t, fun_t>;

using namespace globalVariable;

namespace sphere {
    const std::string example = "Sphere";


    // Computational domain
    const double offset_x = 0.0;
    const double offset_y = 0.0;
    const double offset_z = 0.0;

    const double Lx = 1.3;

    const double X0 = -Lx + offset_x;
    const double Y0 = -Lx + offset_y;
    const double Z0 = -Lx + offset_z;

    const double DX = 2.0 * Lx;
    const double DY = 2.0 * Lx;
    const double DZ = 2.0 * Lx;

    // Mesh settings
    double const h0 = DX / 15.0;
    const int n_refinements = 6; // 18 works for Q1

    // Stabilisation
    #if defined (HEXA_Q1)
    const double tauP[] = {1e-2, 1e-1, 1e0};
    const double tauS[] = {0.0, 1e-3, 1e-2, 1e-1};
    const double tauE[] = {0.0, 1e-3, 1e-2, 1e-1};
    #elif defined(HEXA_Q2) 
    const double tauP[] = {1e1, 1e-1, 1e-2};
    const double tauS[] = {0.0, 1e0, 1e1, 1e2};
    const double tauE[] = {0.0, 1e0, 1e1, 1e2};
    #else
    const double tauP[] = {1e-2, 1e-1, 1e0};
    const double tauS[] = {0.0, 1e-3, 1e-2, 1e-1};
    const double tauE[] = {0.0, 1e-3, 1e-2, 1e-1};
    #endif
    
    
    // Helper functions
    template <typename T> T sq(const T& x) { return x * x; }

    // Templated Level Set Function
    template<typename T>
    T phi(const T& x, const T& y, const T& z){
        using std::sqrt; // Permits ADL for algoim::sqrt
        
        T phi = sqrt(sq(x)+sq(y)+sq(z)) - 1.0;

        return phi;
    }

    // Templated First Derivatives
    template<typename T>
    T phi_x(const T& x, const T& y, const T& z){
        using std::sqrt;
        T factor = phi(x,y,z) +  1.0;

        return x / (factor + 1e-13);
    }

    template<typename T>
    T phi_y(const T& x, const T& y, const T& z){
        using std::sqrt;
        T factor = phi(x,y,z) +  1.0;

        return y / (factor + 1e-13);
    }

    template<typename T>
    T phi_z(const T& x, const T& y, const T& z){
        using std::sqrt;
        T factor = phi(x,y,z) +  1.0;

        return z / (factor + 1e-13);
    }

    // Double-precision Second Derivatives
    double phi_xx(double x, double y, double z){
        double a = y*y + z*z;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        return a / denom;
    }

    double phi_yy(double x, double y, double z){
        double a = x*x + z*z;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        return a / denom;
    }


    double phi_zz(double x, double y, double z){
        double a = y*y + x*x;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        return a / denom;
    }


    double phi_xy(double x, double y, double z){
        double a = - x * y;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        
        return a / denom;
    }

    double phi_xz(double x, double y, double z){
        double a = - x * z;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        
        return a / denom;
    }

    double phi_yz(double x, double y, double z){
        double a = - z * y;
        double denom = pow(phi(x,y,z) + 1.0,3.0);
        
        return a / denom;
    }
}

namespace elliptic_torus {
    const std::string example = "Elliptic_Torus";

    // Geometry
    const double R = 1.0; // Flattened torus major radius
    const double a = 0.4; // Flattened torus minor radius x/y
    const double b = 0.3; // Flattened torus minor radius z

    // Computational domain
    const double offset_x = 0.0;
    const double offset_y = 0.0;
    const double offset_z = 0.0;

    const double Lx = (R + a) * 1.3;
    const double Lz = b * 2.5;

    const double X0 = -Lx + offset_x;
    const double Y0 = -Lx + offset_y;
    const double Z0 = -Lz + offset_z;

    const double DX = 2.0 * Lx;
    const double DY = 2.0 * Lx;
    const double DZ = 2.0 * Lz;

    // Mesh settings
    double const h0 = DX / 15.0;

    // Stabilisation
    #if defined (HEXA_Q1)
    const double tauP[] = {1e-2, 1e-1, 1e0};
    const double tauS[] = {0.0, 1e-3, 1e-2, 1e-1};
    const double tauE[] = {0.0, 1e-3, 1e-2, 1e-1};
    const int n_refinements = 18; // 18 works for Q1
    #elif defined(HEXA_Q2) 
    const double tauP[] = {1e1, 1e-1, 1e-2};
    const double tauS[] = {0.0, 1e0, 1e1, 1e2};
    const double tauE[] = {0.0, 1e0, 1e1, 1e2};
    const int n_refinements = 12; // 18 works for Q1
    #else
    const double tauP[] = {1e-4, 1e-3, 1e-2};
    const double tauS[] = {0.0, 1e-4, 1e-3, 1e-2};
    const double tauE[] = {0.0, 1e-4, 1e-3, 1e-2};
    const int n_refinements = 14; // 18 works for Q1
    #endif
    
    
    // Helper functions
    template <typename T> T sq(const T& x) { return x * x; }

    // Templated Level Set Function
    template<typename T>
    T phi(const T& x, const T& y, const T& z){
        using std::sqrt; // Permits ADL for algoim::sqrt
        
        T r_xy  = sqrt(sq(x) + sq(y));
        T term1 = (r_xy - R) / a;
        T term2 = z / b;

        return sq(term1) + sq(term2) - 1.0;
    }

    // Templated First Derivatives
    template<typename T>
    T phi_x(const T& x, const T& y, const T& z){
        using std::sqrt;
        T r_xy = sqrt(sq(x) + sq(y));
        
        // Note: Assumes evaluation does not occur exactly at the singularity (x=0, y=0)
        T fac = 2.0 / sq(a) * (r_xy - R) / r_xy;
        return fac * x;
    }

    template<typename T>
    T phi_y(const T& x, const T& y, const T& z){
        using std::sqrt;
        T r_xy = sqrt(sq(x) + sq(y));
        
        T fac = 2.0 / sq(a) * (r_xy - R) / r_xy;
        return fac * y;
    }

    template<typename T>
    T phi_z(const T& x, const T& y, const T& z){
        return 2.0 * z / sq(b);
    }

    // Double-precision Second Derivatives
    double phi_xx(double x, double y, double z){
        double r_xy = std::sqrt(x*x + y*y);
        if (r_xy < 1e-10) return 0.0; // Safety for the central singularity
        
        return 2.0 / (a*a) * (1.0 - R/r_xy + (R * x*x) / (r_xy*r_xy*r_xy));
    }

    double phi_yy(double x, double y, double z){
        double r_xy = std::sqrt(x*x + y*y);
        if (r_xy < 1e-10) return 0.0;
        
        return 2.0 / (a*a) * (1.0 - R/r_xy + (R * y*y) / (r_xy*r_xy*r_xy));
    }

    double phi_zz(double x, double y, double z){
        return 2.0 / (b*b);
    }

    double phi_xy(double x, double y, double z){
        double r_xy = std::sqrt(x*x + y*y);
        if (r_xy < 1e-10) return 0.0;
        
        return 2.0 / (a*a) * ((R * x*y) / (r_xy*r_xy*r_xy));
    }

    double phi_xz(double x, double y, double z){
        return 0.0; // Mixed z-plane derivatives vanish for an aligned torus
    }

    double phi_yz(double x, double y, double z){
        return 0.0;
    }
}
namespace decocube {
    const std::string example = "Decocube";

    const double c = 0.86602540378; // sqrt(0.75)
    const double delta = 0.02;

    // Computational Domain
    const double offset_x = 0.0;
    const double offset_y = 0.0;
    const double offset_z = 0.0;

    const double L = 1.4;
    const double X0 = -L + offset_x;
    const double Y0 = -L + offset_y;
    const double Z0 = -L + offset_z;

    const double DX = 2.0 * L;
    const double DY = 2.0 * L;
    const double DZ = 2.0 * L;

    // Stabilisation
    #if defined (HEXA_Q1)
    const double tauP[] = {1e-3, 1e-2, 1e-1};
    const double tauS[] = {0.0, 1e-3, 1e-2, 1e-1};
    const double tauE[] = {0.0, 1e-3, 1e-2, 1e-1};
    // Mesh setting
    const double h0 = DX / (40.0 * sqrt(2)); // Changes n_refinements, needs to be reduced by 2 to have the same final size
    const int n_refinements = 6; // 11 works for Q1
    #elif defined(HEXA_Q2) 
    const double tauP[] = {1e1, 1e-1, 1e-2};
    const double tauS[] = {0.0, 1e0, 1e1, 1e2};
    const double tauE[] = {0.0, 1e0, 1e1, 1e2};
    // Mesh setting
    const double h0 = DX / (40.0 * sqrt(2)); // Changes n_refinements, needs to be reduced by 2 to have the same final size
    const int n_refinements = 6; 
    #else
    const double tauP[] = {1e-5, 1e-4, 1e-3};
    const double tauS[] = {0.0, 1e-6, 1e-5, 1e-4};
    const double tauE[] = {0.0, 1e-6, 1e-5, 1e-4};
    // Mesh setting
    const double h0 = DX / (40.0);
    const int n_refinements = 6; 
    #endif

    // ========================================================================
    // 1. TEMPLATED EVALUATIONS (Compatible with double & algoim::Interval)
    // ========================================================================
    
    template <typename T> T sq(const T& x) { return x * x; }

    template<typename T>
    T phi(const T& x, const T& y, const T& z) {
        T t1 = sq(sq(x) + sq(y) - sq(c)) + sq(sq(z) - 1.0);
        T t2 = sq(sq(y) + sq(z) - sq(c)) + sq(sq(x) - 1.0);
        T t3 = sq(sq(z) + sq(x) - sq(c)) + sq(sq(y) - 1.0);

        return t1 * t2 * t3 - delta;
    }

    template<typename T>
    T phi_x(const T& x, const T& y, const T& z) {
        T t1 = sq(sq(x) + sq(y) - sq(c)) + sq(sq(z) - 1.0);
        T t2 = sq(sq(y) + sq(z) - sq(c)) + sq(sq(x) - 1.0);
        T t3 = sq(sq(z) + sq(x) - sq(c)) + sq(sq(y) - 1.0);

        T dt1 = 4.0 * (sq(x) + sq(y) - sq(c)) * x;
        T dt2 = 4.0 * (sq(x) - 1.0) * x;
        T dt3 = 4.0 * (sq(z) + sq(x) - sq(c)) * x;

        return dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3;
    }

    template<typename T>
    T phi_y(const T& x, const T& y, const T& z) {
        T t1 = sq(sq(x) + sq(y) - sq(c)) + sq(sq(z) - 1.0);
        T t2 = sq(sq(y) + sq(z) - sq(c)) + sq(sq(x) - 1.0);
        T t3 = sq(sq(z) + sq(x) - sq(c)) + sq(sq(y) - 1.0);

        T dt1 = 4.0 * (sq(x) + sq(y) - sq(c)) * y;
        T dt2 = 4.0 * (sq(y) + sq(z) - sq(c)) * y;
        T dt3 = 4.0 * (sq(y) - 1.0) * y;

        return dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3;
    }

    template<typename T>
    T phi_z(const T& x, const T& y, const T& z) {
        T t1 = sq(sq(x) + sq(y) - sq(c)) + sq(sq(z) - 1.0);
        T t2 = sq(sq(y) + sq(z) - sq(c)) + sq(sq(x) - 1.0);
        T t3 = sq(sq(z) + sq(x) - sq(c)) + sq(sq(y) - 1.0);

        T dt1 = 4.0 * (sq(z) - 1.0) * z;
        T dt2 = 4.0 * (sq(y) + sq(z) - sq(c)) * z;
        T dt3 = 4.0 * (sq(z) + sq(x) - sq(c)) * z;

        return dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3;
    }

    // ========================================================================
    // 2. STANDARD DOUBLE EVALUATIONS (For Curvature Hessian)
    // ========================================================================
    
    struct TermProps {
        double val, du, dv, dw, duu, dvv, dww, duv;
    };

    inline TermProps compute_term(double u, double v, double w) {
        TermProps p;
        const double r2 = u*u + v*v;
        const double c2 = c*c;
        const double base1 = r2 - c2;
        const double base2 = w*w - 1.0;

        p.val = base1*base1 + base2*base2;
        p.du  = 4.0 * base1 * u;
        p.dv  = 4.0 * base1 * v;
        p.dw  = 4.0 * base2 * w;
        p.duu = 4.0 * base1 + 8.0 * u * u;
        p.dvv = 4.0 * base1 + 8.0 * v * v;
        p.dww = 4.0 * base2 + 8.0 * w * w;
        p.duv = 8.0 * u * v;

        return p;
    }

    inline double compute_H_diag(double t1, double t2, double t3, 
                                 double d1, double d2, double d3, 
                                 double dd1, double dd2, double dd3) {
        return dd1*t2*t3 + t1*dd2*t3 + t1*t2*dd3 
             + 2.0*(d1*d2*t3 + d1*t2*d3 + t1*d2*d3);
    }

    inline double compute_H_cross(double val1, double val2, double val3,
                                  double d1_a, double d2_a, double d3_a, 
                                  double d1_b, double d2_b, double d3_b, 
                                  double dd1_ab, double dd2_ab, double dd3_ab) {
        double term_ii = dd1_ab*val2*val3 + val1*dd2_ab*val3 + val1*val2*dd3_ab;
        double term_ij = d1_a*d2_b*val3 + d1_a*val2*d3_b +
                         d1_b*d2_a*val3 + val1*d2_a*d3_b +
                         d1_b*val2*d3_a + val1*d2_b*d3_a;
        return term_ii + term_ij;
    }

    // Diagonal Hessian Elements
    double phi_xx(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_diag(T1.val, T2.val, T3.val, T1.du, T2.dw, T3.dv, T1.duu, T2.dww, T3.dvv);
    }

    double phi_yy(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_diag(T1.val, T2.val, T3.val, T1.dv, T2.du, T3.dw, T1.dvv, T2.duu, T3.dww);
    }

    double phi_zz(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_diag(T1.val, T2.val, T3.val, T1.dw, T2.dv, T3.du, T1.dww, T2.dvv, T3.duu);
    }

    // Mixed Hessian Elements
    double phi_xy(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_cross(T1.val, T2.val, T3.val,
                               T1.du, T2.dw, T3.dv,   
                               T1.dv, T2.du, T3.dw,   
                               T1.duv, 0.0,  0.0);    
    }

    double phi_xz(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_cross(T1.val, T2.val, T3.val,
                               T1.du, T2.dw, T3.dv,   
                               T1.dw, T2.dv, T3.du,   
                               0.0,   0.0,   T3.duv); 
    }

    double phi_yz(double x, double y, double z) {
        TermProps T1 = compute_term(x, y, z);
        TermProps T2 = compute_term(y, z, x);
        TermProps T3 = compute_term(z, x, y);
        return compute_H_cross(T1.val, T2.val, T3.val,
                               T1.dv, T2.du, T3.dw,   
                               T1.dw, T2.dv, T3.du,   
                               0.0,   T2.duv, 0.0);   
    }
}
#if defined(ELL_TOR)
    using namespace elliptic_torus;
#elif defined(DECO_CUBE)
    using namespace decocube;
#elif defined(SPHERE)
    using namespace sphere;
#else
    #error "No mesh type defined"
#endif


namespace common {
    // Assuming the library uses R3 for 3D vectors (analogous to R2 in 2D)
    R3 normal(double x, double y, double z) {
        double phix = phi_x(x,y,z);
        double phiy = phi_y(x,y,z);
        double phiz = phi_z(x,y,z);
        double r    = sqrt(phix*phix + phiy*phiy + phiz*phiz);
        
        return R3(phix/r, phiy/r, phiz/r);
    }

    double curvature(double x, double y, double z) {
        // 1. First derivatives
        double phix = phi_x(x,y,z);
        double phiy = phi_y(x,y,z);
        double phiz = phi_z(x,y,z);
        
        // 2. Second derivatives (Diagonal)
        double phixx = phi_xx(x,y,z);
        double phiyy = phi_yy(x,y,z);
        double phizz = phi_zz(x,y,z);
        
        // 3. Second derivatives (Mixed)
        double phixy = phi_xy(x,y,z);
        double phixz = phi_xz(x,y,z);
        double phiyz = phi_yz(x,y,z);

        // 4. Gradient norms
        double gn2 = phix*phix + phiy*phiy + phiz*phiz;
        double g_norm = sqrt(gn2);
        
        // 5. Laplacian (Trace of Hessian)
        double lap = phixx + phiyy + phizz;
        
        // 6. Quadratic form of the Hessian: g^T * H * g
        double gHg = phix * (phixx*phix + phixy*phiy + phixz*phiz) + 
                     phiy * (phixy*phix + phiyy*phiy + phiyz*phiz) + 
                     phiz * (phixz*phix + phiyz*phiy + phizz*phiz);

        // Generalised mean curvature formula: kappa = div(n)
        double kappa = (lap - (gHg / gn2)) / g_norm;
        
        return kappa;
    }

    R3 curvature_vector(double x, double y, double z) {
        return (-curvature(x,y,z)) * normal(x,y,z);
    }

    double fun_levelset(double* P, int i = 0) {
        return phi(P[0], P[1], P[2]);
    }

    double fun_curvature(double* P, int direction) {
        R3 H = curvature_vector(P[0], P[1], P[2]);
        if (direction == 0) return H[0];
        else if (direction == 1) return H[1];
        else if (direction == 2) return H[2];
        else return std::numeric_limits<double>::max();
    }

    double fun_nx(double* P, int domin) {
        return normal(P[0], P[1], P[2])[0];
    }

    double fun_ny(double* P, int domin) {
        return normal(P[0], P[1], P[2])[1];
    }

    double fun_nz(double* P, int domin) {
        return normal(P[0], P[1], P[2])[2];
    }

    double fun_nx_max(double* P, int domin, double t) {
        return fun_nx(P, 0);
    }

    double fun_ny_max(double* P, int domin, double t) {
        return fun_ny(P, 0);
    }

    double fun_nz_max(double* P, int domin, double t) {
        return fun_nz(P, 0);
    }
}

using namespace common;



class CurvatureAlgoim3D {
  public:
    CurvatureAlgoim3D(const space_t &V, const fun_t &phi)
        : phi_h_(phi), interface_(V.Th, phi_h_) {
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_->createSurfaceMesh(interface_);

        Vh_ = std::make_unique<cutspace_t>(*surface_mesh_, V);

        ProblemOption problem_option;
        problem_option.order_space_element_quadrature_ = 8;
        problem_option.algoim_surface_quad_deg_        = 6;
        problem_option.algoim_bernstein_deg_           = 2;
        problem_option.algoim_vol_quad_deg_            = 6;

        problem = std::make_unique<cutfem_t>(*Vh_, phi_h_, problem_option);
    }

    void weak_form() {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        Rnm Id(D_mesh, D_mesh);
        Id = 0.0;
        for (int i = 0; i < D_mesh; ++i) {
            Id(i, i) = 1.0;
        }

        // a_h
        problem->addBilinear(+innerProduct(H, v), interface_);

        // L_h
        problem->addLinear(-contractProduct(Id, gradS(v)), interface_);
    }

    void patch_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();

        problem->addPatchStabilization(
            tau * pow(h, -3) * innerProduct(jump(H), jump(v)),
            *surface_mesh_);
    }

    void face_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();
        Normal n;
        
        auto DnH = grad(H) * n;
        auto Dnv = grad(v) * n;
        #if defined(TET_P2)

        problem->addFaceStabilization(
            tau * innerProduct(jump(DnH), jump(Dnv)),
            *surface_mesh_);
            
        auto Dn2H = grad(DnH) * n;
        auto Dn2v = grad(Dnv) * n;
        problem->addFaceStabilization(
            tau * h*h* innerProduct(jump(Dn2H), jump(Dn2v)),
            *surface_mesh_);
        #endif
    }


    void elem_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();

        auto phi_dx = dx(phi_h_.expr());
        auto phi_dy = dy(phi_h_.expr());
        auto phi_dz = dz(phi_h_.expr());
        auto grad_norm = sqrt(phi_dx * phi_dx + phi_dy * phi_dy + phi_dz * phi_dz);

        auto n_x = phi_dx / grad_norm;
        auto n_y = phi_dy / grad_norm;
        auto n_z = phi_dz / grad_norm;

        auto d_n_H = n_x * dx(H) + n_y * dy(H) + n_z * dz(H);
        auto d_n_v = n_x * dx(v) + n_y * dy(v) + n_z * dz(v);

        problem->BaseFEM::addBilinear(
            tau * pow(h, -1.0) * innerProduct(d_n_H, d_n_v),
            *surface_mesh_);
    }

    fun_t solve() {
        problem->solve();
        data_ = problem->rhs_;
        return fun_t(*Vh_, data_);
    }

    double get_error(const fun_t &approx) {
        return L2_norm_surface(approx, fun_curvature, interface_, phi_h_, 0, 3);
    }

    double normal_error_max() {
        auto phi_dx = dx(phi_h_.expr());
        auto phi_dy = dy(phi_h_.expr());
        auto phi_dz = dz(phi_h_.expr());
        auto grad_norm = sqrt(phi_dx * phi_dx + phi_dy * phi_dy + phi_dz * phi_dz);

        auto nx = phi_dx / grad_norm;
        auto ny = phi_dy / grad_norm;
        auto nz = phi_dz / grad_norm;

        double error1 = max_norm_surface(nx, fun_nx, interface_, phi_h_);
        double error2 = max_norm_surface(ny, fun_ny, interface_, phi_h_);
        double error3 = max_norm_surface(nz, fun_nz, interface_, phi_h_);

        return error1 + error2 + error3;
    }

    void export_global_matrix(const std::string &filename) const {
        matlab::Export(problem->mat_, filename);
    }

    const activemesh_t &mesh() const { return *surface_mesh_; }

  private:
    fun_t phi_h_;
    interface_t interface_;

    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::vector<double> data_;
    std::unique_ptr<cutfem_t> problem;
};

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    std::string example_name = example + "_algoim_" + fe_suffix;
    if (eleming) {
        example_name += "_elem";
    } else {
        example_name += "_patch";
    }

    if (!patching && !facing) {
        example_name += "_only";
    } else if (facing) {
        example_name += "_face";
    }
    example_name += "_" + mesh;

    const std::string base_output_path = "../results/curvature_3d_algoim/" + example_name + "/";
    std::string path_output_data(base_output_path + "data/");
    std::string path_output_figures(base_output_path + "paraview/");

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
        std::cout << "Running example: " << example_name << "\n";
    }

    double h = h0;

    int stab_cases = 1;
    if (eleming) {
        stab_cases = sizeof(tauE) / sizeof(tauE[0]);
    } else if (patching || facing) {
        stab_cases = sizeof(tauP) / sizeof(tauP[0]);
    }

    std::vector<std::vector<double>> error(n_refinements, std::vector<double>(stab_cases, 0.0));
    std::vector<double> hs(n_refinements, 0.0);
    std::vector<double> n_err(n_refinements, 0.0);

    for (int i = 0; i < n_refinements; ++i) {
        hs[i] = h;
        const int nx = static_cast<int>(DX / h);
        const int ny = static_cast<int>(DY / h);
        const int nz = static_cast<int>(DZ / h);
        mesh_t Th(nx, ny, nz, X0, Y0, Z0, DX, DY, DZ);

        bool export_matrix = (i < int(floor(1.5 * n_refinements / 2.0)));

        // Qk curvature approximation on hexahedra.
        lagrange_t FEcurv(fe_order);
        space_t CurvV(Th, FEcurv);

        // Qk interpolation of the level set.
        #if defined(HEXA_Q1)
            space_t Lh(Th, DataFE<mesh_t>::P1);
        #elif defined(HEXA_Q2)
            space_t Lh(Th, DataFE<mesh_t>::P2);
        #elif defined(TET_P2)
            space_t Lh(Th, DataFE<mesh_t>::P2);
        #endif
        fun_t phi_h(Lh, fun_levelset);

        for (int j = 0; j < stab_cases; ++j) {
            {
                CurvatureAlgoim3D curvature(CurvV, phi_h);

                if (eleming && tauE[j] > 1e-15) {
                    curvature.elem_stab(tauE[j]);
                }
                if (patching) {
                    double current_tauP = !eleming ? tauP[j] : tauP[0];
                    if (current_tauP > 1e-15) {
                        curvature.patch_stab(current_tauP);
                    }
                } else if (facing) {
                    double current_tauP = !eleming ? tauP[j] : tauP[0];
                    if (current_tauP > 1e-15) {
                        curvature.face_stab(current_tauP);
                    }
                }

                curvature.weak_form();

                if (MPIcf::IamMaster() && export_matrix) {
                    std::string base_filename = path_output_data + example_name +
                                                "_ref" + std::to_string(i + 1) +
                                                "_stab" + std::to_string(j + 1);
                    curvature.export_global_matrix(base_filename + "_mat.dat");
                }

                fun_t H = curvature.solve();
                error[i][j] = curvature.get_error(H);

                if (j == 0) {
                    n_err[i] = curvature.normal_error_max();
                }
            }
            malloc_trim(0);
        }

        std::cout << "h = " << h << ", error[" << i << "] = " << error[i][0] << "\n";
        h *= sqrt(sqrt(0.5));
    }

    if (MPIcf::IamMaster()) {
        std::ofstream output_data(path_output_data + "output.dat", std::ios::out);
        if (output_data.is_open()) {
            output_data << std::setprecision(16);

            std::string filename_base = example + "_Algoim_3d_" + fe_suffix;
            if (eleming) {
                filename_base += "_elem";
            } else {
                filename_base += "_patch";
            }
            filename_base += "_" + mesh;
            output_data << filename_base << "\n";

            for (int j = 0; j < stab_cases; ++j) {
                double current_tau = eleming ? tauE[j] : tauP[j];
                output_data << current_tau << (j < stab_cases - 1 ? " " : "");
            }
            output_data << "\n";

            for (int i = 0; i < n_refinements; ++i) {
                output_data << hs[i] << " " << n_err[i];
                for (int j = 0; j < stab_cases; ++j) {
                    output_data << " " << error[i][j];
                }
                output_data << "\n";
            }
            output_data.close();
            std::cout << "Data successfully written to: " << path_output_data << "output.dat\n";
        } else {
            std::cerr << "Failed to open output.dat for writing.\n";
        }
    }

    return 0;
}

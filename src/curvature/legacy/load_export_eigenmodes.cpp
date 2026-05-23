#include "FESpace/FESpace.hpp"
#include "FESpace/finiteElement.hpp"
#include "FESpace/paraview.hpp"
#include "common/cut_mesh.hpp"
#include "cutfem.hpp"
#include "num/util.hpp"
#include "num/matlab.hpp"
#include "problem/testFunction.hpp"
#include <string>
#include <vector>
#include <filesystem>

// Types for better readability
// One of these needed
// #define QUADS
#define TRIANGLES

// One of these needed too
#define PEANUT_elem
// #define PEANUT_surf
// #define STAR_elem
// #define STAR_surf
// #define AMOEBA_elem
// #define AMOEBA_surf

const bool patching = true;

const double gamma_a = 0.0; // Gamma >0 destroys convergence

#if defined(QUADS)
    using lagrange_t    = LagrangeQuad2;
    using mesh_t        = MeshQuad2;
    const std::string mesh = "Quads";
#elif defined(TRIANGLES) 
    using lagrange_t    = Lagrange2;
    using mesh_t        = Mesh2;
    const std::string mesh = "Tris";
#else
    #error "No mesh type defined"
#endif

using fun_test_t    = TestFunction<mesh_t>;
using fun_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using space_t       = GFESpace<mesh_t>;
using cutspace_t    = CutFESpace<mesh_t>;
using cutfem_t      = CutFEM<mesh_t>;
using paraview_t    = Paraview<mesh_t>;
using Rd_t          = mesh_t::Rd;

using namespace globalVariable;



namespace peanut {
    const std::string example = "Peanut";

    // Computational domain
    // For testing mesh location robustness
    const double offset_x = 0.0;
    const double offset_y = 0.0;

    const double X0 = -1.65 + offset_x;
    const double Y0 = -0.7 + offset_y;

    const double DX = 1.65 + 1.6;
    const double DY = 1.4;

    // Mesh settings
    double const h0 = DX / 100;
    const int n_refinements = 12;

    // Stabilization
    // const double tauP = 0.1;
    #if defined(QUADS)
        const double tauP = 1e1;
    #else
        const double tauP = 1e-1;
    #endif
    const double tauS[] = {0.0,1e-4, 1e-3, 1e-2, 1e-1, 1e0};
    const double tauE[] = {0.0,1e-3, 1e-2, 1e-1, 1e0, 1e1};


    // Geometry
    const double a = 1.0;
    const double b = 1.02 * a;
    const double c = 0.5;

    
    // Phi and first derivatives need to work both with doubles and algoim::Interval<N>
    template <typename T> T sq(const T& x) { return x * x; }
    template<typename T>
    T phi(const T& x, const T& y){
        T factor1  = sq(x-a)+sq(y);
        T factor2  = sq(x+a)+sq(y);

        T result   = factor1*factor2 - (b*b*b*b) + c*x;
        return result;
    }
    template<typename T>
    T phi_x(const T& x, const T& y){
        T factor1 = sq(x - a) + sq(y);
        T factor2 = sq(x + a) + sq(y);

        T d_factor1_dx = 2.0 * (x - a);
        T d_factor2_dx = 2.0 * (x + a);

        T result = d_factor1_dx * factor2 + factor1 * d_factor2_dx + c;
        return result;
    }
    template<typename T>
    T phi_y(const T& x, const T& y){
        T result = 4.0*y*(sq(x) + a*a + sq(y));
        return result;
    }

    // Higher order derivatives only needed for reference solutions
    double phi_xx(double x, double y){
        // double x = P[0], y = P[1];
        double result = 4.0*(pow(x,2)+pow(a,2)+pow(y,2)) + 8*(pow(x,2)-pow(a,2));

        return result;
    }
    double phi_xy(double x, double y){
        // double x = P[0], y = P[1];
        double result = 8.0*x*y;

        return result;
    }
    double phi_yy(double x, double y){
        // double x = P[0], y = P[1];
        double result = 4.0*(pow(x,2)+pow(a,2)+3*pow(y,2));

        return result;
    }
}

// Gemini made the polar namespace compatible with algoim and interval arithmetic
#include <stdexcept>

namespace polar {
    // ---------------------------------------------------------
    // 1. Exact evaluations for standard doubles (Analytical roots)
    // ---------------------------------------------------------
    inline double r(const double& x, const double& y) { return std::sqrt(x*x + y*y); }
    inline double theta(const double& x, const double& y) { return std::atan2(y, x); }
    
    inline double rx(const double& x, const double& y) { return x / r(x,y); }
    inline double ry(const double& x, const double& y) { return y / r(x,y); }
    
    inline double thetax(const double& x, const double& y) { double r2 = x*x + y*y; return -y / r2; }
    inline double thetay(const double& x, const double& y) { double r2 = x*x + y*y; return x / r2; }

    inline double rxx(const double& x, const double& y) { double r_val = r(x,y); return (y*y) / (r_val*r_val*r_val); }
    inline double ryy(const double& x, const double& y) { double r_val = r(x,y); return (x*x) / (r_val*r_val*r_val); }
    inline double rxy(const double& x, const double& y) { double r_val = r(x,y); return -(x*y) / (r_val*r_val*r_val); }

    inline double thetaxx(const double& x, const double& y) { double r2 = x*x + y*y; return (2.0*x*y) / (r2*r2); }
    inline double thetayy(const double& x, const double& y) { double r2 = x*x + y*y; return -(2.0*x*y) / (r2*r2); }
    inline double thetaxy(const double& x, const double& y) { double r2 = x*x + y*y; return (y*y - x*x) / (r2*r2); }

    // ---------------------------------------------------------
    // 2. Safe bounding evaluations for Algoim intervals
    // ---------------------------------------------------------
    template <int N>
    inline algoim::Interval<N> r(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        try {
            return algoim::sqrt(x*x + y*y);
        } catch (const std::domain_error&) {
            double bx = x.maxDeviation(), by = y.maxDeviation();
            double r_max = std::sqrt((std::abs(x.alpha) + bx)*(std::abs(x.alpha) + bx) +
                                     (std::abs(y.alpha) + by)*(std::abs(y.alpha) + by));
            return algoim::Interval<N>(r_max / 2.0, algoim::uvector<double, N>(0.0), r_max / 2.0);
        }
    }

    template <int N>
    inline algoim::Interval<N> theta(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        algoim::Interval<N> r2 = x*x + y*y;
        
        // Manual check for theta since it relies on analytical double calculations internally
        if (r2.maxDeviation() >= 0.5 * std::abs(r2.alpha) || r2.alpha < 1e-12) {
            return algoim::Interval<N>(0.0, algoim::uvector<double, N>(0.0), M_PI);
        }
        
        double a_x = x.alpha, a_y = y.alpha;
        double r2_val = a_x * a_x + a_y * a_y;
        double val = std::atan2(a_y, a_x);
        
        algoim::uvector<double, N> beta;
        for (int i = 0; i < N; ++i) {
            beta(i) = (a_x * y.beta(i) - a_y * x.beta(i)) / r2_val;
        }
        double eps = (x.eps + y.eps) / std::sqrt(r2_val) + 
                     (x.maxDeviation() * x.maxDeviation() + y.maxDeviation() * y.maxDeviation()) / r2_val;
        return algoim::Interval<N>(val, beta, eps);
    }

    template <int N>
    inline algoim::Interval<N> rx(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        try { return x / algoim::sqrt(x*x + y*y); } 
        catch (const std::domain_error&) { return algoim::Interval<N>(0.0, algoim::uvector<double, N>(0.0), 1.0); }
    }

    template <int N>
    inline algoim::Interval<N> ry(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        try { return y / algoim::sqrt(x*x + y*y); } 
        catch (const std::domain_error&) { return algoim::Interval<N>(0.0, algoim::uvector<double, N>(0.0), 1.0); }
    }

    template <int N>
    inline algoim::Interval<N> thetax(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        try { return -y / (x*x + y*y); } 
        catch (const std::domain_error&) { return algoim::Interval<N>(0.0, algoim::uvector<double, N>(0.0), 1e6); }
    }

    template <int N>
    inline algoim::Interval<N> thetay(const algoim::Interval<N>& x, const algoim::Interval<N>& y) {
        try { return x / (x*x + y*y); } 
        catch (const std::domain_error&) { return algoim::Interval<N>(0.0, algoim::uvector<double, N>(0.0), 1e6); }
    }
}
// End of Gemini contribution

namespace star {
    const std::string example = "Star";

    // Computational domain
    const double offset_x = 3.333e-4;
    const double offset_y = 3.333e-4;

    const double X0 = -1.4 + offset_x;
    const double Y0 = -1.5 + offset_y;

    const double DX = 1.4 * 2.0;
    const double DY = 3.0;

    // Mesh settings
    double const h0 = DX / 70.0; // Changing this to 200 gives a segfault on the first mesh already?
    const int n_refinements = 12;

    // Stabilization
    #if defined(QUADS)
        const double tauP = 1e-1;
    #else
        const double tauP = 1e-2;
    #endif
    const double tauS[] = {0.0,1e-4, 1e-3, 1e-2, 1e-1, 1e0};
    const double tauE[] = {0.0,1e-4, 1e-3, 1e-2, 1e-1, 1e0};

    // Geometry
    using namespace polar;

    const double R0 = 1.0; 
    const double A  = 0.3; 
    const double K  = 7.0; 

    template <typename T> T get_R(const T& t) { 
        using std::cos;
        return R0 + A * cos(K*t); 
    }
    
    template <typename T> T get_dR(const T& t) { 
        using std::sin;
        return -A * K * sin(K*t); 
    }
    
    template <typename T> T get_ddR(const T& t) { 
        using std::cos;
        return -A * K * K * cos(K*t); 
    }

    template <typename T> T phi(const T& x, const T& y) {
        return r(x,y) - get_R(theta(x,y));
    }

    template <typename T> T phi_x(const T& x, const T& y) {
        T t = theta(x,y);
        return rx(x,y) - get_dR(t) * thetax(x,y);
    }

    template <typename T> T phi_y(const T& x, const T& y) {
        T t = theta(x,y);
        return ry(x,y) - get_dR(t) * thetay(x,y);
    }

    template <typename T> T phi_xx(const T& x, const T& y) {
        T t    = theta(x,y);
        T Rp   = get_dR(t);
        T Rpp  = get_ddR(t);
        T th_x = thetax(x,y);
        return rxx(x,y) - (Rpp * th_x * th_x + Rp * thetaxx(x,y));
    }

    template <typename T> T phi_yy(const T& x, const T& y) {
        T t    = theta(x,y);
        T Rp   = get_dR(t);
        T Rpp  = get_ddR(t);
        T th_y = thetay(x,y);
        return ryy(x,y) - (Rpp * th_y * th_y + Rp * thetayy(x,y));
    }

    template <typename T> T phi_xy(const T& x, const T& y) {
        T t    = theta(x,y);
        T Rp   = get_dR(t);
        T Rpp  = get_ddR(t);
        return rxy(x,y) - (Rpp * thetax(x,y) * thetay(x,y) + Rp * thetaxy(x,y));
    }
}

namespace amoeba {
    const std::string example = "Amoeba";

    // Computational domain
    const double offset_x = 0.0;
    const double offset_y = 0.0;

    const double X0 = -1.0 + offset_x;
    const double Y0 = -1.4 + offset_y;

    const double DX = 1.0 + 1.65;
    const double DY = 1.4 + 1.3;

    // Mesh settings
    double const h0 = DX / 70;
    const int n_refinements = 12;

    // Stabilization
    #if defined(QUADS)
        const double tauP = 1e2;
    #else
        const double tauP = 1e-1;
    #endif
    const double tauS[] = {0.0,1e-4, 1e-3, 1e-2, 1e-1, 1e0};
    const double tauE[] = {0.0,1e-4, 1e-3, 1e-2, 1e-1, 1e0};

    // Geometry
    using namespace polar;

    const double R0 = 1.0;
    const double C1 = 0.3;
    const double S2 = 0.2;
    const double C3 = 0.1;

    template <typename T> T get_R(const T& t) {
        using std::cos; 
        using std::sin;
        return R0 + C1*cos(t) + S2*sin(2.0*t) + C3*cos(3.0*t);
    }

    template <typename T> T get_dR(const T& t) {
        using std::cos; 
        using std::sin;
        return -C1*sin(t) + 2.0*S2*cos(2.0*t) - 3.0*C3*sin(3.0*t);
    }

    template <typename T> T get_ddR(const T& t) {
        using std::cos; 
        using std::sin;
        return -C1*cos(t) - 4.0*S2*sin(2.0*t) - 9.0*C3*cos(3.0*t);
    }

    template <typename T> T phi(const T& x, const T& y) {
        return r(x,y) - get_R(theta(x,y));
    }

    template <typename T> T phi_x(const T& x, const T& y) {
        T t = theta(x,y);
        return rx(x,y) - get_dR(t) * thetax(x,y);
    }

    template <typename T> T phi_y(const T& x, const T& y) {
        T t = theta(x,y);
        return ry(x,y) - get_dR(t) * thetay(x,y);
    }

    template <typename T> T phi_xx(const T& x, const T& y) {
        T t   = theta(x,y);
        T Rp  = get_dR(t);
        T Rpp = get_ddR(t);
        T th_x = thetax(x,y);
        return rxx(x,y) - (Rpp * th_x * th_x + Rp * thetaxx(x,y));
    }

    template <typename T> T phi_yy(const T& x, const T& y) {
        T t   = theta(x,y);
        T Rp  = get_dR(t);
        T Rpp = get_ddR(t);
        T th_y = thetay(x,y);
        return ryy(x,y) - (Rpp * th_y * th_y + Rp * thetayy(x,y));
    }

    template <typename T> T phi_xy(const T& x, const T& y) {
        T t   = theta(x,y);
        T Rp  = get_dR(t);
        T Rpp = get_ddR(t);
        return rxy(x,y) - (Rpp * thetax(x,y) * thetay(x,y) + Rp * thetaxy(x,y));
    }
}
#if defined(PEANUT_elem)
    using namespace peanut;
    const bool surfing = false;
    const bool eleming = true;
#elif defined(PEANUT_surf)
    using namespace peanut;
    const bool surfing = true;
    const bool eleming = false;
#elif defined(STAR_elem)
    using namespace star;
    const bool surfing = false;
    const bool eleming = true;
#elif defined(STAR_surf)
    using namespace star;
    const bool surfing = true;
    const bool eleming = false;
#elif defined(AMOEBA_elem)
    using namespace amoeba;
    const bool surfing = false;
    const bool eleming = true;
#elif defined(AMOEBA_surf)
    using namespace amoeba;
    const bool surfing = true;
    const bool eleming = false;
#else
    #error "No example defined"
#endif

namespace common {
    R2 normal(double x, double y) {
        double phix = phi_x(x,y);
        double phiy = phi_y(x,y);
        double r    = sqrt(phix*phix + phiy*phiy);
        return R2(phix/r, phiy/r);
    }

    double curvature(double x, double y) {
        double phix = phi_x(x,y);
        double phiy = phi_y(x,y);
        double phixx = phi_xx(x,y);
        double phixy = phi_xy(x,y);
        double phiyy = phi_yy(x,y);

        double kappa = phiy*phiy*phixx - 2*phix*phiy*phixy + phix*phix*phiyy;
        kappa /= pow(phix*phix + phiy*phiy, 1.5);
        return kappa;
    }

    R2 curvature_vector(double x, double y) {
        return (-curvature(x,y))*normal(x,y);
    }

    double fun_levelset(double* P, int i = 0) {
        return phi(P[0], P[1]);
    }

    double fun_curvature(double* P, int direction) {
        R2 H = curvature_vector(P[0], P[1]);
        if (direction == 0) return H[0];
        else if (direction == 1) return H[1];
        else return std::numeric_limits<double>::max();
    }

    double fun_nx(double* P, int domin) {
        return normal(P[0], P[1])[0];
    }

    double fun_ny(double* P, int domin) {
        return normal(P[0], P[1])[1];
    }

    double fun_nx_max(double* P, int domin, double t) {
        return fun_nx(P, 0);
    }

    double fun_ny_max(double* P, int domin, double t) {
        return fun_ny(P, 0);
    }
}

using namespace common;

class CurvatureP1 {
public:
    CurvatureP1(const space_t &V, const InterfaceLevelSet<mesh_t> &interface, const fun_t &phi_h)
        : interface_(interface), phi_h_(phi_h)
    {
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_->createSurfaceMesh(interface_);
        
        Vh_ = std::make_unique<cutspace_t>(*surface_mesh_, V);

        ProblemOption problem_option;
        problem_option.clear_matrix_ = false;

        problem = std::make_unique<cutfem_t>(*Vh_, problem_option);
    }

    void weak_form(){
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        Rnm Id(D_mesh, D_mesh);
        Id = 0.0;
        for (int i = 0; i < D_mesh; ++i){
            Id(i,i) = 1.0;
        }

        // a_h
        problem->addBilinear(+ innerProduct(H,v), interface_);

        // L_h
        problem->addLinear( - contractProduct(Id, gradS(v)), interface_);
    }

    // Face stabilization scales with h^(d-3)
    // Now we have d=2, so all terms should scale with h^(-1)

    void patch_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();

        problem->addPatchStabilization(
                tau * pow(h, -3 + gamma_a) * innerProduct(jump(H),jump(v)), 
                *surface_mesh_);
    }

    void face_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();
        Normal n;
        
        //TODO: This needs a for loop for higher order terms as well to be proper
        // currently unused as Patch works well
        problem->addFaceStabilization(
            tau * innerProduct(jump(grad(H)*n), jump(grad(v)*n)),
            *surface_mesh_);
    }

    // TODO: Both elem and surf need scaling fixed (and for loops for higher order)

    void elem_stab(const double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();

        auto phi_dx = dx(phi_h_.expr());
        auto phi_dy = dy(phi_h_.expr());
        auto grad_norm = sqrt(phi_dx * phi_dx + phi_dy * phi_dy);

        // Construct the symbolic expressions for the normal components
        auto n_x = phi_dx / grad_norm;
        auto n_y = phi_dy / grad_norm;

        auto d_n_H = n_x * dx(H) + n_y * dy(H);
        auto d_n_v = n_x * dx(v) + n_y * dy(v);

        // Integral in 2D -> h²
        // grad(h)*grad(v) -> h⁻²
        // must scale by h⁻¹
        problem->BaseFEM::addBilinear( 
        tau * pow(h,-1.0 + gamma_a) * innerProduct(d_n_H, d_n_v), 
        *surface_mesh_);
    }

    void surf_stab(double tau) {
        constexpr int D_mesh = mesh_t::Rd::d;
        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        const double h = surface_mesh_->Th.get_mesh_size();

        Normal n;

        auto d_n_H = grad(H) * n;
        auto d_n_v = grad(v) * n;

        // Integral 1D -> h¹
        // grad(H)*grad(v) -> h⁻²
        // no scaling on first term

        problem->addBilinear(
            tau * innerProduct(d_n_H, d_n_v),
            interface_
        );
    }

    fun_t solve() {
        problem->solve();
        return fun_t(*Vh_, problem->rhs_);
    }

    double get_error(const fun_t& approx) {
        // Standard P1 interface integration norm evaluation
        return L2normSurf(approx, fun_curvature, interface_, 0, 2);
    }

    double normal_error_max() {
        auto phi_dx = dx(phi_h_.expr());
        auto phi_dy = dy(phi_h_.expr());
        auto grad_norm = sqrt(phi_dx * phi_dx + phi_dy * phi_dy);

        auto nx = phi_dx / grad_norm;
        auto ny = phi_dy / grad_norm;

        // Note: Relies on MaxnormSurf being defined analogously to L2normSurf in the library framework.
        double error1 = max_norm_surface(nx, fun_nx_max, interface_, 0.0);
        double error2 = max_norm_surface(ny, fun_ny_max, interface_, 0.0);

        return error1 + error2;
    }

    // Replaces the old export method
    void export_active_matrix(const std::string& base_filename) const {
        std::ofstream out_mat(base_filename + "_mat.dat", std::ios::out);
        std::ofstream out_map(base_filename + "_map.dat", std::ios::out);
        if (!out_mat.is_open() || !out_map.is_open()) return;
        
        out_mat << std::setprecision(16);

        std::set<int> unique_dofs;
        int max_dof = 0; // Tracks the true algebraic dimension of the vector block
        
        for (const auto& entry : problem->mat_) {
            unique_dofs.insert(entry.first.first);
            unique_dofs.insert(entry.first.second);
            max_dof = std::max({max_dof, entry.first.first, entry.first.second});
        }

        std::map<int, int> global_to_compressed;
        int active_idx = 1;
        
        // Write the true 2N vector dimension to the header for MATLAB
        out_map << (max_dof + 1) << "\n"; 
        
        for (int dof : unique_dofs) {
            global_to_compressed[dof] = active_idx++;
            out_map << dof + 1 << "\n"; 
        }

        for (const auto& entry : problem->mat_) {
            int i = global_to_compressed[entry.first.first];
            int j = global_to_compressed[entry.first.second];
            out_mat << i << "\t" << j << "\t" << entry.second << "\n";
        }
        
        out_mat.close();
        out_map.close();
    }

    // Add these getters at the bottom of your public methods
    int get_cut_ndof() const { return Vh_->NbDoF(); }
    const cutspace_t& get_cut_space() const { return *Vh_; }

    void clear_matrix() {
        problem->mat_.clear();
    }

    const activemesh_t &mesh() const { return *surface_mesh_; }

private:
    const InterfaceLevelSet<mesh_t> &interface_;
    fun_t phi_h_;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::unique_ptr<cutfem_t> problem;
};

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // 1. Determine the exact example name
    std::string example_name = example;
    if (eleming) { example_name += "_elem"; } 
    else if (surfing) { example_name += "_surf"; } 
    else { example_name += "_patch"; }
    if (!patching) { example_name += "_only"; }
    example_name += ("_"+mesh);

    int num_modes = 3;
    int stab_cases = sizeof(tauE) / sizeof(tauE[0]);

    const std::string base_output_path = "../output_files/curvature_p1_2d/" + example_name + "/";
    std::string path_input_data = base_output_path + "data/";
    std::string path_output_figures = base_output_path + "paraview/";

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_figures);
        std::cout << "Starting batch eigenmode ingestion for: " << example_name << "\n";
    }

    double h = h0;
    
    // Array to store the generated VTK filenames for the PVD playlist
    std::vector<std::vector<std::string>> vtk_collection(stab_cases);

    // 2. Loop through all refinements automatically
    for (int i = 0; i < n_refinements; ++i) {
        
        const int nx = static_cast<int>(DX / h);
        const int ny = static_cast<int>(DY / h);

        mesh_t Th(nx, ny, X0, Y0, DX, DY);
        
        lagrange_t FEcurv(1);
        space_t CurvV(Th, FEcurv);
        
        space_t Lh(Th, DataFE<mesh_t>::P1);
        fun_t phi_h(Lh, fun_levelset);
        InterfaceLevelSet<mesh_t> interface(Th, phi_h);

        CurvatureP1 curvature(CurvV, interface, phi_h);

        // 3. Loop through stabilisation cases
        for (int j = 0; j < stab_cases; ++j) {
            
            std::string check_filepath = path_input_data + "eigmode_1_ref" + std::to_string(i+1) + "_stab" + std::to_string(j+1) + ".dat";
            if (!std::filesystem::exists(check_filepath)) {
                continue; 
            }

            if (MPIcf::IamMaster()) {
                std::string ref_str = (i < 9 ? "0" : "") + std::to_string(i + 1);
                
                // Keep the relative basename for the PVD file
                std::string vtk_basename = example_name + "_stab" + std::to_string(j+1) + "_ref" + ref_str + ".vtk";
                std::string vtk_filename = path_output_figures + vtk_basename;
                
                paraview_t para_eig(curvature.mesh(), vtk_filename);

                int vector_space_size = 2 * curvature.get_cut_ndof();

                for (int mode = 1; mode <= num_modes; ++mode) {
                    std::string eig_filepath = path_input_data + "eigmode_" + std::to_string(mode) + "_ref" + std::to_string(i+1) + "_stab" + std::to_string(j+1) + ".dat";
                    
                    std::vector<double> v_cut_data;
                    matlab::loadVector(v_cut_data, vector_space_size, eig_filepath);

                    fun_t eigenmode_fun(curvature.get_cut_space(), v_cut_data);
                    std::string field_name = "Eigenmode_" + std::to_string(mode);
                    para_eig.add(eigenmode_fun, field_name, 0, 2); 
                }
                
                para_eig.add(phi_h, "levelset", 0, 1);
                
                // Log the successful export for the master playlist
                vtk_collection[j].push_back(vtk_basename);
            }
        }
        h *= sqrt(0.5); 
    }

    // 4. Generate the VTK Series JSON Files
    if (MPIcf::IamMaster()) {
        for (int j = 0; j < stab_cases; ++j) {
            if (vtk_collection[j].empty()) continue;

            // ParaView JSON series files MUST end with the exact extension .vtk.series
            std::string series_filename = path_output_figures + example_name + "_stab" + std::to_string(j+1) + ".vtk.series";
            std::ofstream series_file(series_filename);
            
            if (series_file.is_open()) {
                series_file << "{\n";
                series_file << "  \"file-series-version\" : \"1.0\",\n";
                series_file << "  \"files\" : [\n";
                
                for (size_t t = 0; t < vtk_collection[j].size(); ++t) {
                    // Assign the explicit integer time step to each legacy VTK file
                    series_file << "    { \"name\" : \"" << vtk_collection[j][t] << "\", \"time\" : " << (t + 1) << " }";
                    
                    // JSON requires commas between objects, but not after the last one
                    if (t < vtk_collection[j].size() - 1) {
                        series_file << ",";
                    }
                    series_file << "\n";
                }
                
                series_file << "  ]\n";
                series_file << "}\n";
                series_file.close();
                std::cout << "Successfully generated ParaView series playlist: " << series_filename << "\n";
            }
        }
    }

    return 0;
}
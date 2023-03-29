#ifndef CS_MYJSON_HPP

#include "csuserfuns.hpp"

using namespace std;
using namespace cs_user_funs;

class ud_t :public user_data_t {
public:
  void inituser();
  double issolid(cs_real_t *);
  string fscase;
  cs_real_t zini,gz;
  cs_real_t R1,R2,mass0,Lmom0;
  // For the wall case
  cs_real_t Hwall,twall,gcirc,ucirc,wcirc;
  cs_real_t Dphi,H1,H2,dphi;
  vector<cs_real_t> xsphere,vsphere,force;
  cs_real_t Rsphere;
  cs_lnum_t ncells;
  cs_real_t rho1,rho2,delta,rhosolid,bdymass,visco;
  // Pile data
  cs_real_t wglen,dgz,T,Rpile,xpile,H;
  // Boxwall data
  cs_real_t Lx,Ly,Lz;
  cs_lnum_t Nx,Ny,Nz;
  vector<cs_real_t> eup,edown;
  cs_real_t W;
  // Wavegen data
  cs_real_t Twave, gwave, lagen,wgen,kgen,
    abso_length,Labso;
  vector<cs_real_t> upre;
  vector<cs_real_t> uxstr, dzstrh;
  cs_real_t gx,gzt;
#define USE_PART_STRUCT
#ifdef USE_PART_STRUCT
  // Set the partially strutured case info
  int dummy_case;
  int k0start;
  void z1init();
  int get_gcid(int jlay,int jz);
  vectormi_t<int> zstrucindx;
#endif
  ud_t();
  double Azbdy,wz,Tz;
  double dzfun(double time);
  // Current displacement of the box
  vector<double> dzbox;
  // Stuff for the structure
  Eigen::MatrixXd A,B,C;
  Eigen::VectorXd y;
  double theta,Dt;
  void advance_struct(double fz);
};

// Set visualization to 0 fo non important fields
void set_viz();

extern ud_t user_data;

#endif // CS_MYJSON_HPP

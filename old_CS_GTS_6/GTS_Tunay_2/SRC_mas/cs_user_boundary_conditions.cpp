/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
// #include <vector>
// #include <algorithm>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_elec_model.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_selector.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include "./userdata.hpp"
#include "./csuserfuns.hpp"

using namespace std;
using namespace cs_user_funs;

#define BFNOR(fid,k) bfnor[3*fid+k]
#define BFCOG(fid,k) bfcog[3*fid+k]
#define ICODCL(fid,var) icodcl[((var)-1)*nfaces+(fid)]
#define RCODCL(fid,var,typ) rcodcl[((typ)-1)*nfaces*nvar+((var)-1)*nfaces+(fid)]

extern "C"
void usalcl_(int *itrale, int *nvar , int *nscal, 
             int *icodcl, int *itypfb, int *ialtyb, int *impale , 
             cs_real_t *dt, cs_real_t *rcodcl, cs_real_t *xyzno0,
             cs_real_t *depale) {
  ud_t &d = user_data;
  d.inituser();
  Json::Value &opts = d.opts;

  // If no FSI then do nothing
  GETID(do_fsi);
  if (!do_fsi) return;
  
  string &runcase = d.runcase;

  cs_mesh_t *m = cs_glob_mesh;
  cs_lnum_t nnod = m->n_vertices;
  cs_lnum_t nfaces = m->n_b_faces;

  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *bfcog = fvq->b_face_cog;
  const cs_real_t *bfnor = fvq->b_face_normal;
  const cs_real_t *b_face_surf = (const cs_real_t *)fvq->b_face_surf;

  int step = cs_glob_time_step->nt_cur;
  cs_real_t time = cs_glob_time_step->t_cur;
  cs_real_t wdx = 6.0;
  // From /base/paramx.f90
#define IBFIXE 1
#define IGLISS 2
#define IVIMPO 3

  GETDD(Lx);
  GETDD(Ly);
  GETDD(Lz);
  GETDD(T);
  GETDD(zfs);
  GETDD(hz);
  GETID(dim);
  GETID(has_buoy);
  GETID(use_ysym);

  cs_real_t tol=1e-6;
  // We base the setting of surfaces by normals so that we set that
  // only in the first step
  cs_real_t Aglis=0.0,Aimpo=0.0;
  // We leave IGLISS only on the vertical walls
  for (int j=0; j<nfaces; j++) {
    //   || BFCOG(j,1)<tol || BFCOG(j,1)>Ly-tol;
    int glis = (BFCOG(j,1)<tol || BFCOG(j,1)>Ly-tol);
    if (use_ysym) glis |= BFCOG(j,1)>0.5*Ly-tol;
    ialtyb[j] = (glis ? IGLISS : IVIMPO);
    (glis? Aglis : Aimpo) += b_face_surf[j];
  }
  cs_parall_sum(1,CS_DOUBLE,&Aglis);
  cs_parall_sum(1,CS_DOUBLE,&Aimpo);
  if (step==1 && cs_glob_rank_id<=0) 
    printf("step %d, Aglis %g, Aimpo %g\n",step,Aglis,Aimpo);
  const int ndim=3;
  
#define X0(j,k) xyzno0[ndim*(j)+(k)]
#define X(j,k) xyznod[ndim*(j)+(k)]
#define DX(j,k) depale[ndim*(j)+(k)]

  vector<cs_real_t> xc0(ndim,0.0);

  // Get the position of a reference node. The mesh may have
  // its own reference position, in that case get it from
  // there
  if(has_buoy) {
    Json::Value z = opts["xyz0"];
    for (int j=0; j<ndim; j++) xc0[j] = z[j].asDouble();
  }

  if (d.fscase=="pile") 
    for (cs_lnum_t j=0; j<nnod; j++)
      impale[j]=1;
  
  int count=0;
  cs_real_t dzmin=INFINITY,dzmax=-INFINITY;
  vector<cs_real_t> dxmax(ndim,0.0);
  for (cs_lnum_t j=0; j<nnod; j++) {
    cs_real_t adx = 0.0;
#define ADX(l) fabs(X0(j,l)-xc0[l])
    if (dim==2) adx = ADX(0)+ADX(2);
    else adx = ADX(0)+ADX(1)+ADX(2);
    if (adx<tol) {
      // printf("x %g %g %g\n",X0(j,0),X0(j,1),X0(j,2));
      dzmin = min(dzmin,DX(j,2));
      dzmax = max(dzmax,DX(j,2));
      count++;
    }
    for (int k=0; k<ndim; k++)
      dxmax[k] = max(dxmax[k],fabs(DX(j,k)));
  }
  cs_parall_max(3,CS_DOUBLE,dxmax.data());
  // printf("dxmax %g %g %g\n",dxmax[0],dxmax[1],dxmax[2]);

  cs_parall_sum(1,CS_INT32,&count);
  if (step==0 && cs_glob_rank_id<=0) 
    printf("using xyz0: %g %g %g, count %d\n",xc0[0],xc0[1],xc0[2],count);
  
  // We must find at least one node
  assert(count>=1);
  
  cs_parall_min(1,CS_DOUBLE,&dzmin);
  cs_parall_max(1,CS_DOUBLE,&dzmax);
  assert(fabs(dzmax-dzmin)<1e-10);
  cs_real_t dzd = NAN;
  vector<cs_real_t> &dzstrh = d.dzstrh;
  dzstrh.push_back(dzmin);
  cs_real_t dz = dzmin;
  // Do an extrapolation from last two entries
  if (dzstrh.size()>=2) {
    int n = dzstrh.size();
    dz = 2*dzstrh[n-1]-dzstrh[n-2];
  }

  // if (cs_glob_rank_id<=0) 
  //   printf("step %d, time %g, dz0 %g dze %g, gx %g %g\n",
  //          step,time,dzmin,dz,d.gx,d.gzt);

  // Z0: z-positions defining the intervals in the REFERENCE mesh. 
  // ZNOW: z-positions defining the intervals in the CURRENT mesh.
  // ZFS: the reference position of the FS
  // HZ: the rigid region around the FS is 2*HZ
  // The fitted mesh is in the region ZFS+[-HZ,+HZ]
  vector<cs_real_t> z0(4),znow(4);
  // Current displacement of the buoy
  double w = 2.0*M_PI/T;
  // dz is the current displacement (imposed) of the buoy
  // dz = d.dzfun(time);
  // Take the displacement from the FSI
  dz = d.y(0);
  if (cs_glob_rank_id<=0) 
    printf("step %d, time %g, dz %g\n",step,time,dz);

  z0[0]=0.0; z0[1]=zfs-hz; z0[2]=zfs+hz; z0[3]=Lz; 
  znow=z0;
  znow[0] -= dz;
  znow[3] -= dz;
  for (cs_lnum_t j=0; j<nnod; j++) {
    // Coord of node in the ref mesh
    cs_real_t Z0 = X0(j,2);
    // Set imposed node
    impale[j] = 1;
    // New Z coord of node (interpolate piecewise linearly)
    cs_real_t Znow = lininterp(z0,znow,Z0);
    DX(j,2) = Znow-Z0;
  }
}

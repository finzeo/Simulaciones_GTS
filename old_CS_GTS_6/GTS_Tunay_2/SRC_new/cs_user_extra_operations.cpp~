/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *==========================================================================*/

/* Code_Saturne version 5.1-alpha */

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

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_balance_by_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include <climits>
#include "userdata.hpp"

extern "C" void
cs_user_extra_operations(cs_domain_t *domain) {
  set_viz();
  // CHRONO.elaps("extra ops start");
  // user_data.zrenorm();
  user_data.check_probes();
  // CHRONO.elaps("extra ops end");
  ud_t &d = user_data;
  d.inituser();
  Json::Value &opts = d.opts;

  const int ndim=3;  
  cs_lnum_t ncells = cs_glob_mesh->n_cells;
  cs_real_t *bforces = cs_field_by_name("boundary_forces")->val;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  cs_lnum_t nfaces = cs_glob_mesh->n_b_faces;
  const cs_real_t *bfcog = fvq->b_face_cog;
  const cs_real_t *bfsurf = (const cs_real_t *)fvq->b_face_surf;
  const cs_real_t *cell_vol = fvq->cell_vol;
#define BFCOG(fid,k) bfcog[ndim*(fid)+(k)]
#define BFORCES(fid,k) bforces[3*fid+k]
  GETDD(Lx);
  GETDD(Ly);
  GETDD(Lz);
  double tol = 1e-6;
  int step = cs_glob_time_step->nt_cur;
  cs_real_t time = cs_glob_time_step->t_cur;
/*
  GETID(do_fsi);
  if (do_fsi) {
    double dzbox = +INFINITY;
    for (cs_lnum_t j = 0; j<nfaces; j++) {
      const cs_real_t *cog = &BFCOG(j,0);
      dzbox = min(dzbox,cog[2]);
    }
    cs_parall_min(1,CS_DOUBLE,&dzbox);
    GETDD(Dt);
*/
/*
    double Abdy=0.0,Abox=0.0;
    vector<cs_real_t> force(ndim,0.0);
    for (cs_lnum_t j = 0; j<nfaces; j++) {
      const cs_real_t *cog = &BFCOG(j,0);
      if (cog[0]>tol && cog[0]<Lx-tol
          && cog[1]>tol && cog[1]<Ly-tol
          && cog[2]-dzbox>tol
          && cog[2]-dzbox<Lz-tol) {
        Abdy += bfsurf[j];
        for (int k=0; k<ndim; k++) force[k] += BFORCES(j,k);
      } else {
        Abox += bfsurf[j];
      }
    }

    cs_parall_sum(ndim,CS_DOUBLE,force.data());
    cs_parall_sum(1,CS_DOUBLE,&Abdy);
    cs_parall_sum(1,CS_DOUBLE,&Abox);
    GETDD(mass);
    GETDD(gz);
    // FORCE[2] is the Z force of the fluid, mass*gz is the weight
    double fz=force[2]+mass*gz;
    // Add an external force
    GETDD(zimpulse);
    // Duration of impulse,start
    double t0=5*Dt,dtimp = 0.15;
    double fext = zimpulse*regdelta(time,t0,t0+dtimp);
    fz += fext;
    // We affect the force to be applied to the structure by a
    // factor in order to damp the initial transient
    double tfac = regheavis(time,5*Dt,10*Dt);
    d.advance_struct(tfac*fz);
    double *yp = d.y.data();
    d.dzbox.push_back(yp[0]);
    if (cs_glob_n_ranks>1)
      MPI_Bcast(yp,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    PPRINTF("time %f, f2 %f, force %f, fext %f, dz %f, vz %f\n",
            time,force[2],fz,fext,yp[0],yp[1]);

    if (step<=1)
      PPRINTF("Abdy %f, Abox %f, force %f %f %f\n",
              Abdy,Abox,force[0],force[1],force[2]);
*/
    // Store in the field "mesh_disp_rel" the displacements
    // relative to the box (not to the buoy)
    cs_lnum_t nvrtx = cs_glob_mesh->n_vertices;
    cs_real_t *velx1 = cs_field_by_name("velocity")->val;
    cs_real_t *velx2 = cs_field_by_name("velx2")->val;
#define VELX1(v,k) velx1[3*v+k]
#define VELX2(v,k) velx2[3*v+k]
    for (int j=0; j<nvrtx; j++) {
      VELX2(j,0) = 2.0*VELX1(j,0);
      VELX2(j,1) = 2.0*VELX1(j,1);
      VELX2(j,2) = 2.0*VELX1(j,2);
    }
}

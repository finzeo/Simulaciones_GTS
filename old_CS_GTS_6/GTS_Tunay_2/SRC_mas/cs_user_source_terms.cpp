/* GITVERSION<<undef>> */
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

using namespace std;
using namespace cs_user_funs;

extern "C"
void ustsnv_(int *nvar, int *nscal, int *ncepdp, int *ncesmp,
             int *ivar, int *icepdc, int *icetsm,
             int *itypsm, double *dt, double *ckupdc,
             double *smacel, double *crvexp, double *crvimp) {
  cs_lnum_t ncells = cs_glob_mesh->n_cells;
  cs_real_t *rhovals = cs_field_by_name("density")->val;
  cs_real_t *xcell = cs_glob_mesh_quantities->cell_cen;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *cell_vol = fvq->cell_vol;
  const int ndim=3;
  ud_t &d = user_data;
  d.inituser();
  Json::Value &opts = d.opts;

  int step = cs_glob_time_step->nt_cur;
  cs_real_t time = cs_glob_time_step->t_cur;

  // Compute acceleration of the body by finite differences
  double azbdy = 0.0;
  GETDD (Dt);
  GETID(do_fsi);
  if (do_fsi) {
    if (cs_glob_rank_id<=0) {
      // Add the acceleration of the fluid but not inmediately
      int nstep = d.dzbox.size();
      // Compute the acceleration of the body by finite differences 
      if (nstep>10) {
        azbdy = (d.dzbox[nstep-1]-2*d.dzbox[nstep-2]+d.dzbox[nstep-3])/(Dt*Dt);
      }
    }
    if (cs_glob_n_ranks>1)
      MPI_Bcast(&azbdy,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    PPRINTF("step %d, time %f, azbdy %g\n",step,time,azbdy);
  }

  // Period of the box 
  GETDD(Tx);
  // Amplitude of the box
  GETDD(Axbox);
  // Start and end of the transition period for the box motion
  double GETDP(tbox0,5*Dt);
  double GETDP(tbox1,0.25*Tx);
  
  double
    w = 2*M_PI/Tx,
    // We add the box motion very slowly after a certain
    // number of time steps
    tfac = regheavis(time,tbox0,tbox1),
    aboxt = Axbox*w*w*tfac*sin(w*time);
  // If the bdy is accelerating with Azbdy then the
  // non-inertial system (NIS) (follows the buoy) so that
  // the acceleration of the NIS is azbdy. Then we must add
  // a force term -mass*Azbdy
  double g[] = {aboxt,0,-azbdy};

#define XCELL(cid,k) xcell[(cid)*ndim+k]
#define U(cid,k) uvals[(cid)*ndim+k]
#define CRVEXP(cid,k) crvexp[(cid)*ndim+k]
#define CRVIMP(cid,k,l) crvimp[((cid)*ndim+(k))*ndim+(l)]

  for (cs_lnum_t j=0; j<ncells; j++) {
      for (cs_lnum_t k=0; k<ndim; k++) { 
        CRVEXP(j,k) += g[k]*rhovals[j]*cell_vol[j];
      }
  }
}

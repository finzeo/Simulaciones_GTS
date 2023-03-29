/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

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
#include <vector>

using namespace std;

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

extern "C" void
cs_user_extra_operations(cs_domain_t *domain) {
  const int ndim=3;  
  cs_real_t *forces = cs_field_by_name("boundary_forces")->val;
#define FORCES(fid,k) forces[(fid)*ndim+(k)]

  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *bfcog = fvq->b_face_cog;
#define BFCOG(fid,k) bfcog[ndim*(fid)+(k)]

  cs_lnum_t nfaces = cs_glob_mesh->n_b_faces;
  int step = cs_glob_time_step->nt_cur;
  cs_real_t time = cs_glob_time_step->t_cur;

  vector<cs_lnum_t> face_list(nfaces);
  cs_lnum_t nfcext;
  cs_selector_get_b_face_list("rint",&nfcext,face_list.data()); // colocar label
  vector<cs_real_t> force(ndim,0);
  cs_real_t xmin[ndim], xmax[ndim], xc[ndim];
  for (int j=0; j<ndim; j++) {
    xmin[j] = INFINITY;
    xmax[j] = -INFINITY;
  }
  for (cs_lnum_t j = 0; j<nfcext; j++) {
    cs_lnum_t fid = face_list[j];
    for (int k=0; k<ndim; k++) {
      force[k] += FORCES(fid,k);
      xmin[k] = min(xmin[k],BFCOG(fid,k));
      xmax[k] = max(xmax[k],BFCOG(fid,k));
    }
  }
  cs_parall_min(ndim,CS_DOUBLE,&xmin);
  cs_parall_max(ndim,CS_DOUBLE,&xmax);
  for (int k=0; k<ndim; k++) 
    xc[k] = 0.5*(xmax[k]+xmin[k]);
  cs_parall_sum(ndim,CS_DOUBLE,force.data());
  if (cs_glob_rank_id<=0)
    printf("step %d, time %f, forces %f %f %f, e %f %f %f\n",
           step,time,force[0],force[1],force[2],xc[0],xc[1],xc[2]);
}


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

#include <cassert>
#include <cmath>
#include <cstdlib>

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
#include "cs_stokes_model.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include <climits>
#include <iostream>
#include <vector>

using namespace std;

#include "vecmacros2.hpp"

double max(double x,double y) { return (x>y? x : y); }
double min(double x,double y) { return (x<y? x : y); }

extern "C" void
cs_user_extra_operations(cs_domain_t*) {
  const int ndim=3;
  cs_real_t *uvals = cs_field_by_name("velocity")->val;
#define UVALS(j,k) VEC2(uvals,j,k,ndim)
  cs_real_t time = cs_glob_time_step->t_cur;
  int step = cs_glob_time_step->nt_cur;

  cs_lnum_t ncells = cs_glob_mesh->n_cells;
  cs_real_t *xcell = cs_glob_mesh_quantities->cell_cen;
#define XCELL(j,k) VEC2(xcell,j,k,ndim)
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *cell_vol = fvq->cell_vol;

  double maxvel=2.0;
  // VOLTOT is the total volume, VOL the volume of the cells
  // that have U>MAXVEL. 
  double vol=0.0, voltot=0.0;
  // XMIN, XMAX define the bounding box of the cells (cell
  // centers in fact) whose velocities are > MAXVEL
  vector<double> xmin(ndim,+INFINITY),xmax(ndim,-INFINITY);
  int count=0;
  for (cs_lnum_t j=0; j<ncells; j++) {
    voltot += cell_vol[j];
    double avel = 0.0;
    for (int l=0; l<ndim; l++) avel += UVALS(j,l)*UVALS(j,l);
    avel = sqrt(avel);
    if (avel>=maxvel) {
      count++; vol += cell_vol[j];
      for (int l=0; l<ndim; l++) {
        xmax[l] = max(xmax[l],XCELL(j,l));
        xmin[l] = min(xmin[l],XCELL(j,l));
      }
    }
  }
  cs_parall_sum(1,CS_DOUBLE,&voltot);
  cs_parall_sum(1,CS_DOUBLE,&vol);
  cs_parall_sum(1,CS_INT32,&count);
  if (cs_glob_rank_id<=0) {
    printf("step %d, time %f, maxvel %f, total vol %g, vol(u>maxvel) %g, "
           "cell count %d\n",step,time,maxvel,voltot,vol,count);
    printf("bounding box: xmin %g %g %g, xmax %g %g %g\m \n",
           xmin[0],xmin[1],xmin[2],xmax[0],xmax[1],xmax[2]);
  }
}

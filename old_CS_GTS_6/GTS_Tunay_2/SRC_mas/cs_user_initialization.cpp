/*============================================================================
 * User initialization prior to solving time steps.
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
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_lagr.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"
#include "cs_post.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include "./userdata.hpp"

extern "C" void
cs_user_initialization(cs_domain_t *) {

  set_viz();
  user_data.inituser();

  cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t ncells = m->n_cells;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *xcell = fvq->cell_cen;
  const cs_real_t *cell_vol = fvq->cell_vol;
  const int ndim=3;
#define XCELL(cid,k) xcell[(cid)*ndim+k]

  cs_real_t *vfvals = cs_field_by_name("void_fraction")->val;
  ud_t &d = user_data;
  Json::Value opts = d.opts;

  GETDD(zfs);
  GETDD(Lx);
  GETDD(Ly);
  // Slope of free surface in Y direction
  double GETDP(zfsy,0.0);
  double GETDP(zfsx,0.0);
  for (cs_lnum_t j=0; j < ncells; j++) {
    // Coords of cell COG
    const double *x = &XCELL(j,0);
    // Free surface level here 
    double zfsh = zfs+zfsx*(x[0]-0.5*Lx)+zfsy*(x[1]-0.5*Ly);
    vfvals[j] = (x[2]>zfsh);
  }
}

// Wrapper that defines the value in the fortran variable
extern "C"
void vof_init_vals_(double *rho1p,double *rho2p,double *mu1p,double *mu2p) {
  ud_t &d = user_data;
  d.inituser();
  Json::Value opts = d.opts;
  // Gets the value in n aux variable and then copy to the fortran variable
#define G(name) GETDD(name); *name##p=name
  G(rho1);
  G(rho2);
  G(mu1);
  G(mu2);
}

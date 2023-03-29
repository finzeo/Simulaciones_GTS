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

// ud_t::ud_t() : dummy_case(-1), k0start(-1), dzbox(NAN) { }

void set_viz() {
  int kpost=cs_field_key_id("post_vis");
  int nf = cs_field_n_fields();
  for (int j=0; j<nf; j++) {
    cs_field_t *f = cs_field_by_id (j);
    int pv0 = cs_field_get_key_int(f,kpost);
    if (strcmp(f->name,"velx2")) {
      // PPRINTF("setting post_vis=0 for field %s\n",f->name);
      cs_field_set_key_int(f,kpost,1);
    }
    int pv1 = cs_field_get_key_int(f,kpost);
    if (pv1!=pv0) PPRINTF("%d %s post_vis before %d, after %d\n",j,f->name,pv0,pv1);
  }
}


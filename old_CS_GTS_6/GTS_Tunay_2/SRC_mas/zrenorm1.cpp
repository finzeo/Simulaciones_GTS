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

#ifdef USE_PART_STRUCT
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ud_t::z1init() {
  Json::Value z1info = opts["z1info"];
  // DUMMY_CASE==1 implies we are dealing with the dummy
  // case that is a full structured mesh and we treat as
  // partially strutured just to 
  dummy_case = z1info.isNull();
  if (dummy_case) {
    user_data_t::z1init();
    GETI(Nz);
    GETI(z1nz);
    GETI(k0start);
    GETD(Lz0);
    GETD(Lz1);
    PPRINTF("ud_t::z1init: Nz %d, z1nz %d, k0start %d, Lz0 %f, Lz1 %f\n",
            Nz,z1nz,k0start,Lz0,Lz1);
  } else {
    // We set the reference position of the FS as the ZINI
    // value. Perhaps we could have used also the mean
    // position of the actual FS (using the zcuts). 
    GETDD(zini);
    zref0 = zini;
    {
      Json::Value &opts=z1info;
      GETID(nlay);
      GETID(nz);
      GETD(Lz0);
      GETD(Lz1);
      z1nlay=nlay;
      z1nz=nz;
    }
    PPRINTF("partially structured buoy mesh params: "
            "z1nlay %d, z1nz %d, Lz0 %f, Lz1 %f\n",
            z1nlay,z1nz,Lz0,Lz1);
    string h5 = opts["mesh"].asString();
    h5 = "../../MESH/meshes/" + h5 + "/zstruc.h5";
    PPRINTF("Reading z-struct indices from %s\n",h5.c_str());
    zstrucindx.read(h5.c_str(),"/indxz/value");
    zstrucindx.reshape(z1nlay,z1nz);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int ud_t::get_gcid(int jlay,int k) {
  int rv;
  if (dummy_case) rv=jlay*Nz+k0start+k;
  else rv = zstrucindx(jlay,k);
  return rv;
}
#endif

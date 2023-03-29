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

extern "C"
void usstr1_(cs_lnum_t *idfstr,cs_real_t *aexxst,
             cs_real_t *bexxst,cs_real_t *cfopre,
             cs_real_t *xstr0 ,cs_real_t *vstr0,
             cs_real_t *xstreq) { }

extern "C"
void usstr2_(cs_lnum_t *nbstru,cs_lnum_t *idfstr,
             cs_real_t *dtcel, cs_real_t * xmstru ,
             cs_real_t *xcstru ,cs_real_t *xkstru ,
             cs_real_t *xstreq , cs_real_t *xstr ,
             cs_real_t *vstr, cs_real_t *forstr ,
             cs_real_t *dtstr) { }

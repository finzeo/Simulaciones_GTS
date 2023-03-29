/* GITVERSION<<undef>> */
/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

// DEPENDE DE QUE H=0.4508 Y QUE EL SISTEMA DE REFERENCIA SEA EL USADO POR GUILLERMO

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

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in]     nvar          total number of variable BC's
 * \param[in]     bc_type       boundary face types
 * \param[in]     icodcl        boundary face code
 *                                - 1  -> Dirichlet
 *                                - 2  -> convective outlet
 *                                - 3  -> flux density
 *                                - 4  -> sliding wall and u.n=0 (velocity)
 *                                - 5  -> friction and u.n=0 (velocity)
 *                                - 6  -> roughness and u.n=0 (velocity)
 *                                - 9  -> free inlet/outlet (velocity)
 *                                inflowing possibly blocked
 * \param[in]     rcodcl        boundary condition values
 *                                rcodcl(3) = flux density value
 *                                (negative for gain) in W/m2
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         bc_type[],
                            int         icodcl[],
                            cs_real_t   rcodcl[]) {

///// cs_real_t ux_inf = cs_notebook_parameter_value_by_name("Vel. del flujo libre Re2.7e4");

  cs_field_t *f;
  const int keyvar = cs_field_key_id("variable_id");
/////  cs_field_t *scalar = cs_field_by_name_try("velocity_x");
/////  int iscal = cs_field_get_key_int(scalar, keyvar) - 1; // access to u_x

 // f = cs_field_by_name("temperature");
 // int tvar = cs_field_get_key_int(f, keyvar);
  f = cs_field_by_name("velocity");
  int uvar = cs_field_get_key_int(f, keyvar);

  cs_lnum_t nfaces = cs_glob_mesh->n_b_faces;

#define BFNOR(fid,k) bfnor[3*fid+k]
#define BFCOG(fid,k) bfcog[3*fid+k]
#define ICODCL(fid,var) icodcl[((var)-1)*nfaces+(fid)]
#define RCODCL(fid,var,typ) rcodcl[((typ)-1)*nfaces*nvar+((var)-1)*nfaces+(fid)]

  cs_real_t *bfcog = cs_glob_mesh_quantities->b_face_cog;
  cs_real_t *bfnor = cs_glob_mesh_quantities->b_face_normal;
  double tol=1e-5;
  const unsigned int ndim=3;
  double h=0.4508;

  // Here we loop over all boundary faces. This includes
  // the external faces and the boundary faces that have been created in
  // the split in cs_internal_coupling_add_volume()
  for (cs_lnum_t fid = 0; fid < nfaces; fid++) {
/////    CORREGIR if ((BFCOG(fid,0)>=-7.5*h && BFCOG(fid,0)<=7.99*h) && (BFCOG(fid,1)<=3.0*h && BFCOG(fid,1)>=-3.0*h) && (BFCOG(fid,2)>=-0.01*h && BFCOG(fid,2)<=0.01*h))
/////	// Para ver que identificadores usar, ver https://www.code-saturne.org/documentation/7.1/doxygen/src/cs__parameters_8h.html
/////	bc_type[fid] = CS_SMOOTHWALL;
/////	icodcl[iscal*nfaces + fid] = 1;
/////   rcodcl[iscal*nfaces + fid] = ux_inf;
    if ((BFCOG(fid,2)>=-0.01*h && BFCOG(fid,2)<=0.01*h) && ((BFCOG(fid,0)<-7.5*h || BFCOG(fid,0)>7.99*h) || (BFCOG(fid,1)>3.0*h || BFCOG(fid,1)<-3.0*h)) )
	bc_type[fid] = CS_SYMMETRY;

/*
    // Set smooth wall condition in all faces 
    bc_type[fid] = CS_SMOOTHWALL;
    const double *x = &BFCOG(fid,0);
    // Set null flux Neumann condition by default
    // on all faces
    ICODCL(fid,tvar) = 3;
    RCODCL(fid,tvar,3) = 0.0;
    if (x[0]>1.0-tol || x[0]<tol) {
      // Set Dirichlet boundary conditions in the left and right boundaries
      // Set T=1 in the lower right part (x==, y<0.5), and T=0 on the rest
      // of those two faces.
      ICODCL(fid,tvar) = 1;
      RCODCL(fid,tvar,1) = x[0]*(x[1]<0.5);
    }
*/
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

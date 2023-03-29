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

ud_t::ud_t() : dummy_case(-1), k0start(-1), dzbox(NAN) { }

void ud_t::inituser() {
  // DATA initially is null so this can be checked as
  // a flag of whether the problem was initialized or not. 
  if (!opts) {
    // reads the JSON and  loads some internal fields
    init();
    GETD(gz);
    GETD(zini);

    // Read fields for this special case
    GETDP(ucirc,0);
    GETDP(gcirc,0);
    GETD(rho1);
    GETD(rho2);
    GETD(visco);
    GETDP(rhosolid,NAN);
    GETDP(delta,NAN);
    if (runcase == "circumwall") {
      GETD(R1);
      GETD(R2);
      GETD(Hwall);
      GETD(twall);
      cs_real_t Rmean = 0.5*(R1+R2);
      wcirc = ucirc/Rmean;
    }
    if (runcase == "dambreak") {
      GETD(H1);
      GETD(H2);
      GETD(Dphi);
      GETD(dphi);
    }
    if (runcase == "boxwall" || runcase == "wavegen") {
      GETD(Lx);
      GETD(Ly);
      GETD(Lz);
      GETD(Hwall);
      GETD(twall);
    }
    if (runcase=="wavegen") {
      GETD(Twave); 
      GETD(gwave);
      GETI(Nx);
      GETI(Nz);
      GETD(Lx);
      GETD(abso_length);

      wgen = 2.0*M_PI/Twave;
      kgen = sqr(wgen)/fabs(gz);
      // Wave length for the wave generator
      lagen = 2*M_PI/kgen;
      Labso = abso_length*lagen;
      if (cs_glob_rank_id <=0) {
        printf("wavegen init: Twave %f, gwave %f, wgen %f, kgen %f, "
               "lagen %f, abso_length %f, Labso %f\n",
               Twave,gwave,wgen,kgen,lagen,abso_length,Labso);
      }
    }
    if (runcase=="fsi-sphere") {
      int ndim=3;
      uxstr.resize(ndim,0.0);
      gx = 0.0;
    }

    int ndim=3;
    xsphere.resize(ndim,0.0);
    vsphere.resize(ndim,0.0);
    force.resize(ndim,0.0);
    if (runcase == "fallsphere") {
      GETD(Rsphere);
      for (int k=0; k<ndim; k++)
        xsphere[k] = user_data.opts.get("xsphere",NAN)[k].asDouble(); 
      cs_real_t hz;
      GETD(hz);
      bdymass = M_PI*sqr(Rsphere)*hz*rhosolid;
    }
    if (runcase == "pile") {
      GETD(wglen);
      GETD(dgz);
      GETD(T);
      GETD(Rpile);
      GETD(xpile);
      GETD(Ly);
      H = user_data.opts.get("zini",NAN).asDouble();
      if (cs_glob_rank_id<=0)
        printf("PILE CASE: wglen %f, dgz %f, T %f, Rpile %f, xpile %f\n",
               wglen,dgz,T,Rpile,xpile);
    }
    int GETIP(do_fsi,1);
    if (do_fsi) {
      GETDP(Tz,NAN);
      wz = 2.0*M_PI/Tz;
      GETD(Azbdy);
  
      GETSNR(fscase);

      A.resize(2,2);
      B.resize(2,2);
      GETDD(mass);
      GETDD(Cdamp);
      GETDD(Kspring);
      GETD(Dt);
      A(0,0) = 1;
      A(1,1) = mass;
      B(0,1) = -1;
      B(1,0) = Kspring;
      B(1,1) = Cdamp;
      GETD(theta);

      C = A+theta*Dt*B;
      y.resize(2);
      y.fill(0.0);
    }      
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ud_t::advance_struct(double fz) {
  Eigen::VectorXd rhs(2);
  rhs.fill(0.0);
  rhs(1) = fz;
  rhs -= B*y;
  y += Dt*C.lu().solve(rhs);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
cs_real_t ud_t
::issolid(cs_real_t *xx) {
  if (runcase=="square") {
    cs_real_t x = xx[0]-0.5, z=xx[2]-0.5;
    return fabs(x)<0.05 && fabs(z)<0.3;
  } else if (runcase=="circumwall") {
    return (xx[2]<Hwall) && (fabs(xx[1])<0.5*twall) && (xx[0]>0);
  } else if (runcase=="fallsphere") {
    cs_real_t r2 = sqr(xx[0]-xsphere[0])+sqr(xx[2]-xsphere[2]);
    return r2<sqr(Rsphere);
  } else if (runcase=="pile") {
    cs_real_t r2 = sqr(xx[0]-xpile)+sqr(xx[1]-0.5*Ly);
    return r2<sqr(Rpile);
  } else if (runcase=="boxwall") {
    return (xx[2]<Hwall) 
      && (fabs(xx[0]-0.5*Lx)<0.5*twall) 
      && (xx[1]>0.5*Ly);
    return 0;
  } else if (runcase=="wavegen") {
    // A full wave-length at the exit plane
    // with a ramp from 0 at Lx-lambda to 1 at the exit
    // XI goes from 0 to 1 at the outlet.
    
    cs_real_t 
      xi = (xx[0]-Lx+Labso)/Labso,
      Kfac=0.0;
    // if (xi>0.0) Kfac = 0.5*(1-cos(M_PI*xi));

    // The 3.0 factor is such that we have a normalized
    // integral of intensity. We will keep the mean of Kfac
    // to 1 and then try scaling with Kdarcy. So we have two
    // "variables", the "shape" Kfac(x) and the scale "Kdarcy". 
    if (xi>0.0) Kfac = 3.0*sqr(xi);
    return Kfac;
  }
  return 0;
}

cs_real_t ud_t::dzfun(double time) {
  return Azbdy*regheavis(time,0.5*Tz)*sin(wz*time);
}

void set_viz() {
  int kpost=cs_field_key_id("post_vis");
  int nf = cs_field_n_fields();
  for (int j=0; j<nf; j++) {
    cs_field_t *f = cs_field_by_id (j);
    int pv0 = cs_field_get_key_int(f,kpost);
    if (strcmp(f->name,"velocity")
        && strcmp(f->name,"void_fraction")
        && strcmp(f->name,"mesh_disp_rel")) {
      // PPRINTF("setting post_vis=0 for field %s\n",f->name);
      cs_field_set_key_int(f,kpost,0);
    }
    int pv1 = cs_field_get_key_int(f,kpost);
    if (pv1!=pv0) PPRINTF("%d %s post_vis before %d, after %d\n",j,f->name,pv0,pv1);
  }
}


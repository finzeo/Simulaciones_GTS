
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
  cs_parall_min(3,CS_DOUBLE,xmin.data());
  cs_parall_max(3,CS_DOUBLE,xmax.data());
  if (cs_glob_rank_id<=0) {
    printf("step %d, time %f, maxvel %f, total vol %g, vol(u>maxvel) %g, "
           "cell count %d\n",step,time,maxvel,voltot,vol,count);
    printf("bounding box (containing cells with vel>maxvel): xmin %g %g %g, xmax %g %g %g\m \n",
           xmin[0],xmin[1],xmin[2],xmax[0],xmax[1],xmax[2]);
  }
  //
  // Check yplus over cells adjacent to vehicle surface
  //
  // Puntero a cs_real_t (double), que apunta al conjunto de valores de la variable yplus en todas las caras (yplus es un campo en caras)
  cs_real_t *yplusvals = cs_field_by_name("yplus")->val;
  // No es necesario definir un macro para manipular mejor el conjunto de valores de yplus, dado que es un campo escalar (a diferencia de, por ejemplo, la velocidad o la fuerza)
  // cs_lnum_t es un tipo asociado a enteros con signo, específicamente usado para numerar entidades de la malla (en este caso, para almacenar el número de celdas en la frontera)
  cs_lnum_t nbfaces = cs_glob_mesh->n_b_faces;
  // (-15 -2 -0.2) (1 2 4) // bbox asociada al vehículo
  // Definir un puntero a un vector que contiene las coordenadas de los centros de gravedad (Center Of Gravity; COG) de cada cara en la frontera (se manipula mejor llevándolo a una matriz, con el macro BFCOG, de forma que cada fila contiene las coordenadas de cierta cara)
  const cs_real_t *bfcog = fvq->b_face_cog;
  #define BFCOG(fid,k) bfcog[ndim*(fid)+(k)]
  // Puntero a vector que contiene las áreas de cada cara de la frontera (usado para calcular yplusprom)
  const cs_real_t *bfsurf = fvq->b_face_surf;
  // Vector de enteros con signo de tamaño igual al número de caras totales, que contendrá útilmente solo a los ids de las caras en la frontera deseada (el resto del vector quedará con valores basura)
  vector<cs_lnum_t> face_list(nbfaces);
  // Se guardará el número de caras sobre la frontera deseada
  cs_lnum_t nfcext;
  // Se guarda la lista de las caras sobre la frontera deseada (indicar el nombre de la frontera deseada, según se cargó en setup.xml de la carpeta DATA de CS; NO el label)
  // Además, se guarda el número de caras sobre la frontera deseada en nfcext
  cs_selector_get_b_face_list("gtsModel",&nfcext,face_list.data());

  // Chequear coordenadas mínimas y máximas de las caras sobre la frontera deseada (para ver si las toma correctamente)
  cs_real_t xmin2[ndim], xmax2[ndim], xc2[ndim];
  cs_real_t yplusmin, yplusmax, yplusprom, gts_area;
  cs_real_t yplusmax_cog[ndim];
  yplusmin = INFINITY;
  yplusmax = -INFINITY;
  yplusprom = 0;
  gts_area = 0;
  for (int i=0;i<ndim;i++) yplusmax_cog[i]=-INFINITY;
  for (int j=0; j<ndim; j++) {
    xmin2[j] = INFINITY;
    xmax2[j] = -INFINITY;
  }
  for (cs_lnum_t j = 0; j<nfcext; j++) {
    cs_lnum_t fid = face_list[j];
    for (int k=0; k<ndim; k++) {
      xmin2[k] = min(xmin2[k],BFCOG(fid,k));
      xmax2[k] = max(xmax2[k],BFCOG(fid,k));
    }
    // Puedo acceder al yplus asociado a cada cara con yplusvals[fid]
    yplusmin = min(yplusmin,yplusvals[fid]);
    //yplusmax = max(yplusmax,yplusvals[fid]);
    if(yplusvals[fid]>yplusmax)
    {
	yplusmax = yplusvals[fid];
        for(int j=0;j<ndim;j++) yplusmax_cog[j]=BFCOG(fid,j);
    }
    yplusprom += yplusvals[fid]*bfsurf[fid];
    gts_area += bfsurf[fid];
  }
  cs_parall_min(ndim,CS_DOUBLE,&xmin2);
  cs_parall_max(ndim,CS_DOUBLE,&xmax2);
  cs_parall_min(1,CS_DOUBLE,&yplusmin);
  cs_parall_max_loc_vals(ndim,&yplusmax,&yplusmax_cog[0]);
  cs_parall_sum(1,CS_DOUBLE,&yplusprom);
  cs_parall_sum(1,CS_DOUBLE,&gts_area);
  yplusprom = yplusprom/gts_area;
  if (cs_glob_rank_id<=0)
  {
     printf("bounding box of vehicle: xmin %g %g %g, xmax %g %g %g\m \n",
        xmin2[0],xmin2[1],xmin2[2],xmax2[0],xmax2[1],xmax2[2]);
     printf("yplus max: %g, yplus min: %g \n",yplusmax,yplusmin);
     printf("position of yplus max: %g %g %g \n",yplusmax_cog[0],yplusmax_cog[1],yplusmax_cog[2]);
     printf("yplus prom: %g \n", yplusprom);
  }

  //
  // Calculate net force applied over the vehicle surface
  // (There is a field named "boundary_forces" that does not appear in the GUI, but it exists if stress calculation on postprocessing actions is activated)
  //
  // Se define un puntero a cs_real_t (double), el cual apuntará a un arreglo que contiene a los valores de las fuerzas en cada cara de toda la frontera
  // (a priori, es un arreglo unidimensional; se usa un macro para tratarlo como matriz, convenientemente).
  cs_real_t *forces = cs_field_by_name("boundary_forces")->val;
  #define FORCES(fid,k) forces[(fid)*ndim+(k)]
  // Se usa el vector face_list para tener la lista de ids de las caras sobre la frontera deseada (superficie del vehículo), ya definido en el bloque de código
  // anterior para cálculos de yplus. Tampoco se vuelve a definir nfcext (número de caras sobre la frontera deseada)
  // Vector en el que se acumularán todas las fuerzas ejercidas sobre cada cara de la frontera deseada, para obtener la fuerza neta sobre toda la frontera deseada
  vector<cs_real_t> force(ndim,0);
  for (cs_lnum_t j = 0; j<nfcext; j++) {
    cs_lnum_t fid2 = face_list[j];
    for (int k=0; k<ndim; k++) {
      force[k] += FORCES(fid2,k);
    }
  }
  cs_parall_sum(ndim,CS_DOUBLE,force.data());
  if (cs_glob_rank_id<=0)
  {
    printf("step %d, time %f, forces %f %f %f \n",
           step,time,force[0],force[1],force[2]);
    printf("----------------------------------------------------- \n");
  }
}

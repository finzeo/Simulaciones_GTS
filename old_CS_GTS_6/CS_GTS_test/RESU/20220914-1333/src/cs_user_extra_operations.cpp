
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
#include <cs_physical_constants.h> // For fluid properties
#include <cs_turbulence_model.h> // For turbulence properties
#include <fstream>
#include "userdata.hpp"

using namespace std;

#include "vecmacros2.hpp"

using namespace cs_user_funs;

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

  double maxvel=3.0;
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
  // Puntero a cs_real_t (double), que apunta al conjunto de valores de la variable yplus en todas las caras de toda la frontera (yplus es un campo en caras)
  cs_real_t *yplusvals = cs_field_by_name("yplus")->val;
  // No es necesario definir un macro para manipular mejor el conjunto de valores de yplus, dado que es un campo escalar (a diferencia de, por ejemplo, la velocidad o la fuerza)
  // cs_lnum_t es un tipo asociado a enteros sin signo, específicamente usado para numerar entidades de la malla (en este caso, para almacenar el número de celdas en la frontera)
  cs_lnum_t nbfaces = cs_glob_mesh->n_b_faces;
  // (-15 -2 -0.2) (1 2 4) // bbox asociada al vehículo
  // Definir un puntero a un vector que contiene las coordenadas de los centros de gravedad (Center Of Gravity; COG) de cada cara en la frontera (se manipula mejor llevándolo a una matriz, con el macro BFCOG, de forma que cada fila contiene las coordenadas de cierta cara)
  const cs_real_t *bfcog = fvq->b_face_cog;
  #define BFCOG(fid,k) bfcog[ndim*(fid)+(k)]
  // Puntero a vector que contiene las áreas de cada cara de la frontera (usado para calcular yplusprom)
  const cs_real_t *bfsurf = fvq->b_face_surf;
  // Vector de enteros con signo de tamaño igual al número de caras totales en toda la frontera, que contendrá útilmente solo a los ids de las caras en la frontera deseada (el resto del vector quedará con valores basura)
  vector<cs_lnum_t> face_list(nbfaces);
  // Se guardará el número de caras sobre la frontera deseada
  cs_lnum_t nfcext;
  // Se guarda la lista de las caras sobre la frontera deseada (indicar el nombre de la frontera deseada, según se cargó en setup.xml de la carpeta DATA de CS; NO el label)
  // Además, se guarda el número de caras sobre la frontera deseada en nfcext
  cs_selector_get_b_face_list("GTS",&nfcext,face_list.data());
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
     printf("gts area: %g \n",gts_area); // Para verificar que coincide con el área del GTS que calcula el módulo de geometría de Salome
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
    printf("forces: %f %f %f \n",
           force[0],force[1],force[2]);
  }

  // Cálculo de coeficiente de arrastre
  if (cs_glob_rank_id<=0)
  {
    // Coeficiente de arrastre
    cs_real_t c_drag;
    // Densidad  
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
    cs_real_t density_fluid = fp->ro0;
    // Velocidad de referencia considerada en modelo de turbulencia (es la velocidad a la entrada)
    cs_turb_ref_values_t *turb_ref_values = cs_get_glob_turb_ref_values();
    cs_real_t u_ref = turb_ref_values->uref;
    // Área frontal (calculada a partir de la bbox del vehículo)
    cs_real_t frontal_area = (xmax2[1]-xmin2[1])*(xmax2[2]-xmin2[2]);
    c_drag = force[0]/(density_fluid*u_ref*u_ref*frontal_area/2);
    printf("Drag coefficient: %f \n", c_drag);
  }
  
  
  // Cálculo del torque neto asociado al conjunto de las fuerzas aplicadas sobre cada cara de toda la frontera, respecto a un punto de referencia
  // (a priori, el CG)
  vector<cs_real_t> torque(ndim,0);
  vector<cs_real_t> ref_point(ndim);
  ref_point[0] = -1.2374;
  ref_point[1] =  0.0;
  ref_point[2] =  0.2384; // Valor colocado provisoriamente en prtclsys
  // Se usa macro BFCOG(fid,k) para acceder a los centros de gravedad de las caras de la frontera, ya definido en el bloque de código asociado al cálculo de yplus
  // Producto vectorial: lo hago sin eigen por dificultades asociadas al uso de la librería en el programa
  // r = (x  y  z )
  // F = (fx fy fz)
  // r^F = (y*fz-fy*z)i + (fx*z-x*fz)j + (x*fy-fx*y)k
  for (cs_lnum_t j = 0; j<nfcext; j++) {
    cs_lnum_t fid3 = face_list[j];
    cs_real_t x  = ref_point[0];
    cs_real_t y  = ref_point[1];
    cs_real_t z  = ref_point[2];
    cs_real_t fx = FORCES(fid3,0);
    cs_real_t fy = FORCES(fid3,1);
    cs_real_t fz = FORCES(fid3,2);
    torque[0] += y*fz-fy*z;
    torque[1] += fx*z-x*fz;
    torque[2] += x*fy-fx*y;
  }
  cs_parall_sum(ndim,CS_DOUBLE,torque.data());
  if (cs_glob_rank_id<=0)
  {
    printf("torque (ref. CG): %f (roll) %f (pitch) %f (yaw) \n",
           torque[0],torque[1],torque[2]);
    printf("----------------------------------------------------- \n");
  }
/*
  // Test torque calculation
    if (cs_glob_rank_id<=0)
  {
    vector<cs_real_t> torque_test(ndim,0);
    cs_real_t x  = 1;
    cs_real_t y  = 0;
    cs_real_t z  = 0;
    cs_real_t fx = 0;
    cs_real_t fy = 0;
    cs_real_t fz = 1;
    torque_test[0] += y*fz-fy*z;
    torque_test[1] += fx*z-x*fz;
    torque_test[2] += x*fy-fx*y;
    printf("torque_test %f %f %f \n",
           torque_test[0],torque_test[1],torque_test[2]);
    printf("----------------------------------------------------- \n");  
  }
*/

  // Check mesh symmetry
  // 1. Check number of cells on each side of the vertical midplane
  if (cs_glob_time_step->nt_cur==1)
  {
     ofstream out("checkMeshSymmetry.txt");
     cs_real_t *xcell = cs_glob_mesh_quantities->cell_cen;
     #define BCOG(fid,k) xcell[ndim*(fid)+(k)]
     cs_lnum_t nc_izq=0;
     cs_lnum_t nc_der=0;
     cs_lnum_t num_cells = cs_glob_mesh->n_cells;
     double tol_ccells = 1e-7;
     for (cs_lnum_t i=0; i<num_cells; i++)
     {
        if (BCOG(i,1)>tol_ccells) {nc_der++;}
        else if (BCOG(i,1)<(tol_ccells*(-1))) {nc_izq++;}
     }
     cs_parall_sum(1,CS_INT32,&nc_der);
     cs_parall_sum(1,CS_INT32,&nc_izq);
     cs_parall_sum(1,CS_INT32,&num_cells);
     if (cs_glob_rank_id<=0)
     {
        out << "Número total de celdas: " << num_cells << "\n";
        out << "Celdas sobre semiespacio izquierdo: " << nc_izq << "\n";
        out << "Celdas sobre semiespacio derecho: " << nc_der << "\n";
        if (nc_izq==nc_der) out << "OK \n";
     }
  }

  //
  // Manipulate an user-defined field
  // (store values in an user-defined field, as written in cs_user_parameters.cpp)
  set_viz();
  // CHRONO.elaps("extra ops start");
  // user_data.zrenorm();
  user_data.check_probes();
  // CHRONO.elaps("extra ops end");
  ud_t &d = user_data;
  d.inituser();
  Json::Value &opts = d.opts;
  cs_lnum_t nvrtx = cs_glob_mesh->n_vertices;
  cs_real_t *velx1 = cs_field_by_name("velocity")->val;
  cs_real_t *velx2 = cs_field_by_name("velx2")->val;
#define VELX1(v,k) velx1[3*v+k]
#define VELX2(v,k) velx2[3*v+k]
  for (int j=0; j<nvrtx; j++) {
    VELX2(j,0) = 2.0*VELX1(j,0);
    VELX2(j,1) = 2.0*VELX1(j,1);
    VELX2(j,2) = 2.0*VELX1(j,2);
  }
}

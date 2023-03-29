!-------------------------------------------------------------------------------

!                      Code_Saturne version 5.0.8-patch
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization-base.f90
!>
!> \brief Basic examples
!>
!> See \subpage cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use vof

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, iutile
integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:), pointer :: voidf
double precision x,z,zfs,dz

!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [alloc]
allocate(lstelt(ncel)) ! temporary array for cells selection
!< [alloc]

!===============================================================================
! Variables initialization:
!
!   isca(1) is the number related to the first user-defined scalar variable.
!   cvar_scal(iel) is the value of this variable in cell number iel.
!
!   ONLY done if there is no restart computation
!===============================================================================

!< [init]

call field_get_val_s(ivarfl(ivolf2), voidf)
if (.false.) then
   do iel = 1, ncel
      voidf(iel) = 1
      if( xyzcen(1,iel).lt. 0.1461) then
         if(xyzcen(3,iel).lt.0.292) then
            voidf(iel) = 0
         endif
      endif
   enddo
else if (.false.) then
   !! FS with a slope
   !! dz = 0.2
   dz = 0.0
   do iel = 1, ncel
      voidf(iel) = 1
      x = xyzcen(1,iel)
      z = xyzcen(3,iel)
      zfs = 0.3 + dz*(x/0.3-1.0)
      if(z .lt. zfs) then
         voidf(iel) = 0
      endif
   enddo
endif

!! This allows to set the values
rho1 = 1.d3
rho2 = 20.0
mu1 = 1.d-3
mu2 = 1.d-5
call vof_init_vals(rho1,rho2,mu1,mu2)

!< [init]

!--------
! Formats
!--------

!----
! End
!----

!< [finalize]
deallocate(lstelt)  ! temporary array for cells selection
!< [finalize]

return
end subroutine cs_user_f_initialization

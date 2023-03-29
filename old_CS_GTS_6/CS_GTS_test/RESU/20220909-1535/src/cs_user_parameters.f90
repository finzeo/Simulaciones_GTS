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

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------


!===============================================================================

!> \brief User subroutine for input of model selection parameters.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      ixmlpu       indicates if the XML file from the GUI is used
!>                              used (1: yes, 0: no
!> \param[in, out] iturb        turbulence model
!> \param[in, out] itherm       thermal model
!> \param[in, out] iale         ale module
!> \param[in, out] ivofmt       vof method
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

subroutine usipph &
 ( ixmlpu, iturb , itherm, iale , ivofmt, icavit )

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, ivofmt, icavit

! Local variables

!===============================================================================

!>    In this subroutine, only the parameters which already appear may
!>    be set, to the exclusion of any other.
!>
!>    If we are not using the Code_Saturne GUI:
!>    All the parameters which appear in this subroutine must be set.
!>
!>    If we are using the Code_Saturne GUI:
!>    parameters protected by a test of the form:
!>
!>      if (ixmlpu.eq.0) then
!>         ...
!>      endif
!>
!>    should already have been defined using the GUI, so only
!>    experts should consider removing the test and adapting them here.

!===============================================================================

!< [usipph]

! --- Turbulence
!       0: Laminar
!      10: Mixing length
!      20: k-epsilon
!      21: k-epsilon (linear production)
!      30: Rij-epsilon, (standard LRR)
!      31: Rij-epsilon (SSG)
!      32: Rij-epsilon (EBRSM)
!      40: LES (Smagorinsky)
!      41: LES (Dynamic)
!      42: LES (WALE)
!      50: v2f (phi-model)
!      51: v2f (BL-v2/k)
!      60: k-omega SST
!      70: Spalart Allmaras
!  For 10, contact the development team before use

if (ixmlpu.eq.0) then

  iturb = 0

endif

!IVOFMT =              0 ( -1: inactive                )
!                        (  0: active                  )

ivofmt = 0

!< [usipph]

!----
! Formats
!----


return
end subroutine usipph


!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

logical       inoprv
integer       ii, jj, ivar, kscmin, kscmax, keydri, kbfid, kccmin, kccmax
integer       klimiter
integer       f_id, idim1, itycat, ityloc, iscdri, iscal, ifcvsl, b_f_id

type(var_cal_opt) :: vcopt

!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================

!< [usipsu]

! Calculation options (optcal)
! ============================


! --- Convective scheme for user (and non-user) scalars

! ischcv is the type of convective scheme:
!   - 0: second order linear upwind
!   - 1: centered
!   - 2: pure upwind gradient in SOLU

! isstpc is the slope test, Min/Max limiter or Roe and Sweby limiters
!   - 0: swich on the slope test
!   - 1: swich off the slope test (default)
!   - 2: continuous limiter ensuring boundedness (beta limiter)
!   - 3: NVD/TVD Scheme
!        Then "limiter_choice" keyword must be set:
!        * 0: Gamma
!        * 1: SMART
!        * 2: CUBISTA
!        * 3: SUPERBEE
!        * 4: MUSCL
!        * 5: MINMOD
!        * 6: CLAM
!        * 7: STOIC
!        * 8: OSHER
!        * 9: WASEB
!        * --- VOF scheme ---
!        * 10: M-HRIC
!        * 11: M-CICSAM

! Get the Key for the Sup and Inf for the convective scheme
call field_get_key_id("min_scalar", kccmin)
call field_get_key_id("max_scalar", kccmax)

if(.true.) then

call field_get_key_struct_var_cal_opt(ivarfl(ivolf2), vcopt)
vcopt%ischcv = 0
vcopt%isstpc = 3
vcopt%blencv = 1

call field_set_key_struct_var_cal_opt(ivarfl(ivolf2), vcopt)

! Get the Key for the limiter choice of the studied scalar
call field_get_key_id("limiter_choice", klimiter)
call field_set_key_int(ivarfl(ivolf2), klimiter, 11)
  
! Set the Value for the Sup and Inf of the studied scalar
call field_set_key_double(ivarfl(ivolf2), kccmin, 0.d0)
call field_set_key_double(ivarfl(ivolf2), kccmax, 1.d0)

endif
!< [usipsu]

!----
! Formats
!----

return
end subroutine usipsu



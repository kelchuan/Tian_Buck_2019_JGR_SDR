!-*- F90 -*-

! --------- Flac ------------------------- 

subroutine flac

use arrays
include 'precision.inc' 
include 'params.inc'
include 'arrays.inc'

! Update Thermal State
! Skip the therm calculations if itherm = 3
if (itherm .eq. 3) goto 222
call fl_therm
222 continue

if (itherm .eq.2) goto 500  ! Thermal calculation only

! Calculation of strain rates from velocity
call fl_srate

! Changing marker phases
! XXX: change_phase is slow, don't call it every loop
!if( mod(nloop, 10).eq.0 ) call change_phase

! Update stresses by constitutive law (and mix isotropic stresses)
call fl_rheol

! update stress boundary conditions
if (ynstressbc.eq.1.) call bc_update

! Calculations in a node: forces, balance, velocities, new coordinates
call fl_node

! New coordinates
call fl_move

! Heat Dike Injection     Tian 2017 adapting from Behn and Ito 2008
!if (ny_inject.gt.0 .and. iynts.ne.1) then
if (ny_inject.gt.0) then
  call fl_injectheat
endif

! Adjust real masses due to temperature
if( mod(nloop,ifreq_rmasses).eq.0 ) call rmasses

! Adjust inertial masses or time step due to deformations
if( mod(nloop,ifreq_imasses) .eq. 0 ) call dt_mass

! Adjust time Step 
call dt_adjust

500 continue


return
end subroutine flac

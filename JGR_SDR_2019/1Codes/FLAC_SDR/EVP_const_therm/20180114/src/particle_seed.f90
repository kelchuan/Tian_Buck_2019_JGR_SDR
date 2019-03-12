!======================================================================
  subroutine particle_seed
! Initial locations of surface particles
! G. Ito 3/07
!======================================================================
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
!integer SDR_ind
!SDR_ind = mod(int(time*0.5*3.171d-8*1.d-6),10)+1  !Tian every 2Myr increment the index by 1

dxp=rxbo
do i = 1,nzonx    !Tian find the minimum distance (dxp) between particles 
  dxzone=rxbo*sizez_x(i)/nelz_x(i)/npelem
  if (dxzone.lt.dxp) dxp=dxzone
enddo
np=rxbo/dxp   !Tian total number of particles (more than npelem * (nx - 1) if variable horizontal sizes of elements
  
if (mod(np,2).ne.0) then !Tian make sure np is an even number
  np=np+1
endif

if (np.gt.mnp) then
  write(*,*) 'Too many particles, np > mnp :', np, mnp
endif

write(*,*) '****Surface Particles for RIDGE models: np=',np

!dxp=rxbo/dble(np)
dxp=rxbo/dble(np)/2 !Tian only first half has particles to save calc

!do 100 i=1,np/2
!   xp(i)=x0+dxp*(i-1)
!!   xp(np-i+1)=x0+rxbo-dxp*(i-1)
!   xp(np-i+1)=x0+rxbo/2.-dxp*(i-1) !Tian only first half has particles to save calc
!100 continue
do 100 i=1,np
   xp(i,:)=dble(x0)+dble(npelem)*dble(dxp)*2*2+dble(dxp)*dble(i-1)   !(+npelem*dxp to shift off axis &
   !  0.5 element width, *2*2, is to shift off axis 2 elements width, &
   ! so that viscous widening dike will not affect the leftmost particles)
   !because fortran store multi-dim arrays in column-major order, different period of
   !SDR layers should be store as a column not a row to improve efficiency
   !meaing i should be at row
   !ref: https://en.wikipedia.org/wiki/Row-_and_column-major_order
!   xp(np-i+1)=x0+rxbo-dxp*(i-1)
!   xp(np-i+1)=x0+rxbo/2.-dxp*(i-1) !Tian only first half has particles to save calc
   !   zp(:,:) = 0 - 0.5 * npelem*dxp ! all zp start at mid poit of the surface element
!   zp(i,:) = -200.0d0 ! all zp start at mid poit of the surface element
   zp(i,:) = -1.0d0 ! 
   !   write(*,*) 'zp  ', i, '=', zp(i,1)
!   write(*,*) 'xp(', i,') = ', xp(i,1)
   !   write(*,*) 'zp(', i,') = ', zp(i,1)
   !   write(*,*) xp(i,1), zp(i,1)
!   write(*,*) "dble(npelem)*dble(dxp)=", dble(npelem)*dble(dxp)
!   write(*,*) "dble(dxp)=", dble(dxp)
100 continue


!do 210 i=1,np
!  do 200 ii=1,nx-1
!    xl=cord(1,ii,1)
!    xr=cord(1,ii+1,1)
!    if (xp(i,:).ge.xl.and.xp(i,:).le.xr) then
!       !      zl=cord(1,ii,2) 
!       !      zr=cord(1,ii+1,2)
!       zl=cord(1,ii,2) - 0.5 * npelem*dxp  ! to 
!       zr=cord(1,ii+1,2) - 0.5 * npelem*dxp
!       !      zp(i)=zl + ((zr-zl)/(xr-xl)) * (xp(i)-xl)  !linearly interpolating zp(i) based on its xp(i) distance from xl
!       zp(i,:)=zl + ((zr-zl)/(xr-xl)) * (xp(i,:)-xl)  !linearly interpolating zp(i) based on its xp(i) distance from xl 
!       exit
!    endif
!  200 continue
!210  continue

return
end

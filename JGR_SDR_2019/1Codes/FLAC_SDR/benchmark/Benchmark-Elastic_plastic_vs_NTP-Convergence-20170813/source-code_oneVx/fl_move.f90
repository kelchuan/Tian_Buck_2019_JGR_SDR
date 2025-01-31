
! Move grid and adjust stresses due to rotation

subroutine fl_move
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

! one stretch test begins Tian 20180103
sec_year = 3.1558e+7
if (time/sec_year/1.e+6 > 0.1) then  ! after 0.1 Myr vx = 0
   vel(:,nx,1) = 0.0
endif
! one stretch test ends Tian 20180103

! Move Grid
if (movegrid .eq. 0) return

! UPDATING COORDINATES

!$OMP parallel
!$OMP do
do i = 1,nx
!    write(*,*) cord(j,i,1),cord(j,i,2),vel(j,i,1),vel(j,i,2),dt
    cord(:,i,1) = cord(:,i,1) + vel(:,i,1)*dt
    cord(:,i,2) = cord(:,i,2) + vel(:,i,2)*dt
!    write(*,*) cord(j,i,1),cord(j,i,2)
enddo
!$OMP end do
!$OMP end parallel

! Diffuse topography
if( topo_kappa.gt.0.) call diff_topo


!$OMP parallel private(i,j,x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  det,dw12,s11,s22,s12)
!$OMP do
!--- Adjusting Stresses And Updating Areas Of Elements
do  i = 1,nx-1
    do  j = 1,nz-1

        ! Coordinates
        x1 = cord (j  ,i  ,1)
        y1 = cord (j  ,i  ,2)
        x2 = cord (j+1,i  ,1)
        y2 = cord (j+1,i  ,2)
        x3 = cord (j  ,i+1,1)
        y3 = cord (j  ,i+1,2)
        x4 = cord (j+1,i+1,1)
        y4 = cord (j+1,i+1,2)

        ! Velocities
        vx1 = vel (j  ,i  ,1)
        vy1 = vel (j  ,i  ,2)
        vx2 = vel (j+1,i  ,1)
        vy2 = vel (j+1,i  ,2)
        vx3 = vel (j  ,i+1,1)
        vy3 = vel (j  ,i+1,2)
        vx4 = vel (j+1,i+1,1)
        vy4 = vel (j+1,i+1,2)

        ! (1) Element A:
        det=((x2*y3-y2*x3)-(x1*y3-y1*x3)+(x1*y2-y1*x2))
        dvol(j,i,1) = det*area(j,i,1) - 1
        area(j,i,1) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x2)+vx2*(x1-x3)+vx3*(x2-x1) - &
            vy1*(y2-y3)-vy2*(y3-y1)-vy3*(y1-y2))/det*dt
        s11 = stress0(j,i,1,1)
        s22 = stress0(j,i,2,1)
        s12 = stress0(j,i,3,1)
        stress0(j,i,1,1) = s11 + s12*2.*dw12
        stress0(j,i,2,1) = s22 - s12*2.*dw12
        stress0(j,i,3,1) = s12 + dw12*(s22-s11)

        ! rotate strains 
        s11 = strain(j,i,1)
        s22 = strain(j,i,2)
        s12 = strain(j,i,3)
        strain(j,i,1) = s11 + s12*2.*dw12
        strain(j,i,2) = s22 - s12*2.*dw12
        strain(j,i,3) = s12 + dw12*(s22-s11)

        ! (2) Element B:
        det=((x2*y4-y2*x4)-(x3*y4-y3*x4)+(x3*y2-y3*x2))
        dvol(j,i,2) = det*area(j,i,2) - 1
        area(j,i,2) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx3*(x4-x2)+vx2*(x3-x4)+vx4*(x2-x3) - &
           vy3*(y2-y4)-vy2*(y4-y3)-vy4*(y3-y2))/det*dt
        s11 = stress0(j,i,1,2)
        s22 = stress0(j,i,2,2)
        s12 = stress0(j,i,3,2)
        stress0(j,i,1,2) = s11 + s12*2.*dw12
        stress0(j,i,2,2) = s22 - s12*2.*dw12
        stress0(j,i,3,2) = s12 + dw12*(s22-s11)

        ! (3) Element C:
        det=((x2*y4-y2*x4)-(x1*y4-y1*x4)+(x1*y2-y1*x2))
        dvol(j,i,3) = det*area(j,i,3) - 1
        area(j,i,3) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x4-x2)+vx2*(x1-x4)+vx4*(x2-x1) - &
           vy1*(y2-y4)-vy2*(y4-y1)-vy4*(y1-y2))/det*dt
        s11 = stress0(j,i,1,3)
        s22 = stress0(j,i,2,3)
        s12 = stress0(j,i,3,3)
        stress0(j,i,1,3) = s11 + s12*2.*dw12
        stress0(j,i,2,3) = s22 - s12*2.*dw12
        stress0(j,i,3,3) = s12 + dw12*(s22-s11)

        ! (4) Element D:
        det=((x4*y3-y4*x3)-(x1*y3-y1*x3)+(x1*y4-y1*x4))
        dvol(j,i,4) = det*area(j,i,4) - 1
        area(j,i,4) = 1./det

        ! Adjusting stresses due to rotation
        dw12 = 0.5*(vx1*(x3-x4)+vx4*(x1-x3)+vx3*(x4-x1) - &
            vy1*(y4-y3)-vy4*(y3-y1)-vy3*(y1-y4))/det*dt
        s11 = stress0(j,i,1,4)
        s22 = stress0(j,i,2,4)
        s12 = stress0(j,i,3,4)
        stress0(j,i,1,4) = s11 + s12*2.*dw12
        stress0(j,i,2,4) = s22 - s12*2.*dw12
        stress0(j,i,3,4) = s12 + dw12*(s22-s11)
    enddo
enddo
!$OMP end do
!$OMP end parallel
return
end subroutine fl_move


!============================================================
! Diffuse topography
!============================================================
subroutine diff_topo
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dimension dh(mnx+1)
real(8) :: infill_level

!EROSION PROCESSES
if( topo_kappa .gt. 0. ) then             !Tian 20170731 topo_kappa has to be larger than zero to allow sedimentation
    do i = 2, nx-1
        water_depth = 0.5*(cord(1,i+1,2)+cord(1,i,2))
        if (water_depth.lt.0) then
!          topo_kappa2 = topo_kappa/10
          topo_kappa2 = topo_kappa/100 !(Tian201607: change from 10 to 100)
        else
          topo_kappa2 = topo_kappa
        endif
        snder = ( (cord(1,i+1,2)-cord(1,i  ,2))/(cord(1,i+1,1)-cord(1,i  ,1)) - &
            (cord(1,i  ,2)-cord(1,i-1,2))/(cord(1,i  ,1)-cord(1,i-1,1)) ) / &
            (cord(1,i+1,1)-cord(1,i-1,1))
        dh(i) = topo_kappa2 * dt * snder
 !       print *, 'erosion activated dh = ', dh(i)
    end do

    dh(1) = dh(2)
    dh(nx) = dh(nx-1)
    cord(1,1:nx,2) = cord(1,1:nx,2) + dh(1:nx)

    ! accumulated topo change since last resurface
    dhacc(1:nx-1) = dhacc(1:nx-1) + 0.5 * (dh(1:nx-1) + dh(2:nx))
!    print *, 'dhacc = ', dhacc(1:nx-1)
    ! Forced sedimentation to a prescribed infill level
!    infill_level = -1250.0d0

    infill_level = 0.0d0 !(Tian: try an infill level)
    if (time*3.171d-8*1.d-6 > 6 .and. time*3.171d-8*1.d-6 < 8) then !(Tian: stop infill after 5Myr
       infill_level = 0.0d0 !(Tian: try an infill level)
    else
       infill_level = 0.0d0
!       infill_level = -10000.0d0   !Tian for analytic benchmark
    endif !(Tian: try an infill level)
    do i = 1, nx-1
       if (cord(1,i,2) < infill_level)then
          dhacc(i) = dhacc(i) + (infill_level - cord(1,i,2))
          cord(1,i,2) = infill_level    ! Tian 20170731 This line causes hardwired flat surface
       endif
    end do

    ! adjust markers
    if(mod(nloop, 100) .eq. 0) then
        !print *, 'max sed/erosion rate (m/yr):' & !(Tian: uncommented)
        !     , maxval(dh(1:nx)) * 3.16e7 / dt & !(Tian: uncommented) 
        !     , minval(dh(1:nx)) * 3.16e7 / dt !(Tian: uncommented) 
        call resurface
 !       print *, 'resurface called' !(Tian)
    end if
endif


return
end subroutine diff_topo




subroutine resurface
  use marker_data
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  include 'phases.inc'

  dimension shp2(2,3,2)

  do i = 1, nx-1
      call shape_functions(1,i,shp2)

      ! add/remove markers if topo changed too much
      surface = 0.5 * (cord(1,i,2) + cord(1,i+1,2))
      elz = surface - 0.5 * (cord(2,i,2) + cord(2,i+1,2))
      diff = dhacc(i)
!      print *, 'add sediment diff, kinc, diff*kinc, elz', diff, kinc, diff*kinc, elz !(Tian1607)
      kinc = sum(nphase_counter(:,1,i))
!      if (diff*kinc .ge. elz) then
      if (diff .ge. 100) then !Tian 20170731 whenever deflection >= 100 meters / (0.1 vertical element size, add marker(
          ! sedimentation, add a sediment marker
!          print *, 'add sediment', i, diff, elz !(Tian1607 uncommented)
          do while (.true.)
              call random_number(rx)
              xx = cord(1,i,1) + rx * (cord(1,i+1,1) - cord(1,i,1))
              yy = min(cord(1,i,2), cord(1,i+1,2)) - 0.05 * elz
              if (xx .le. 30000) then !Tian 20170718 added for restricted lava in pure elastic case
              !! time*3.171d-8*1d-5 = time in 0.1 Myrs: 1 = 0.1 Myrs, 2 = 0.2 Myrs, 10 = 1.0 Myrs, etc.
 !             call add_marker(xx, yy, ksed1+mod(int(time*3.171d-8*1.d-5),10), time, nmarkers, 1, i, inc) !(Tian CM)
              call add_marker(xx, yy, 8, time, nmarkers, 1, i, inc)!(Tian1607)
!              call add_marker(xx, yy, 7 + mod(int(time*3.171d-8*1.d-6),10), time, nmarkers, 1, i, inc)!(Tian1607)
!              mod(int(time*3.171d-8*1.d-6),10) each Myr will increase phase by one (Tian1607)
!             call add_marker(xx, yy, 8 + mod(int(time*0.5*3.171d-8*1.d-6),2), time, nmarkers, 1, i, inc)!(Tian1607)
!              mod(int(time*3.171d-8*1.d-6),10) each two Myr alternate phases between 8 and 9 (Tian1607)
!!                 if (time*3.171d-8*1.d-6 > 6 .and. time*3.171d-8*1.d-6 < 10) then !(Tian: stop infill after 5Myr
!!                    call add_marker(xx, yy, 11 + mod(int(time*0.1*3.171d-8*1.d-5),2), time, nmarkers, 1, i, inc)!(Tian1607)
!!                 else
!!                    call add_marker(xx, yy, 8 + mod(int(time*0.1*3.171d-8*1.d-5),2), time, nmarkers, 1, i, inc)!(Tian1607)
!!                endif !(Tian: try an infill level)
!!              endif
!             call add_marker(xx, yy, 8 + mod(int(time*0.1*3.171d-8*1.d-5),2), time, nmarkers, 1, i, inc)!(Tian1607)
!              mod(int(time*3.171d-8*1.d-6),10) each 1 Myr alternate phases between 8 and 9 (Tian1607)
              !Tian 20170801 initiate remesh when new infill is added
                    ! Some calculations was not performed every time step, need to
              ! perform these calculation before remeshing
              !      if(topo_kappa .gt. 0 .and. mod(nloop, 100) .ne. 0) call resurface !(Tian:commented)
              
              ! If there are markers recalculate their x,y global coordinate and assign them aps, eII, press, temp
              if(iint_marker.eq.1) then
                 call bar2euler
              endif
!              call remesh
              ! If markers are present recalculate a1,a2 local coordinates and assign elements with phase ratio vector
              if (iint_marker.eq.1) then
                 call lpeuler2bar
                 call marker2elem
              endif
              !Tian20170801 end of remesh due to adding marker
           endif
              if(inc.ne.0) exit
!              write(*,*) cord(1,i:i+1,1), cord(1:2,i,2), xx, yy
          enddo
!          print *, 'debug end' !(Tian 1607 for finding where the problem comes from)
          dhacc(i) = 0

          ! recalculate phase ratio
          kinc = sum(nphase_counter(:,1,i))
          phase_ratio(1:nphase,1,i) = nphase_counter(1:nphase,1,i) / float(kinc)

      else if(-diff*kinc .ge. elz) then
          ! erosion, remove the top marker
          !print *, 'erosion', i, diff, elz
          ymax = -1e30
          nmax = 0
          kmax = 0
          do k = 1, ntopmarker(i)
              n = itopmarker(k, i)
              ntriag = mark(n)%ntriag
              m = mod(ntriag,2) + 1
              call bar2xy(mark(n)%a1, mark(n)%a2, shp2(:,:,m), x, y)
              if(ymax < y) then
                  ymax = y
                  nmax = n
                  kmax = k
              endif
          end do
          mark(nmax)%dead = 0
          ! replace topmarker k with last topmarker
          itopmarker(k,i) = itopmarker(ntopmarker(i),i)
          ntopmarker(i) = ntopmarker(i) - 1

          dhacc(i) = 0

          ! recalculate phase ratio
          kinc = sum(nphase_counter(:,1,i))
          phase_ratio(1:nphase,1,i) = nphase_counter(1:nphase,1,i) / float(kinc)
      end if
  end do

end subroutine resurface


!======================================================================
  subroutine particle_move
! Move particles every time step
! G. Ito 3/07
!======================================================================
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'
!implicit none
real*8 x1,x2,x3,x4,y1,y2,y3,y4,vx1,vx2,vx3,vx4,vy1,vy2,vy3,vy4,xl,xr,zu,zd ! add xl,xr,zu,zd as real8
real*8 area_triangle, areaPBC,areaAPC,areaABP
!real*8 area, areaPBC,areaAPC,areaABP    !This area thing cause rank 0 and rank 3 mismatch
!which cause an hour to find out, area was already declared in arrays.f90
integer SDR_ind
integer ind
!ind = mod(int(time*0.5*3.171d-8*1.d-6),10)+1  !Tian every 2Myr increment the index by 1
!ind = mod(int(time*3.171d-8*1.d-6),10)+1  !Tian every 1Myr increment the index by 1
ind = mod(int(time*3.171d-8*1.d-6),(num_SDR+1))+1  !Tian every 1Myr increment the index by 1
!write(*,*) "ind=",ind
!integer iorder(np)
!real*8 xpord(np)
!--------------------------------------------------------------------------
!Interpolate velocities to particles and move them
!--------------------------------------------------------------------------
dxp=rxbo/dble(np)
!xpmax=x0
!xpmin=x0+rxbo
!---------------------------
!Tian Paralellize the code 20170417
!---------------------------
!$OMP parallel private(SDR_ind,i,vxp,vzp,jj,ii,xl,xr,zu,zd, &
!$OMP                  x1,y1,x2,y2,x3,y3,x4,y4, &
!$OMP                  vx1,vy1,vx2,vy2,vx3,vy3,vx4,vy4, &
!$OMP                  in_out,in_out123,in_out134)


do SDR_ind = 1,ind !SDR_ind is the index of # of reflector
!$OMP do
   do i=1,np    !i is the index of each particle within one reflector
      vxp=999.d0
      vzp=999.d0
      if (xp(i,SDR_ind).lt.cord(1,1,1)) then		!particle outside left side of nodes
         vxp=0.d0
         vzp=0.d0
      elseif (xp(i,SDR_ind).gt.cord(1,nx,1)) then		!particle outside right side of nodes
         vxp=vel(1,nx,1)
         vzp=vel(1,nx,2)
!      elseif (zp(i,SDR_ind).gt.cord(1,1,2)) then		!particle outside top side of nodes
!         vxp=0.d0
!         vzp=0.d0
      elseif (zp(i,SDR_ind).lt.cord(nz,1,2)) then		!particle outside bottom side of nodes
         vxp=0.d0
         vzp=0.d0
      else						!particle within bounds of nodes
         ! find out corresponding element for each particle (Tian)
!         do ii=1,nx-1
!            do jj=1,nz-1    Exchange the order of loop to inprove speed
         do jj=1,nz-1
            do ii=1,nx-1
               !xl=cord(1,ii,1)
               !xr=cord(1,ii+1,1)
               !ref about continuation of a line in fortarn using ampersand "&":
               !https://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap01/continue.html
               xl=cord(jj,ii,1) + (cord(jj,ii,2)-zp(i,SDR_ind))/(cord(jj,ii,2)-cord(jj+1,ii,2))*&
                    (cord(jj+1,ii,1)-cord(jj,ii,1))  !real left bound of xp
               xr=cord(jj,ii+1,1) + (cord(jj,ii+1,2)-zp(i,SDR_ind))/(cord(jj,ii+1,2)-cord(jj+1,ii+1,2))*&
                    (cord(jj+1,ii+1,1)-cord(jj,ii+1,1))  !real right bound of xp

               !zu=cord(jj,ii,2)  !upper bound of z
               !zd=cord(jj+1,ii,2) !lower/down bound of z
               zu=cord(jj,ii,2) + (xp(i,SDR_ind)-xl)/(xr-xl)*(cord(jj,ii+1,2)-cord(jj,ii,2))  !real upper bound of zp
               zd=cord(jj+1,ii,2)+ (xp(i,SDR_ind)-xl)/(xr-xl)*(cord(jj+1,ii+1,2)-cord(jj+1,ii,2)) !lower/down bound of zp
               !for particle moves above the top layer, bring it down into the element(Tian20170515)
               if (xp(i,SDR_ind).ge.xl.and.xp(i,SDR_ind).le.xr) then
                  if (zp(i,SDR_ind).gt.zu .and. zp(i,SDR_ind).le.(zu+200.0d0)) then
                     zp(i,SDR_ind) = zu
                  endif
               endif
               
               if (zp(i,SDR_ind).ge.zd.and.zp(i,SDR_ind).le.zu) then
                  if (xp(i,SDR_ind).ge.xl.and.xp(i,SDR_ind).le.xr) then
                     ! Coordinates of nodes in the element that surround the particle
                     !                 write(*,*) 'got in inner loop'
                     !                 write(*,*) 'ii=',ii,'jj=',jj
                     ! Figurative representation of node points
                     !    1-------3
                     !    |       |
                     !    |       |
                     !    2-------4
                     x1 = cord(jj  ,ii  ,1)
                     y1 = cord(jj  ,ii  ,2)
                     x2 = cord(jj+1,ii  ,1)
                     y2 = cord(jj+1,ii  ,2)
                     x3 = cord(jj  ,ii+1,1)
                     y3 = cord(jj  ,ii+1,2)
                     x4 = cord(jj+1,ii+1,1)
                     y4 = cord(jj+1,ii+1,2)
                     ! Velocities of nodes in the element that surround the particle
                     vx1 = vel(jj  ,ii  ,1)
                     vy1 = vel(jj  ,ii  ,2)
                     vx2 = vel(jj+1,ii  ,1)
                     vy2 = vel(jj+1,ii  ,2)
                     vx3 = vel(jj  ,ii+1,1)
                     vy3 = vel(jj  ,ii+1,2)
                     vx4 = vel(jj+1,ii+1,1)
                     vy4 = vel(jj+1,ii+1,2)
                     goto 222  ! element that surround the particle found, jump out of the double for loop
                  endif
               endif
            end do
         end do
222      continue
         ! four regimes
         !    1-------------3
         !    |  .  III  .  |
         !    | I   . .  IV |
         !    |    .   .    |
         !    |  .  II    . |
         !    2-------------4
         !test if (xp,zp) in triangle 123 or 234 of the found element
         in_out123 = 0  ! zero means not in triangle 123
         in_out134 = 0  ! zero means not in triangle 134
         in_out = 0 ! value for passing in and out of function if_in_triangle
         ! if neither in triangle 123 nor in triangle 124 then it would be in regime II
         ! four regimes:
         ! I :   in_out123 == 1 .and. in_out134 == 0   123+124
         ! II :  in_out123 == 0 .and. in_out134 == 0   124_234
         ! III : in_out123 == 1 .and. in_out134 == 1   123+134
         ! IV :  in_out123 == 0 .and. in_out134 == 1   134+234
         ! test if the particle is in triangle 123
         call if_in_triangle(xp(i,SDR_ind),zp(i,SDR_ind),x1,x2,x3,y1,y2,y3,in_out)
         in_out123 = in_out
         in_out = 0
         ! test if the particle is in triangle 134, replace point B(2) with piont D(4)
         call if_in_triangle(xp(i,SDR_ind),zp(i,SDR_ind),x1,x4,x3,y1,y4,y3,in_out)
         in_out134 = in_out
         ! 
         !     write(*,*) "calling if_in_triangle"
         !     call if_in_triangle(kkx,kky,x1,x2,x3,x4,y1,y2,y3,y4,in_out)
         ! four regimes:
         ! I :   in_out123 == 1 .and. in_out134 == 0   123+124
         if (in_out123 == 1 .and. in_out134 == 0) then
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x2,x3,y1,y2,y3,vx1,vx2,vx3,vy1,vy2,vy3,vxp,vzp)
            tempx = vxp
            tempz = vzp
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x2,x4,y1,y2,y4,vx1,vx2,vx4,vy1,vy2,vy4,vxp,vzp)
            vxp = 0.50d0 * (tempx + vxp)
            vzp = 0.50d0 * (tempz + vzp)
         endif
         ! II :  in_out123 == 0 .and. in_out134 == 0   124_234
         if (in_out123 == 0 .and. in_out134 == 0) then
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x2,x4,y1,y2,y4,vx1,vx2,vx4,vy1,vy2,vy4,vxp,vzp)
            tempx = vxp
            tempz = vzp
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x2,x3,x4,y2,y3,y4,vx2,vx3,vx4,vy2,vy3,vy4,vxp,vzp)
            vxp = 0.50d0 * (tempx + vxp)
            vzp = 0.50d0 * (tempz + vzp)
         endif
         ! III : in_out123 == 1 .and. in_out134 == 1   123+134
         if (in_out123 == 1 .and. in_out134 == 1) then
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x2,x3,y1,y2,y3,vx1,vx2,vx3,vy1,vy2,vy3,vxp,vzp)
            tempx = vxp
            tempz = vzp
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x3,x4,y1,y3,y4,vx1,vx3,vx4,vy1,vy3,vy4,vxp,vzp)
            vxp = 0.50d0 * (tempx + vxp)
            vzp = 0.50d0 * (tempz + vzp)
         endif
         ! IV :  in_out123 == 0 .and. in_out134 == 1   134+234
         if (in_out123 == 0 .and. in_out134 == 1) then
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x1,x3,x4,y1,y3,y4,vx1,vx3,vx4,vy1,vy3,vy4,vxp,vzp)
            tempx = vxp
            tempz = vzp
            call velocity_interpolation(xp(i,SDR_ind),zp(i,SDR_ind),x2,x3,x4,y2,y3,y4,vx2,vx3,vx4,vy2,vy3,vy4,vxp,vzp)
            vxp = 0.50d0 * (tempx + vxp)
            vzp = 0.50d0 * (tempz + vzp)
         endif
!         write(*,*) "i_np=", i
!         if (in_out == 1) then !in triangle 123
            !ref: calc area of triangle
            !https://people.richland.edu/james/lecture/m116/matrices/area.html
!            area_triangle = abs(x1*y2-x1*y3+y1*x3-y1*x2+x2*y3-y2*x3)
!            areaPBC = abs(xp(i,SDR_ind)*y2-xp(i,SDR_ind)*y3+zp(i,SDR_ind)*x3-zp(i,SDR_ind)*x2+x2*y3-y2*x3)
!            areaAPC = abs(x1*zp(i,SDR_ind)-x1*y3+y1*x3-y1*xp(i,SDR_ind)+xp(i,SDR_ind)*y3-zp(i,SDR_ind)*x3)
!            areaABP = abs(x1*y2-x1*zp(i,SDR_ind)+y1*xp(i,SDR_ind)-y1*x2+x2*zp(i,SDR_ind)-y2*xp(i,SDR_ind))
!            vxp = 1.0/area_triangle*(areaPBC*vx1+areaAPC*vx2+areaABP*vx3)
!            vzp = 1.0/area_triangle*(areaPBC*vy1+areaAPC*vy2+areaABP*vy3)
!        else !in triangle 234
!            area_triangle = abs(x4*y2-x4*y3+y4*x3-y4*x2+x2*y3-y2*x3)
!            areaPBC = abs(x4*zp(i,SDR_ind)-x4*y3+y4*x3-y4*xp(i,SDR_ind)+xp(i,SDR_ind)*y3-zp(i,SDR_ind)*x3)
!            areaAPC = abs(x4*y2-x4*zp(i,SDR_ind)+y4*xp(i,SDR_ind)-y4*x2+x2*zp(i,SDR_ind)-y2*xp(i,SDR_ind))
!            areaABP = abs(xp(i,SDR_ind)*y2-xp(i,SDR_ind)*y3+zp(i,SDR_ind)*x3-zp(i,SDR_ind)*x2+x2*y3-y2*x3)
!            vxp = 1.0/area_triangle*(areaPBC*vx2+areaAPC*vx3+areaABP*vx4)
!            vzp = 1.0/area_triangle*(areaPBC*vy2+areaAPC*vy3+areaABP*vy4)
!         endif
         !     write(*,*) 'vxp=',vxp,'vzp=',vzp
         !     write(*,*) 'area_triangle=',area_triangle,'areaPBC=',areaPBC
         xp(i,SDR_ind)=xp(i,SDR_ind) + vxp*dt
         zp(i,SDR_ind)=zp(i,SDR_ind) + vzp*dt
      endif
      if (vxp.eq.999.0d0.or.vzp.eq.999.0d0) then
         write(*,*) 'Problems in particle_move, setting particle velocity'
         write(*,*) "at location of i=", i
         stop
      endif
      !  write(*,*) xp(i), zp(i)
      !100  continue
   end do
!$OMP end do   !Tian parallelized
end do


!$OMP end parallel

return
end subroutine particle_move

!subroutine if_in_triangle(xp,yp,x1,x2,x3,x4,y1,y2,y3,y4,in_out123,in_out134)
subroutine if_in_triangle(xp,yp,x1,x2,x3,y1,y2,y3,in_out)
  real(kind=8), intent(in) :: xp,yp,x1,x2,x3,y1,y2,y3
!  real*8 :: xp,yp,x1,x2,x3,y1,y2,y3
  integer(kind=4), intent(inout) :: in_out
  real*8 :: dot00,dot01,dot02,dot11,dot12,invDenom,u,v   !add precision
!  real*8 :: ypc, y1c, y2c, y3c
  !Tian applied baricentric coordinate to check which triangle xp,zp is in (1,2,3) or (2,3,4)
  !if not in (1,2,3) or (A,B,C) then it must be in (2,3,4), no need to test again
  !ref:http://blackpawn.com/texts/pointinpoly/
  !v0 = C - A or 3 - 1
  !v1 = B - A or 2 - 1
  !v2 = P - A or xp,zp - 1
!  ypc = dabs(yp)
!  y1c = dabs(y1)
!  y2c = dabs(y2)
!  y3c = dabs(y3)
  dot00 = x3*x3+x1*x1-2*x1*x3+y3*y3+y1*y1-2*y1*y3
  dot01 = x3*x2-x3*x1-x1*x2+x1*x1+y3*y2-y3*y1-y1*y2+y1*y1
  dot02 = x3*xp-x3*x1-x1*xp+x1*x1+y3*yp-y3*x1-y1*yp+y1*x1
  dot11 = x2*x2+x1*x1-2*x1*x2+y2*y2+y1*y1-2*y1*y2
  dot12 = x2*xp-x1*x2-x1*xp+x1*x1+y2*yp-y2*x1-y1*yp+x1*y1
!  dot00 = x3*x3+x1*x1-2*x1*x3+y3c*y3c+y1c*y1c-2*y1c*y3c
!  dot01 = x3*x2-x3*x1-x1*x2+x1*x1+y3c*y2c-y3c*y1c-y1c*y2c+y1c*y1c
!  dot02 = x3*xp-x3*x1-x1*xp+x1*x1+y3c*ypc-y3c*x1-y1c*ypc+y1c*x1
!  dot11 = x2*x2+x1*x1-2*x1*x2+y2c*y2c+y1c*y1c-2*y1c*y2c
!  dot12 = x2*xp-x1*x2-x1*xp+x1*x1+y2c*ypc-y2c*x1-y1c*ypc+x1*y1c
  ! Compute barycentric coordinates
  !  invDenom = 1 / (dot00 * dot11 - dot01 * dot01)  ! This line could be the root of problem Tian!!
  invDenom = 1.0d0 / (dot00 * dot11 - dot01 * dot01)
  u = (dot11 * dot02 - dot01 * dot12) * invDenom
  v = (dot00 * dot12 - dot01 * dot02) * invDenom

  ! Check if point is in triangle
  if (u .ge. 0.0d0 .and. v .ge. 0.0d0 .and. (u + v) .le. 1.0d0) then !in the triangle 123
     in_out = 1
  else
     in_out = 0
  endif
!  write(*,*) 'in_out = ', in_out
  return
end subroutine if_in_triangle

subroutine velocity_interpolation(xp,yp,x1,x2,x3,y1,y2,y3,vx1,vx2,vx3,vy1,vy2,vy3,vxp,vzp)
!  real(kind=8), intent(in) :: xp,yp,x1,x2,x3,y1,y2,y3,vx1,vx2,vx3,vy1,vy2,vy3
!  real(kind=8), intent(inout) :: vxp,vzp
  real*8, intent(in) :: xp,yp,x1,x2,x3,y1,y2,y3,vx1,vx2,vx3,vy1,vy2,vy3
  real*8, intent(inout) :: vxp,vzp
  !ref: calc area of triangle
  !https://people.richland.edu/james/lecture/m116/matrices/area.html
!  area_triangle = dabs(x1*y2-x1*y3+y1*x3-y1*x2+x2*y3-y2*x3)
!  areaPBC = dabs(xp*y2-xp*y3+yp*x3-yp*x2+x2*y3-y2*x3)
!  areaAPC = dabs(x1*yp-x1*y3+y1*x3-y1*xp+xp*y3-yp*x3)
!  areaABP = dabs(x1*y2-x1*yp+y1*xp-y1*x2+x2*yp-y2*xp)
  area_triangle = (x1*y2-x1*y3+y1*x3-y1*x2+x2*y3-y2*x3)
  areaPBC = (xp*y2-xp*y3+yp*x3-yp*x2+x2*y3-y2*x3)
  areaAPC = (x1*yp-x1*y3+y1*x3-y1*xp+xp*y3-yp*x3)
  areaABP = (x1*y2-x1*yp+y1*xp-y1*x2+x2*yp-y2*xp)
  vxp = 1.0d0/area_triangle*(areaPBC*vx1+areaAPC*vx2+areaABP*vx3)
  vzp = 1.0d0/area_triangle*(areaPBC*vy1+areaAPC*vy2+areaABP*vy3)
  !  write(*,*) "1.0/area_triangle - 1.0d0/area_triangle=", (1.0/area_triangle - 1.0/area_triangle)

!  write(*,*) "(vxpD0 - vxp)%=", ((1.0/area_triangle*(areaPBC*vx1+areaAPC*vx2+areaABP*vx3) - vxp))/vxp*100.0d0
!  write(*,*) "(vzpD0 - vzp)%=", ((1.0/area_triangle*(areaPBC*vy1+areaAPC*vy2+areaABP*vy3) - vzp))/vzp*100.0d0
  return
end subroutine velocity_interpolation

! Ref for barycentric tecnique
!Barycentric Technique

!The advantage of the method above is that it's very simple to understand so that once you read it you should be able to remember it forever and code it up at any time without having to refer back to anything. It's just - hey the point has to be on the same side of each line as the triangle point that's not in the line. Cake.

!Well, there's another method that is also as easy conceptually but executes faster. The downside is there's a little more math involved, but once you see it worked out it should be no problem.

!So remember that the three points of the triangle define a plane in space. Pick one of the points and we can consider all other locations on the plane as relative to that point. Let's go with A -- it'll be our origin on the plane. Now what we need are basis vectors so we can give coordinate values to all the locations on the plane. We'll pick the two edges of the triangle that touch A, (C - A) and (B - A). Now we can get to any point on the plane just by starting at A and walking some distance along (C - A) and then from there walking some more in the direction (B - A).

!With that in mind we can now describe any point on the plane as

!P = A + u * (C - A) + v * (B - A)
!Notice now that if u or v < 0 then we've walked in the wrong direction and must be outside the triangle. Also if u or v > 1 then we've walked too far in a direction and are outside the triangle. Finally if u + v > 1 then we've crossed the edge BC again leaving the triangle.

!Given u and v we can easily calculate the point P with the above equation, but how can we go in the reverse direction and calculate u and v from a given point P? Time for some math!

!    P = A + u * (C - A) + v * (B - A)       // Original equation
!    (P - A) = u * (C - A) + v * (B - A)     // Subtract A from both sides
!    v2 = u * v0 + v * v1                    // Substitute v0, v1, v2 for less writing
    
!    // We have two unknowns (u and v) so we need two equations to solve
!    // for them.  Dot both sides by v0 to get one and dot both sides by
!    // v1 to get a second.
!    (v2) . v0 = (u * v0 + v * v1) . v0
!    (v2) . v1 = (u * v0 + v * v1) . v1

!    // Distribute v0 and v1
!    v2 . v0 = u * (v0 . v0) + v * (v1 . v0)
!    v2 . v1 = u * (v0 . v1) + v * (v1 . v1)

!   // Now we have two equations and two unknowns and can solve one 
!   // equation for one variable and substitute into the other.  Or
!    // if you're lazy like me, fire up Mathematica and save yourself
!    // some handwriting.
!    Solve[v2.v0 == {u(v0.v0) + v(v1.v0), v2.v1 == u(v0.v1) + v(v1.v1)}, {u, v}]
!    u = ((v1.v1)(v2.v0)-(v1.v0)(v2.v1)) / ((v0.v0)(v1.v1) - (v0.v1)(v1.v0))
!    v = ((v0.v0)(v2.v1)-(v0.v1)(v2.v0)) / ((v0.v0)(v1.v1) - (v0.v1)(v1.v0))
!    // Compute vectors
!    v0 = C - A
!    v1 = B - A
!    v2 = P - A

!    // Compute dot products
!    dot00 = dot(v0, v0)
!    dot01 = dot(v0, v1)
!    dot02 = dot(v0, v2)
!    dot11 = dot(v1, v1)
!    dot12 = dot(v1, v2)

!    // Compute barycentric coordinates
!    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
!    u = (dot11 * dot02 - dot01 * dot12) * invDenom
!    v = (dot00 * dot12 - dot01 * dot02) * invDenom
!
!    // Check if point is in triangle
!    return (u >= 0) && (v >= 0) && (u + v < 1)

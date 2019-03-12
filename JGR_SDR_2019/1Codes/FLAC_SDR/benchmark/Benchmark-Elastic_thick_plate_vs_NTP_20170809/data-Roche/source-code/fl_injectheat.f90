!==============================================================================
!  Adding magma and latent heat in injection zone
!  G. Ito 8/11/06 
!==============================================================================
subroutine fl_injectheat    
use arrays                       
include 'precision.inc' ! vs BI(Behn and Ito 2008): the same
include 'params.inc'   ! vs BI: not the same
include 'arrays.inc'   ! vs BI: not the same
                          ! (J-1,I-1)------(J-1,I)----(J-1,I+1)
                          !     |            |            |
                          !     |   a11      |   a12      |
                          !     | (j-1,i-1)  | (j,i-1)    |
                          !     |            |            |
                          !  (J,I-1)-------(J,I)--------(J,I+1)
                          !     |            |            |
                          !     |   a21      |    a22     |
                          !     | (j-1,i)    |   (j,i)    |
                          !     |            |            |
                          ! (J+1,I-1)------(J+1,I)----(J+1,I+1)
 
dimension njTinj(nz) 

ninj=iinj2-iinj1+1   !new ninj(horizontal #nodes), iinj1(left bound of dike), iinj2(right bound)
!fdum=xlatheat/( cp(1)*(Tliq-Tsol) )   !L/(Cp*(T_liq - T_sol))
iph = 6   ! phase of the dike
fdum=xlatheat/( cp(iph)*(tl(iph)-ts(iph)) )   !L/(Cp*(T_liq - T_sol))
ninjbot=jinj2     !new ninjbot, jinj2(bottom # of dike)
!print *, 'iinj2=', iinj2, 'ninj=',ninj,'xlatheat=',xlatheat,'fdum=',fdum,'ninjbot=',ninjbot !Tian debug
!----------------------------------------
! Tian: don't understand this block (begin)
!----------------------------------------
do 110 i=iinj1,iinj2+1
goto 222 !Tian Comment this block, it make ninjbot become weird
   jcnt=1
   njTinj(1) = ninjbot				!Start: same code as in fl_rheol
   do j = 1,nz-1
     dcord = cord(1,i,2) - cord(j+1,i,2)
     if(dcord.gt.Tcinj) then
       njTinj(jcnt) = j
       jcnt = jcnt+1
     endif
   end do
   ninjbot = min(ninjbot,njTinj(1))		!End: same code as in fl_rheol
222 continue
!----------------------------------------
! Tian: don't understand this block (end)
!----------------------------------------

   do 100 j=jinj1,ninjbot+1
!-------------------------------------------------------------------------------
! Compute areas to centers of adjacent cells
!-------------------------------------------------------------------------------
   i1=max0(i-1,1)
   j1=max0(j-1,1)
   i2=min0(i+1,nx)
   j2=min0(j+1,nz)
     
   dx11=0.5*(cord(j1,i ,1)-cord(j1,i1,1)+cord(j ,i ,1)-cord(j ,i1,1))
   dy11=0.5*(cord(j1,i1,2)-cord(j ,i1,2)+cord(j1,i ,2)-cord(j ,i ,2))
   
   dx12=0.5*(cord(j1,i2,1)-cord(j1,i ,1)+cord(j ,i2,1)-cord(j ,i ,1))
   dy12=0.5*(cord(j1,i ,2)-cord(j ,i ,2)+cord(j1,i2,2)-cord(j ,i2,2))
   
   dx21=0.5*(cord(j ,i ,1)-cord(j ,i1,1)+cord(j2,i ,1)-cord(j2,i1,1))
   dy21=0.5*(cord(j ,i1,2)-cord(j2,i1,2)+cord(j ,i ,2)-cord(j2,i  ,2))
   
   dx22=0.5*(cord(j ,i2,1)-cord(j ,i ,1)+cord(j2,i2,1)-cord(j2,i  ,1))
   dy22=0.5*(cord(j ,i ,2)-cord(j2,i ,2)+cord(j ,i2,2)-cord(j2,i2,2))

   a11=dabs(dx11*dy11)/4.  !no real need for 1/4 but it reminds us areas are to center of element
   a12=dabs(dx12*dy12)/4.
   a21=dabs(dx21*dy21)/4.
   a22=dabs(dx22*dy22)/4.
   atot=a11+a12+a21+a22
!-------------------------------------------------------------------------------
! Compute dT's
!-------------------------------------------------------------------------------

! different M values in the brittle and ductile zones - Olive 4/08
!         if (temp(j,i).lt.600) then
!          rate_inject=rate_inject_brittle
!         else
!          rate_inject=rate_inject_ductile
!         endif


   rfac=ratfac*rate_inject*dt/(dx11*dble(ninj))   ! no need to have ratfac
   T0=0.25*(temp(j1,i1)+temp(j1,i)+temp(j,i1)+temp(j,i))
   dtemp=dtemp_inj(T0,tl(iph),ts(iph),rfac,fdum)
   dtemp=dmin1(dtemp,tl(iph)-T0)
   dtemp11=dmax1(dtemp,0.0)

   rfac=ratfac*rate_inject*dt/(dx12*dble(ninj))
   T0=0.25*(temp(j1,i)+temp(j1,i2)+temp(j,i)+temp(j,i2))
   dtemp=dtemp_inj(T0,tl(iph),ts(iph),rfac,fdum)
   dtemp=dmin1(dtemp,tl(iph)-T0)
   dtemp12=dmax1(dtemp,0.0)

   rfac=ratfac*rate_inject*dt/(dx21*dble(ninj))
   T0=0.25*(temp(j,i1)+temp(j,i)+temp(j2,i1)+temp(j2,i))
   dtemp=dtemp_inj(T0,tl(iph),ts(iph),rfac,fdum)
   dtemp=dmin1(dtemp,tl(iph)-T0)
   dtemp21=dmax1(dtemp,0.0)

   rfac=ratfac*rate_inject*dt/(dx22*dble(ninj))
   T0=0.25*(temp(j,i)+temp(j,i2)+temp(j2,i)+temp(j2,i2))
   dtemp=dtemp_inj(T0,tl(iph),ts(iph),rfac,fdum)
   dtemp=dmin1(dtemp,tl(iph)-T0)
   dtemp22=dmax1(dtemp,0.0)
!-------------------------------------------------------------------------------
! Eliminate contribution from elements outside diking zone
!-------------------------------------------------------------------------------
   f11=1.0		
   f12=1.0
   f21=1.0
   f22=1.0 
   if     (i.eq.iinj1) then
     f11=0.
     f21=0.
   elseif (i.eq.iinj2+1) then
     f12=0.
     f22=0.
   endif
   
   if (j.eq.jinj1) then  
     f11=0.                               
     f12=0. 
   elseif (j.eq.ninjbot+1) then
     f21=0.
     f22=0.
   endif
   
!-------------------------------------------------------------------------------
! Update Nodal temperatures
!-------------------------------------------------------------------------------
  dtemp_ave=(dtemp11*f11*a22+dtemp12*f12*a21+dtemp21*f21*a12+dtemp22*f22*a11)/atot

  temp(j,i)=temp(j,i)+dtemp_ave
  if(temp(j,i)>2000 .or. isnan(temp(j,i)) .or. temp(j,i)<0) then
     print *, 'j=',j,'i=',i,'dtemp_ave=',dtemp_ave ,'temp(',j,',',i,')=',temp(j  ,i  )
  endif
!------------------------------------------------------------------------------------
! There is/was a bug that make dtemp_ave=NAN, with iinj1=1,
! but it goes away with this if statement.  Keep and eye out!
!------------------------------------------------------------------------------------
  if (dabs(dtemp_ave).gt.10.*tl(iph)) then
    write(*,'(8ES11.3,3i4)') dtemp11, dtemp12, dtemp21, dtemp22, dtemp_ave, dble(ninj), &
    rfac, atot, i,j,nloop
    stop
  endif

100 continue
110 continue

! Boundary conditions (top)
!
if (jinj1.eq.1) then
  do i = iinj1,iinj2+1
    temp(1,i) = t_top
  end do
endif

! Boundary conditions: dt/dx =0 on left and right  
!$DIR PREFER_PARALLEL
if (iinj1.eq.1) then
  do j = jinj1,jinj2
    temp(j ,1)  = temp(j,2)
  end do
endif


return
end  

!==============================================================================
! Compute dT including effects of latent heat
!==============================================================================
function dtemp_inj(T0,Tliq,Tsol,rfac,fdum)

include 'precision.inc'   

! Garrett's original implementation (Page 4 of his notes)
!Tdum=(T0*(1.-rfac) + Tliq*(1.+fdum)*rfac)/(1.0+fdum*rfac)

!if (Tdum .lt. Tsol) then
!  dtemp_inj=( (Tliq-T0) + fdum*(Tliq-Tsol) )*rfac   
!else
!  dtemp_inj=(Tliq-T0)*(1+fdum)*rfac/(1.0+fdum*rfac)
!endif

!Revised by M. Behn Nov 2007
Tdum=(T0 + Tliq*rfac)/(1+rfac)  !<?> What is Tdum, 

if (Tdum .lt. Tsol) then
   dtemp_inj=( (Tliq-T0) + fdum*(Tliq-Tsol) ) * (rfac/(1+rfac))
   !<?> why rfac/(1+rfac) isn't rfac is melt fraction in dike zone?
!else
!  dtemp_inj=(Tliq-T0) * (rfac/(1+rfac))
elseif (Tdum.lt.Tliq .and. Tdum.ge.Tsol) then
   dtemp_inj=(Tliq-T0) * (rfac/(1+rfac))
elseif (Tdum.ge.Tliq) then
   dtemp_inj = 0.
endif

return
end




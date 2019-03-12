!=======================================================
! Surface particles to Keep track of SDRs (adapted from Ito2007) Tian20170405
!=======================================================

subroutine ReadParticles
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

call AdvanceToNextInputLine( 4 ) ! add to the last line of the input file and param.inc
read(4,*) npelem,num_SDR
write(*,*) 'npelem=',npelem,'num_SDR=',num_SDR
!Surface particles, number of particles per element;
! num_SDR is total number of SDRs, should be less than 20, which is set in arrays.inc
write(*,*) 'mnp=',mnp
!if(npelem*mnx.gt.mnp) then
!    write(*,*) '# of particles exceed maximum mnp. Increase mnp in "arrays.inc".'
!    stop 1
!endif
write(*,*) 'max_num_SDR=', max_num_SDR
if(num_SDR.gt.max_num_SDR) then
    write(*,*) '# of number of SDRs exceed max_num_SDR. Increase max_num_SDR in "arrays.inc".'
    stop 1
endif

return
end

!=======================================================
! Heat injection due to diking (adapted from Behn and Ito 2008) Tian 2017
!=======================================================

subroutine ReadHeatinject
include 'precision.inc'
include 'params.inc'

call AdvanceToNextInputLine( 4 ) ! add to the last line of the input file and param.inc
read(4,*) iinj1, iinj2, jinj1, jinj2, xlatheat, ratfac

return
end

!=======================================================
! HYDROTHERMAL ALTERATION OF THERMAL DIFFUSIVITY
! Modified from Luc's code, G. Ito 8/7/06
!=======================================================
!-------------------------------------------------------------
subroutine ReadHydro()
!-------------------------------------------------------------
include 'precision.inc'
include 'params.inc'


!open( 9, file='hydrother.inp',status='old',err=2001 )
if_hydro = 1
call AdvanceToNextInputLine(4)
read (4,*) xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2
write(*,*) '>>Hydrothermal effects on Diffusivity<<<'
if (xenhc2.lt.xenhc1) then
  write(*,*) 'WARNING:  xenhc2 corrected to be >/= xenhc1'
endif
write(*,*) ' xmaxdepth,     xmaxt,   xmaxstr,    xenhc1,    xenhc2'
write(*,'(5f10.2)') xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2
write(*,*) '>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<'
!close( 9 )


return

!2001 if_hydro = 0
!write(*,*) '>>NO Hydrothermal effects on Diffusivity<<<'
!return

end

!-------------------------------------------------------------
function HydroCond(j,i)
!-------------------------------------------------------------
  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
!  common /hydroth/ xmaxdepth,xmaxt,xmaxstr,xenhc1,xenhc2

!  iph = iphase(i,j,phasez(j,i))
  iph = iphase(j,i)
  tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
  yc = 0.25*(cord(j,i,2)+cord(j+1,i,2)+cord(j,i+1,2)+cord(j+1,i+1,2))
  
  if( tmpr.le.xmaxt .and. yc.ge.xmaxdepth) then
    cdum=dmin1(aps(j,i)/xmaxstr, 1.0)
    HydroCond=(xenhc1 + cdum*(xenhc2-xenhc1))*conduct(iph) 
  else
    HydroCond=conduct(iph)
  endif
  
return

end


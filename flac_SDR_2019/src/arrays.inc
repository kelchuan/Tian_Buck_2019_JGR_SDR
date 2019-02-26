! -*- F90 -*-

!!! maximum number of ELEMENTS !!!
parameter(mnz=200, mnx=700, mnp=20*mnx, max_num_SDR=20) !Tian20170405 adding mnp for particle SDRs
parameter(max_markers_per_elem=32)

common /phase_arrays/ iphase(mnz,mnx), nphase_counter(20,mnz,mnx), phase_ratio(20,mnz,mnx)

common /topo_arrays/ dhacc(mnx+1), &
     ntopmarker(mnx), itopmarker(max_markers_per_elem,mnx)

common /all_arrays/ strainr(3,4,mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),e2sr(mnz,mnx),irheol_fl(mnz,mnx)

common /temp_arrays/ temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bcstress((mnz+mnx)*2,3), &
     nopbou((mnz+mnx)*2,4), ncodbou((mnz+mnx)*2,3)

!!common /particles/ xp(mnp), zp(mnp) !Tian20170405 for particles to keep track of SDR
!!real, dimension (1:mnp,1:max_num_SDR) :: xp, zp
!!real, dimension (mnp,max_num_SDR) :: xp, zp
!common /particles_arrays/ xp(mnp,max_num_SDR), zp(mnp,max_num_SDR) !Tian20170405 for particles to keep track of SDR

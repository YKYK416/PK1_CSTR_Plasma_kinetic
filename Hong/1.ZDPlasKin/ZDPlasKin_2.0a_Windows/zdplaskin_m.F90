!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Mon Jul  6 15:13:42 2026
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E N H S M 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 65, species_electrons = 1, species_length = 7, reactions_max = 158, reactions_length = 24
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 2, bolsig_species_length = 2, bolsig_rates_max = 24 
  character(*), parameter, private          :: bolsigfile = "BOLSIGDB.DAT"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  /-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,-1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0/
  data species_name(1:species_max) &
  /"E      ","N2     ","H2     ","N      ","H      ","NH     ","NH2    ","NH3    ","N2H    ","M      ","N2(V1) ","N2(V2) ",&
   "N2(V3) ","N2(V4) ","N2(V5) ","N2(V6) ","N2(V7) ","N2(V8) ","N2(V9) ","N2(V10)","N2(V11)","N2(V12)","N2(V13)","N2(V14)",&
   "N2(V15)","N2(V16)","N2(V17)","N2(V18)","N2(V19)","N2(V20)","N2(V21)","N2(V22)","N2(V23)","N2(V24)","H2(V1) ","H2(V2) ",&
   "H2(V3) ","N2(A3) ","N2(B3) ","N2(A1) ","N2(C3) ","H2(B3) ","H2(B1) ","H2(C3) ","H2(A3) ","N(2D)  ","N(2P)  ","N2^+   ",&
   "N3^+   ","N4^+   ","N^+    ","H^+    ","H2^+   ","H3^+   ","H^-    ","NH^+   ","NH2^+  ","NH3^+  ","NH4^+  ","N2H^+  ",&
   "S      ","NS     ","HS     ","NHS    ","NH2S   "/
  data reaction_sign(1:72) &
  /"N+NH=>H+N2              ","H+NH=>N+H2              ","NH+NH=>H2+N2            ","NH+NH=>N+NH2            ",&
   "NH+NH=>N2+H+H           ","H+NH2=>H2+NH            ","N+NH2=>N2+H+H           ","N+NH2=>N2+H2            ",&
   "NH+NH2=>NH3+N           ","H2+N=>NH+H              ","H2+NH2=>NH3+H           ","H+NH3=>NH2+H2           ",&
   "N2(A3)+H=>N2+H          ","N2(A3)+H2=>N2+H+H       ","N2(A3)+NH3=>N2+NH3      ","N2(B3)+H2=>N2(A3)+H2    ",&
   "N2(A1)+H=>N2+H          ","N2(A1)+H2=>N2+H+H       ","N+H2(B3)=>H+NH          ","N(2D)+H2=>H+NH          ",&
   "N(2D)+NH3=>NH+NH2       ","N(2P)+H2=>H+NH          ","N+N+M=>N2+M             ","N+N+N=>N2(A3)+N         ",&
   "N+N+N=>N2(B3)+N         ","N+N+N2=>N2(A3)+N2       ","N+N+N2=>N2(B3)+N2       ","N+N+H2=>N2+H2           ",&
   "H+H+N2=>H2+N2           ","H+N+M=>NH+M             ","N+H2+M=>NH2+M           ","H+NH+M=>NH2+M           ",&
   "H+NH2+M=>NH3+M          ","NH+H2+M=>NH3+M          ","N+N+H2=>N2(A3)+H2       ","N+N+H=>N2(A3)+H         ",&
   "N+N+H2=>N2(B3)+H2       ","N+N+H=>N2(B3)+H         ","H+H+H2=>H2+H2           ","N+N=>N2^++E             ",&
   "N2(A1)+N2(A1)=>N2^++N2+E","N2(A3)+N2(A1)=>N2^++N2+E","N2(A3)=>N2              ","N2(B3)=>N2(A3)          ",&
   "N2(A1)=>N2              ","N2(C3)=>N2(B3)          ","N^++H2=>NH^++H          ","N^++NH3=>NH2^++NH       ",&
   "N^++NH3=>NH3^++N        ","N^++NH3=>N2H^++H2       ","N2^++N=>N^++N2          ","N2^++H2=>N2H^++H        ",&
   "N2^++N2(A3)=>N3^++N     ","N2^++NH3=>NH3^++N2      ","N3^++N=>N2^++N2         ","N4^++N=>N^++N2+N2       ",&
   "N4^++N2=>N2^++N2+N2     ","H^++NH3=>NH3^++H        ","H2^++H=>H^++H2          ","H2^++H2=>H3^++H         ",&
   "H2^++N2=>N2H^++H        ","H2^++NH3=>NH3^++H2      ","NH^++H2=>H3^++N         ","NH^++H2=>NH2^++H        ",&
   "NH^++NH3=>NH3^++NH      ","NH^++NH3=>NH4^++N       ","NH^++N2=>N2H^++N        ","NH2^++H2=>NH3^++H       ",&
   "NH2^++NH3=>NH3^++NH2    ","NH2^++NH3=>NH4^++NH     ","NH3^++NH3=>NH4^++NH2    ","N2H^++NH3=>NH4^++N2     "/
  data reaction_sign(73:144) &
  /"N2^++N+N2=>N3^++N2      ","N^++N2+N2=>N3^++N2      ","N2^++N2+N2=>N4^++N2     ","N^++N+N2=>N2^++N2       ",&
   "H^-+H^+=>H+H            ","H^-+H2^+=>H+H+H         ","H^-+H3^+=>H2+H+H        ","H^-+N2^+=>N2+H          ",&
   "H^-+N4^+=>N2+N2+H       ","H^-+N2H^+=>H2+N2        ","H^-+H^++M=>H2+M         ","H^-+H2^++M=>H2+H+M      ",&
   "H^-+H3^++M=>H2+H2+M     ","H^-+N2^++M=>N2+H+M      ","H^-+N4^++M=>N2+N2+H+M   ","H^-+N2H^++M=>H2+N2+M    ",&
   "N2(V1)+N2=>N2+N2        ","N2(V2)+N2=>N2(V1)+N2    ","N2(V3)+N2=>N2(V2)+N2    ","N2(V4)+N2=>N2(V3)+N2    ",&
   "N2(V5)+N2=>N2(V4)+N2    ","N2(V6)+N2=>N2(V5)+N2    ","N2(V7)+N2=>N2(V6)+N2    ","N2(V8)+N2=>N2(V7)+N2    ",&
   "H2(V1)+H2=>H2+H2        ","H2(V2)+H2=>H2(V1)+H2    ","H2(V3)+H2=>H2(V2)+H2    ","N2(V1)+H=>N2+H          ",&
   "N2(V2)+H=>N2(V1)+H      ","N2(V3)+H=>N2(V2)+H      ","N2(V4)+H=>N2(V3)+H      ","N2(V5)+H=>N2(V4)+H      ",&
   "N2(V6)+H=>N2(V5)+H      ","N2(V7)+H=>N2(V6)+H      ","N2(V8)+H=>N2(V7)+H      ","H2(V1)+N=>H2+N          ",&
   "H2(V2)+N=>H2(V1)+N      ","H2(V3)+N=>H2(V2)+N      ","N2(A3)+S=>N2+S          ","H+S=>HS                 ",&
   "N2(B3)+S=>N2+S          ","H2+S=>H2+S              ","N+S=>NS                 ","NH+S=>NHS               ",&
   "NH2+S=>NH2S             ","N+HS=>NHS               ","H+HS=>H2+S              ","N2+HS=>N2H+S            ",&
   "NH+HS=>NH2S             ","H+NS=>NHS               ","H+NHS=>NH2S             ","H+NH2S=>NH3+S           ",&
   "NH2+HS=>NH3+S           ","NS+HS=>NHS+S            ","NHS+HS=>NH2S+S          ","NH2S+HS=>NH3+S+S        ",&
   "N2+S+S=>NS+NS           ","N2(A3)+S+S=>NS+NS       ","H2+S+S=>HS+HS           ","H2(V1)+S+S=>HS+HS       ",&
   "H2(V2)+S+S=>HS+HS       ","H2(V3)+S+S=>HS+HS       ","bolsig:H2->H2(B3)       ","bolsig:H2->H2(B1)       ",&
   "bolsig:H2->H2(C3)       ","bolsig:H2->H2(A3)       ","bolsig:N2->N2(A3)       ","bolsig:N2->N2(B3)       ",&
   "bolsig:N2->N2(A1)       ","bolsig:N2->N2(C3)       ","bolsig:N2->N2^+         ","bolsig:H2->H2^+         "/
  data reaction_sign(145:158) &
  /"bolsig:N->N^+           ","bolsig:H->H^+           ","bolsig:NH->NH^+         ","bolsig:NH2->NH2^+       ",&
   "bolsig:NH3->NH3^+       ","bolsig:N2->N^+          ","bolsig:H2->H^+          ","bolsig:NH->H^+          ",&
   "bolsig:NH2->H^+         ","bolsig:NH3->H^+         ","bolsig:H2->H+H          ","bolsig:N2->N+N          ",&
   "bolsig:N2->N2(V1)       ","bolsig:H2->H2(V1)       "/
  data bolsig_species(1:bolsig_species_max) &
  /"N2","H2"/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A7)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A24)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(291x,65(1x,i9))
312 format(A3,1x,A24,1x,65(1x,A9))
313 format(i3,1x,A24,1x,65(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(158(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,65(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,158(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,65(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(159(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(02,001) = + reac_rate_local(001) 
  reac_source_local(04,001) = - reac_rate_local(001) 
  reac_source_local(05,001) = + reac_rate_local(001) 
  reac_source_local(06,001) = - reac_rate_local(001) 
  reac_source_local(03,002) = + reac_rate_local(002) 
  reac_source_local(04,002) = + reac_rate_local(002) 
  reac_source_local(05,002) = - reac_rate_local(002) 
  reac_source_local(06,002) = - reac_rate_local(002) 
  reac_source_local(02,003) = + reac_rate_local(003) 
  reac_source_local(03,003) = + reac_rate_local(003) 
  reac_source_local(06,003) = - reac_rate_local(003) * 2.d0
  reac_source_local(04,004) = + reac_rate_local(004) 
  reac_source_local(06,004) = - reac_rate_local(004) * 2.d0
  reac_source_local(07,004) = + reac_rate_local(004) 
  reac_source_local(02,005) = + reac_rate_local(005) 
  reac_source_local(05,005) = + reac_rate_local(005) * 2.d0
  reac_source_local(06,005) = - reac_rate_local(005) * 2.d0
  reac_source_local(03,006) = + reac_rate_local(006) 
  reac_source_local(05,006) = - reac_rate_local(006) 
  reac_source_local(06,006) = + reac_rate_local(006) 
  reac_source_local(07,006) = - reac_rate_local(006) 
  reac_source_local(02,007) = + reac_rate_local(007) 
  reac_source_local(04,007) = - reac_rate_local(007) 
  reac_source_local(05,007) = + reac_rate_local(007) * 2.d0
  reac_source_local(07,007) = - reac_rate_local(007) 
  reac_source_local(02,008) = + reac_rate_local(008) 
  reac_source_local(03,008) = + reac_rate_local(008) 
  reac_source_local(04,008) = - reac_rate_local(008) 
  reac_source_local(07,008) = - reac_rate_local(008) 
  reac_source_local(04,009) = + reac_rate_local(009) 
  reac_source_local(06,009) = - reac_rate_local(009) 
  reac_source_local(07,009) = - reac_rate_local(009) 
  reac_source_local(08,009) = + reac_rate_local(009) 
  reac_source_local(03,010) = - reac_rate_local(010) 
  reac_source_local(04,010) = - reac_rate_local(010) 
  reac_source_local(05,010) = + reac_rate_local(010) 
  reac_source_local(06,010) = + reac_rate_local(010) 
  reac_source_local(03,011) = - reac_rate_local(011) 
  reac_source_local(05,011) = + reac_rate_local(011) 
  reac_source_local(07,011) = - reac_rate_local(011) 
  reac_source_local(08,011) = + reac_rate_local(011) 
  reac_source_local(03,012) = + reac_rate_local(012) 
  reac_source_local(05,012) = - reac_rate_local(012) 
  reac_source_local(07,012) = + reac_rate_local(012) 
  reac_source_local(08,012) = - reac_rate_local(012) 
  reac_source_local(02,013) = + reac_rate_local(013) 
  reac_source_local(38,013) = - reac_rate_local(013) 
  reac_source_local(02,014) = + reac_rate_local(014) 
  reac_source_local(03,014) = - reac_rate_local(014) 
  reac_source_local(05,014) = + reac_rate_local(014) * 2.d0
  reac_source_local(38,014) = - reac_rate_local(014) 
  reac_source_local(02,015) = + reac_rate_local(015) 
  reac_source_local(38,015) = - reac_rate_local(015) 
  reac_source_local(38,016) = + reac_rate_local(016) 
  reac_source_local(39,016) = - reac_rate_local(016) 
  reac_source_local(02,017) = + reac_rate_local(017) 
  reac_source_local(40,017) = - reac_rate_local(017) 
  reac_source_local(02,018) = + reac_rate_local(018) 
  reac_source_local(03,018) = - reac_rate_local(018) 
  reac_source_local(05,018) = + reac_rate_local(018) * 2.d0
  reac_source_local(40,018) = - reac_rate_local(018) 
  reac_source_local(04,019) = - reac_rate_local(019) 
  reac_source_local(05,019) = + reac_rate_local(019) 
  reac_source_local(06,019) = + reac_rate_local(019) 
  reac_source_local(42,019) = - reac_rate_local(019) 
  reac_source_local(03,020) = - reac_rate_local(020) 
  reac_source_local(05,020) = + reac_rate_local(020) 
  reac_source_local(06,020) = + reac_rate_local(020) 
  reac_source_local(46,020) = - reac_rate_local(020) 
  reac_source_local(06,021) = + reac_rate_local(021) 
  reac_source_local(07,021) = + reac_rate_local(021) 
  reac_source_local(08,021) = - reac_rate_local(021) 
  reac_source_local(46,021) = - reac_rate_local(021) 
  reac_source_local(03,022) = - reac_rate_local(022) 
  reac_source_local(05,022) = + reac_rate_local(022) 
  reac_source_local(06,022) = + reac_rate_local(022) 
  reac_source_local(47,022) = - reac_rate_local(022) 
  reac_source_local(02,023) = + reac_rate_local(023) 
  reac_source_local(04,023) = - reac_rate_local(023) * 2.d0
  reac_source_local(04,024) = - reac_rate_local(024) * 2.d0
  reac_source_local(38,024) = + reac_rate_local(024) 
  reac_source_local(04,025) = - reac_rate_local(025) * 2.d0
  reac_source_local(39,025) = + reac_rate_local(025) 
  reac_source_local(04,026) = - reac_rate_local(026) * 2.d0
  reac_source_local(38,026) = + reac_rate_local(026) 
  reac_source_local(04,027) = - reac_rate_local(027) * 2.d0
  reac_source_local(39,027) = + reac_rate_local(027) 
  reac_source_local(02,028) = + reac_rate_local(028) 
  reac_source_local(04,028) = - reac_rate_local(028) * 2.d0
  reac_source_local(03,029) = + reac_rate_local(029) 
  reac_source_local(05,029) = - reac_rate_local(029) * 2.d0
  reac_source_local(04,030) = - reac_rate_local(030) 
  reac_source_local(05,030) = - reac_rate_local(030) 
  reac_source_local(06,030) = + reac_rate_local(030) 
  reac_source_local(03,031) = - reac_rate_local(031) 
  reac_source_local(04,031) = - reac_rate_local(031) 
  reac_source_local(07,031) = + reac_rate_local(031) 
  reac_source_local(05,032) = - reac_rate_local(032) 
  reac_source_local(06,032) = - reac_rate_local(032) 
  reac_source_local(07,032) = + reac_rate_local(032) 
  reac_source_local(05,033) = - reac_rate_local(033) 
  reac_source_local(07,033) = - reac_rate_local(033) 
  reac_source_local(08,033) = + reac_rate_local(033) 
  reac_source_local(03,034) = - reac_rate_local(034) 
  reac_source_local(06,034) = - reac_rate_local(034) 
  reac_source_local(08,034) = + reac_rate_local(034) 
  reac_source_local(04,035) = - reac_rate_local(035) * 2.d0
  reac_source_local(38,035) = + reac_rate_local(035) 
  reac_source_local(04,036) = - reac_rate_local(036) * 2.d0
  reac_source_local(38,036) = + reac_rate_local(036) 
  reac_source_local(04,037) = - reac_rate_local(037) * 2.d0
  reac_source_local(39,037) = + reac_rate_local(037) 
  reac_source_local(04,038) = - reac_rate_local(038) * 2.d0
  reac_source_local(39,038) = + reac_rate_local(038) 
  reac_source_local(03,039) = + reac_rate_local(039) 
  reac_source_local(05,039) = - reac_rate_local(039) * 2.d0
  reac_source_local(01,040) = + reac_rate_local(040) 
  reac_source_local(04,040) = - reac_rate_local(040) * 2.d0
  reac_source_local(48,040) = + reac_rate_local(040) 
  reac_source_local(01,041) = + reac_rate_local(041) 
  reac_source_local(02,041) = + reac_rate_local(041) 
  reac_source_local(40,041) = - reac_rate_local(041) * 2.d0
  reac_source_local(48,041) = + reac_rate_local(041) 
  reac_source_local(01,042) = + reac_rate_local(042) 
  reac_source_local(02,042) = + reac_rate_local(042) 
  reac_source_local(38,042) = - reac_rate_local(042) 
  reac_source_local(40,042) = - reac_rate_local(042) 
  reac_source_local(48,042) = + reac_rate_local(042) 
  reac_source_local(02,043) = + reac_rate_local(043) 
  reac_source_local(38,043) = - reac_rate_local(043) 
  reac_source_local(38,044) = + reac_rate_local(044) 
  reac_source_local(39,044) = - reac_rate_local(044) 
  reac_source_local(02,045) = + reac_rate_local(045) 
  reac_source_local(40,045) = - reac_rate_local(045) 
  reac_source_local(39,046) = + reac_rate_local(046) 
  reac_source_local(41,046) = - reac_rate_local(046) 
  reac_source_local(03,047) = - reac_rate_local(047) 
  reac_source_local(05,047) = + reac_rate_local(047) 
  reac_source_local(51,047) = - reac_rate_local(047) 
  reac_source_local(56,047) = + reac_rate_local(047) 
  reac_source_local(06,048) = + reac_rate_local(048) 
  reac_source_local(08,048) = - reac_rate_local(048) 
  reac_source_local(51,048) = - reac_rate_local(048) 
  reac_source_local(57,048) = + reac_rate_local(048) 
  reac_source_local(04,049) = + reac_rate_local(049) 
  reac_source_local(08,049) = - reac_rate_local(049) 
  reac_source_local(51,049) = - reac_rate_local(049) 
  reac_source_local(58,049) = + reac_rate_local(049) 
  reac_source_local(03,050) = + reac_rate_local(050) 
  reac_source_local(08,050) = - reac_rate_local(050) 
  reac_source_local(51,050) = - reac_rate_local(050) 
  reac_source_local(60,050) = + reac_rate_local(050) 
  reac_source_local(02,051) = + reac_rate_local(051) 
  reac_source_local(04,051) = - reac_rate_local(051) 
  reac_source_local(48,051) = - reac_rate_local(051) 
  reac_source_local(51,051) = + reac_rate_local(051) 
  reac_source_local(03,052) = - reac_rate_local(052) 
  reac_source_local(05,052) = + reac_rate_local(052) 
  reac_source_local(48,052) = - reac_rate_local(052) 
  reac_source_local(60,052) = + reac_rate_local(052) 
  reac_source_local(04,053) = + reac_rate_local(053) 
  reac_source_local(38,053) = - reac_rate_local(053) 
  reac_source_local(48,053) = - reac_rate_local(053) 
  reac_source_local(49,053) = + reac_rate_local(053) 
  reac_source_local(02,054) = + reac_rate_local(054) 
  reac_source_local(08,054) = - reac_rate_local(054) 
  reac_source_local(48,054) = - reac_rate_local(054) 
  reac_source_local(58,054) = + reac_rate_local(054) 
  reac_source_local(02,055) = + reac_rate_local(055) 
  reac_source_local(04,055) = - reac_rate_local(055) 
  reac_source_local(48,055) = + reac_rate_local(055) 
  reac_source_local(49,055) = - reac_rate_local(055) 
  reac_source_local(02,056) = + reac_rate_local(056) * 2.d0
  reac_source_local(04,056) = - reac_rate_local(056) 
  reac_source_local(50,056) = - reac_rate_local(056) 
  reac_source_local(51,056) = + reac_rate_local(056) 
  reac_source_local(02,057) = + reac_rate_local(057) 
  reac_source_local(48,057) = + reac_rate_local(057) 
  reac_source_local(50,057) = - reac_rate_local(057) 
  reac_source_local(05,058) = + reac_rate_local(058) 
  reac_source_local(08,058) = - reac_rate_local(058) 
  reac_source_local(52,058) = - reac_rate_local(058) 
  reac_source_local(58,058) = + reac_rate_local(058) 
  reac_source_local(03,059) = + reac_rate_local(059) 
  reac_source_local(05,059) = - reac_rate_local(059) 
  reac_source_local(52,059) = + reac_rate_local(059) 
  reac_source_local(53,059) = - reac_rate_local(059) 
  reac_source_local(03,060) = - reac_rate_local(060) 
  reac_source_local(05,060) = + reac_rate_local(060) 
  reac_source_local(53,060) = - reac_rate_local(060) 
  reac_source_local(54,060) = + reac_rate_local(060) 
  reac_source_local(02,061) = - reac_rate_local(061) 
  reac_source_local(05,061) = + reac_rate_local(061) 
  reac_source_local(53,061) = - reac_rate_local(061) 
  reac_source_local(60,061) = + reac_rate_local(061) 
  reac_source_local(03,062) = + reac_rate_local(062) 
  reac_source_local(08,062) = - reac_rate_local(062) 
  reac_source_local(53,062) = - reac_rate_local(062) 
  reac_source_local(58,062) = + reac_rate_local(062) 
  reac_source_local(03,063) = - reac_rate_local(063) 
  reac_source_local(04,063) = + reac_rate_local(063) 
  reac_source_local(54,063) = + reac_rate_local(063) 
  reac_source_local(56,063) = - reac_rate_local(063) 
  reac_source_local(03,064) = - reac_rate_local(064) 
  reac_source_local(05,064) = + reac_rate_local(064) 
  reac_source_local(56,064) = - reac_rate_local(064) 
  reac_source_local(57,064) = + reac_rate_local(064) 
  reac_source_local(06,065) = + reac_rate_local(065) 
  reac_source_local(08,065) = - reac_rate_local(065) 
  reac_source_local(56,065) = - reac_rate_local(065) 
  reac_source_local(58,065) = + reac_rate_local(065) 
  reac_source_local(04,066) = + reac_rate_local(066) 
  reac_source_local(08,066) = - reac_rate_local(066) 
  reac_source_local(56,066) = - reac_rate_local(066) 
  reac_source_local(59,066) = + reac_rate_local(066) 
  reac_source_local(02,067) = - reac_rate_local(067) 
  reac_source_local(04,067) = + reac_rate_local(067) 
  reac_source_local(56,067) = - reac_rate_local(067) 
  reac_source_local(60,067) = + reac_rate_local(067) 
  reac_source_local(03,068) = - reac_rate_local(068) 
  reac_source_local(05,068) = + reac_rate_local(068) 
  reac_source_local(57,068) = - reac_rate_local(068) 
  reac_source_local(58,068) = + reac_rate_local(068) 
  reac_source_local(07,069) = + reac_rate_local(069) 
  reac_source_local(08,069) = - reac_rate_local(069) 
  reac_source_local(57,069) = - reac_rate_local(069) 
  reac_source_local(58,069) = + reac_rate_local(069) 
  reac_source_local(06,070) = + reac_rate_local(070) 
  reac_source_local(08,070) = - reac_rate_local(070) 
  reac_source_local(57,070) = - reac_rate_local(070) 
  reac_source_local(59,070) = + reac_rate_local(070) 
  reac_source_local(07,071) = + reac_rate_local(071) 
  reac_source_local(08,071) = - reac_rate_local(071) 
  reac_source_local(58,071) = - reac_rate_local(071) 
  reac_source_local(59,071) = + reac_rate_local(071) 
  reac_source_local(02,072) = + reac_rate_local(072) 
  reac_source_local(08,072) = - reac_rate_local(072) 
  reac_source_local(59,072) = + reac_rate_local(072) 
  reac_source_local(60,072) = - reac_rate_local(072) 
  reac_source_local(04,073) = - reac_rate_local(073) 
  reac_source_local(48,073) = - reac_rate_local(073) 
  reac_source_local(49,073) = + reac_rate_local(073) 
  reac_source_local(02,074) = - reac_rate_local(074) 
  reac_source_local(49,074) = + reac_rate_local(074) 
  reac_source_local(51,074) = - reac_rate_local(074) 
  reac_source_local(02,075) = - reac_rate_local(075) 
  reac_source_local(48,075) = - reac_rate_local(075) 
  reac_source_local(50,075) = + reac_rate_local(075) 
  reac_source_local(04,076) = - reac_rate_local(076) 
  reac_source_local(48,076) = + reac_rate_local(076) 
  reac_source_local(51,076) = - reac_rate_local(076) 
  reac_source_local(05,077) = + reac_rate_local(077) * 2.d0
  reac_source_local(52,077) = - reac_rate_local(077) 
  reac_source_local(55,077) = - reac_rate_local(077) 
  reac_source_local(05,078) = + reac_rate_local(078) * 3.d0
  reac_source_local(53,078) = - reac_rate_local(078) 
  reac_source_local(55,078) = - reac_rate_local(078) 
  reac_source_local(03,079) = + reac_rate_local(079) 
  reac_source_local(05,079) = + reac_rate_local(079) * 2.d0
  reac_source_local(54,079) = - reac_rate_local(079) 
  reac_source_local(55,079) = - reac_rate_local(079) 
  reac_source_local(02,080) = + reac_rate_local(080) 
  reac_source_local(05,080) = + reac_rate_local(080) 
  reac_source_local(48,080) = - reac_rate_local(080) 
  reac_source_local(55,080) = - reac_rate_local(080) 
  reac_source_local(02,081) = + reac_rate_local(081) * 2.d0
  reac_source_local(05,081) = + reac_rate_local(081) 
  reac_source_local(50,081) = - reac_rate_local(081) 
  reac_source_local(55,081) = - reac_rate_local(081) 
  reac_source_local(02,082) = + reac_rate_local(082) 
  reac_source_local(03,082) = + reac_rate_local(082) 
  reac_source_local(55,082) = - reac_rate_local(082) 
  reac_source_local(60,082) = - reac_rate_local(082) 
  reac_source_local(03,083) = + reac_rate_local(083) 
  reac_source_local(52,083) = - reac_rate_local(083) 
  reac_source_local(55,083) = - reac_rate_local(083) 
  reac_source_local(03,084) = + reac_rate_local(084) 
  reac_source_local(05,084) = + reac_rate_local(084) 
  reac_source_local(53,084) = - reac_rate_local(084) 
  reac_source_local(55,084) = - reac_rate_local(084) 
  reac_source_local(03,085) = + reac_rate_local(085) * 2.d0
  reac_source_local(54,085) = - reac_rate_local(085) 
  reac_source_local(55,085) = - reac_rate_local(085) 
  reac_source_local(02,086) = + reac_rate_local(086) 
  reac_source_local(05,086) = + reac_rate_local(086) 
  reac_source_local(48,086) = - reac_rate_local(086) 
  reac_source_local(55,086) = - reac_rate_local(086) 
  reac_source_local(02,087) = + reac_rate_local(087) * 2.d0
  reac_source_local(05,087) = + reac_rate_local(087) 
  reac_source_local(50,087) = - reac_rate_local(087) 
  reac_source_local(55,087) = - reac_rate_local(087) 
  reac_source_local(02,088) = + reac_rate_local(088) 
  reac_source_local(03,088) = + reac_rate_local(088) 
  reac_source_local(55,088) = - reac_rate_local(088) 
  reac_source_local(60,088) = - reac_rate_local(088) 
  reac_source_local(02,089) = + reac_rate_local(089) 
  reac_source_local(11,089) = - reac_rate_local(089) 
  reac_source_local(11,090) = + reac_rate_local(090) 
  reac_source_local(12,090) = - reac_rate_local(090) 
  reac_source_local(12,091) = + reac_rate_local(091) 
  reac_source_local(13,091) = - reac_rate_local(091) 
  reac_source_local(13,092) = + reac_rate_local(092) 
  reac_source_local(14,092) = - reac_rate_local(092) 
  reac_source_local(14,093) = + reac_rate_local(093) 
  reac_source_local(15,093) = - reac_rate_local(093) 
  reac_source_local(15,094) = + reac_rate_local(094) 
  reac_source_local(16,094) = - reac_rate_local(094) 
  reac_source_local(16,095) = + reac_rate_local(095) 
  reac_source_local(17,095) = - reac_rate_local(095) 
  reac_source_local(17,096) = + reac_rate_local(096) 
  reac_source_local(18,096) = - reac_rate_local(096) 
  reac_source_local(03,097) = + reac_rate_local(097) 
  reac_source_local(35,097) = - reac_rate_local(097) 
  reac_source_local(35,098) = + reac_rate_local(098) 
  reac_source_local(36,098) = - reac_rate_local(098) 
  reac_source_local(36,099) = + reac_rate_local(099) 
  reac_source_local(37,099) = - reac_rate_local(099) 
  reac_source_local(02,100) = + reac_rate_local(100) 
  reac_source_local(11,100) = - reac_rate_local(100) 
  reac_source_local(11,101) = + reac_rate_local(101) 
  reac_source_local(12,101) = - reac_rate_local(101) 
  reac_source_local(12,102) = + reac_rate_local(102) 
  reac_source_local(13,102) = - reac_rate_local(102) 
  reac_source_local(13,103) = + reac_rate_local(103) 
  reac_source_local(14,103) = - reac_rate_local(103) 
  reac_source_local(14,104) = + reac_rate_local(104) 
  reac_source_local(15,104) = - reac_rate_local(104) 
  reac_source_local(15,105) = + reac_rate_local(105) 
  reac_source_local(16,105) = - reac_rate_local(105) 
  reac_source_local(16,106) = + reac_rate_local(106) 
  reac_source_local(17,106) = - reac_rate_local(106) 
  reac_source_local(17,107) = + reac_rate_local(107) 
  reac_source_local(18,107) = - reac_rate_local(107) 
  reac_source_local(03,108) = + reac_rate_local(108) 
  reac_source_local(35,108) = - reac_rate_local(108) 
  reac_source_local(35,109) = + reac_rate_local(109) 
  reac_source_local(36,109) = - reac_rate_local(109) 
  reac_source_local(36,110) = + reac_rate_local(110) 
  reac_source_local(37,110) = - reac_rate_local(110) 
  reac_source_local(02,111) = + reac_rate_local(111) 
  reac_source_local(38,111) = - reac_rate_local(111) 
  reac_source_local(05,112) = - reac_rate_local(112) 
  reac_source_local(61,112) = - reac_rate_local(112) 
  reac_source_local(63,112) = + reac_rate_local(112) 
  reac_source_local(02,113) = + reac_rate_local(113) 
  reac_source_local(39,113) = - reac_rate_local(113) 
  reac_source_local(04,115) = - reac_rate_local(115) 
  reac_source_local(61,115) = - reac_rate_local(115) 
  reac_source_local(62,115) = + reac_rate_local(115) 
  reac_source_local(06,116) = - reac_rate_local(116) 
  reac_source_local(61,116) = - reac_rate_local(116) 
  reac_source_local(64,116) = + reac_rate_local(116) 
  reac_source_local(07,117) = - reac_rate_local(117) 
  reac_source_local(61,117) = - reac_rate_local(117) 
  reac_source_local(65,117) = + reac_rate_local(117) 
  reac_source_local(04,118) = - reac_rate_local(118) 
  reac_source_local(63,118) = - reac_rate_local(118) 
  reac_source_local(64,118) = + reac_rate_local(118) 
  reac_source_local(03,119) = + reac_rate_local(119) 
  reac_source_local(05,119) = - reac_rate_local(119) 
  reac_source_local(61,119) = + reac_rate_local(119) 
  reac_source_local(63,119) = - reac_rate_local(119) 
  reac_source_local(02,120) = - reac_rate_local(120) 
  reac_source_local(09,120) = + reac_rate_local(120) 
  reac_source_local(61,120) = + reac_rate_local(120) 
  reac_source_local(63,120) = - reac_rate_local(120) 
  reac_source_local(06,121) = - reac_rate_local(121) 
  reac_source_local(63,121) = - reac_rate_local(121) 
  reac_source_local(65,121) = + reac_rate_local(121) 
  reac_source_local(05,122) = - reac_rate_local(122) 
  reac_source_local(62,122) = - reac_rate_local(122) 
  reac_source_local(64,122) = + reac_rate_local(122) 
  reac_source_local(05,123) = - reac_rate_local(123) 
  reac_source_local(64,123) = - reac_rate_local(123) 
  reac_source_local(65,123) = + reac_rate_local(123) 
  reac_source_local(05,124) = - reac_rate_local(124) 
  reac_source_local(08,124) = + reac_rate_local(124) 
  reac_source_local(61,124) = + reac_rate_local(124) 
  reac_source_local(65,124) = - reac_rate_local(124) 
  reac_source_local(07,125) = - reac_rate_local(125) 
  reac_source_local(08,125) = + reac_rate_local(125) 
  reac_source_local(61,125) = + reac_rate_local(125) 
  reac_source_local(63,125) = - reac_rate_local(125) 
  reac_source_local(61,126) = + reac_rate_local(126) 
  reac_source_local(62,126) = - reac_rate_local(126) 
  reac_source_local(63,126) = - reac_rate_local(126) 
  reac_source_local(64,126) = + reac_rate_local(126) 
  reac_source_local(61,127) = + reac_rate_local(127) 
  reac_source_local(63,127) = - reac_rate_local(127) 
  reac_source_local(64,127) = - reac_rate_local(127) 
  reac_source_local(65,127) = + reac_rate_local(127) 
  reac_source_local(08,128) = + reac_rate_local(128) 
  reac_source_local(61,128) = + reac_rate_local(128) * 2.d0
  reac_source_local(63,128) = - reac_rate_local(128) 
  reac_source_local(65,128) = - reac_rate_local(128) 
  reac_source_local(02,129) = - reac_rate_local(129) 
  reac_source_local(61,129) = - reac_rate_local(129) * 2.d0
  reac_source_local(62,129) = + reac_rate_local(129) * 2.d0
  reac_source_local(38,130) = - reac_rate_local(130) 
  reac_source_local(61,130) = - reac_rate_local(130) * 2.d0
  reac_source_local(62,130) = + reac_rate_local(130) * 2.d0
  reac_source_local(03,131) = - reac_rate_local(131) 
  reac_source_local(61,131) = - reac_rate_local(131) * 2.d0
  reac_source_local(63,131) = + reac_rate_local(131) * 2.d0
  reac_source_local(35,132) = - reac_rate_local(132) 
  reac_source_local(61,132) = - reac_rate_local(132) * 2.d0
  reac_source_local(63,132) = + reac_rate_local(132) * 2.d0
  reac_source_local(36,133) = - reac_rate_local(133) 
  reac_source_local(61,133) = - reac_rate_local(133) * 2.d0
  reac_source_local(63,133) = + reac_rate_local(133) * 2.d0
  reac_source_local(37,134) = - reac_rate_local(134) 
  reac_source_local(61,134) = - reac_rate_local(134) * 2.d0
  reac_source_local(63,134) = + reac_rate_local(134) * 2.d0
  reac_source_local(03,135) = - reac_rate_local(135) 
  reac_source_local(42,135) = + reac_rate_local(135) 
  reac_source_local(03,136) = - reac_rate_local(136) 
  reac_source_local(43,136) = + reac_rate_local(136) 
  reac_source_local(03,137) = - reac_rate_local(137) 
  reac_source_local(44,137) = + reac_rate_local(137) 
  reac_source_local(03,138) = - reac_rate_local(138) 
  reac_source_local(45,138) = + reac_rate_local(138) 
  reac_source_local(02,139) = - reac_rate_local(139) 
  reac_source_local(38,139) = + reac_rate_local(139) 
  reac_source_local(02,140) = - reac_rate_local(140) 
  reac_source_local(39,140) = + reac_rate_local(140) 
  reac_source_local(02,141) = - reac_rate_local(141) 
  reac_source_local(40,141) = + reac_rate_local(141) 
  reac_source_local(02,142) = - reac_rate_local(142) 
  reac_source_local(41,142) = + reac_rate_local(142) 
  reac_source_local(01,143) = + reac_rate_local(143) 
  reac_source_local(02,143) = - reac_rate_local(143) 
  reac_source_local(48,143) = + reac_rate_local(143) 
  reac_source_local(01,144) = + reac_rate_local(144) 
  reac_source_local(03,144) = - reac_rate_local(144) 
  reac_source_local(53,144) = + reac_rate_local(144) 
  reac_source_local(01,145) = + reac_rate_local(145) 
  reac_source_local(04,145) = - reac_rate_local(145) 
  reac_source_local(51,145) = + reac_rate_local(145) 
  reac_source_local(01,146) = + reac_rate_local(146) 
  reac_source_local(05,146) = - reac_rate_local(146) 
  reac_source_local(52,146) = + reac_rate_local(146) 
  reac_source_local(01,147) = + reac_rate_local(147) 
  reac_source_local(06,147) = - reac_rate_local(147) 
  reac_source_local(56,147) = + reac_rate_local(147) 
  reac_source_local(01,148) = + reac_rate_local(148) 
  reac_source_local(07,148) = - reac_rate_local(148) 
  reac_source_local(57,148) = + reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(08,149) = - reac_rate_local(149) 
  reac_source_local(58,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(02,150) = - reac_rate_local(150) 
  reac_source_local(04,150) = + reac_rate_local(150) 
  reac_source_local(51,150) = + reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(03,151) = - reac_rate_local(151) 
  reac_source_local(05,151) = + reac_rate_local(151) 
  reac_source_local(52,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(04,152) = + reac_rate_local(152) 
  reac_source_local(06,152) = - reac_rate_local(152) 
  reac_source_local(52,152) = + reac_rate_local(152) 
  reac_source_local(01,153) = + reac_rate_local(153) 
  reac_source_local(06,153) = + reac_rate_local(153) 
  reac_source_local(07,153) = - reac_rate_local(153) 
  reac_source_local(52,153) = + reac_rate_local(153) 
  reac_source_local(01,154) = + reac_rate_local(154) 
  reac_source_local(07,154) = + reac_rate_local(154) 
  reac_source_local(08,154) = - reac_rate_local(154) 
  reac_source_local(52,154) = + reac_rate_local(154) 
  reac_source_local(03,155) = - reac_rate_local(155) 
  reac_source_local(05,155) = + reac_rate_local(155) * 2.d0
  reac_source_local(02,156) = - reac_rate_local(156) 
  reac_source_local(04,156) = + reac_rate_local(156) * 2.d0
  reac_source_local(02,157) = - reac_rate_local(157) 
  reac_source_local(11,157) = + reac_rate_local(157) 
  reac_source_local(03,158) = - reac_rate_local(158) 
  reac_source_local(35,158) = + reac_rate_local(158) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(66)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(04) * density(06) 
  rrt(002) = rrt(002) * density(05) * density(06) 
  rrt(003) = rrt(003) * density(06)**2 
  rrt(004) = rrt(004) * density(06)**2 
  rrt(005) = rrt(005) * density(06)**2 
  rrt(006) = rrt(006) * density(05) * density(07) 
  rrt(007) = rrt(007) * density(04) * density(07) 
  rrt(008) = rrt(008) * density(04) * density(07) 
  rrt(009) = rrt(009) * density(06) * density(07) 
  rrt(010) = rrt(010) * density(03) * density(04) 
  rrt(011) = rrt(011) * density(03) * density(07) 
  rrt(012) = rrt(012) * density(05) * density(08) 
  rrt(013) = rrt(013) * density(05) * density(38) 
  rrt(014) = rrt(014) * density(03) * density(38) 
  rrt(015) = rrt(015) * density(08) * density(38) 
  rrt(016) = rrt(016) * density(03) * density(39) 
  rrt(017) = rrt(017) * density(05) * density(40) 
  rrt(018) = rrt(018) * density(03) * density(40) 
  rrt(019) = rrt(019) * density(04) * density(42) 
  rrt(020) = rrt(020) * density(03) * density(46) 
  rrt(021) = rrt(021) * density(08) * density(46) 
  rrt(022) = rrt(022) * density(03) * density(47) 
  rrt(023) = rrt(023) * density(04)**2 * density(10) 
  rrt(024) = rrt(024) * density(04)**3 
  rrt(025) = rrt(025) * density(04)**3 
  rrt(026) = rrt(026) * density(02) * density(04)**2 
  rrt(027) = rrt(027) * density(02) * density(04)**2 
  rrt(028) = rrt(028) * density(03) * density(04)**2 
  rrt(029) = rrt(029) * density(02) * density(05)**2 
  rrt(030) = rrt(030) * density(04) * density(05) * density(10) 
  rrt(031) = rrt(031) * density(03) * density(04) * density(10) 
  rrt(032) = rrt(032) * density(05) * density(06) * density(10) 
  rrt(033) = rrt(033) * density(05) * density(07) * density(10) 
  rrt(034) = rrt(034) * density(03) * density(06) * density(10) 
  rrt(035) = rrt(035) * density(03) * density(04)**2 
  rrt(036) = rrt(036) * density(04)**2 * density(05) 
  rrt(037) = rrt(037) * density(03) * density(04)**2 
  rrt(038) = rrt(038) * density(04)**2 * density(05) 
  rrt(039) = rrt(039) * density(03) * density(05)**2 
  rrt(040) = rrt(040) * density(04)**2 
  rrt(041) = rrt(041) * density(40)**2 
  rrt(042) = rrt(042) * density(38) * density(40) 
  rrt(043) = rrt(043) * density(38) 
  rrt(044) = rrt(044) * density(39) 
  rrt(045) = rrt(045) * density(40) 
  rrt(046) = rrt(046) * density(41) 
  rrt(047) = rrt(047) * density(03) * density(51) 
  rrt(048) = rrt(048) * density(08) * density(51) 
  rrt(049) = rrt(049) * density(08) * density(51) 
  rrt(050) = rrt(050) * density(08) * density(51) 
  rrt(051) = rrt(051) * density(04) * density(48) 
  rrt(052) = rrt(052) * density(03) * density(48) 
  rrt(053) = rrt(053) * density(38) * density(48) 
  rrt(054) = rrt(054) * density(08) * density(48) 
  rrt(055) = rrt(055) * density(04) * density(49) 
  rrt(056) = rrt(056) * density(04) * density(50) 
  rrt(057) = rrt(057) * density(02) * density(50) 
  rrt(058) = rrt(058) * density(08) * density(52) 
  rrt(059) = rrt(059) * density(05) * density(53) 
  rrt(060) = rrt(060) * density(03) * density(53) 
  rrt(061) = rrt(061) * density(02) * density(53) 
  rrt(062) = rrt(062) * density(08) * density(53) 
  rrt(063) = rrt(063) * density(03) * density(56) 
  rrt(064) = rrt(064) * density(03) * density(56) 
  rrt(065) = rrt(065) * density(08) * density(56) 
  rrt(066) = rrt(066) * density(08) * density(56) 
  rrt(067) = rrt(067) * density(02) * density(56) 
  rrt(068) = rrt(068) * density(03) * density(57) 
  rrt(069) = rrt(069) * density(08) * density(57) 
  rrt(070) = rrt(070) * density(08) * density(57) 
  rrt(071) = rrt(071) * density(08) * density(58) 
  rrt(072) = rrt(072) * density(08) * density(60) 
  rrt(073) = rrt(073) * density(02) * density(04) * density(48) 
  rrt(074) = rrt(074) * density(02)**2 * density(51) 
  rrt(075) = rrt(075) * density(02)**2 * density(48) 
  rrt(076) = rrt(076) * density(02) * density(04) * density(51) 
  rrt(077) = rrt(077) * density(52) * density(55) 
  rrt(078) = rrt(078) * density(53) * density(55) 
  rrt(079) = rrt(079) * density(54) * density(55) 
  rrt(080) = rrt(080) * density(48) * density(55) 
  rrt(081) = rrt(081) * density(50) * density(55) 
  rrt(082) = rrt(082) * density(55) * density(60) 
  rrt(083) = rrt(083) * density(10) * density(52) * density(55) 
  rrt(084) = rrt(084) * density(10) * density(53) * density(55) 
  rrt(085) = rrt(085) * density(10) * density(54) * density(55) 
  rrt(086) = rrt(086) * density(10) * density(48) * density(55) 
  rrt(087) = rrt(087) * density(10) * density(50) * density(55) 
  rrt(088) = rrt(088) * density(10) * density(55) * density(60) 
  rrt(089) = rrt(089) * density(02) * density(11) 
  rrt(090) = rrt(090) * density(02) * density(12) 
  rrt(091) = rrt(091) * density(02) * density(13) 
  rrt(092) = rrt(092) * density(02) * density(14) 
  rrt(093) = rrt(093) * density(02) * density(15) 
  rrt(094) = rrt(094) * density(02) * density(16) 
  rrt(095) = rrt(095) * density(02) * density(17) 
  rrt(096) = rrt(096) * density(02) * density(18) 
  rrt(097) = rrt(097) * density(03) * density(35) 
  rrt(098) = rrt(098) * density(03) * density(36) 
  rrt(099) = rrt(099) * density(03) * density(37) 
  rrt(100) = rrt(100) * density(05) * density(11) 
  rrt(101) = rrt(101) * density(05) * density(12) 
  rrt(102) = rrt(102) * density(05) * density(13) 
  rrt(103) = rrt(103) * density(05) * density(14) 
  rrt(104) = rrt(104) * density(05) * density(15) 
  rrt(105) = rrt(105) * density(05) * density(16) 
  rrt(106) = rrt(106) * density(05) * density(17) 
  rrt(107) = rrt(107) * density(05) * density(18) 
  rrt(108) = rrt(108) * density(04) * density(35) 
  rrt(109) = rrt(109) * density(04) * density(36) 
  rrt(110) = rrt(110) * density(04) * density(37) 
  rrt(111) = rrt(111) * density(38) * density(61) 
  rrt(112) = rrt(112) * density(05) * density(61) 
  rrt(113) = rrt(113) * density(39) * density(61) 
  rrt(114) = rrt(114) * density(03) * density(61) 
  rrt(115) = rrt(115) * density(04) * density(61) 
  rrt(116) = rrt(116) * density(06) * density(61) 
  rrt(117) = rrt(117) * density(07) * density(61) 
  rrt(118) = rrt(118) * density(04) * density(63) 
  rrt(119) = rrt(119) * density(05) * density(63) 
  rrt(120) = rrt(120) * density(02) * density(63) 
  rrt(121) = rrt(121) * density(06) * density(63) 
  rrt(122) = rrt(122) * density(05) * density(62) 
  rrt(123) = rrt(123) * density(05) * density(64) 
  rrt(124) = rrt(124) * density(05) * density(65) 
  rrt(125) = rrt(125) * density(07) * density(63) 
  rrt(126) = rrt(126) * density(62) * density(63) 
  rrt(127) = rrt(127) * density(63) * density(64) 
  rrt(128) = rrt(128) * density(63) * density(65) 
  rrt(129) = rrt(129) * density(02) * density(61)**2 
  rrt(130) = rrt(130) * density(38) * density(61)**2 
  rrt(131) = rrt(131) * density(03) * density(61)**2 
  rrt(132) = rrt(132) * density(35) * density(61)**2 
  rrt(133) = rrt(133) * density(36) * density(61)**2 
  rrt(134) = rrt(134) * density(37) * density(61)**2 
  rrt(135) = rrt(135) * density(01) * density(03) 
  rrt(136) = rrt(136) * density(01) * density(03) 
  rrt(137) = rrt(137) * density(01) * density(03) 
  rrt(138) = rrt(138) * density(01) * density(03) 
  rrt(139) = rrt(139) * density(01) * density(02) 
  rrt(140) = rrt(140) * density(01) * density(02) 
  rrt(141) = rrt(141) * density(01) * density(02) 
  rrt(142) = rrt(142) * density(01) * density(02) 
  rrt(143) = rrt(143) * density(01) * density(02) 
  rrt(144) = rrt(144) * density(01) * density(03) 
  rrt(145) = rrt(145) * density(01) * density(04) 
  rrt(146) = rrt(146) * density(01) * density(05) 
  rrt(147) = rrt(147) * density(01) * density(06) 
  rrt(148) = rrt(148) * density(01) * density(07) 
  rrt(149) = rrt(149) * density(01) * density(08) 
  rrt(150) = rrt(150) * density(01) * density(02) 
  rrt(151) = rrt(151) * density(01) * density(03) 
  rrt(152) = rrt(152) * density(01) * density(06) 
  rrt(153) = rrt(153) * density(01) * density(07) 
  rrt(154) = rrt(154) * density(01) * density(08) 
  rrt(155) = rrt(155) * density(01) * density(03) 
  rrt(156) = rrt(156) * density(01) * density(02) 
  rrt(157) = rrt(157) * density(01) * density(02) 
  rrt(158) = rrt(158) * density(01) * density(03) 
  ydot(01) = +rrt(040)+rrt(041)+rrt(042)+rrt(143)+rrt(144)+rrt(145)+rrt(146)+rrt(147)+rrt(148)+rrt(149)+rrt(150)+rrt(151)+rrt(152)&
             +rrt(153)+rrt(154) 
  ydot(02) = +rrt(001)+rrt(003)+rrt(005)+rrt(007)+rrt(008)+rrt(013)+rrt(014)+rrt(015)+rrt(017)+rrt(018)+rrt(023)+rrt(028)+rrt(041)&
             +rrt(042)+rrt(043)+rrt(045)+rrt(051)+rrt(054)+rrt(055)+  2.d0 * rrt(056)+rrt(057)-rrt(061)-rrt(067)+rrt(072)-rrt(074)&
             -rrt(075)+rrt(080)+  2.d0 * rrt(081)+rrt(082)+rrt(086)+  2.d0 * rrt(087)+rrt(088)+rrt(089)+rrt(100)+rrt(111)+rrt(113)&
             -rrt(120)-rrt(129)-rrt(139)-rrt(140)-rrt(141)-rrt(142)-rrt(143)-rrt(150)-rrt(156)-rrt(157) 
  ydot(03) = +rrt(002)+rrt(003)+rrt(006)+rrt(008)-rrt(010)-rrt(011)+rrt(012)-rrt(014)-rrt(018)-rrt(020)-rrt(022)+rrt(029)-rrt(031)&
             -rrt(034)+rrt(039)-rrt(047)+rrt(050)-rrt(052)+rrt(059)-rrt(060)+rrt(062)-rrt(063)-rrt(064)-rrt(068)+rrt(079)+rrt(082)&
             +rrt(083)+rrt(084)+  2.d0 * rrt(085)+rrt(088)+rrt(097)+rrt(108)+rrt(119)-rrt(131)-rrt(135)-rrt(136)-rrt(137)-rrt(138)&
             -rrt(144)-rrt(151)-rrt(155)-rrt(158) 
  ydot(04) = -rrt(001)+rrt(002)+rrt(004)-rrt(007)-rrt(008)+rrt(009)-rrt(010)-rrt(019)-  2.d0 * rrt(023)-  2.d0 * rrt(024)&
             -  2.d0 * rrt(025)-  2.d0 * rrt(026)-  2.d0 * rrt(027)-  2.d0 * rrt(028)-rrt(030)-rrt(031)-  2.d0 * rrt(035)&
             -  2.d0 * rrt(036)-  2.d0 * rrt(037)-  2.d0 * rrt(038)-  2.d0 * rrt(040)+rrt(049)-rrt(051)+rrt(053)-rrt(055)-rrt(056)&
             +rrt(063)+rrt(066)+rrt(067)-rrt(073)-rrt(076)-rrt(115)-rrt(118)-rrt(145)+rrt(150)+rrt(152)+  2.d0 * rrt(156) 
  ydot(05) = +rrt(001)-rrt(002)+  2.d0 * rrt(005)-rrt(006)+  2.d0 * rrt(007)+rrt(010)+rrt(011)-rrt(012)+  2.d0 * rrt(014)&
             +  2.d0 * rrt(018)+rrt(019)+rrt(020)+rrt(022)-  2.d0 * rrt(029)-rrt(030)-rrt(032)-rrt(033)-  2.d0 * rrt(039)+rrt(047)&
             +rrt(052)+rrt(058)-rrt(059)+rrt(060)+rrt(061)+rrt(064)+rrt(068)+  2.d0 * rrt(077)+  3.d0 * rrt(078)+  2.d0 * rrt(079)&
             +rrt(080)+rrt(081)+rrt(084)+rrt(086)+rrt(087)-rrt(112)-rrt(119)-rrt(122)-rrt(123)-rrt(124)-rrt(146)+rrt(151)&
             +  2.d0 * rrt(155) 
  ydot(06) = -rrt(001)-rrt(002)-  2.d0 * rrt(003)-  2.d0 * rrt(004)-  2.d0 * rrt(005)+rrt(006)-rrt(009)+rrt(010)+rrt(019)+rrt(020)&
             +rrt(021)+rrt(022)+rrt(030)-rrt(032)-rrt(034)+rrt(048)+rrt(065)+rrt(070)-rrt(116)-rrt(121)-rrt(147)-rrt(152)+rrt(153) 
  ydot(07) = +rrt(004)-rrt(006)-rrt(007)-rrt(008)-rrt(009)-rrt(011)+rrt(012)+rrt(021)+rrt(031)+rrt(032)-rrt(033)+rrt(069)+rrt(071)&
             -rrt(117)-rrt(125)-rrt(148)-rrt(153)+rrt(154) 
  ydot(08) = +rrt(009)+rrt(011)-rrt(012)-rrt(021)+rrt(033)+rrt(034)-rrt(048)-rrt(049)-rrt(050)-rrt(054)-rrt(058)-rrt(062)-rrt(065)&
             -rrt(066)-rrt(069)-rrt(070)-rrt(071)-rrt(072)+rrt(124)+rrt(125)+rrt(128)-rrt(149)-rrt(154) 
  ydot(09) = +rrt(120) 
  ydot(10) = 0.0d0
  ydot(11) = -rrt(089)+rrt(090)-rrt(100)+rrt(101)+rrt(157) 
  ydot(12) = -rrt(090)+rrt(091)-rrt(101)+rrt(102) 
  ydot(13) = -rrt(091)+rrt(092)-rrt(102)+rrt(103) 
  ydot(14) = -rrt(092)+rrt(093)-rrt(103)+rrt(104) 
  ydot(15) = -rrt(093)+rrt(094)-rrt(104)+rrt(105) 
  ydot(16) = -rrt(094)+rrt(095)-rrt(105)+rrt(106) 
  ydot(17) = -rrt(095)+rrt(096)-rrt(106)+rrt(107) 
  ydot(18) = -rrt(096)-rrt(107) 
  ydot(19) = 0.0d0
  ydot(20) = 0.0d0
  ydot(21) = 0.0d0
  ydot(22) = 0.0d0
  ydot(23) = 0.0d0
  ydot(24) = 0.0d0
  ydot(25) = 0.0d0
  ydot(26) = 0.0d0
  ydot(27) = 0.0d0
  ydot(28) = 0.0d0
  ydot(29) = 0.0d0
  ydot(30) = 0.0d0
  ydot(31) = 0.0d0
  ydot(32) = 0.0d0
  ydot(33) = 0.0d0
  ydot(34) = 0.0d0
  ydot(35) = -rrt(097)+rrt(098)-rrt(108)+rrt(109)-rrt(132)+rrt(158) 
  ydot(36) = -rrt(098)+rrt(099)-rrt(109)+rrt(110)-rrt(133) 
  ydot(37) = -rrt(099)-rrt(110)-rrt(134) 
  ydot(38) = -rrt(013)-rrt(014)-rrt(015)+rrt(016)+rrt(024)+rrt(026)+rrt(035)+rrt(036)-rrt(042)-rrt(043)+rrt(044)-rrt(053)-rrt(111)&
             -rrt(130)+rrt(139) 
  ydot(39) = -rrt(016)+rrt(025)+rrt(027)+rrt(037)+rrt(038)-rrt(044)+rrt(046)-rrt(113)+rrt(140) 
  ydot(40) = -rrt(017)-rrt(018)-  2.d0 * rrt(041)-rrt(042)-rrt(045)+rrt(141) 
  ydot(41) = -rrt(046)+rrt(142) 
  ydot(42) = -rrt(019)+rrt(135) 
  ydot(43) = +rrt(136) 
  ydot(44) = +rrt(137) 
  ydot(45) = +rrt(138) 
  ydot(46) = -rrt(020)-rrt(021) 
  ydot(47) = -rrt(022) 
  ydot(48) = +rrt(040)+rrt(041)+rrt(042)-rrt(051)-rrt(052)-rrt(053)-rrt(054)+rrt(055)+rrt(057)-rrt(073)-rrt(075)+rrt(076)-rrt(080)&
             -rrt(086)+rrt(143) 
  ydot(49) = +rrt(053)-rrt(055)+rrt(073)+rrt(074) 
  ydot(50) = -rrt(056)-rrt(057)+rrt(075)-rrt(081)-rrt(087) 
  ydot(51) = -rrt(047)-rrt(048)-rrt(049)-rrt(050)+rrt(051)+rrt(056)-rrt(074)-rrt(076)+rrt(145)+rrt(150) 
  ydot(52) = -rrt(058)+rrt(059)-rrt(077)-rrt(083)+rrt(146)+rrt(151)+rrt(152)+rrt(153)+rrt(154) 
  ydot(53) = -rrt(059)-rrt(060)-rrt(061)-rrt(062)-rrt(078)-rrt(084)+rrt(144) 
  ydot(54) = +rrt(060)+rrt(063)-rrt(079)-rrt(085) 
  ydot(55) = -rrt(077)-rrt(078)-rrt(079)-rrt(080)-rrt(081)-rrt(082)-rrt(083)-rrt(084)-rrt(085)-rrt(086)-rrt(087)-rrt(088) 
  ydot(56) = +rrt(047)-rrt(063)-rrt(064)-rrt(065)-rrt(066)-rrt(067)+rrt(147) 
  ydot(57) = +rrt(048)+rrt(064)-rrt(068)-rrt(069)-rrt(070)+rrt(148) 
  ydot(58) = +rrt(049)+rrt(054)+rrt(058)+rrt(062)+rrt(065)+rrt(068)+rrt(069)-rrt(071)+rrt(149) 
  ydot(59) = +rrt(066)+rrt(070)+rrt(071)+rrt(072) 
  ydot(60) = +rrt(050)+rrt(052)+rrt(061)+rrt(067)-rrt(072)-rrt(082)-rrt(088) 
  ydot(61) = -rrt(112)-rrt(115)-rrt(116)-rrt(117)+rrt(119)+rrt(120)+rrt(124)+rrt(125)+rrt(126)+rrt(127)+  2.d0 * rrt(128)&
             -  2.d0 * rrt(129)-  2.d0 * rrt(130)-  2.d0 * rrt(131)-  2.d0 * rrt(132)-  2.d0 * rrt(133)-  2.d0 * rrt(134) 
  ydot(62) = +rrt(115)-rrt(122)-rrt(126)+  2.d0 * rrt(129)+  2.d0 * rrt(130) 
  ydot(63) = +rrt(112)-rrt(118)-rrt(119)-rrt(120)-rrt(121)-rrt(125)-rrt(126)-rrt(127)-rrt(128)+  2.d0 * rrt(131)+  2.d0 * rrt(132)&
             +  2.d0 * rrt(133)+  2.d0 * rrt(134) 
  ydot(64) = +rrt(116)+rrt(118)+rrt(122)-rrt(123)+rrt(126)-rrt(127) 
  ydot(65) = +rrt(117)+rrt(121)+rrt(123)-rrt(124)+rrt(127)-rrt(128) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(66) = 0.0d0
  if( lgas_heating ) then
    ydot(66) = ( ZDPlasKin_cfg(14)/k_B + ydot(66) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(66) = ydot(66) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(66)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(02,04) = pd(02,04) + rrt(001) * density(06) 
  pd(02,06) = pd(02,06) + rrt(001) * density(04) 
  pd(04,04) = pd(04,04) - rrt(001) * density(06) 
  pd(04,06) = pd(04,06) - rrt(001) * density(04) 
  pd(05,04) = pd(05,04) + rrt(001) * density(06) 
  pd(05,06) = pd(05,06) + rrt(001) * density(04) 
  pd(06,04) = pd(06,04) - rrt(001) * density(06) 
  pd(06,06) = pd(06,06) - rrt(001) * density(04) 
  pd(03,05) = pd(03,05) + rrt(002) * density(06) 
  pd(03,06) = pd(03,06) + rrt(002) * density(05) 
  pd(04,05) = pd(04,05) + rrt(002) * density(06) 
  pd(04,06) = pd(04,06) + rrt(002) * density(05) 
  pd(05,05) = pd(05,05) - rrt(002) * density(06) 
  pd(05,06) = pd(05,06) - rrt(002) * density(05) 
  pd(06,05) = pd(06,05) - rrt(002) * density(06) 
  pd(06,06) = pd(06,06) - rrt(002) * density(05) 
  pd(02,06) = pd(02,06) + rrt(003) * density(06) * 2.0d0
  pd(03,06) = pd(03,06) + rrt(003) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(003) * density(06) * 4.0d0
  pd(04,06) = pd(04,06) + rrt(004) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(004) * density(06) * 4.0d0
  pd(07,06) = pd(07,06) + rrt(004) * density(06) * 2.0d0
  pd(02,06) = pd(02,06) + rrt(005) * density(06) * 2.0d0
  pd(05,06) = pd(05,06) + rrt(005) * density(06) * 4.0d0
  pd(06,06) = pd(06,06) - rrt(005) * density(06) * 4.0d0
  pd(03,05) = pd(03,05) + rrt(006) * density(07) 
  pd(03,07) = pd(03,07) + rrt(006) * density(05) 
  pd(05,05) = pd(05,05) - rrt(006) * density(07) 
  pd(05,07) = pd(05,07) - rrt(006) * density(05) 
  pd(06,05) = pd(06,05) + rrt(006) * density(07) 
  pd(06,07) = pd(06,07) + rrt(006) * density(05) 
  pd(07,05) = pd(07,05) - rrt(006) * density(07) 
  pd(07,07) = pd(07,07) - rrt(006) * density(05) 
  pd(02,04) = pd(02,04) + rrt(007) * density(07) 
  pd(02,07) = pd(02,07) + rrt(007) * density(04) 
  pd(04,04) = pd(04,04) - rrt(007) * density(07) 
  pd(04,07) = pd(04,07) - rrt(007) * density(04) 
  pd(05,04) = pd(05,04) + rrt(007) * density(07) * 2.0d0
  pd(05,07) = pd(05,07) + rrt(007) * density(04) * 2.0d0
  pd(07,04) = pd(07,04) - rrt(007) * density(07) 
  pd(07,07) = pd(07,07) - rrt(007) * density(04) 
  pd(02,04) = pd(02,04) + rrt(008) * density(07) 
  pd(02,07) = pd(02,07) + rrt(008) * density(04) 
  pd(03,04) = pd(03,04) + rrt(008) * density(07) 
  pd(03,07) = pd(03,07) + rrt(008) * density(04) 
  pd(04,04) = pd(04,04) - rrt(008) * density(07) 
  pd(04,07) = pd(04,07) - rrt(008) * density(04) 
  pd(07,04) = pd(07,04) - rrt(008) * density(07) 
  pd(07,07) = pd(07,07) - rrt(008) * density(04) 
  pd(04,06) = pd(04,06) + rrt(009) * density(07) 
  pd(04,07) = pd(04,07) + rrt(009) * density(06) 
  pd(06,06) = pd(06,06) - rrt(009) * density(07) 
  pd(06,07) = pd(06,07) - rrt(009) * density(06) 
  pd(07,06) = pd(07,06) - rrt(009) * density(07) 
  pd(07,07) = pd(07,07) - rrt(009) * density(06) 
  pd(08,06) = pd(08,06) + rrt(009) * density(07) 
  pd(08,07) = pd(08,07) + rrt(009) * density(06) 
  pd(03,03) = pd(03,03) - rrt(010) * density(04) 
  pd(03,04) = pd(03,04) - rrt(010) * density(03) 
  pd(04,03) = pd(04,03) - rrt(010) * density(04) 
  pd(04,04) = pd(04,04) - rrt(010) * density(03) 
  pd(05,03) = pd(05,03) + rrt(010) * density(04) 
  pd(05,04) = pd(05,04) + rrt(010) * density(03) 
  pd(06,03) = pd(06,03) + rrt(010) * density(04) 
  pd(06,04) = pd(06,04) + rrt(010) * density(03) 
  pd(03,03) = pd(03,03) - rrt(011) * density(07) 
  pd(03,07) = pd(03,07) - rrt(011) * density(03) 
  pd(05,03) = pd(05,03) + rrt(011) * density(07) 
  pd(05,07) = pd(05,07) + rrt(011) * density(03) 
  pd(07,03) = pd(07,03) - rrt(011) * density(07) 
  pd(07,07) = pd(07,07) - rrt(011) * density(03) 
  pd(08,03) = pd(08,03) + rrt(011) * density(07) 
  pd(08,07) = pd(08,07) + rrt(011) * density(03) 
  pd(03,05) = pd(03,05) + rrt(012) * density(08) 
  pd(03,08) = pd(03,08) + rrt(012) * density(05) 
  pd(05,05) = pd(05,05) - rrt(012) * density(08) 
  pd(05,08) = pd(05,08) - rrt(012) * density(05) 
  pd(07,05) = pd(07,05) + rrt(012) * density(08) 
  pd(07,08) = pd(07,08) + rrt(012) * density(05) 
  pd(08,05) = pd(08,05) - rrt(012) * density(08) 
  pd(08,08) = pd(08,08) - rrt(012) * density(05) 
  pd(02,05) = pd(02,05) + rrt(013) * density(38) 
  pd(02,38) = pd(02,38) + rrt(013) * density(05) 
  pd(38,05) = pd(38,05) - rrt(013) * density(38) 
  pd(38,38) = pd(38,38) - rrt(013) * density(05) 
  pd(02,03) = pd(02,03) + rrt(014) * density(38) 
  pd(02,38) = pd(02,38) + rrt(014) * density(03) 
  pd(03,03) = pd(03,03) - rrt(014) * density(38) 
  pd(03,38) = pd(03,38) - rrt(014) * density(03) 
  pd(05,03) = pd(05,03) + rrt(014) * density(38) * 2.0d0
  pd(05,38) = pd(05,38) + rrt(014) * density(03) * 2.0d0
  pd(38,03) = pd(38,03) - rrt(014) * density(38) 
  pd(38,38) = pd(38,38) - rrt(014) * density(03) 
  pd(02,08) = pd(02,08) + rrt(015) * density(38) 
  pd(02,38) = pd(02,38) + rrt(015) * density(08) 
  pd(38,08) = pd(38,08) - rrt(015) * density(38) 
  pd(38,38) = pd(38,38) - rrt(015) * density(08) 
  pd(38,03) = pd(38,03) + rrt(016) * density(39) 
  pd(38,39) = pd(38,39) + rrt(016) * density(03) 
  pd(39,03) = pd(39,03) - rrt(016) * density(39) 
  pd(39,39) = pd(39,39) - rrt(016) * density(03) 
  pd(02,05) = pd(02,05) + rrt(017) * density(40) 
  pd(02,40) = pd(02,40) + rrt(017) * density(05) 
  pd(40,05) = pd(40,05) - rrt(017) * density(40) 
  pd(40,40) = pd(40,40) - rrt(017) * density(05) 
  pd(02,03) = pd(02,03) + rrt(018) * density(40) 
  pd(02,40) = pd(02,40) + rrt(018) * density(03) 
  pd(03,03) = pd(03,03) - rrt(018) * density(40) 
  pd(03,40) = pd(03,40) - rrt(018) * density(03) 
  pd(05,03) = pd(05,03) + rrt(018) * density(40) * 2.0d0
  pd(05,40) = pd(05,40) + rrt(018) * density(03) * 2.0d0
  pd(40,03) = pd(40,03) - rrt(018) * density(40) 
  pd(40,40) = pd(40,40) - rrt(018) * density(03) 
  pd(04,04) = pd(04,04) - rrt(019) * density(42) 
  pd(04,42) = pd(04,42) - rrt(019) * density(04) 
  pd(05,04) = pd(05,04) + rrt(019) * density(42) 
  pd(05,42) = pd(05,42) + rrt(019) * density(04) 
  pd(06,04) = pd(06,04) + rrt(019) * density(42) 
  pd(06,42) = pd(06,42) + rrt(019) * density(04) 
  pd(42,04) = pd(42,04) - rrt(019) * density(42) 
  pd(42,42) = pd(42,42) - rrt(019) * density(04) 
  pd(03,03) = pd(03,03) - rrt(020) * density(46) 
  pd(03,46) = pd(03,46) - rrt(020) * density(03) 
  pd(05,03) = pd(05,03) + rrt(020) * density(46) 
  pd(05,46) = pd(05,46) + rrt(020) * density(03) 
  pd(06,03) = pd(06,03) + rrt(020) * density(46) 
  pd(06,46) = pd(06,46) + rrt(020) * density(03) 
  pd(46,03) = pd(46,03) - rrt(020) * density(46) 
  pd(46,46) = pd(46,46) - rrt(020) * density(03) 
  pd(06,08) = pd(06,08) + rrt(021) * density(46) 
  pd(06,46) = pd(06,46) + rrt(021) * density(08) 
  pd(07,08) = pd(07,08) + rrt(021) * density(46) 
  pd(07,46) = pd(07,46) + rrt(021) * density(08) 
  pd(08,08) = pd(08,08) - rrt(021) * density(46) 
  pd(08,46) = pd(08,46) - rrt(021) * density(08) 
  pd(46,08) = pd(46,08) - rrt(021) * density(46) 
  pd(46,46) = pd(46,46) - rrt(021) * density(08) 
  pd(03,03) = pd(03,03) - rrt(022) * density(47) 
  pd(03,47) = pd(03,47) - rrt(022) * density(03) 
  pd(05,03) = pd(05,03) + rrt(022) * density(47) 
  pd(05,47) = pd(05,47) + rrt(022) * density(03) 
  pd(06,03) = pd(06,03) + rrt(022) * density(47) 
  pd(06,47) = pd(06,47) + rrt(022) * density(03) 
  pd(47,03) = pd(47,03) - rrt(022) * density(47) 
  pd(47,47) = pd(47,47) - rrt(022) * density(03) 
  pd(02,04) = pd(02,04) + rrt(023) * density(04) * density(10) * 2.0d0
  pd(02,10) = pd(02,10) + rrt(023) * density(04)**2 
  pd(04,04) = pd(04,04) - rrt(023) * density(04) * density(10) * 4.0d0
  pd(04,10) = pd(04,10) - rrt(023) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(024) * density(04)**2 * 6.0d0
  pd(38,04) = pd(38,04) + rrt(024) * density(04)**2 * 3.0d0
  pd(04,04) = pd(04,04) - rrt(025) * density(04)**2 * 6.0d0
  pd(39,04) = pd(39,04) + rrt(025) * density(04)**2 * 3.0d0
  pd(04,02) = pd(04,02) - rrt(026) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(026) * density(02) * density(04) * 4.0d0
  pd(38,02) = pd(38,02) + rrt(026) * density(04)**2 
  pd(38,04) = pd(38,04) + rrt(026) * density(02) * density(04) * 2.0d0
  pd(04,02) = pd(04,02) - rrt(027) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(027) * density(02) * density(04) * 4.0d0
  pd(39,02) = pd(39,02) + rrt(027) * density(04)**2 
  pd(39,04) = pd(39,04) + rrt(027) * density(02) * density(04) * 2.0d0
  pd(02,03) = pd(02,03) + rrt(028) * density(04)**2 
  pd(02,04) = pd(02,04) + rrt(028) * density(03) * density(04) * 2.0d0
  pd(04,03) = pd(04,03) - rrt(028) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(028) * density(03) * density(04) * 4.0d0
  pd(03,02) = pd(03,02) + rrt(029) * density(05)**2 
  pd(03,05) = pd(03,05) + rrt(029) * density(02) * density(05) * 2.0d0
  pd(05,02) = pd(05,02) - rrt(029) * density(05)**2 * 2.0d0
  pd(05,05) = pd(05,05) - rrt(029) * density(02) * density(05) * 4.0d0
  pd(04,04) = pd(04,04) - rrt(030) * density(05) * density(10) 
  pd(04,05) = pd(04,05) - rrt(030) * density(04) * density(10) 
  pd(04,10) = pd(04,10) - rrt(030) * density(04) * density(05) 
  pd(05,04) = pd(05,04) - rrt(030) * density(05) * density(10) 
  pd(05,05) = pd(05,05) - rrt(030) * density(04) * density(10) 
  pd(05,10) = pd(05,10) - rrt(030) * density(04) * density(05) 
  pd(06,04) = pd(06,04) + rrt(030) * density(05) * density(10) 
  pd(06,05) = pd(06,05) + rrt(030) * density(04) * density(10) 
  pd(06,10) = pd(06,10) + rrt(030) * density(04) * density(05) 
  pd(03,03) = pd(03,03) - rrt(031) * density(04) * density(10) 
  pd(03,04) = pd(03,04) - rrt(031) * density(03) * density(10) 
  pd(03,10) = pd(03,10) - rrt(031) * density(03) * density(04) 
  pd(04,03) = pd(04,03) - rrt(031) * density(04) * density(10) 
  pd(04,04) = pd(04,04) - rrt(031) * density(03) * density(10) 
  pd(04,10) = pd(04,10) - rrt(031) * density(03) * density(04) 
  pd(07,03) = pd(07,03) + rrt(031) * density(04) * density(10) 
  pd(07,04) = pd(07,04) + rrt(031) * density(03) * density(10) 
  pd(07,10) = pd(07,10) + rrt(031) * density(03) * density(04) 
  pd(05,05) = pd(05,05) - rrt(032) * density(06) * density(10) 
  pd(05,06) = pd(05,06) - rrt(032) * density(05) * density(10) 
  pd(05,10) = pd(05,10) - rrt(032) * density(05) * density(06) 
  pd(06,05) = pd(06,05) - rrt(032) * density(06) * density(10) 
  pd(06,06) = pd(06,06) - rrt(032) * density(05) * density(10) 
  pd(06,10) = pd(06,10) - rrt(032) * density(05) * density(06) 
  pd(07,05) = pd(07,05) + rrt(032) * density(06) * density(10) 
  pd(07,06) = pd(07,06) + rrt(032) * density(05) * density(10) 
  pd(07,10) = pd(07,10) + rrt(032) * density(05) * density(06) 
  pd(05,05) = pd(05,05) - rrt(033) * density(07) * density(10) 
  pd(05,07) = pd(05,07) - rrt(033) * density(05) * density(10) 
  pd(05,10) = pd(05,10) - rrt(033) * density(05) * density(07) 
  pd(07,05) = pd(07,05) - rrt(033) * density(07) * density(10) 
  pd(07,07) = pd(07,07) - rrt(033) * density(05) * density(10) 
  pd(07,10) = pd(07,10) - rrt(033) * density(05) * density(07) 
  pd(08,05) = pd(08,05) + rrt(033) * density(07) * density(10) 
  pd(08,07) = pd(08,07) + rrt(033) * density(05) * density(10) 
  pd(08,10) = pd(08,10) + rrt(033) * density(05) * density(07) 
  pd(03,03) = pd(03,03) - rrt(034) * density(06) * density(10) 
  pd(03,06) = pd(03,06) - rrt(034) * density(03) * density(10) 
  pd(03,10) = pd(03,10) - rrt(034) * density(03) * density(06) 
  pd(06,03) = pd(06,03) - rrt(034) * density(06) * density(10) 
  pd(06,06) = pd(06,06) - rrt(034) * density(03) * density(10) 
  pd(06,10) = pd(06,10) - rrt(034) * density(03) * density(06) 
  pd(08,03) = pd(08,03) + rrt(034) * density(06) * density(10) 
  pd(08,06) = pd(08,06) + rrt(034) * density(03) * density(10) 
  pd(08,10) = pd(08,10) + rrt(034) * density(03) * density(06) 
  pd(04,03) = pd(04,03) - rrt(035) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(035) * density(03) * density(04) * 4.0d0
  pd(38,03) = pd(38,03) + rrt(035) * density(04)**2 
  pd(38,04) = pd(38,04) + rrt(035) * density(03) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(036) * density(04) * density(05) * 4.0d0
  pd(04,05) = pd(04,05) - rrt(036) * density(04)**2 * 2.0d0
  pd(38,04) = pd(38,04) + rrt(036) * density(04) * density(05) * 2.0d0
  pd(38,05) = pd(38,05) + rrt(036) * density(04)**2 
  pd(04,03) = pd(04,03) - rrt(037) * density(04)**2 * 2.0d0
  pd(04,04) = pd(04,04) - rrt(037) * density(03) * density(04) * 4.0d0
  pd(39,03) = pd(39,03) + rrt(037) * density(04)**2 
  pd(39,04) = pd(39,04) + rrt(037) * density(03) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(038) * density(04) * density(05) * 4.0d0
  pd(04,05) = pd(04,05) - rrt(038) * density(04)**2 * 2.0d0
  pd(39,04) = pd(39,04) + rrt(038) * density(04) * density(05) * 2.0d0
  pd(39,05) = pd(39,05) + rrt(038) * density(04)**2 
  pd(03,03) = pd(03,03) + rrt(039) * density(05)**2 
  pd(03,05) = pd(03,05) + rrt(039) * density(03) * density(05) * 2.0d0
  pd(05,03) = pd(05,03) - rrt(039) * density(05)**2 * 2.0d0
  pd(05,05) = pd(05,05) - rrt(039) * density(03) * density(05) * 4.0d0
  pd(01,04) = pd(01,04) + rrt(040) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(040) * density(04) * 4.0d0
  pd(48,04) = pd(48,04) + rrt(040) * density(04) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(041) * density(40) * 2.0d0
  pd(02,40) = pd(02,40) + rrt(041) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(041) * density(40) * 4.0d0
  pd(48,40) = pd(48,40) + rrt(041) * density(40) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(042) * density(40) 
  pd(01,40) = pd(01,40) + rrt(042) * density(38) 
  pd(02,38) = pd(02,38) + rrt(042) * density(40) 
  pd(02,40) = pd(02,40) + rrt(042) * density(38) 
  pd(38,38) = pd(38,38) - rrt(042) * density(40) 
  pd(38,40) = pd(38,40) - rrt(042) * density(38) 
  pd(40,38) = pd(40,38) - rrt(042) * density(40) 
  pd(40,40) = pd(40,40) - rrt(042) * density(38) 
  pd(48,38) = pd(48,38) + rrt(042) * density(40) 
  pd(48,40) = pd(48,40) + rrt(042) * density(38) 
  pd(02,38) = pd(02,38) + rrt(043) 
  pd(38,38) = pd(38,38) - rrt(043) 
  pd(38,39) = pd(38,39) + rrt(044) 
  pd(39,39) = pd(39,39) - rrt(044) 
  pd(02,40) = pd(02,40) + rrt(045) 
  pd(40,40) = pd(40,40) - rrt(045) 
  pd(39,41) = pd(39,41) + rrt(046) 
  pd(41,41) = pd(41,41) - rrt(046) 
  pd(03,03) = pd(03,03) - rrt(047) * density(51) 
  pd(03,51) = pd(03,51) - rrt(047) * density(03) 
  pd(05,03) = pd(05,03) + rrt(047) * density(51) 
  pd(05,51) = pd(05,51) + rrt(047) * density(03) 
  pd(51,03) = pd(51,03) - rrt(047) * density(51) 
  pd(51,51) = pd(51,51) - rrt(047) * density(03) 
  pd(56,03) = pd(56,03) + rrt(047) * density(51) 
  pd(56,51) = pd(56,51) + rrt(047) * density(03) 
  pd(06,08) = pd(06,08) + rrt(048) * density(51) 
  pd(06,51) = pd(06,51) + rrt(048) * density(08) 
  pd(08,08) = pd(08,08) - rrt(048) * density(51) 
  pd(08,51) = pd(08,51) - rrt(048) * density(08) 
  pd(51,08) = pd(51,08) - rrt(048) * density(51) 
  pd(51,51) = pd(51,51) - rrt(048) * density(08) 
  pd(57,08) = pd(57,08) + rrt(048) * density(51) 
  pd(57,51) = pd(57,51) + rrt(048) * density(08) 
  pd(04,08) = pd(04,08) + rrt(049) * density(51) 
  pd(04,51) = pd(04,51) + rrt(049) * density(08) 
  pd(08,08) = pd(08,08) - rrt(049) * density(51) 
  pd(08,51) = pd(08,51) - rrt(049) * density(08) 
  pd(51,08) = pd(51,08) - rrt(049) * density(51) 
  pd(51,51) = pd(51,51) - rrt(049) * density(08) 
  pd(58,08) = pd(58,08) + rrt(049) * density(51) 
  pd(58,51) = pd(58,51) + rrt(049) * density(08) 
  pd(03,08) = pd(03,08) + rrt(050) * density(51) 
  pd(03,51) = pd(03,51) + rrt(050) * density(08) 
  pd(08,08) = pd(08,08) - rrt(050) * density(51) 
  pd(08,51) = pd(08,51) - rrt(050) * density(08) 
  pd(51,08) = pd(51,08) - rrt(050) * density(51) 
  pd(51,51) = pd(51,51) - rrt(050) * density(08) 
  pd(60,08) = pd(60,08) + rrt(050) * density(51) 
  pd(60,51) = pd(60,51) + rrt(050) * density(08) 
  pd(02,04) = pd(02,04) + rrt(051) * density(48) 
  pd(02,48) = pd(02,48) + rrt(051) * density(04) 
  pd(04,04) = pd(04,04) - rrt(051) * density(48) 
  pd(04,48) = pd(04,48) - rrt(051) * density(04) 
  pd(48,04) = pd(48,04) - rrt(051) * density(48) 
  pd(48,48) = pd(48,48) - rrt(051) * density(04) 
  pd(51,04) = pd(51,04) + rrt(051) * density(48) 
  pd(51,48) = pd(51,48) + rrt(051) * density(04) 
  pd(03,03) = pd(03,03) - rrt(052) * density(48) 
  pd(03,48) = pd(03,48) - rrt(052) * density(03) 
  pd(05,03) = pd(05,03) + rrt(052) * density(48) 
  pd(05,48) = pd(05,48) + rrt(052) * density(03) 
  pd(48,03) = pd(48,03) - rrt(052) * density(48) 
  pd(48,48) = pd(48,48) - rrt(052) * density(03) 
  pd(60,03) = pd(60,03) + rrt(052) * density(48) 
  pd(60,48) = pd(60,48) + rrt(052) * density(03) 
  pd(04,38) = pd(04,38) + rrt(053) * density(48) 
  pd(04,48) = pd(04,48) + rrt(053) * density(38) 
  pd(38,38) = pd(38,38) - rrt(053) * density(48) 
  pd(38,48) = pd(38,48) - rrt(053) * density(38) 
  pd(48,38) = pd(48,38) - rrt(053) * density(48) 
  pd(48,48) = pd(48,48) - rrt(053) * density(38) 
  pd(49,38) = pd(49,38) + rrt(053) * density(48) 
  pd(49,48) = pd(49,48) + rrt(053) * density(38) 
  pd(02,08) = pd(02,08) + rrt(054) * density(48) 
  pd(02,48) = pd(02,48) + rrt(054) * density(08) 
  pd(08,08) = pd(08,08) - rrt(054) * density(48) 
  pd(08,48) = pd(08,48) - rrt(054) * density(08) 
  pd(48,08) = pd(48,08) - rrt(054) * density(48) 
  pd(48,48) = pd(48,48) - rrt(054) * density(08) 
  pd(58,08) = pd(58,08) + rrt(054) * density(48) 
  pd(58,48) = pd(58,48) + rrt(054) * density(08) 
  pd(02,04) = pd(02,04) + rrt(055) * density(49) 
  pd(02,49) = pd(02,49) + rrt(055) * density(04) 
  pd(04,04) = pd(04,04) - rrt(055) * density(49) 
  pd(04,49) = pd(04,49) - rrt(055) * density(04) 
  pd(48,04) = pd(48,04) + rrt(055) * density(49) 
  pd(48,49) = pd(48,49) + rrt(055) * density(04) 
  pd(49,04) = pd(49,04) - rrt(055) * density(49) 
  pd(49,49) = pd(49,49) - rrt(055) * density(04) 
  pd(02,04) = pd(02,04) + rrt(056) * density(50) * 2.0d0
  pd(02,50) = pd(02,50) + rrt(056) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(056) * density(50) 
  pd(04,50) = pd(04,50) - rrt(056) * density(04) 
  pd(50,04) = pd(50,04) - rrt(056) * density(50) 
  pd(50,50) = pd(50,50) - rrt(056) * density(04) 
  pd(51,04) = pd(51,04) + rrt(056) * density(50) 
  pd(51,50) = pd(51,50) + rrt(056) * density(04) 
  pd(02,02) = pd(02,02) + rrt(057) * density(50) 
  pd(02,50) = pd(02,50) + rrt(057) * density(02) 
  pd(48,02) = pd(48,02) + rrt(057) * density(50) 
  pd(48,50) = pd(48,50) + rrt(057) * density(02) 
  pd(50,02) = pd(50,02) - rrt(057) * density(50) 
  pd(50,50) = pd(50,50) - rrt(057) * density(02) 
  pd(05,08) = pd(05,08) + rrt(058) * density(52) 
  pd(05,52) = pd(05,52) + rrt(058) * density(08) 
  pd(08,08) = pd(08,08) - rrt(058) * density(52) 
  pd(08,52) = pd(08,52) - rrt(058) * density(08) 
  pd(52,08) = pd(52,08) - rrt(058) * density(52) 
  pd(52,52) = pd(52,52) - rrt(058) * density(08) 
  pd(58,08) = pd(58,08) + rrt(058) * density(52) 
  pd(58,52) = pd(58,52) + rrt(058) * density(08) 
  pd(03,05) = pd(03,05) + rrt(059) * density(53) 
  pd(03,53) = pd(03,53) + rrt(059) * density(05) 
  pd(05,05) = pd(05,05) - rrt(059) * density(53) 
  pd(05,53) = pd(05,53) - rrt(059) * density(05) 
  pd(52,05) = pd(52,05) + rrt(059) * density(53) 
  pd(52,53) = pd(52,53) + rrt(059) * density(05) 
  pd(53,05) = pd(53,05) - rrt(059) * density(53) 
  pd(53,53) = pd(53,53) - rrt(059) * density(05) 
  pd(03,03) = pd(03,03) - rrt(060) * density(53) 
  pd(03,53) = pd(03,53) - rrt(060) * density(03) 
  pd(05,03) = pd(05,03) + rrt(060) * density(53) 
  pd(05,53) = pd(05,53) + rrt(060) * density(03) 
  pd(53,03) = pd(53,03) - rrt(060) * density(53) 
  pd(53,53) = pd(53,53) - rrt(060) * density(03) 
  pd(54,03) = pd(54,03) + rrt(060) * density(53) 
  pd(54,53) = pd(54,53) + rrt(060) * density(03) 
  pd(02,02) = pd(02,02) - rrt(061) * density(53) 
  pd(02,53) = pd(02,53) - rrt(061) * density(02) 
  pd(05,02) = pd(05,02) + rrt(061) * density(53) 
  pd(05,53) = pd(05,53) + rrt(061) * density(02) 
  pd(53,02) = pd(53,02) - rrt(061) * density(53) 
  pd(53,53) = pd(53,53) - rrt(061) * density(02) 
  pd(60,02) = pd(60,02) + rrt(061) * density(53) 
  pd(60,53) = pd(60,53) + rrt(061) * density(02) 
  pd(03,08) = pd(03,08) + rrt(062) * density(53) 
  pd(03,53) = pd(03,53) + rrt(062) * density(08) 
  pd(08,08) = pd(08,08) - rrt(062) * density(53) 
  pd(08,53) = pd(08,53) - rrt(062) * density(08) 
  pd(53,08) = pd(53,08) - rrt(062) * density(53) 
  pd(53,53) = pd(53,53) - rrt(062) * density(08) 
  pd(58,08) = pd(58,08) + rrt(062) * density(53) 
  pd(58,53) = pd(58,53) + rrt(062) * density(08) 
  pd(03,03) = pd(03,03) - rrt(063) * density(56) 
  pd(03,56) = pd(03,56) - rrt(063) * density(03) 
  pd(04,03) = pd(04,03) + rrt(063) * density(56) 
  pd(04,56) = pd(04,56) + rrt(063) * density(03) 
  pd(54,03) = pd(54,03) + rrt(063) * density(56) 
  pd(54,56) = pd(54,56) + rrt(063) * density(03) 
  pd(56,03) = pd(56,03) - rrt(063) * density(56) 
  pd(56,56) = pd(56,56) - rrt(063) * density(03) 
  pd(03,03) = pd(03,03) - rrt(064) * density(56) 
  pd(03,56) = pd(03,56) - rrt(064) * density(03) 
  pd(05,03) = pd(05,03) + rrt(064) * density(56) 
  pd(05,56) = pd(05,56) + rrt(064) * density(03) 
  pd(56,03) = pd(56,03) - rrt(064) * density(56) 
  pd(56,56) = pd(56,56) - rrt(064) * density(03) 
  pd(57,03) = pd(57,03) + rrt(064) * density(56) 
  pd(57,56) = pd(57,56) + rrt(064) * density(03) 
  pd(06,08) = pd(06,08) + rrt(065) * density(56) 
  pd(06,56) = pd(06,56) + rrt(065) * density(08) 
  pd(08,08) = pd(08,08) - rrt(065) * density(56) 
  pd(08,56) = pd(08,56) - rrt(065) * density(08) 
  pd(56,08) = pd(56,08) - rrt(065) * density(56) 
  pd(56,56) = pd(56,56) - rrt(065) * density(08) 
  pd(58,08) = pd(58,08) + rrt(065) * density(56) 
  pd(58,56) = pd(58,56) + rrt(065) * density(08) 
  pd(04,08) = pd(04,08) + rrt(066) * density(56) 
  pd(04,56) = pd(04,56) + rrt(066) * density(08) 
  pd(08,08) = pd(08,08) - rrt(066) * density(56) 
  pd(08,56) = pd(08,56) - rrt(066) * density(08) 
  pd(56,08) = pd(56,08) - rrt(066) * density(56) 
  pd(56,56) = pd(56,56) - rrt(066) * density(08) 
  pd(59,08) = pd(59,08) + rrt(066) * density(56) 
  pd(59,56) = pd(59,56) + rrt(066) * density(08) 
  pd(02,02) = pd(02,02) - rrt(067) * density(56) 
  pd(02,56) = pd(02,56) - rrt(067) * density(02) 
  pd(04,02) = pd(04,02) + rrt(067) * density(56) 
  pd(04,56) = pd(04,56) + rrt(067) * density(02) 
  pd(56,02) = pd(56,02) - rrt(067) * density(56) 
  pd(56,56) = pd(56,56) - rrt(067) * density(02) 
  pd(60,02) = pd(60,02) + rrt(067) * density(56) 
  pd(60,56) = pd(60,56) + rrt(067) * density(02) 
  pd(03,03) = pd(03,03) - rrt(068) * density(57) 
  pd(03,57) = pd(03,57) - rrt(068) * density(03) 
  pd(05,03) = pd(05,03) + rrt(068) * density(57) 
  pd(05,57) = pd(05,57) + rrt(068) * density(03) 
  pd(57,03) = pd(57,03) - rrt(068) * density(57) 
  pd(57,57) = pd(57,57) - rrt(068) * density(03) 
  pd(58,03) = pd(58,03) + rrt(068) * density(57) 
  pd(58,57) = pd(58,57) + rrt(068) * density(03) 
  pd(07,08) = pd(07,08) + rrt(069) * density(57) 
  pd(07,57) = pd(07,57) + rrt(069) * density(08) 
  pd(08,08) = pd(08,08) - rrt(069) * density(57) 
  pd(08,57) = pd(08,57) - rrt(069) * density(08) 
  pd(57,08) = pd(57,08) - rrt(069) * density(57) 
  pd(57,57) = pd(57,57) - rrt(069) * density(08) 
  pd(58,08) = pd(58,08) + rrt(069) * density(57) 
  pd(58,57) = pd(58,57) + rrt(069) * density(08) 
  pd(06,08) = pd(06,08) + rrt(070) * density(57) 
  pd(06,57) = pd(06,57) + rrt(070) * density(08) 
  pd(08,08) = pd(08,08) - rrt(070) * density(57) 
  pd(08,57) = pd(08,57) - rrt(070) * density(08) 
  pd(57,08) = pd(57,08) - rrt(070) * density(57) 
  pd(57,57) = pd(57,57) - rrt(070) * density(08) 
  pd(59,08) = pd(59,08) + rrt(070) * density(57) 
  pd(59,57) = pd(59,57) + rrt(070) * density(08) 
  pd(07,08) = pd(07,08) + rrt(071) * density(58) 
  pd(07,58) = pd(07,58) + rrt(071) * density(08) 
  pd(08,08) = pd(08,08) - rrt(071) * density(58) 
  pd(08,58) = pd(08,58) - rrt(071) * density(08) 
  pd(58,08) = pd(58,08) - rrt(071) * density(58) 
  pd(58,58) = pd(58,58) - rrt(071) * density(08) 
  pd(59,08) = pd(59,08) + rrt(071) * density(58) 
  pd(59,58) = pd(59,58) + rrt(071) * density(08) 
  pd(02,08) = pd(02,08) + rrt(072) * density(60) 
  pd(02,60) = pd(02,60) + rrt(072) * density(08) 
  pd(08,08) = pd(08,08) - rrt(072) * density(60) 
  pd(08,60) = pd(08,60) - rrt(072) * density(08) 
  pd(59,08) = pd(59,08) + rrt(072) * density(60) 
  pd(59,60) = pd(59,60) + rrt(072) * density(08) 
  pd(60,08) = pd(60,08) - rrt(072) * density(60) 
  pd(60,60) = pd(60,60) - rrt(072) * density(08) 
  pd(04,02) = pd(04,02) - rrt(073) * density(04) * density(48) 
  pd(04,04) = pd(04,04) - rrt(073) * density(02) * density(48) 
  pd(04,48) = pd(04,48) - rrt(073) * density(02) * density(04) 
  pd(48,02) = pd(48,02) - rrt(073) * density(04) * density(48) 
  pd(48,04) = pd(48,04) - rrt(073) * density(02) * density(48) 
  pd(48,48) = pd(48,48) - rrt(073) * density(02) * density(04) 
  pd(49,02) = pd(49,02) + rrt(073) * density(04) * density(48) 
  pd(49,04) = pd(49,04) + rrt(073) * density(02) * density(48) 
  pd(49,48) = pd(49,48) + rrt(073) * density(02) * density(04) 
  pd(02,02) = pd(02,02) - rrt(074) * density(02) * density(51) * 2.0d0
  pd(02,51) = pd(02,51) - rrt(074) * density(02)**2 
  pd(49,02) = pd(49,02) + rrt(074) * density(02) * density(51) * 2.0d0
  pd(49,51) = pd(49,51) + rrt(074) * density(02)**2 
  pd(51,02) = pd(51,02) - rrt(074) * density(02) * density(51) * 2.0d0
  pd(51,51) = pd(51,51) - rrt(074) * density(02)**2 
  pd(02,02) = pd(02,02) - rrt(075) * density(02) * density(48) * 2.0d0
  pd(02,48) = pd(02,48) - rrt(075) * density(02)**2 
  pd(48,02) = pd(48,02) - rrt(075) * density(02) * density(48) * 2.0d0
  pd(48,48) = pd(48,48) - rrt(075) * density(02)**2 
  pd(50,02) = pd(50,02) + rrt(075) * density(02) * density(48) * 2.0d0
  pd(50,48) = pd(50,48) + rrt(075) * density(02)**2 
  pd(04,02) = pd(04,02) - rrt(076) * density(04) * density(51) 
  pd(04,04) = pd(04,04) - rrt(076) * density(02) * density(51) 
  pd(04,51) = pd(04,51) - rrt(076) * density(02) * density(04) 
  pd(48,02) = pd(48,02) + rrt(076) * density(04) * density(51) 
  pd(48,04) = pd(48,04) + rrt(076) * density(02) * density(51) 
  pd(48,51) = pd(48,51) + rrt(076) * density(02) * density(04) 
  pd(51,02) = pd(51,02) - rrt(076) * density(04) * density(51) 
  pd(51,04) = pd(51,04) - rrt(076) * density(02) * density(51) 
  pd(51,51) = pd(51,51) - rrt(076) * density(02) * density(04) 
  pd(05,52) = pd(05,52) + rrt(077) * density(55) * 2.0d0
  pd(05,55) = pd(05,55) + rrt(077) * density(52) * 2.0d0
  pd(52,52) = pd(52,52) - rrt(077) * density(55) 
  pd(52,55) = pd(52,55) - rrt(077) * density(52) 
  pd(55,52) = pd(55,52) - rrt(077) * density(55) 
  pd(55,55) = pd(55,55) - rrt(077) * density(52) 
  pd(05,53) = pd(05,53) + rrt(078) * density(55) * 3.0d0
  pd(05,55) = pd(05,55) + rrt(078) * density(53) * 3.0d0
  pd(53,53) = pd(53,53) - rrt(078) * density(55) 
  pd(53,55) = pd(53,55) - rrt(078) * density(53) 
  pd(55,53) = pd(55,53) - rrt(078) * density(55) 
  pd(55,55) = pd(55,55) - rrt(078) * density(53) 
  pd(03,54) = pd(03,54) + rrt(079) * density(55) 
  pd(03,55) = pd(03,55) + rrt(079) * density(54) 
  pd(05,54) = pd(05,54) + rrt(079) * density(55) * 2.0d0
  pd(05,55) = pd(05,55) + rrt(079) * density(54) * 2.0d0
  pd(54,54) = pd(54,54) - rrt(079) * density(55) 
  pd(54,55) = pd(54,55) - rrt(079) * density(54) 
  pd(55,54) = pd(55,54) - rrt(079) * density(55) 
  pd(55,55) = pd(55,55) - rrt(079) * density(54) 
  pd(02,48) = pd(02,48) + rrt(080) * density(55) 
  pd(02,55) = pd(02,55) + rrt(080) * density(48) 
  pd(05,48) = pd(05,48) + rrt(080) * density(55) 
  pd(05,55) = pd(05,55) + rrt(080) * density(48) 
  pd(48,48) = pd(48,48) - rrt(080) * density(55) 
  pd(48,55) = pd(48,55) - rrt(080) * density(48) 
  pd(55,48) = pd(55,48) - rrt(080) * density(55) 
  pd(55,55) = pd(55,55) - rrt(080) * density(48) 
  pd(02,50) = pd(02,50) + rrt(081) * density(55) * 2.0d0
  pd(02,55) = pd(02,55) + rrt(081) * density(50) * 2.0d0
  pd(05,50) = pd(05,50) + rrt(081) * density(55) 
  pd(05,55) = pd(05,55) + rrt(081) * density(50) 
  pd(50,50) = pd(50,50) - rrt(081) * density(55) 
  pd(50,55) = pd(50,55) - rrt(081) * density(50) 
  pd(55,50) = pd(55,50) - rrt(081) * density(55) 
  pd(55,55) = pd(55,55) - rrt(081) * density(50) 
  pd(02,55) = pd(02,55) + rrt(082) * density(60) 
  pd(02,60) = pd(02,60) + rrt(082) * density(55) 
  pd(03,55) = pd(03,55) + rrt(082) * density(60) 
  pd(03,60) = pd(03,60) + rrt(082) * density(55) 
  pd(55,55) = pd(55,55) - rrt(082) * density(60) 
  pd(55,60) = pd(55,60) - rrt(082) * density(55) 
  pd(60,55) = pd(60,55) - rrt(082) * density(60) 
  pd(60,60) = pd(60,60) - rrt(082) * density(55) 
  pd(03,10) = pd(03,10) + rrt(083) * density(52) * density(55) 
  pd(03,52) = pd(03,52) + rrt(083) * density(10) * density(55) 
  pd(03,55) = pd(03,55) + rrt(083) * density(10) * density(52) 
  pd(52,10) = pd(52,10) - rrt(083) * density(52) * density(55) 
  pd(52,52) = pd(52,52) - rrt(083) * density(10) * density(55) 
  pd(52,55) = pd(52,55) - rrt(083) * density(10) * density(52) 
  pd(55,10) = pd(55,10) - rrt(083) * density(52) * density(55) 
  pd(55,52) = pd(55,52) - rrt(083) * density(10) * density(55) 
  pd(55,55) = pd(55,55) - rrt(083) * density(10) * density(52) 
  pd(03,10) = pd(03,10) + rrt(084) * density(53) * density(55) 
  pd(03,53) = pd(03,53) + rrt(084) * density(10) * density(55) 
  pd(03,55) = pd(03,55) + rrt(084) * density(10) * density(53) 
  pd(05,10) = pd(05,10) + rrt(084) * density(53) * density(55) 
  pd(05,53) = pd(05,53) + rrt(084) * density(10) * density(55) 
  pd(05,55) = pd(05,55) + rrt(084) * density(10) * density(53) 
  pd(53,10) = pd(53,10) - rrt(084) * density(53) * density(55) 
  pd(53,53) = pd(53,53) - rrt(084) * density(10) * density(55) 
  pd(53,55) = pd(53,55) - rrt(084) * density(10) * density(53) 
  pd(55,10) = pd(55,10) - rrt(084) * density(53) * density(55) 
  pd(55,53) = pd(55,53) - rrt(084) * density(10) * density(55) 
  pd(55,55) = pd(55,55) - rrt(084) * density(10) * density(53) 
  pd(03,10) = pd(03,10) + rrt(085) * density(54) * density(55) * 2.0d0
  pd(03,54) = pd(03,54) + rrt(085) * density(10) * density(55) * 2.0d0
  pd(03,55) = pd(03,55) + rrt(085) * density(10) * density(54) * 2.0d0
  pd(54,10) = pd(54,10) - rrt(085) * density(54) * density(55) 
  pd(54,54) = pd(54,54) - rrt(085) * density(10) * density(55) 
  pd(54,55) = pd(54,55) - rrt(085) * density(10) * density(54) 
  pd(55,10) = pd(55,10) - rrt(085) * density(54) * density(55) 
  pd(55,54) = pd(55,54) - rrt(085) * density(10) * density(55) 
  pd(55,55) = pd(55,55) - rrt(085) * density(10) * density(54) 
  pd(02,10) = pd(02,10) + rrt(086) * density(48) * density(55) 
  pd(02,48) = pd(02,48) + rrt(086) * density(10) * density(55) 
  pd(02,55) = pd(02,55) + rrt(086) * density(10) * density(48) 
  pd(05,10) = pd(05,10) + rrt(086) * density(48) * density(55) 
  pd(05,48) = pd(05,48) + rrt(086) * density(10) * density(55) 
  pd(05,55) = pd(05,55) + rrt(086) * density(10) * density(48) 
  pd(48,10) = pd(48,10) - rrt(086) * density(48) * density(55) 
  pd(48,48) = pd(48,48) - rrt(086) * density(10) * density(55) 
  pd(48,55) = pd(48,55) - rrt(086) * density(10) * density(48) 
  pd(55,10) = pd(55,10) - rrt(086) * density(48) * density(55) 
  pd(55,48) = pd(55,48) - rrt(086) * density(10) * density(55) 
  pd(55,55) = pd(55,55) - rrt(086) * density(10) * density(48) 
  pd(02,10) = pd(02,10) + rrt(087) * density(50) * density(55) * 2.0d0
  pd(02,50) = pd(02,50) + rrt(087) * density(10) * density(55) * 2.0d0
  pd(02,55) = pd(02,55) + rrt(087) * density(10) * density(50) * 2.0d0
  pd(05,10) = pd(05,10) + rrt(087) * density(50) * density(55) 
  pd(05,50) = pd(05,50) + rrt(087) * density(10) * density(55) 
  pd(05,55) = pd(05,55) + rrt(087) * density(10) * density(50) 
  pd(50,10) = pd(50,10) - rrt(087) * density(50) * density(55) 
  pd(50,50) = pd(50,50) - rrt(087) * density(10) * density(55) 
  pd(50,55) = pd(50,55) - rrt(087) * density(10) * density(50) 
  pd(55,10) = pd(55,10) - rrt(087) * density(50) * density(55) 
  pd(55,50) = pd(55,50) - rrt(087) * density(10) * density(55) 
  pd(55,55) = pd(55,55) - rrt(087) * density(10) * density(50) 
  pd(02,10) = pd(02,10) + rrt(088) * density(55) * density(60) 
  pd(02,55) = pd(02,55) + rrt(088) * density(10) * density(60) 
  pd(02,60) = pd(02,60) + rrt(088) * density(10) * density(55) 
  pd(03,10) = pd(03,10) + rrt(088) * density(55) * density(60) 
  pd(03,55) = pd(03,55) + rrt(088) * density(10) * density(60) 
  pd(03,60) = pd(03,60) + rrt(088) * density(10) * density(55) 
  pd(55,10) = pd(55,10) - rrt(088) * density(55) * density(60) 
  pd(55,55) = pd(55,55) - rrt(088) * density(10) * density(60) 
  pd(55,60) = pd(55,60) - rrt(088) * density(10) * density(55) 
  pd(60,10) = pd(60,10) - rrt(088) * density(55) * density(60) 
  pd(60,55) = pd(60,55) - rrt(088) * density(10) * density(60) 
  pd(60,60) = pd(60,60) - rrt(088) * density(10) * density(55) 
  pd(02,02) = pd(02,02) + rrt(089) * density(11) 
  pd(02,11) = pd(02,11) + rrt(089) * density(02) 
  pd(11,02) = pd(11,02) - rrt(089) * density(11) 
  pd(11,11) = pd(11,11) - rrt(089) * density(02) 
  pd(11,02) = pd(11,02) + rrt(090) * density(12) 
  pd(11,12) = pd(11,12) + rrt(090) * density(02) 
  pd(12,02) = pd(12,02) - rrt(090) * density(12) 
  pd(12,12) = pd(12,12) - rrt(090) * density(02) 
  pd(12,02) = pd(12,02) + rrt(091) * density(13) 
  pd(12,13) = pd(12,13) + rrt(091) * density(02) 
  pd(13,02) = pd(13,02) - rrt(091) * density(13) 
  pd(13,13) = pd(13,13) - rrt(091) * density(02) 
  pd(13,02) = pd(13,02) + rrt(092) * density(14) 
  pd(13,14) = pd(13,14) + rrt(092) * density(02) 
  pd(14,02) = pd(14,02) - rrt(092) * density(14) 
  pd(14,14) = pd(14,14) - rrt(092) * density(02) 
  pd(14,02) = pd(14,02) + rrt(093) * density(15) 
  pd(14,15) = pd(14,15) + rrt(093) * density(02) 
  pd(15,02) = pd(15,02) - rrt(093) * density(15) 
  pd(15,15) = pd(15,15) - rrt(093) * density(02) 
  pd(15,02) = pd(15,02) + rrt(094) * density(16) 
  pd(15,16) = pd(15,16) + rrt(094) * density(02) 
  pd(16,02) = pd(16,02) - rrt(094) * density(16) 
  pd(16,16) = pd(16,16) - rrt(094) * density(02) 
  pd(16,02) = pd(16,02) + rrt(095) * density(17) 
  pd(16,17) = pd(16,17) + rrt(095) * density(02) 
  pd(17,02) = pd(17,02) - rrt(095) * density(17) 
  pd(17,17) = pd(17,17) - rrt(095) * density(02) 
  pd(17,02) = pd(17,02) + rrt(096) * density(18) 
  pd(17,18) = pd(17,18) + rrt(096) * density(02) 
  pd(18,02) = pd(18,02) - rrt(096) * density(18) 
  pd(18,18) = pd(18,18) - rrt(096) * density(02) 
  pd(03,03) = pd(03,03) + rrt(097) * density(35) 
  pd(03,35) = pd(03,35) + rrt(097) * density(03) 
  pd(35,03) = pd(35,03) - rrt(097) * density(35) 
  pd(35,35) = pd(35,35) - rrt(097) * density(03) 
  pd(35,03) = pd(35,03) + rrt(098) * density(36) 
  pd(35,36) = pd(35,36) + rrt(098) * density(03) 
  pd(36,03) = pd(36,03) - rrt(098) * density(36) 
  pd(36,36) = pd(36,36) - rrt(098) * density(03) 
  pd(36,03) = pd(36,03) + rrt(099) * density(37) 
  pd(36,37) = pd(36,37) + rrt(099) * density(03) 
  pd(37,03) = pd(37,03) - rrt(099) * density(37) 
  pd(37,37) = pd(37,37) - rrt(099) * density(03) 
  pd(02,05) = pd(02,05) + rrt(100) * density(11) 
  pd(02,11) = pd(02,11) + rrt(100) * density(05) 
  pd(11,05) = pd(11,05) - rrt(100) * density(11) 
  pd(11,11) = pd(11,11) - rrt(100) * density(05) 
  pd(11,05) = pd(11,05) + rrt(101) * density(12) 
  pd(11,12) = pd(11,12) + rrt(101) * density(05) 
  pd(12,05) = pd(12,05) - rrt(101) * density(12) 
  pd(12,12) = pd(12,12) - rrt(101) * density(05) 
  pd(12,05) = pd(12,05) + rrt(102) * density(13) 
  pd(12,13) = pd(12,13) + rrt(102) * density(05) 
  pd(13,05) = pd(13,05) - rrt(102) * density(13) 
  pd(13,13) = pd(13,13) - rrt(102) * density(05) 
  pd(13,05) = pd(13,05) + rrt(103) * density(14) 
  pd(13,14) = pd(13,14) + rrt(103) * density(05) 
  pd(14,05) = pd(14,05) - rrt(103) * density(14) 
  pd(14,14) = pd(14,14) - rrt(103) * density(05) 
  pd(14,05) = pd(14,05) + rrt(104) * density(15) 
  pd(14,15) = pd(14,15) + rrt(104) * density(05) 
  pd(15,05) = pd(15,05) - rrt(104) * density(15) 
  pd(15,15) = pd(15,15) - rrt(104) * density(05) 
  pd(15,05) = pd(15,05) + rrt(105) * density(16) 
  pd(15,16) = pd(15,16) + rrt(105) * density(05) 
  pd(16,05) = pd(16,05) - rrt(105) * density(16) 
  pd(16,16) = pd(16,16) - rrt(105) * density(05) 
  pd(16,05) = pd(16,05) + rrt(106) * density(17) 
  pd(16,17) = pd(16,17) + rrt(106) * density(05) 
  pd(17,05) = pd(17,05) - rrt(106) * density(17) 
  pd(17,17) = pd(17,17) - rrt(106) * density(05) 
  pd(17,05) = pd(17,05) + rrt(107) * density(18) 
  pd(17,18) = pd(17,18) + rrt(107) * density(05) 
  pd(18,05) = pd(18,05) - rrt(107) * density(18) 
  pd(18,18) = pd(18,18) - rrt(107) * density(05) 
  pd(03,04) = pd(03,04) + rrt(108) * density(35) 
  pd(03,35) = pd(03,35) + rrt(108) * density(04) 
  pd(35,04) = pd(35,04) - rrt(108) * density(35) 
  pd(35,35) = pd(35,35) - rrt(108) * density(04) 
  pd(35,04) = pd(35,04) + rrt(109) * density(36) 
  pd(35,36) = pd(35,36) + rrt(109) * density(04) 
  pd(36,04) = pd(36,04) - rrt(109) * density(36) 
  pd(36,36) = pd(36,36) - rrt(109) * density(04) 
  pd(36,04) = pd(36,04) + rrt(110) * density(37) 
  pd(36,37) = pd(36,37) + rrt(110) * density(04) 
  pd(37,04) = pd(37,04) - rrt(110) * density(37) 
  pd(37,37) = pd(37,37) - rrt(110) * density(04) 
  pd(02,38) = pd(02,38) + rrt(111) * density(61) 
  pd(02,61) = pd(02,61) + rrt(111) * density(38) 
  pd(38,38) = pd(38,38) - rrt(111) * density(61) 
  pd(38,61) = pd(38,61) - rrt(111) * density(38) 
  pd(05,05) = pd(05,05) - rrt(112) * density(61) 
  pd(05,61) = pd(05,61) - rrt(112) * density(05) 
  pd(61,05) = pd(61,05) - rrt(112) * density(61) 
  pd(61,61) = pd(61,61) - rrt(112) * density(05) 
  pd(63,05) = pd(63,05) + rrt(112) * density(61) 
  pd(63,61) = pd(63,61) + rrt(112) * density(05) 
  pd(02,39) = pd(02,39) + rrt(113) * density(61) 
  pd(02,61) = pd(02,61) + rrt(113) * density(39) 
  pd(39,39) = pd(39,39) - rrt(113) * density(61) 
  pd(39,61) = pd(39,61) - rrt(113) * density(39) 
  pd(04,04) = pd(04,04) - rrt(115) * density(61) 
  pd(04,61) = pd(04,61) - rrt(115) * density(04) 
  pd(61,04) = pd(61,04) - rrt(115) * density(61) 
  pd(61,61) = pd(61,61) - rrt(115) * density(04) 
  pd(62,04) = pd(62,04) + rrt(115) * density(61) 
  pd(62,61) = pd(62,61) + rrt(115) * density(04) 
  pd(06,06) = pd(06,06) - rrt(116) * density(61) 
  pd(06,61) = pd(06,61) - rrt(116) * density(06) 
  pd(61,06) = pd(61,06) - rrt(116) * density(61) 
  pd(61,61) = pd(61,61) - rrt(116) * density(06) 
  pd(64,06) = pd(64,06) + rrt(116) * density(61) 
  pd(64,61) = pd(64,61) + rrt(116) * density(06) 
  pd(07,07) = pd(07,07) - rrt(117) * density(61) 
  pd(07,61) = pd(07,61) - rrt(117) * density(07) 
  pd(61,07) = pd(61,07) - rrt(117) * density(61) 
  pd(61,61) = pd(61,61) - rrt(117) * density(07) 
  pd(65,07) = pd(65,07) + rrt(117) * density(61) 
  pd(65,61) = pd(65,61) + rrt(117) * density(07) 
  pd(04,04) = pd(04,04) - rrt(118) * density(63) 
  pd(04,63) = pd(04,63) - rrt(118) * density(04) 
  pd(63,04) = pd(63,04) - rrt(118) * density(63) 
  pd(63,63) = pd(63,63) - rrt(118) * density(04) 
  pd(64,04) = pd(64,04) + rrt(118) * density(63) 
  pd(64,63) = pd(64,63) + rrt(118) * density(04) 
  pd(03,05) = pd(03,05) + rrt(119) * density(63) 
  pd(03,63) = pd(03,63) + rrt(119) * density(05) 
  pd(05,05) = pd(05,05) - rrt(119) * density(63) 
  pd(05,63) = pd(05,63) - rrt(119) * density(05) 
  pd(61,05) = pd(61,05) + rrt(119) * density(63) 
  pd(61,63) = pd(61,63) + rrt(119) * density(05) 
  pd(63,05) = pd(63,05) - rrt(119) * density(63) 
  pd(63,63) = pd(63,63) - rrt(119) * density(05) 
  pd(02,02) = pd(02,02) - rrt(120) * density(63) 
  pd(02,63) = pd(02,63) - rrt(120) * density(02) 
  pd(09,02) = pd(09,02) + rrt(120) * density(63) 
  pd(09,63) = pd(09,63) + rrt(120) * density(02) 
  pd(61,02) = pd(61,02) + rrt(120) * density(63) 
  pd(61,63) = pd(61,63) + rrt(120) * density(02) 
  pd(63,02) = pd(63,02) - rrt(120) * density(63) 
  pd(63,63) = pd(63,63) - rrt(120) * density(02) 
  pd(06,06) = pd(06,06) - rrt(121) * density(63) 
  pd(06,63) = pd(06,63) - rrt(121) * density(06) 
  pd(63,06) = pd(63,06) - rrt(121) * density(63) 
  pd(63,63) = pd(63,63) - rrt(121) * density(06) 
  pd(65,06) = pd(65,06) + rrt(121) * density(63) 
  pd(65,63) = pd(65,63) + rrt(121) * density(06) 
  pd(05,05) = pd(05,05) - rrt(122) * density(62) 
  pd(05,62) = pd(05,62) - rrt(122) * density(05) 
  pd(62,05) = pd(62,05) - rrt(122) * density(62) 
  pd(62,62) = pd(62,62) - rrt(122) * density(05) 
  pd(64,05) = pd(64,05) + rrt(122) * density(62) 
  pd(64,62) = pd(64,62) + rrt(122) * density(05) 
  pd(05,05) = pd(05,05) - rrt(123) * density(64) 
  pd(05,64) = pd(05,64) - rrt(123) * density(05) 
  pd(64,05) = pd(64,05) - rrt(123) * density(64) 
  pd(64,64) = pd(64,64) - rrt(123) * density(05) 
  pd(65,05) = pd(65,05) + rrt(123) * density(64) 
  pd(65,64) = pd(65,64) + rrt(123) * density(05) 
  pd(05,05) = pd(05,05) - rrt(124) * density(65) 
  pd(05,65) = pd(05,65) - rrt(124) * density(05) 
  pd(08,05) = pd(08,05) + rrt(124) * density(65) 
  pd(08,65) = pd(08,65) + rrt(124) * density(05) 
  pd(61,05) = pd(61,05) + rrt(124) * density(65) 
  pd(61,65) = pd(61,65) + rrt(124) * density(05) 
  pd(65,05) = pd(65,05) - rrt(124) * density(65) 
  pd(65,65) = pd(65,65) - rrt(124) * density(05) 
  pd(07,07) = pd(07,07) - rrt(125) * density(63) 
  pd(07,63) = pd(07,63) - rrt(125) * density(07) 
  pd(08,07) = pd(08,07) + rrt(125) * density(63) 
  pd(08,63) = pd(08,63) + rrt(125) * density(07) 
  pd(61,07) = pd(61,07) + rrt(125) * density(63) 
  pd(61,63) = pd(61,63) + rrt(125) * density(07) 
  pd(63,07) = pd(63,07) - rrt(125) * density(63) 
  pd(63,63) = pd(63,63) - rrt(125) * density(07) 
  pd(61,62) = pd(61,62) + rrt(126) * density(63) 
  pd(61,63) = pd(61,63) + rrt(126) * density(62) 
  pd(62,62) = pd(62,62) - rrt(126) * density(63) 
  pd(62,63) = pd(62,63) - rrt(126) * density(62) 
  pd(63,62) = pd(63,62) - rrt(126) * density(63) 
  pd(63,63) = pd(63,63) - rrt(126) * density(62) 
  pd(64,62) = pd(64,62) + rrt(126) * density(63) 
  pd(64,63) = pd(64,63) + rrt(126) * density(62) 
  pd(61,63) = pd(61,63) + rrt(127) * density(64) 
  pd(61,64) = pd(61,64) + rrt(127) * density(63) 
  pd(63,63) = pd(63,63) - rrt(127) * density(64) 
  pd(63,64) = pd(63,64) - rrt(127) * density(63) 
  pd(64,63) = pd(64,63) - rrt(127) * density(64) 
  pd(64,64) = pd(64,64) - rrt(127) * density(63) 
  pd(65,63) = pd(65,63) + rrt(127) * density(64) 
  pd(65,64) = pd(65,64) + rrt(127) * density(63) 
  pd(08,63) = pd(08,63) + rrt(128) * density(65) 
  pd(08,65) = pd(08,65) + rrt(128) * density(63) 
  pd(61,63) = pd(61,63) + rrt(128) * density(65) * 2.0d0
  pd(61,65) = pd(61,65) + rrt(128) * density(63) * 2.0d0
  pd(63,63) = pd(63,63) - rrt(128) * density(65) 
  pd(63,65) = pd(63,65) - rrt(128) * density(63) 
  pd(65,63) = pd(65,63) - rrt(128) * density(65) 
  pd(65,65) = pd(65,65) - rrt(128) * density(63) 
  pd(02,02) = pd(02,02) - rrt(129) * density(61)**2 
  pd(02,61) = pd(02,61) - rrt(129) * density(02) * density(61) * 2.0d0
  pd(61,02) = pd(61,02) - rrt(129) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(129) * density(02) * density(61) * 4.0d0
  pd(62,02) = pd(62,02) + rrt(129) * density(61)**2 * 2.0d0
  pd(62,61) = pd(62,61) + rrt(129) * density(02) * density(61) * 4.0d0
  pd(38,38) = pd(38,38) - rrt(130) * density(61)**2 
  pd(38,61) = pd(38,61) - rrt(130) * density(38) * density(61) * 2.0d0
  pd(61,38) = pd(61,38) - rrt(130) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(130) * density(38) * density(61) * 4.0d0
  pd(62,38) = pd(62,38) + rrt(130) * density(61)**2 * 2.0d0
  pd(62,61) = pd(62,61) + rrt(130) * density(38) * density(61) * 4.0d0
  pd(03,03) = pd(03,03) - rrt(131) * density(61)**2 
  pd(03,61) = pd(03,61) - rrt(131) * density(03) * density(61) * 2.0d0
  pd(61,03) = pd(61,03) - rrt(131) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(131) * density(03) * density(61) * 4.0d0
  pd(63,03) = pd(63,03) + rrt(131) * density(61)**2 * 2.0d0
  pd(63,61) = pd(63,61) + rrt(131) * density(03) * density(61) * 4.0d0
  pd(35,35) = pd(35,35) - rrt(132) * density(61)**2 
  pd(35,61) = pd(35,61) - rrt(132) * density(35) * density(61) * 2.0d0
  pd(61,35) = pd(61,35) - rrt(132) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(132) * density(35) * density(61) * 4.0d0
  pd(63,35) = pd(63,35) + rrt(132) * density(61)**2 * 2.0d0
  pd(63,61) = pd(63,61) + rrt(132) * density(35) * density(61) * 4.0d0
  pd(36,36) = pd(36,36) - rrt(133) * density(61)**2 
  pd(36,61) = pd(36,61) - rrt(133) * density(36) * density(61) * 2.0d0
  pd(61,36) = pd(61,36) - rrt(133) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(133) * density(36) * density(61) * 4.0d0
  pd(63,36) = pd(63,36) + rrt(133) * density(61)**2 * 2.0d0
  pd(63,61) = pd(63,61) + rrt(133) * density(36) * density(61) * 4.0d0
  pd(37,37) = pd(37,37) - rrt(134) * density(61)**2 
  pd(37,61) = pd(37,61) - rrt(134) * density(37) * density(61) * 2.0d0
  pd(61,37) = pd(61,37) - rrt(134) * density(61)**2 * 2.0d0
  pd(61,61) = pd(61,61) - rrt(134) * density(37) * density(61) * 4.0d0
  pd(63,37) = pd(63,37) + rrt(134) * density(61)**2 * 2.0d0
  pd(63,61) = pd(63,61) + rrt(134) * density(37) * density(61) * 4.0d0
  pd(03,01) = pd(03,01) - rrt(135) * density(03) 
  pd(03,03) = pd(03,03) - rrt(135) * density(01) 
  pd(42,01) = pd(42,01) + rrt(135) * density(03) 
  pd(42,03) = pd(42,03) + rrt(135) * density(01) 
  pd(03,01) = pd(03,01) - rrt(136) * density(03) 
  pd(03,03) = pd(03,03) - rrt(136) * density(01) 
  pd(43,01) = pd(43,01) + rrt(136) * density(03) 
  pd(43,03) = pd(43,03) + rrt(136) * density(01) 
  pd(03,01) = pd(03,01) - rrt(137) * density(03) 
  pd(03,03) = pd(03,03) - rrt(137) * density(01) 
  pd(44,01) = pd(44,01) + rrt(137) * density(03) 
  pd(44,03) = pd(44,03) + rrt(137) * density(01) 
  pd(03,01) = pd(03,01) - rrt(138) * density(03) 
  pd(03,03) = pd(03,03) - rrt(138) * density(01) 
  pd(45,01) = pd(45,01) + rrt(138) * density(03) 
  pd(45,03) = pd(45,03) + rrt(138) * density(01) 
  pd(02,01) = pd(02,01) - rrt(139) * density(02) 
  pd(02,02) = pd(02,02) - rrt(139) * density(01) 
  pd(38,01) = pd(38,01) + rrt(139) * density(02) 
  pd(38,02) = pd(38,02) + rrt(139) * density(01) 
  pd(02,01) = pd(02,01) - rrt(140) * density(02) 
  pd(02,02) = pd(02,02) - rrt(140) * density(01) 
  pd(39,01) = pd(39,01) + rrt(140) * density(02) 
  pd(39,02) = pd(39,02) + rrt(140) * density(01) 
  pd(02,01) = pd(02,01) - rrt(141) * density(02) 
  pd(02,02) = pd(02,02) - rrt(141) * density(01) 
  pd(40,01) = pd(40,01) + rrt(141) * density(02) 
  pd(40,02) = pd(40,02) + rrt(141) * density(01) 
  pd(02,01) = pd(02,01) - rrt(142) * density(02) 
  pd(02,02) = pd(02,02) - rrt(142) * density(01) 
  pd(41,01) = pd(41,01) + rrt(142) * density(02) 
  pd(41,02) = pd(41,02) + rrt(142) * density(01) 
  pd(01,01) = pd(01,01) + rrt(143) * density(02) 
  pd(01,02) = pd(01,02) + rrt(143) * density(01) 
  pd(02,01) = pd(02,01) - rrt(143) * density(02) 
  pd(02,02) = pd(02,02) - rrt(143) * density(01) 
  pd(48,01) = pd(48,01) + rrt(143) * density(02) 
  pd(48,02) = pd(48,02) + rrt(143) * density(01) 
  pd(01,01) = pd(01,01) + rrt(144) * density(03) 
  pd(01,03) = pd(01,03) + rrt(144) * density(01) 
  pd(03,01) = pd(03,01) - rrt(144) * density(03) 
  pd(03,03) = pd(03,03) - rrt(144) * density(01) 
  pd(53,01) = pd(53,01) + rrt(144) * density(03) 
  pd(53,03) = pd(53,03) + rrt(144) * density(01) 
  pd(01,01) = pd(01,01) + rrt(145) * density(04) 
  pd(01,04) = pd(01,04) + rrt(145) * density(01) 
  pd(04,01) = pd(04,01) - rrt(145) * density(04) 
  pd(04,04) = pd(04,04) - rrt(145) * density(01) 
  pd(51,01) = pd(51,01) + rrt(145) * density(04) 
  pd(51,04) = pd(51,04) + rrt(145) * density(01) 
  pd(01,01) = pd(01,01) + rrt(146) * density(05) 
  pd(01,05) = pd(01,05) + rrt(146) * density(01) 
  pd(05,01) = pd(05,01) - rrt(146) * density(05) 
  pd(05,05) = pd(05,05) - rrt(146) * density(01) 
  pd(52,01) = pd(52,01) + rrt(146) * density(05) 
  pd(52,05) = pd(52,05) + rrt(146) * density(01) 
  pd(01,01) = pd(01,01) + rrt(147) * density(06) 
  pd(01,06) = pd(01,06) + rrt(147) * density(01) 
  pd(06,01) = pd(06,01) - rrt(147) * density(06) 
  pd(06,06) = pd(06,06) - rrt(147) * density(01) 
  pd(56,01) = pd(56,01) + rrt(147) * density(06) 
  pd(56,06) = pd(56,06) + rrt(147) * density(01) 
  pd(01,01) = pd(01,01) + rrt(148) * density(07) 
  pd(01,07) = pd(01,07) + rrt(148) * density(01) 
  pd(07,01) = pd(07,01) - rrt(148) * density(07) 
  pd(07,07) = pd(07,07) - rrt(148) * density(01) 
  pd(57,01) = pd(57,01) + rrt(148) * density(07) 
  pd(57,07) = pd(57,07) + rrt(148) * density(01) 
  pd(01,01) = pd(01,01) + rrt(149) * density(08) 
  pd(01,08) = pd(01,08) + rrt(149) * density(01) 
  pd(08,01) = pd(08,01) - rrt(149) * density(08) 
  pd(08,08) = pd(08,08) - rrt(149) * density(01) 
  pd(58,01) = pd(58,01) + rrt(149) * density(08) 
  pd(58,08) = pd(58,08) + rrt(149) * density(01) 
  pd(01,01) = pd(01,01) + rrt(150) * density(02) 
  pd(01,02) = pd(01,02) + rrt(150) * density(01) 
  pd(02,01) = pd(02,01) - rrt(150) * density(02) 
  pd(02,02) = pd(02,02) - rrt(150) * density(01) 
  pd(04,01) = pd(04,01) + rrt(150) * density(02) 
  pd(04,02) = pd(04,02) + rrt(150) * density(01) 
  pd(51,01) = pd(51,01) + rrt(150) * density(02) 
  pd(51,02) = pd(51,02) + rrt(150) * density(01) 
  pd(01,01) = pd(01,01) + rrt(151) * density(03) 
  pd(01,03) = pd(01,03) + rrt(151) * density(01) 
  pd(03,01) = pd(03,01) - rrt(151) * density(03) 
  pd(03,03) = pd(03,03) - rrt(151) * density(01) 
  pd(05,01) = pd(05,01) + rrt(151) * density(03) 
  pd(05,03) = pd(05,03) + rrt(151) * density(01) 
  pd(52,01) = pd(52,01) + rrt(151) * density(03) 
  pd(52,03) = pd(52,03) + rrt(151) * density(01) 
  pd(01,01) = pd(01,01) + rrt(152) * density(06) 
  pd(01,06) = pd(01,06) + rrt(152) * density(01) 
  pd(04,01) = pd(04,01) + rrt(152) * density(06) 
  pd(04,06) = pd(04,06) + rrt(152) * density(01) 
  pd(06,01) = pd(06,01) - rrt(152) * density(06) 
  pd(06,06) = pd(06,06) - rrt(152) * density(01) 
  pd(52,01) = pd(52,01) + rrt(152) * density(06) 
  pd(52,06) = pd(52,06) + rrt(152) * density(01) 
  pd(01,01) = pd(01,01) + rrt(153) * density(07) 
  pd(01,07) = pd(01,07) + rrt(153) * density(01) 
  pd(06,01) = pd(06,01) + rrt(153) * density(07) 
  pd(06,07) = pd(06,07) + rrt(153) * density(01) 
  pd(07,01) = pd(07,01) - rrt(153) * density(07) 
  pd(07,07) = pd(07,07) - rrt(153) * density(01) 
  pd(52,01) = pd(52,01) + rrt(153) * density(07) 
  pd(52,07) = pd(52,07) + rrt(153) * density(01) 
  pd(01,01) = pd(01,01) + rrt(154) * density(08) 
  pd(01,08) = pd(01,08) + rrt(154) * density(01) 
  pd(07,01) = pd(07,01) + rrt(154) * density(08) 
  pd(07,08) = pd(07,08) + rrt(154) * density(01) 
  pd(08,01) = pd(08,01) - rrt(154) * density(08) 
  pd(08,08) = pd(08,08) - rrt(154) * density(01) 
  pd(52,01) = pd(52,01) + rrt(154) * density(08) 
  pd(52,08) = pd(52,08) + rrt(154) * density(01) 
  pd(03,01) = pd(03,01) - rrt(155) * density(03) 
  pd(03,03) = pd(03,03) - rrt(155) * density(01) 
  pd(05,01) = pd(05,01) + rrt(155) * density(03) * 2.0d0
  pd(05,03) = pd(05,03) + rrt(155) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(156) * density(02) 
  pd(02,02) = pd(02,02) - rrt(156) * density(01) 
  pd(04,01) = pd(04,01) + rrt(156) * density(02) * 2.0d0
  pd(04,02) = pd(04,02) + rrt(156) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(157) * density(02) 
  pd(02,02) = pd(02,02) - rrt(157) * density(01) 
  pd(11,01) = pd(11,01) + rrt(157) * density(02) 
  pd(11,02) = pd(11,02) + rrt(157) * density(01) 
  pd(03,01) = pd(03,01) - rrt(158) * density(03) 
  pd(03,03) = pd(03,03) - rrt(158) * density(01) 
  pd(35,01) = pd(35,01) + rrt(158) * density(03) 
  pd(35,03) = pd(35,03) + rrt(158) * density(01) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(66,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(66,:) = pd(66,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  rrt(001) = 5.0D-11
  rrt(002) = 5.4D-11*EXP(-165.0D0/TGAS)
  rrt(003) = 5.0D-14*(TGAS/300.0D0)
  rrt(004) = 1.7D-12*(TGAS/300.0D0)**1.5D0
  rrt(005) = 8.5D-11
  rrt(006) = 6.6D-11*EXP(-1840.0D0/TGAS)
  rrt(007) = 1.2D-10
  rrt(008) = 1.2D-10
  rrt(009) = 1.66D-12
  rrt(010) = 4.0D-10*(TGAS/300.0D0)**0.5D0*EXP(-16600.0D0/TGAS)
  rrt(011) = 5.4D-11*EXP(-6492.0D0/TGAS)
  rrt(012) = 8.4D-14*(TGAS/300.0D0)**4.1D0*EXP(-4760.0D0/TGAS)
  rrt(013) = 5.0D-11
  rrt(014) = 2.0D-10*EXP(-3500.0D0/TGAS)
  rrt(015) = 1.6D-10
  rrt(016) = 2.5D-11
  rrt(017) = 1.5D-11
  rrt(018) = 2.6D-11
  rrt(019) = 4.0D-10*(TGAS/300.0D0)**0.5D0
  rrt(020) = 2.3D-12
  rrt(021) = 1.1D-10
  rrt(022) = 2.5D-14
  rrt(023) = 1.38D-33*EXP(502.978D0/TGAS)
  rrt(024) = 1.0D-32
  rrt(025) = 1.4D-32
  rrt(026) = 1.7D-33
  rrt(027) = 2.4D-33
  rrt(028) = 8.3D-33*EXP(500.0D0/TGAS)
  rrt(029) = 8.3D-33*(300.0D0/TGAS)
  rrt(030) = 1.0D-33
  rrt(031) = 1.0D-34
  rrt(032) = 1.0D-32
  rrt(033) = 5.5D-30
  rrt(034) = 2.5D-35*(TGAS/300.0D0)*EXP(1700.0D0/TGAS)
  rrt(035) = 1.7D-33
  rrt(036) = 1.0D-32
  rrt(037) = 2.4D-33
  rrt(038) = 1.4D-32
  rrt(039) = 8.8D-33*(300.0D0/TGAS)**0.6D0
  rrt(040) = 2.7D-11*EXP(-67400.0D0/TGAS)
  rrt(041) = 5.0D-13
  rrt(042) = 1.0D-12
  rrt(043) = 0.5D0
  rrt(044) = 1.34D5
  rrt(045) = 1.0D2
  rrt(046) = 2.45D7
  rrt(047) = 5.0D-10
  rrt(048) = 0.20D0*2.35D-9
  rrt(049) = 0.71D0*2.35D-9
  rrt(050) = 0.09D0*2.35D-9
  rrt(051) = 7.2D-13*(TION/300.0D0)
  rrt(052) = 2.00D-9
  rrt(053) = 3.0D-10
  rrt(054) = 1.95D-9
  rrt(055) = 6.6D-11
  rrt(056) = 1.0D-11
  rrt(057) = 2.1D-16*EXP(TION/121.0D0)
  rrt(058) = 5.20D-9
  rrt(059) = 6.4D-10
  rrt(060) = 2.0D-9
  rrt(061) = 2.00D-9
  rrt(062) = 5.70D-9
  rrt(063) = 0.15D0*1.23D-9
  rrt(064) = 0.85D0*1.23D-9
  rrt(065) = 0.75D0*2.40D-9
  rrt(066) = 0.25D0*2.40D-9
  rrt(067) = 6.50D-10
  rrt(068) = 1.95D-10
  rrt(069) = 0.5D0*2.30D-9
  rrt(070) = rrt(69)
  rrt(071) = 2.10D-9
  rrt(072) = 2.3D-9
  rrt(073) = 9.0D-30*EXP(400.0D0/TION)
  rrt(074) = 1.7D-29*(300.0D0/TION)**2.1D0
  rrt(075) = 5.2D-29*(300.0D0/TION)**2.2D0
  rrt(076) = 1.0D-29
  rrt(077) = 2.0D-7*(300.0D0/TGAS)
  rrt(078) = rrt(77)
  rrt(079) = rrt(77)
  rrt(080) = rrt(77)
  rrt(081) = rrt(77)
  rrt(082) = rrt(77)
  rrt(083) = 2.0D-25*(300.0D0/TGAS)**2.5D0
  rrt(084) = rrt(83)
  rrt(085) = rrt(83)
  rrt(086) = rrt(83)
  rrt(087) = rrt(83)
  rrt(088) = rrt(83)
  rrt(089) = 1.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(090) = 2.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(091) = 3.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(092) = 4.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(093) = 5.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(094) = 6.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(095) = 7.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(096) = 8.00D-17*(TGAS/300.0D0)**0.5D0
  rrt(097) = 1.00D-16*(TGAS/300.0D0)**0.5D0
  rrt(098) = 2.00D-16*(TGAS/300.0D0)**0.5D0
  rrt(099) = 3.00D-16*(TGAS/300.0D0)**0.5D0
  rrt(100) = 5.00D-18*(TGAS/300.0D0)**0.5D0
  rrt(101) = rrt(89)
  rrt(102) = 1.50D-17*(TGAS/300.0D0)**0.5D0
  rrt(103) = rrt(90)
  rrt(104) = 2.50D-17*(TGAS/300.0D0)**0.5D0
  rrt(105) = rrt(91)
  rrt(106) = 3.50D-17*(TGAS/300.0D0)**0.5D0
  rrt(107) = rrt(92)
  rrt(108) = rrt(90)
  rrt(109) = rrt(92)
  rrt(110) = rrt(94)
  rrt(111) = 1.0D-3
  rrt(112) = 4.5D-4
  rrt(113) = 1.0D-3
  rrt(114) = 1.0D-4
  rrt(115) = 1.0D0
  rrt(116) = 1.0D0
  rrt(117) = 1.0D0
  rrt(118) = 6.0D-3
  rrt(119) = 1.0D-3
  rrt(120) = 6.0D-3
  rrt(121) = 1.5D-2
  rrt(122) = 1.0D-2
  rrt(123) = 1.0D-3
  rrt(124) = 8.0D-3
  rrt(125) = 1.0D-2
  rrt(126) = 8.0D-4
  rrt(127) = 8.0D-4
  rrt(128) = 8.0D-4
  rrt(129) = 1.0D-2
  rrt(130) = 1.0D-2
  rrt(131) = 1.0D-1
  rrt(132) = 1.0D-2
  rrt(133) = 5.0D-2
  rrt(134) = 1.0D-1
  rrt(135) = bolsig_rates(bolsig_pointer(1))
  rrt(136) = bolsig_rates(bolsig_pointer(2))
  rrt(137) = bolsig_rates(bolsig_pointer(3))
  rrt(138) = bolsig_rates(bolsig_pointer(4))
  rrt(139) = bolsig_rates(bolsig_pointer(5))
  rrt(140) = bolsig_rates(bolsig_pointer(6))
  rrt(141) = bolsig_rates(bolsig_pointer(7))
  rrt(142) = bolsig_rates(bolsig_pointer(8))
  rrt(143) = bolsig_rates(bolsig_pointer(9))
  rrt(144) = bolsig_rates(bolsig_pointer(10))
  rrt(145) = bolsig_rates(bolsig_pointer(11))
  rrt(146) = bolsig_rates(bolsig_pointer(12))
  rrt(147) = bolsig_rates(bolsig_pointer(13))
  rrt(148) = bolsig_rates(bolsig_pointer(14))
  rrt(149) = bolsig_rates(bolsig_pointer(15))
  rrt(150) = bolsig_rates(bolsig_pointer(16))
  rrt(151) = bolsig_rates(bolsig_pointer(17))
  rrt(152) = bolsig_rates(bolsig_pointer(18))
  rrt(153) = bolsig_rates(bolsig_pointer(19))
  rrt(154) = bolsig_rates(bolsig_pointer(20))
  rrt(155) = bolsig_rates(bolsig_pointer(21))
  rrt(156) = bolsig_rates(bolsig_pointer(22))
  rrt(157) = bolsig_rates(bolsig_pointer(23))
  rrt(158) = bolsig_rates(bolsig_pointer(24))
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------

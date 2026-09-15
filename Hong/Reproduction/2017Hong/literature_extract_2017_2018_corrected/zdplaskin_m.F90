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
! Tue Aug 25 12:29:14 2026
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS:    E    N    H SURF 
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
  integer, parameter :: species_max = 48, species_electrons = 43, species_length = 15, reactions_max = 504, reactions_length = 32
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
  integer, parameter, private               :: bolsig_species_max = 18, bolsig_species_length = 15, bolsig_rates_max = 41, &
                                               bolsig_addsect_max = 8 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0, &
                                               bolsig_addsect(2,bolsig_addsect_max) 
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
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1, 0, 0, 0, 1, 1, 1, 1, 1,&
   -1, 0, 0, 0, 0, 0/
  data species_name(1:species_max) &
  /"N2             ","N2(V1)         ","N2(V2)         ","N2(V3)         ","N2(V4)         ","N2(V5)         ","N2(V6)         ",&
   "N2(V7)         ","N2(V8)         ","N2(A3)         ","N2(B3)         ","N2(A`1)        ","N2(C3)         ","N              ",&
   "N(2D)          ","N(2P)          ","N^+            ","N2^+           ","N3^+           ","N4^+           ","H2             ",&
   "H2(B3SIG)      ","H2(B1SIG)      ","H2(C3PI)       ","H2(A3SIG)      ","H2(RYDBERG_SUM)","H2(V1)         ","H2(V2)         ",&
   "H2(V3)         ","H              ","H2^+           ","H3^+           ","H^+            ","H^-            ","NH             ",&
   "NH2            ","NH3            ","NH^+           ","N2H^+          ","NH2^+          ","NH3^+          ","NH4^+          ",&
   "E              ","SURF           ","HSURF          ","NSURF          ","NHSURF         ","NH2SURF        "/
  data reaction_sign(1:54) &
  /"bolsig:N2->N2(A3,V0-4)          ","bolsig:N2->N2(A3,V5-9)          ","bolsig:N2->N2(A3,V10-)          ",&
   "bolsig:N2->N2(B3)               ","bolsig:N2->N2(A`1)              ","bolsig:N2->N2(C3)               ",&
   "bolsig:H2->H2(B3SIG)            ","bolsig:H2->H2(B1SIG)            ","bolsig:H2->H2(C3PI)             ",&
   "bolsig:H2->H2(A3SIG)            ","bolsig:H2->H2(RYDBERG_SUM)      ","bolsig:N2->N2^+                 ",&
   "bolsig:H2->H2^+                 ","E+N=>E+E+N^+                    ","E+H=>H^++E+E                    ",&
   "E+NH3=>E+E+NH3^+                ","bolsig:N2->N2(V1)               ","bolsig:N2->N2(V2)               ",&
   "bolsig:N2->N2(V3)               ","bolsig:N2->N2(V4)               ","bolsig:N2->N2(V5)               ",&
   "bolsig:N2->N2(V6)               ","bolsig:N2->N2(V7)               ","bolsig:N2->N2(V8)               ",&
   "bolsig:N2(V1)->N2               ","bolsig:N2(V2)->N2               ","bolsig:N2(V3)->N2               ",&
   "bolsig:N2(V4)->N2               ","bolsig:N2(V5)->N2               ","bolsig:N2(V6)->N2               ",&
   "bolsig:N2(V7)->N2               ","bolsig:N2(V8)->N2               ","bolsig:H2->H2(V1)               ",&
   "bolsig:H2->H2(V2)               ","bolsig:H2->H2(V3)               ","bolsig:H2(V1)->H2               ",&
   "bolsig:H2(V2)->H2               ","bolsig:H2(V3)->H2               ","bolsig:H2(B3SIG)->H2            ",&
   "bolsig:H2(B1SIG)->H2            ","bolsig:H2(C3PI)->H2             ","bolsig:H2(A3SIG)->H2            ",&
   "N2(V1)+N2=>N2+N2                ","N2+N2=>N2(V1)+N2                ","N2(V2)+N2=>N2(V1)+N2            ",&
   "N2(V1)+N2=>N2(V2)+N2            ","N2(V3)+N2=>N2(V2)+N2            ","N2(V2)+N2=>N2(V3)+N2            ",&
   "N2(V4)+N2=>N2(V3)+N2            ","N2(V3)+N2=>N2(V4)+N2            ","N2(V5)+N2=>N2(V4)+N2            ",&
   "N2(V4)+N2=>N2(V5)+N2            ","N2(V6)+N2=>N2(V5)+N2            ","N2(V5)+N2=>N2(V6)+N2            "/
  data reaction_sign(55:108) &
  /"N2(V7)+N2=>N2(V6)+N2            ","N2(V6)+N2=>N2(V7)+N2            ","N2(V8)+N2=>N2(V7)+N2            ",&
   "N2(V7)+N2=>N2(V8)+N2            ","N2(V1)+H2=>N2+H2                ","N2+H2=>N2(V1)+H2                ",&
   "N2(V2)+H2=>N2(V1)+H2            ","N2(V1)+H2=>N2(V2)+H2            ","N2(V3)+H2=>N2(V2)+H2            ",&
   "N2(V2)+H2=>N2(V3)+H2            ","N2(V4)+H2=>N2(V3)+H2            ","N2(V3)+H2=>N2(V4)+H2            ",&
   "N2(V5)+H2=>N2(V4)+H2            ","N2(V4)+H2=>N2(V5)+H2            ","N2(V6)+H2=>N2(V5)+H2            ",&
   "N2(V5)+H2=>N2(V6)+H2            ","N2(V7)+H2=>N2(V6)+H2            ","N2(V6)+H2=>N2(V7)+H2            ",&
   "N2(V8)+H2=>N2(V7)+H2            ","N2(V7)+H2=>N2(V8)+H2            ","N2(V1)+N=>N2+N                  ",&
   "N2+N=>N2(V1)+N                  ","N2(V2)+N=>N2(V1)+N              ","N2(V1)+N=>N2(V2)+N              ",&
   "N2(V3)+N=>N2(V2)+N              ","N2(V2)+N=>N2(V3)+N              ","N2(V4)+N=>N2(V3)+N              ",&
   "N2(V3)+N=>N2(V4)+N              ","N2(V5)+N=>N2(V4)+N              ","N2(V4)+N=>N2(V5)+N              ",&
   "N2(V6)+N=>N2(V5)+N              ","N2(V5)+N=>N2(V6)+N              ","N2(V7)+N=>N2(V6)+N              ",&
   "N2(V6)+N=>N2(V7)+N              ","N2(V8)+N=>N2(V7)+N              ","N2(V7)+N=>N2(V8)+N              ",&
   "N2(V1)+H=>N2+H                  ","N2+H=>N2(V1)+H                  ","N2(V2)+H=>N2(V1)+H              ",&
   "N2(V1)+H=>N2(V2)+H              ","N2(V3)+H=>N2(V2)+H              ","N2(V2)+H=>N2(V3)+H              ",&
   "N2(V4)+H=>N2(V3)+H              ","N2(V3)+H=>N2(V4)+H              ","N2(V5)+H=>N2(V4)+H              ",&
   "N2(V4)+H=>N2(V5)+H              ","N2(V6)+H=>N2(V5)+H              ","N2(V5)+H=>N2(V6)+H              ",&
   "N2(V7)+H=>N2(V6)+H              ","N2(V6)+H=>N2(V7)+H              ","N2(V8)+H=>N2(V7)+H              ",&
   "N2(V7)+H=>N2(V8)+H              ","H2(V1)+H2=>H2+H2                ","H2+H2=>H2(V1)+H2                "/
  data reaction_sign(109:162) &
  /"H2(V2)+H2=>H2(V1)+H2            ","H2(V1)+H2=>H2(V2)+H2            ","H2(V3)+H2=>H2(V2)+H2            ",&
   "H2(V2)+H2=>H2(V3)+H2            ","H2(V1)+H=>H2+H                  ","H2+H=>H2(V1)+H                  ",&
   "H2(V2)+H=>H2(V1)+H              ","H2(V1)+H=>H2(V2)+H              ","H2(V3)+H=>H2(V2)+H              ",&
   "H2(V2)+H=>H2(V3)+H              ","N2(V1)+N2(V1)=>N2+N2(V2)        ","N2+N2(V2)=>N2(V1)+N2(V1)        ",&
   "N2(V1)+N2(V2)=>N2+N2(V3)        ","N2+N2(V3)=>N2(V1)+N2(V2)        ","N2(V1)+N2(V3)=>N2+N2(V4)        ",&
   "N2+N2(V4)=>N2(V1)+N2(V3)        ","N2(V1)+N2(V4)=>N2+N2(V5)        ","N2+N2(V5)=>N2(V1)+N2(V4)        ",&
   "N2(V1)+N2(V5)=>N2+N2(V6)        ","N2+N2(V6)=>N2(V1)+N2(V5)        ","N2(V1)+N2(V6)=>N2+N2(V7)        ",&
   "N2+N2(V7)=>N2(V1)+N2(V6)        ","N2(V1)+N2(V7)=>N2+N2(V8)        ","N2+N2(V8)=>N2(V1)+N2(V7)        ",&
   "N2(V2)+N2(V2)=>N2(V1)+N2(V3)    ","N2(V1)+N2(V3)=>N2(V2)+N2(V2)    ","N2(V2)+N2(V3)=>N2(V1)+N2(V4)    ",&
   "N2(V1)+N2(V4)=>N2(V2)+N2(V3)    ","N2(V2)+N2(V4)=>N2(V1)+N2(V5)    ","N2(V1)+N2(V5)=>N2(V2)+N2(V4)    ",&
   "N2(V2)+N2(V5)=>N2(V1)+N2(V6)    ","N2(V1)+N2(V6)=>N2(V2)+N2(V5)    ","N2(V2)+N2(V6)=>N2(V1)+N2(V7)    ",&
   "N2(V1)+N2(V7)=>N2(V2)+N2(V6)    ","N2(V2)+N2(V7)=>N2(V1)+N2(V8)    ","N2(V1)+N2(V8)=>N2(V2)+N2(V7)    ",&
   "N2(V3)+N2(V3)=>N2(V2)+N2(V4)    ","N2(V2)+N2(V4)=>N2(V3)+N2(V3)    ","N2(V3)+N2(V4)=>N2(V2)+N2(V5)    ",&
   "N2(V2)+N2(V5)=>N2(V3)+N2(V4)    ","N2(V3)+N2(V5)=>N2(V2)+N2(V6)    ","N2(V2)+N2(V6)=>N2(V3)+N2(V5)    ",&
   "N2(V3)+N2(V6)=>N2(V2)+N2(V7)    ","N2(V2)+N2(V7)=>N2(V3)+N2(V6)    ","N2(V3)+N2(V7)=>N2(V2)+N2(V8)    ",&
   "N2(V2)+N2(V8)=>N2(V3)+N2(V7)    ","N2(V4)+N2(V4)=>N2(V3)+N2(V5)    ","N2(V3)+N2(V5)=>N2(V4)+N2(V4)    ",&
   "N2(V4)+N2(V5)=>N2(V3)+N2(V6)    ","N2(V3)+N2(V6)=>N2(V4)+N2(V5)    ","N2(V4)+N2(V6)=>N2(V3)+N2(V7)    ",&
   "N2(V3)+N2(V7)=>N2(V4)+N2(V6)    ","N2(V4)+N2(V7)=>N2(V3)+N2(V8)    ","N2(V3)+N2(V8)=>N2(V4)+N2(V7)    "/
  data reaction_sign(163:216) &
  /"N2(V5)+N2(V5)=>N2(V4)+N2(V6)    ","N2(V4)+N2(V6)=>N2(V5)+N2(V5)    ","N2(V5)+N2(V6)=>N2(V4)+N2(V7)    ",&
   "N2(V4)+N2(V7)=>N2(V5)+N2(V6)    ","N2(V5)+N2(V7)=>N2(V4)+N2(V8)    ","N2(V4)+N2(V8)=>N2(V5)+N2(V7)    ",&
   "N2(V6)+N2(V6)=>N2(V5)+N2(V7)    ","N2(V5)+N2(V7)=>N2(V6)+N2(V6)    ","N2(V6)+N2(V7)=>N2(V5)+N2(V8)    ",&
   "N2(V5)+N2(V8)=>N2(V6)+N2(V7)    ","N2(V7)+N2(V7)=>N2(V6)+N2(V8)    ","N2(V6)+N2(V8)=>N2(V7)+N2(V7)    ",&
   "H2(V1)+N2=>H2+N2(V1)            ","H2+N2(V1)=>H2(V1)+N2            ","H2(V1)+N2(V1)=>H2+N2(V2)        ",&
   "H2+N2(V2)=>H2(V1)+N2(V1)        ","H2(V1)+N2(V2)=>H2+N2(V3)        ","H2+N2(V3)=>H2(V1)+N2(V2)        ",&
   "H2(V1)+N2(V3)=>H2+N2(V4)        ","H2+N2(V4)=>H2(V1)+N2(V3)        ","H2(V1)+N2(V4)=>H2+N2(V5)        ",&
   "H2+N2(V5)=>H2(V1)+N2(V4)        ","H2(V1)+N2(V5)=>H2+N2(V6)        ","H2+N2(V6)=>H2(V1)+N2(V5)        ",&
   "H2(V1)+N2(V6)=>H2+N2(V7)        ","H2+N2(V7)=>H2(V1)+N2(V6)        ","H2(V1)+N2(V7)=>H2+N2(V8)        ",&
   "H2+N2(V8)=>H2(V1)+N2(V7)        ","H2(V2)+N2=>H2(V1)+N2(V1)        ","H2(V1)+N2(V1)=>H2(V2)+N2        ",&
   "H2(V2)+N2(V1)=>H2(V1)+N2(V2)    ","H2(V1)+N2(V2)=>H2(V2)+N2(V1)    ","H2(V2)+N2(V2)=>H2(V1)+N2(V3)    ",&
   "H2(V1)+N2(V3)=>H2(V2)+N2(V2)    ","H2(V2)+N2(V3)=>H2(V1)+N2(V4)    ","H2(V1)+N2(V4)=>H2(V2)+N2(V3)    ",&
   "H2(V2)+N2(V4)=>H2(V1)+N2(V5)    ","H2(V1)+N2(V5)=>H2(V2)+N2(V4)    ","H2(V2)+N2(V5)=>H2(V1)+N2(V6)    ",&
   "H2(V1)+N2(V6)=>H2(V2)+N2(V5)    ","H2(V2)+N2(V6)=>H2(V1)+N2(V7)    ","H2(V1)+N2(V7)=>H2(V2)+N2(V6)    ",&
   "H2(V2)+N2(V7)=>H2(V1)+N2(V8)    ","H2(V1)+N2(V8)=>H2(V2)+N2(V7)    ","H2(V3)+N2=>H2(V2)+N2(V1)        ",&
   "H2(V2)+N2(V1)=>H2(V3)+N2        ","H2(V3)+N2(V1)=>H2(V2)+N2(V2)    ","H2(V2)+N2(V2)=>H2(V3)+N2(V1)    ",&
   "H2(V3)+N2(V2)=>H2(V2)+N2(V3)    ","H2(V2)+N2(V3)=>H2(V3)+N2(V2)    ","H2(V3)+N2(V3)=>H2(V2)+N2(V4)    ",&
   "H2(V2)+N2(V4)=>H2(V3)+N2(V3)    ","H2(V3)+N2(V4)=>H2(V2)+N2(V5)    ","H2(V2)+N2(V5)=>H2(V3)+N2(V4)    "/
  data reaction_sign(217:270) &
  /"H2(V3)+N2(V5)=>H2(V2)+N2(V6)    ","H2(V2)+N2(V6)=>H2(V3)+N2(V5)    ","H2(V3)+N2(V6)=>H2(V2)+N2(V7)    ",&
   "H2(V2)+N2(V7)=>H2(V3)+N2(V6)    ","H2(V3)+N2(V7)=>H2(V2)+N2(V8)    ","H2(V2)+N2(V8)=>H2(V3)+N2(V7)    ",&
   "N2(V2)+H2=>N2+H2(V1)            ","N2+H2(V1)=>N2(V2)+H2            ","N2(V3)+H2=>N2(V1)+H2(V1)        ",&
   "N2(V1)+H2(V1)=>N2(V3)+H2        ","N2(V4)+H2=>N2(V2)+H2(V1)        ","N2(V2)+H2(V1)=>N2(V4)+H2        ",&
   "N2(V5)+H2=>N2(V3)+H2(V1)        ","N2(V3)+H2(V1)=>N2(V5)+H2        ","N2(V6)+H2=>N2(V4)+H2(V1)        ",&
   "N2(V4)+H2(V1)=>N2(V6)+H2        ","N2(V7)+H2=>N2(V5)+H2(V1)        ","N2(V5)+H2(V1)=>N2(V7)+H2        ",&
   "N2(V8)+H2=>N2(V6)+H2(V1)        ","N2(V6)+H2(V1)=>N2(V8)+H2        ","N2(V2)+H2(V1)=>N2+H2(V2)        ",&
   "N2+H2(V2)=>N2(V2)+H2(V1)        ","N2(V3)+H2(V1)=>N2(V1)+H2(V2)    ","N2(V1)+H2(V2)=>N2(V3)+H2(V1)    ",&
   "N2(V4)+H2(V1)=>N2(V2)+H2(V2)    ","N2(V2)+H2(V2)=>N2(V4)+H2(V1)    ","N2(V5)+H2(V1)=>N2(V3)+H2(V2)    ",&
   "N2(V3)+H2(V2)=>N2(V5)+H2(V1)    ","N2(V6)+H2(V1)=>N2(V4)+H2(V2)    ","N2(V4)+H2(V2)=>N2(V6)+H2(V1)    ",&
   "N2(V7)+H2(V1)=>N2(V5)+H2(V2)    ","N2(V5)+H2(V2)=>N2(V7)+H2(V1)    ","N2(V8)+H2(V1)=>N2(V6)+H2(V2)    ",&
   "N2(V6)+H2(V2)=>N2(V8)+H2(V1)    ","N2(V2)+H2(V2)=>N2+H2(V3)        ","N2+H2(V3)=>N2(V2)+H2(V2)        ",&
   "N2(V3)+H2(V2)=>N2(V1)+H2(V3)    ","N2(V1)+H2(V3)=>N2(V3)+H2(V2)    ","N2(V4)+H2(V2)=>N2(V2)+H2(V3)    ",&
   "N2(V2)+H2(V3)=>N2(V4)+H2(V2)    ","N2(V5)+H2(V2)=>N2(V3)+H2(V3)    ","N2(V3)+H2(V3)=>N2(V5)+H2(V2)    ",&
   "N2(V6)+H2(V2)=>N2(V4)+H2(V3)    ","N2(V4)+H2(V3)=>N2(V6)+H2(V2)    ","N2(V7)+H2(V2)=>N2(V5)+H2(V3)    ",&
   "N2(V5)+H2(V3)=>N2(V7)+H2(V2)    ","N2(V8)+H2(V2)=>N2(V6)+H2(V3)    ","N2(V6)+H2(V3)=>N2(V8)+H2(V2)    ",&
   "H2(V2)+H2(V2)=>H2(V3)+H2(V1)    ","H2(V3)+H2(V1)=>H2(V2)+H2(V2)    ","H2(V2)+H2(V1)=>H2(V3)+H2        ",&
   "H2(V3)+H2=>H2(V2)+H2(V1)        ","H2(V1)+H2(V1)=>H2(V2)+H2        ","H2(V2)+H2=>H2(V1)+H2(V1)        "/
  data reaction_sign(271:324) &
  /"N2(A3)+N2(V6)=>N2(B3)+N2        ","N2(A3)+N2(V7)=>N2(B3)+N2(V1)    ","N2(A3)+N2(V8)=>N2(B3)+N2(V2)    ",&
   "N2(B3)+N2=>N2(A3)+N2(V6)        ","N2(B3)+N2(V1)=>N2(A3)+N2(V7)    ","N2(B3)+N2(V2)=>N2(A3)+N2(V8)    ",&
   "N2(B3)+N2(V3)=>N2(A3)+N2(V8)    ","N2(B3)+N2(V4)=>N2(A3)+N2(V8)    ","N2(B3)+N2(V5)=>N2(A3)+N2(V8)    ",&
   "N2(B3)+N2(V6)=>N2(A3)+N2(V8)    ","N2(B3)+N2(V7)=>N2(A3)+N2(V8)    ","N2(B3)+N2(V8)=>N2(A3)+N2(V8)    ",&
   "N2(B3)=>N2(A3)                  ","N2(A`1)=>N2                     ","N2(C3)=>N2(B3)                  ",&
   "E+H2=>2H+E                      ","bolsig:N2->N2(SUM)              ","E+NH=>E+N+H                     ",&
   "E+NH2=>E+N+H2                   ","E+NH2=>E+NH+H                   ","E+NH3=>E+NH2+H                  ",&
   "E+NH3=>E+NH+H2                  ","E+N2^+=>N+N                     ","E+N2^+=>N+N(2D)                 ",&
   "E+N2^+=>N+N(2P)                 ","E+N3^+=>N2+N                    ","E+N4^+=>N2+N2                   ",&
   "E+H2=>E+H+H^++E                 ","E+H2^+=>2H                      ","E+H3^+=>3H                      ",&
   "E+H3^+=>H2+H                    ","E+NH^+=>N+H                     ","E+NH2^+=>NH+H                   ",&
   "E+NH2^+=>N+2H                   ","E+NH3^+=>NH+2H                  ","E+NH3^+=>NH2+H                  ",&
   "E+NH4^+=>NH3+H                  ","E+NH4^+=>NH2+2H                 ","E+N2H^+=>N2+H                   ",&
   "N2^++H2=>N2H^++H                ","N2^++N2(A3)=>N3^++N             ","N2^++N=>N^++N2                  ",&
   "N2^++N+N2=>N3^++N2              ","N2^++N2+N2=>N4^++N2             ","N2^++NH3=>NH3^++N2              ",&
   "N3^++N=>N2^++N2                 ","N4^++N2=>N2^++N2+N2             ","N4^++N=>N^++N2+N2               ",&
   "N^++H2=>NH^++H                  ","N^++NH3=>NH2^++NH               ","N^++NH3=>NH3^++N                ",&
   "N^++NH3=>N2H^++H2               ","H2^++H=>H2+H^+                  ","H2^++H2=>H3^++H                 "/
  data reaction_sign(325:378) &
  /"H2^++NH3=>NH3^++H2              ","H2^++N2=>N2H^++H                ","H^++NH3=>NH3^++H                ",&
   "NH^++H2=>H3^++N                 ","NH^++H2=>NH2^++H                ","NH^++NH3=>NH3^++NH              ",&
   "NH^++NH3=>NH4^++N               ","NH^++N2=>N2H^++N                ","NH2^++H2=>NH3^++H               ",&
   "NH2^++NH3=>NH3^++NH2            ","NH3+NH2^+=>NH+NH4^+             ","NH3^++NH3=>NH4^++NH2            ",&
   "N2(A3)+N=>N2+N                  ","N2(A3)+N=>N2+N(2P)              ","N2(A3)+N2=>N2+N2                ",&
   "N2(A3)+N2(A3)=>N2+N2(B3)        ","N2(A3)+N2(A3)=>N2+N2(C3)        ","N2(A3)+N2(A3)=>N2+2N            ",&
   "N2(B3)+N2=>N2(A3)+N2            ","N2(B3)+N2=>N2+N2                ","N2(A`1)+N2=>N2(B3)+N2           ",&
   "N2(A`1)+N2(A`1)=>N2^++N2+E      ","N2(A`1)+N2(A`1)=>N4^++E         ","N2(A`1)+N2(A3)=>N4^++E          ",&
   "N(2D)+N2=>N+N2                  ","N(2P)+N=>N+N                    ","N(2P)+N=>N(2D)+N                ",&
   "N(2P)+N2=>N+N2                  ","N2(A3)+H=>N2+H                  ","N2(A3)+H2=>N2+2H                ",&
   "N2(A3)+NH3=>N2+NH3              ","N2(B3)+H2=>N2(A3)+H2            ","N2(A`1)+H=>N2+H                 ",&
   "N2(A`1)+H2=>N2+2H               ","N+H2(V1)=>H+NH                  ","N+H2(V2)=>H+NH                  ",&
   "N+H2(V3)=>H+NH                  ","N+H2(RYDBERG_SUM)=>H+NH         ","N+H2(B3SIG)=>H+NH               ",&
   "N+H2(B1SIG)=>H+NH               ","N+H2(C3PI)=>H+NH                ","N+H2(A3SIG)=>H+NH               ",&
   "N(2D)+H2=>H+NH                  ","N(2D)+NH3=>NH+NH2               ","N(2P)+H2=>H+NH                  ",&
   "N+NH=>H+N2                      ","H+NH=>N+H2                      ","NH+NH=>H2+N2                    ",&
   "NH+NH=>N+NH2                    ","NH+NH=>N2+2H                    ","H+NH2=>H2+NH                    ",&
   "N+NH2=>N2+2H                    ","N+NH2=>N2+H2                    ","NH+NH2=>NH3+N                   "/
  data reaction_sign(379:432) &
  /"H2+NH2=>NH3+H                   ","H+NH3=>NH2+H2                   ","N+N+N2=>N2(A3)+N2               ",&
   "N+N+H2=>N2(A3)+H2               ","N+N+N2(V1)=>N2(A3)+N2(V1)       ","N+N+N2(V2)=>N2(A3)+N2(V2)       ",&
   "N+N+N2(V3)=>N2(A3)+N2(V3)       ","N+N+N2(V4)=>N2(A3)+N2(V4)       ","N+N+N2(V5)=>N2(A3)+N2(V5)       ",&
   "N+N+N2(V6)=>N2(A3)+N2(V6)       ","N+N+N2(V7)=>N2(A3)+N2(V7)       ","N+N+N2(V8)=>N2(A3)+N2(V8)       ",&
   "N+N+H2(V1)=>N2(A3)+H2(V1)       ","N+N+H2(V2)=>N2(A3)+H2(V2)       ","N+N+H2(V3)=>N2(A3)+H2(V3)       ",&
   "N+N+N=>N2(A3)+N                 ","N+N+H=>N2(A3)+H                 ","N+N+N2=>N2(B3)+N2               ",&
   "N+N+H2=>N2(B3)+H2               ","N+N+N2(V1)=>N2(B3)+N2(V1)       ","N+N+N2(V2)=>N2(B3)+N2(V2)       ",&
   "N+N+N2(V3)=>N2(B3)+N2(V3)       ","N+N+N2(V4)=>N2(B3)+N2(V4)       ","N+N+N2(V5)=>N2(B3)+N2(V5)       ",&
   "N+N+N2(V6)=>N2(B3)+N2(V6)       ","N+N+N2(V7)=>N2(B3)+N2(V7)       ","N+N+N2(V8)=>N2(B3)+N2(V8)       ",&
   "N+N+H2(V1)=>N2(B3)+H2(V1)       ","N+N+H2(V2)=>N2(B3)+H2(V2)       ","N+N+H2(V3)=>N2(B3)+H2(V3)       ",&
   "N+N+N=>N2(B3)+N                 ","N+N+H=>N2(B3)+H                 ","N+N+N2=>N2+N2                   ",&
   "N+N+H2=>N2+H2                   ","H+H+H2=>H2+H2                   ","H+H+N2=>H2+N2                   ",&
   "H+N+N2=>NH+N2                   ","H+N+H2=>NH+H2                   ","N+H2+N2=>NH2+N2                 ",&
   "N+H2+H2=>NH2+H2                 ","H+NH+N2=>NH2+N2                 ","H+NH+H2=>NH2+H2                 ",&
   "H+NH2+N2=>NH3+N2                ","H+NH2+H2=>NH3+H2                ","NH+H2+N2=>NH3+N2                ",&
   "NH+H2+H2=>NH3+H2                ","N2(A3)=>N2                      ","N2(A`1)=>N2(B3)                 ",&
   "N2(V1)=>N2                      ","N2(V2)=>N2(V1)                  ","N2(V3)=>N2(V2)                  ",&
   "N2(V4)=>N2(V3)                  ","N2(V5)=>N2(V4)                  ","N2(V6)=>N2(V5)                  "/
  data reaction_sign(433:486) &
  /"N2(V7)=>N2(V6)                  ","N2(V8)=>N2(V7)                  ","H2(B3SIG)=>H2                   ",&
   "H2(B1SIG)=>H2                   ","H2(C3PI)=>H2                    ","H2(A3SIG)=>H2                   ",&
   "H2(V1)=>H2                      ","H2(V2)=>H2(V1)                  ","H2(V3)=>H2(V2)                  ",&
   "bolsig:H2->H2(RYDBERG_SUM)      ","H^-+H2^+=>3H                    ","H^-+H3^+=>H2+2H                 ",&
   "H^-+N2^+=>N2+H                  ","H^-+N4^+=>2N2+H                 ","H^-+N2H^+=>H2+N2                ",&
   "H^-+H2^++N2=>H2+H+N2            ","H^-+H2^++H2=>H2+H+H2            ","H^-+H2^++N=>H2+H+N              ",&
   "H^-+H2^++H=>H2+H+H              ","H^-+H3^++N2=>2H2+N2             ","H^-+H3^++H2=>2H2+H2             ",&
   "H^-+H3^++N=>2H2+N               ","H^-+H3^++H=>2H2+H               ","H^-+N2^++N2=>N2+H+N2            ",&
   "H^-+N2^++H2=>N2+H+H2            ","H^-+N2^++N=>N2+H+N              ","H^-+N2^++H=>N2+H+H              ",&
   "H^-+N4^++N2=>2N2+H+N2           ","H^-+N4^++H2=>2N2+H+H2           ","H^-+N4^++N=>2N2+H+N             ",&
   "H^-+N4^++H=>2N2+H+H             ","H^-+N2H^++N2=>H2+N2+N2          ","H^-+N2H^++H2=>H2+N2+H2          ",&
   "H^-+N2H^++N=>H2+N2+N            ","H^-+N2H^++H=>H2+N2+H            ","N+SURF=>NSURF                   ",&
   "H+SURF=>HSURF                   ","NH+SURF=>NHSURF                 ","NH2+SURF=>NH2SURF               ",&
   "N+NSURF=>N2+SURF                ","H+HSURF=>H2+SURF                ","N+HSURF=>NHSURF                 ",&
   "NH+HSURF=>NH2SURF               ","NH2+HSURF=>NH3+SURF             ","H+NSURF=>NHSURF                 ",&
   "H+NHSURF=>NH2SURF               ","H+NH2SURF=>NH3+SURF             ","H2+NHSURF=>NH3+SURF             ",&
   "NSURF+HSURF=>NHSURF+SURF        ","NHSURF+HSURF=>NH2SURF+SURF      ","NH2SURF+HSURF=>NH3+2SURF        ",&
   "N2+SURF+SURF=>NSURF+NSURF       ","N2(V1)+SURF+SURF=>NSURF+NSURF   ","N2(V2)+SURF+SURF=>NSURF+NSURF   "/
  data reaction_sign(487:504) &
  /"N2(V3)+SURF+SURF=>NSURF+NSURF   ","N2(V4)+SURF+SURF=>NSURF+NSURF   ","N2(V5)+SURF+SURF=>NSURF+NSURF   ",&
   "N2(V6)+SURF+SURF=>NSURF+NSURF   ","N2(V7)+SURF+SURF=>NSURF+NSURF   ","N2(V8)+SURF+SURF=>NSURF+NSURF   ",&
   "N2(A3)+SURF+SURF=>NSURF+NSURF   ","N2(B3)+SURF+SURF=>NSURF+NSURF   ","N2(A`1)+SURF+SURF=>NSURF+NSURF  ",&
   "N2(C3)+SURF+SURF=>NSURF+NSURF   ","H2+SURF+SURF=>HSURF+HSURF       ","H2(V1)+SURF+SURF=>HSURF+HSURF   ",&
   "H2(V2)+SURF+SURF=>HSURF+HSURF   ","H2(V3)+SURF+SURF=>HSURF+HSURF   ","H2(B3SIG)+SURF+SURF=>HSURF+HSURF",&
   "H2(B1SIG)+SURF+SURF=>HSURF+HSURF","H2(C3PI)+SURF+SURF=>HSURF+HSURF ","H2(A3SIG)+SURF+SURF=>HSURF+HSURF"/
  data bolsig_species(1:bolsig_species_max) &
  /"N2             ","N2(V1)         ","N2(V2)         ","N2(V3)         ","N2(V4)         ","N2(V5)         ","N2(V6)         ",&
   "N2(V7)         ","N2(V8)         ","H2             ","H2(V1)         ","H2(V2)         ","H2(V3)         ","H2(B3SIG)      ",&
   "H2(B1SIG)      ","H2(C3PI)       ","H2(A3SIG)      ","H2(RYDBERG_SUM)"/
  data bolsig_addsect(1:2,1:bolsig_addsect_max) &
  / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9/
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
    write(*,"(A,1pd11.4)") "    EN =", ZDPlasKin_cfg(3)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    write(*,"(A,1pd11.4)") "   Dif =", ZDPlasKin_cfg(6)
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
111 format(i2,1x,A15)
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
211 format(i3,1x,A32)
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
311 format(371x,48(1x,i15))
312 format(A3,1x,A32,1x,48(1x,A15))
313 format(i3,1x,A32,1x,48(1x,1pd15.2))
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
      write(ifile_unit,"(504(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,48(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,504(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,48(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(505(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(01,001) = - reac_rate_local(001) 
  reac_source_local(10,001) = + reac_rate_local(001) 
  reac_source_local(01,002) = - reac_rate_local(002) 
  reac_source_local(10,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = - reac_rate_local(003) 
  reac_source_local(10,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = - reac_rate_local(004) 
  reac_source_local(11,004) = + reac_rate_local(004) 
  reac_source_local(01,005) = - reac_rate_local(005) 
  reac_source_local(12,005) = + reac_rate_local(005) 
  reac_source_local(01,006) = - reac_rate_local(006) 
  reac_source_local(13,006) = + reac_rate_local(006) 
  reac_source_local(21,007) = - reac_rate_local(007) 
  reac_source_local(22,007) = + reac_rate_local(007) 
  reac_source_local(21,008) = - reac_rate_local(008) 
  reac_source_local(23,008) = + reac_rate_local(008) 
  reac_source_local(21,009) = - reac_rate_local(009) 
  reac_source_local(24,009) = + reac_rate_local(009) 
  reac_source_local(21,010) = - reac_rate_local(010) 
  reac_source_local(25,010) = + reac_rate_local(010) 
  reac_source_local(21,011) = - reac_rate_local(011) 
  reac_source_local(26,011) = + reac_rate_local(011) 
  reac_source_local(01,012) = - reac_rate_local(012) 
  reac_source_local(18,012) = + reac_rate_local(012) 
  reac_source_local(43,012) = + reac_rate_local(012) 
  reac_source_local(21,013) = - reac_rate_local(013) 
  reac_source_local(31,013) = + reac_rate_local(013) 
  reac_source_local(43,013) = + reac_rate_local(013) 
  reac_source_local(14,014) = - reac_rate_local(014) 
  reac_source_local(17,014) = + reac_rate_local(014) 
  reac_source_local(43,014) = + reac_rate_local(014) 
  reac_source_local(30,015) = - reac_rate_local(015) 
  reac_source_local(33,015) = + reac_rate_local(015) 
  reac_source_local(43,015) = + reac_rate_local(015) 
  reac_source_local(37,016) = - reac_rate_local(016) 
  reac_source_local(41,016) = + reac_rate_local(016) 
  reac_source_local(43,016) = + reac_rate_local(016) 
  reac_source_local(01,017) = - reac_rate_local(017) 
  reac_source_local(02,017) = + reac_rate_local(017) 
  reac_source_local(01,018) = - reac_rate_local(018) 
  reac_source_local(03,018) = + reac_rate_local(018) 
  reac_source_local(01,019) = - reac_rate_local(019) 
  reac_source_local(04,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = - reac_rate_local(020) 
  reac_source_local(05,020) = + reac_rate_local(020) 
  reac_source_local(01,021) = - reac_rate_local(021) 
  reac_source_local(06,021) = + reac_rate_local(021) 
  reac_source_local(01,022) = - reac_rate_local(022) 
  reac_source_local(07,022) = + reac_rate_local(022) 
  reac_source_local(01,023) = - reac_rate_local(023) 
  reac_source_local(08,023) = + reac_rate_local(023) 
  reac_source_local(01,024) = - reac_rate_local(024) 
  reac_source_local(09,024) = + reac_rate_local(024) 
  reac_source_local(01,025) = + reac_rate_local(025) 
  reac_source_local(02,025) = - reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(03,026) = - reac_rate_local(026) 
  reac_source_local(01,027) = + reac_rate_local(027) 
  reac_source_local(04,027) = - reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(05,028) = - reac_rate_local(028) 
  reac_source_local(01,029) = + reac_rate_local(029) 
  reac_source_local(06,029) = - reac_rate_local(029) 
  reac_source_local(01,030) = + reac_rate_local(030) 
  reac_source_local(07,030) = - reac_rate_local(030) 
  reac_source_local(01,031) = + reac_rate_local(031) 
  reac_source_local(08,031) = - reac_rate_local(031) 
  reac_source_local(01,032) = + reac_rate_local(032) 
  reac_source_local(09,032) = - reac_rate_local(032) 
  reac_source_local(21,033) = - reac_rate_local(033) 
  reac_source_local(27,033) = + reac_rate_local(033) 
  reac_source_local(21,034) = - reac_rate_local(034) 
  reac_source_local(28,034) = + reac_rate_local(034) 
  reac_source_local(21,035) = - reac_rate_local(035) 
  reac_source_local(29,035) = + reac_rate_local(035) 
  reac_source_local(21,036) = + reac_rate_local(036) 
  reac_source_local(27,036) = - reac_rate_local(036) 
  reac_source_local(21,037) = + reac_rate_local(037) 
  reac_source_local(28,037) = - reac_rate_local(037) 
  reac_source_local(21,038) = + reac_rate_local(038) 
  reac_source_local(29,038) = - reac_rate_local(038) 
  reac_source_local(21,039) = + reac_rate_local(039) 
  reac_source_local(22,039) = - reac_rate_local(039) 
  reac_source_local(21,040) = + reac_rate_local(040) 
  reac_source_local(23,040) = - reac_rate_local(040) 
  reac_source_local(21,041) = + reac_rate_local(041) 
  reac_source_local(24,041) = - reac_rate_local(041) 
  reac_source_local(21,042) = + reac_rate_local(042) 
  reac_source_local(25,042) = - reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) 
  reac_source_local(02,043) = - reac_rate_local(043) 
  reac_source_local(01,044) = - reac_rate_local(044) 
  reac_source_local(02,044) = + reac_rate_local(044) 
  reac_source_local(02,045) = + reac_rate_local(045) 
  reac_source_local(03,045) = - reac_rate_local(045) 
  reac_source_local(02,046) = - reac_rate_local(046) 
  reac_source_local(03,046) = + reac_rate_local(046) 
  reac_source_local(03,047) = + reac_rate_local(047) 
  reac_source_local(04,047) = - reac_rate_local(047) 
  reac_source_local(03,048) = - reac_rate_local(048) 
  reac_source_local(04,048) = + reac_rate_local(048) 
  reac_source_local(04,049) = + reac_rate_local(049) 
  reac_source_local(05,049) = - reac_rate_local(049) 
  reac_source_local(04,050) = - reac_rate_local(050) 
  reac_source_local(05,050) = + reac_rate_local(050) 
  reac_source_local(05,051) = + reac_rate_local(051) 
  reac_source_local(06,051) = - reac_rate_local(051) 
  reac_source_local(05,052) = - reac_rate_local(052) 
  reac_source_local(06,052) = + reac_rate_local(052) 
  reac_source_local(06,053) = + reac_rate_local(053) 
  reac_source_local(07,053) = - reac_rate_local(053) 
  reac_source_local(06,054) = - reac_rate_local(054) 
  reac_source_local(07,054) = + reac_rate_local(054) 
  reac_source_local(07,055) = + reac_rate_local(055) 
  reac_source_local(08,055) = - reac_rate_local(055) 
  reac_source_local(07,056) = - reac_rate_local(056) 
  reac_source_local(08,056) = + reac_rate_local(056) 
  reac_source_local(08,057) = + reac_rate_local(057) 
  reac_source_local(09,057) = - reac_rate_local(057) 
  reac_source_local(08,058) = - reac_rate_local(058) 
  reac_source_local(09,058) = + reac_rate_local(058) 
  reac_source_local(01,059) = + reac_rate_local(059) 
  reac_source_local(02,059) = - reac_rate_local(059) 
  reac_source_local(01,060) = - reac_rate_local(060) 
  reac_source_local(02,060) = + reac_rate_local(060) 
  reac_source_local(02,061) = + reac_rate_local(061) 
  reac_source_local(03,061) = - reac_rate_local(061) 
  reac_source_local(02,062) = - reac_rate_local(062) 
  reac_source_local(03,062) = + reac_rate_local(062) 
  reac_source_local(03,063) = + reac_rate_local(063) 
  reac_source_local(04,063) = - reac_rate_local(063) 
  reac_source_local(03,064) = - reac_rate_local(064) 
  reac_source_local(04,064) = + reac_rate_local(064) 
  reac_source_local(04,065) = + reac_rate_local(065) 
  reac_source_local(05,065) = - reac_rate_local(065) 
  reac_source_local(04,066) = - reac_rate_local(066) 
  reac_source_local(05,066) = + reac_rate_local(066) 
  reac_source_local(05,067) = + reac_rate_local(067) 
  reac_source_local(06,067) = - reac_rate_local(067) 
  reac_source_local(05,068) = - reac_rate_local(068) 
  reac_source_local(06,068) = + reac_rate_local(068) 
  reac_source_local(06,069) = + reac_rate_local(069) 
  reac_source_local(07,069) = - reac_rate_local(069) 
  reac_source_local(06,070) = - reac_rate_local(070) 
  reac_source_local(07,070) = + reac_rate_local(070) 
  reac_source_local(07,071) = + reac_rate_local(071) 
  reac_source_local(08,071) = - reac_rate_local(071) 
  reac_source_local(07,072) = - reac_rate_local(072) 
  reac_source_local(08,072) = + reac_rate_local(072) 
  reac_source_local(08,073) = + reac_rate_local(073) 
  reac_source_local(09,073) = - reac_rate_local(073) 
  reac_source_local(08,074) = - reac_rate_local(074) 
  reac_source_local(09,074) = + reac_rate_local(074) 
  reac_source_local(01,075) = + reac_rate_local(075) 
  reac_source_local(02,075) = - reac_rate_local(075) 
  reac_source_local(01,076) = - reac_rate_local(076) 
  reac_source_local(02,076) = + reac_rate_local(076) 
  reac_source_local(02,077) = + reac_rate_local(077) 
  reac_source_local(03,077) = - reac_rate_local(077) 
  reac_source_local(02,078) = - reac_rate_local(078) 
  reac_source_local(03,078) = + reac_rate_local(078) 
  reac_source_local(03,079) = + reac_rate_local(079) 
  reac_source_local(04,079) = - reac_rate_local(079) 
  reac_source_local(03,080) = - reac_rate_local(080) 
  reac_source_local(04,080) = + reac_rate_local(080) 
  reac_source_local(04,081) = + reac_rate_local(081) 
  reac_source_local(05,081) = - reac_rate_local(081) 
  reac_source_local(04,082) = - reac_rate_local(082) 
  reac_source_local(05,082) = + reac_rate_local(082) 
  reac_source_local(05,083) = + reac_rate_local(083) 
  reac_source_local(06,083) = - reac_rate_local(083) 
  reac_source_local(05,084) = - reac_rate_local(084) 
  reac_source_local(06,084) = + reac_rate_local(084) 
  reac_source_local(06,085) = + reac_rate_local(085) 
  reac_source_local(07,085) = - reac_rate_local(085) 
  reac_source_local(06,086) = - reac_rate_local(086) 
  reac_source_local(07,086) = + reac_rate_local(086) 
  reac_source_local(07,087) = + reac_rate_local(087) 
  reac_source_local(08,087) = - reac_rate_local(087) 
  reac_source_local(07,088) = - reac_rate_local(088) 
  reac_source_local(08,088) = + reac_rate_local(088) 
  reac_source_local(08,089) = + reac_rate_local(089) 
  reac_source_local(09,089) = - reac_rate_local(089) 
  reac_source_local(08,090) = - reac_rate_local(090) 
  reac_source_local(09,090) = + reac_rate_local(090) 
  reac_source_local(01,091) = + reac_rate_local(091) 
  reac_source_local(02,091) = - reac_rate_local(091) 
  reac_source_local(01,092) = - reac_rate_local(092) 
  reac_source_local(02,092) = + reac_rate_local(092) 
  reac_source_local(02,093) = + reac_rate_local(093) 
  reac_source_local(03,093) = - reac_rate_local(093) 
  reac_source_local(02,094) = - reac_rate_local(094) 
  reac_source_local(03,094) = + reac_rate_local(094) 
  reac_source_local(03,095) = + reac_rate_local(095) 
  reac_source_local(04,095) = - reac_rate_local(095) 
  reac_source_local(03,096) = - reac_rate_local(096) 
  reac_source_local(04,096) = + reac_rate_local(096) 
  reac_source_local(04,097) = + reac_rate_local(097) 
  reac_source_local(05,097) = - reac_rate_local(097) 
  reac_source_local(04,098) = - reac_rate_local(098) 
  reac_source_local(05,098) = + reac_rate_local(098) 
  reac_source_local(05,099) = + reac_rate_local(099) 
  reac_source_local(06,099) = - reac_rate_local(099) 
  reac_source_local(05,100) = - reac_rate_local(100) 
  reac_source_local(06,100) = + reac_rate_local(100) 
  reac_source_local(06,101) = + reac_rate_local(101) 
  reac_source_local(07,101) = - reac_rate_local(101) 
  reac_source_local(06,102) = - reac_rate_local(102) 
  reac_source_local(07,102) = + reac_rate_local(102) 
  reac_source_local(07,103) = + reac_rate_local(103) 
  reac_source_local(08,103) = - reac_rate_local(103) 
  reac_source_local(07,104) = - reac_rate_local(104) 
  reac_source_local(08,104) = + reac_rate_local(104) 
  reac_source_local(08,105) = + reac_rate_local(105) 
  reac_source_local(09,105) = - reac_rate_local(105) 
  reac_source_local(08,106) = - reac_rate_local(106) 
  reac_source_local(09,106) = + reac_rate_local(106) 
  reac_source_local(21,107) = + reac_rate_local(107) 
  reac_source_local(27,107) = - reac_rate_local(107) 
  reac_source_local(21,108) = - reac_rate_local(108) 
  reac_source_local(27,108) = + reac_rate_local(108) 
  reac_source_local(27,109) = + reac_rate_local(109) 
  reac_source_local(28,109) = - reac_rate_local(109) 
  reac_source_local(27,110) = - reac_rate_local(110) 
  reac_source_local(28,110) = + reac_rate_local(110) 
  reac_source_local(28,111) = + reac_rate_local(111) 
  reac_source_local(29,111) = - reac_rate_local(111) 
  reac_source_local(28,112) = - reac_rate_local(112) 
  reac_source_local(29,112) = + reac_rate_local(112) 
  reac_source_local(21,113) = + reac_rate_local(113) 
  reac_source_local(27,113) = - reac_rate_local(113) 
  reac_source_local(21,114) = - reac_rate_local(114) 
  reac_source_local(27,114) = + reac_rate_local(114) 
  reac_source_local(27,115) = + reac_rate_local(115) 
  reac_source_local(28,115) = - reac_rate_local(115) 
  reac_source_local(27,116) = - reac_rate_local(116) 
  reac_source_local(28,116) = + reac_rate_local(116) 
  reac_source_local(28,117) = + reac_rate_local(117) 
  reac_source_local(29,117) = - reac_rate_local(117) 
  reac_source_local(28,118) = - reac_rate_local(118) 
  reac_source_local(29,118) = + reac_rate_local(118) 
  reac_source_local(01,119) = + reac_rate_local(119) 
  reac_source_local(02,119) = - reac_rate_local(119) * 2.d0
  reac_source_local(03,119) = + reac_rate_local(119) 
  reac_source_local(01,120) = - reac_rate_local(120) 
  reac_source_local(02,120) = + reac_rate_local(120) * 2.d0
  reac_source_local(03,120) = - reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(02,121) = - reac_rate_local(121) 
  reac_source_local(03,121) = - reac_rate_local(121) 
  reac_source_local(04,121) = + reac_rate_local(121) 
  reac_source_local(01,122) = - reac_rate_local(122) 
  reac_source_local(02,122) = + reac_rate_local(122) 
  reac_source_local(03,122) = + reac_rate_local(122) 
  reac_source_local(04,122) = - reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(02,123) = - reac_rate_local(123) 
  reac_source_local(04,123) = - reac_rate_local(123) 
  reac_source_local(05,123) = + reac_rate_local(123) 
  reac_source_local(01,124) = - reac_rate_local(124) 
  reac_source_local(02,124) = + reac_rate_local(124) 
  reac_source_local(04,124) = + reac_rate_local(124) 
  reac_source_local(05,124) = - reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(02,125) = - reac_rate_local(125) 
  reac_source_local(05,125) = - reac_rate_local(125) 
  reac_source_local(06,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = - reac_rate_local(126) 
  reac_source_local(02,126) = + reac_rate_local(126) 
  reac_source_local(05,126) = + reac_rate_local(126) 
  reac_source_local(06,126) = - reac_rate_local(126) 
  reac_source_local(01,127) = + reac_rate_local(127) 
  reac_source_local(02,127) = - reac_rate_local(127) 
  reac_source_local(06,127) = - reac_rate_local(127) 
  reac_source_local(07,127) = + reac_rate_local(127) 
  reac_source_local(01,128) = - reac_rate_local(128) 
  reac_source_local(02,128) = + reac_rate_local(128) 
  reac_source_local(06,128) = + reac_rate_local(128) 
  reac_source_local(07,128) = - reac_rate_local(128) 
  reac_source_local(01,129) = + reac_rate_local(129) 
  reac_source_local(02,129) = - reac_rate_local(129) 
  reac_source_local(07,129) = - reac_rate_local(129) 
  reac_source_local(08,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = - reac_rate_local(130) 
  reac_source_local(02,130) = + reac_rate_local(130) 
  reac_source_local(07,130) = + reac_rate_local(130) 
  reac_source_local(08,130) = - reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(02,131) = - reac_rate_local(131) 
  reac_source_local(08,131) = - reac_rate_local(131) 
  reac_source_local(09,131) = + reac_rate_local(131) 
  reac_source_local(01,132) = - reac_rate_local(132) 
  reac_source_local(02,132) = + reac_rate_local(132) 
  reac_source_local(08,132) = + reac_rate_local(132) 
  reac_source_local(09,132) = - reac_rate_local(132) 
  reac_source_local(02,133) = + reac_rate_local(133) 
  reac_source_local(03,133) = - reac_rate_local(133) * 2.d0
  reac_source_local(04,133) = + reac_rate_local(133) 
  reac_source_local(02,134) = - reac_rate_local(134) 
  reac_source_local(03,134) = + reac_rate_local(134) * 2.d0
  reac_source_local(04,134) = - reac_rate_local(134) 
  reac_source_local(02,135) = + reac_rate_local(135) 
  reac_source_local(03,135) = - reac_rate_local(135) 
  reac_source_local(04,135) = - reac_rate_local(135) 
  reac_source_local(05,135) = + reac_rate_local(135) 
  reac_source_local(02,136) = - reac_rate_local(136) 
  reac_source_local(03,136) = + reac_rate_local(136) 
  reac_source_local(04,136) = + reac_rate_local(136) 
  reac_source_local(05,136) = - reac_rate_local(136) 
  reac_source_local(02,137) = + reac_rate_local(137) 
  reac_source_local(03,137) = - reac_rate_local(137) 
  reac_source_local(05,137) = - reac_rate_local(137) 
  reac_source_local(06,137) = + reac_rate_local(137) 
  reac_source_local(02,138) = - reac_rate_local(138) 
  reac_source_local(03,138) = + reac_rate_local(138) 
  reac_source_local(05,138) = + reac_rate_local(138) 
  reac_source_local(06,138) = - reac_rate_local(138) 
  reac_source_local(02,139) = + reac_rate_local(139) 
  reac_source_local(03,139) = - reac_rate_local(139) 
  reac_source_local(06,139) = - reac_rate_local(139) 
  reac_source_local(07,139) = + reac_rate_local(139) 
  reac_source_local(02,140) = - reac_rate_local(140) 
  reac_source_local(03,140) = + reac_rate_local(140) 
  reac_source_local(06,140) = + reac_rate_local(140) 
  reac_source_local(07,140) = - reac_rate_local(140) 
  reac_source_local(02,141) = + reac_rate_local(141) 
  reac_source_local(03,141) = - reac_rate_local(141) 
  reac_source_local(07,141) = - reac_rate_local(141) 
  reac_source_local(08,141) = + reac_rate_local(141) 
  reac_source_local(02,142) = - reac_rate_local(142) 
  reac_source_local(03,142) = + reac_rate_local(142) 
  reac_source_local(07,142) = + reac_rate_local(142) 
  reac_source_local(08,142) = - reac_rate_local(142) 
  reac_source_local(02,143) = + reac_rate_local(143) 
  reac_source_local(03,143) = - reac_rate_local(143) 
  reac_source_local(08,143) = - reac_rate_local(143) 
  reac_source_local(09,143) = + reac_rate_local(143) 
  reac_source_local(02,144) = - reac_rate_local(144) 
  reac_source_local(03,144) = + reac_rate_local(144) 
  reac_source_local(08,144) = + reac_rate_local(144) 
  reac_source_local(09,144) = - reac_rate_local(144) 
  reac_source_local(03,145) = + reac_rate_local(145) 
  reac_source_local(04,145) = - reac_rate_local(145) * 2.d0
  reac_source_local(05,145) = + reac_rate_local(145) 
  reac_source_local(03,146) = - reac_rate_local(146) 
  reac_source_local(04,146) = + reac_rate_local(146) * 2.d0
  reac_source_local(05,146) = - reac_rate_local(146) 
  reac_source_local(03,147) = + reac_rate_local(147) 
  reac_source_local(04,147) = - reac_rate_local(147) 
  reac_source_local(05,147) = - reac_rate_local(147) 
  reac_source_local(06,147) = + reac_rate_local(147) 
  reac_source_local(03,148) = - reac_rate_local(148) 
  reac_source_local(04,148) = + reac_rate_local(148) 
  reac_source_local(05,148) = + reac_rate_local(148) 
  reac_source_local(06,148) = - reac_rate_local(148) 
  reac_source_local(03,149) = + reac_rate_local(149) 
  reac_source_local(04,149) = - reac_rate_local(149) 
  reac_source_local(06,149) = - reac_rate_local(149) 
  reac_source_local(07,149) = + reac_rate_local(149) 
  reac_source_local(03,150) = - reac_rate_local(150) 
  reac_source_local(04,150) = + reac_rate_local(150) 
  reac_source_local(06,150) = + reac_rate_local(150) 
  reac_source_local(07,150) = - reac_rate_local(150) 
  reac_source_local(03,151) = + reac_rate_local(151) 
  reac_source_local(04,151) = - reac_rate_local(151) 
  reac_source_local(07,151) = - reac_rate_local(151) 
  reac_source_local(08,151) = + reac_rate_local(151) 
  reac_source_local(03,152) = - reac_rate_local(152) 
  reac_source_local(04,152) = + reac_rate_local(152) 
  reac_source_local(07,152) = + reac_rate_local(152) 
  reac_source_local(08,152) = - reac_rate_local(152) 
  reac_source_local(03,153) = + reac_rate_local(153) 
  reac_source_local(04,153) = - reac_rate_local(153) 
  reac_source_local(08,153) = - reac_rate_local(153) 
  reac_source_local(09,153) = + reac_rate_local(153) 
  reac_source_local(03,154) = - reac_rate_local(154) 
  reac_source_local(04,154) = + reac_rate_local(154) 
  reac_source_local(08,154) = + reac_rate_local(154) 
  reac_source_local(09,154) = - reac_rate_local(154) 
  reac_source_local(04,155) = + reac_rate_local(155) 
  reac_source_local(05,155) = - reac_rate_local(155) * 2.d0
  reac_source_local(06,155) = + reac_rate_local(155) 
  reac_source_local(04,156) = - reac_rate_local(156) 
  reac_source_local(05,156) = + reac_rate_local(156) * 2.d0
  reac_source_local(06,156) = - reac_rate_local(156) 
  reac_source_local(04,157) = + reac_rate_local(157) 
  reac_source_local(05,157) = - reac_rate_local(157) 
  reac_source_local(06,157) = - reac_rate_local(157) 
  reac_source_local(07,157) = + reac_rate_local(157) 
  reac_source_local(04,158) = - reac_rate_local(158) 
  reac_source_local(05,158) = + reac_rate_local(158) 
  reac_source_local(06,158) = + reac_rate_local(158) 
  reac_source_local(07,158) = - reac_rate_local(158) 
  reac_source_local(04,159) = + reac_rate_local(159) 
  reac_source_local(05,159) = - reac_rate_local(159) 
  reac_source_local(07,159) = - reac_rate_local(159) 
  reac_source_local(08,159) = + reac_rate_local(159) 
  reac_source_local(04,160) = - reac_rate_local(160) 
  reac_source_local(05,160) = + reac_rate_local(160) 
  reac_source_local(07,160) = + reac_rate_local(160) 
  reac_source_local(08,160) = - reac_rate_local(160) 
  reac_source_local(04,161) = + reac_rate_local(161) 
  reac_source_local(05,161) = - reac_rate_local(161) 
  reac_source_local(08,161) = - reac_rate_local(161) 
  reac_source_local(09,161) = + reac_rate_local(161) 
  reac_source_local(04,162) = - reac_rate_local(162) 
  reac_source_local(05,162) = + reac_rate_local(162) 
  reac_source_local(08,162) = + reac_rate_local(162) 
  reac_source_local(09,162) = - reac_rate_local(162) 
  reac_source_local(05,163) = + reac_rate_local(163) 
  reac_source_local(06,163) = - reac_rate_local(163) * 2.d0
  reac_source_local(07,163) = + reac_rate_local(163) 
  reac_source_local(05,164) = - reac_rate_local(164) 
  reac_source_local(06,164) = + reac_rate_local(164) * 2.d0
  reac_source_local(07,164) = - reac_rate_local(164) 
  reac_source_local(05,165) = + reac_rate_local(165) 
  reac_source_local(06,165) = - reac_rate_local(165) 
  reac_source_local(07,165) = - reac_rate_local(165) 
  reac_source_local(08,165) = + reac_rate_local(165) 
  reac_source_local(05,166) = - reac_rate_local(166) 
  reac_source_local(06,166) = + reac_rate_local(166) 
  reac_source_local(07,166) = + reac_rate_local(166) 
  reac_source_local(08,166) = - reac_rate_local(166) 
  reac_source_local(05,167) = + reac_rate_local(167) 
  reac_source_local(06,167) = - reac_rate_local(167) 
  reac_source_local(08,167) = - reac_rate_local(167) 
  reac_source_local(09,167) = + reac_rate_local(167) 
  reac_source_local(05,168) = - reac_rate_local(168) 
  reac_source_local(06,168) = + reac_rate_local(168) 
  reac_source_local(08,168) = + reac_rate_local(168) 
  reac_source_local(09,168) = - reac_rate_local(168) 
  reac_source_local(06,169) = + reac_rate_local(169) 
  reac_source_local(07,169) = - reac_rate_local(169) * 2.d0
  reac_source_local(08,169) = + reac_rate_local(169) 
  reac_source_local(06,170) = - reac_rate_local(170) 
  reac_source_local(07,170) = + reac_rate_local(170) * 2.d0
  reac_source_local(08,170) = - reac_rate_local(170) 
  reac_source_local(06,171) = + reac_rate_local(171) 
  reac_source_local(07,171) = - reac_rate_local(171) 
  reac_source_local(08,171) = - reac_rate_local(171) 
  reac_source_local(09,171) = + reac_rate_local(171) 
  reac_source_local(06,172) = - reac_rate_local(172) 
  reac_source_local(07,172) = + reac_rate_local(172) 
  reac_source_local(08,172) = + reac_rate_local(172) 
  reac_source_local(09,172) = - reac_rate_local(172) 
  reac_source_local(07,173) = + reac_rate_local(173) 
  reac_source_local(08,173) = - reac_rate_local(173) * 2.d0
  reac_source_local(09,173) = + reac_rate_local(173) 
  reac_source_local(07,174) = - reac_rate_local(174) 
  reac_source_local(08,174) = + reac_rate_local(174) * 2.d0
  reac_source_local(09,174) = - reac_rate_local(174) 
  reac_source_local(01,175) = - reac_rate_local(175) 
  reac_source_local(02,175) = + reac_rate_local(175) 
  reac_source_local(21,175) = + reac_rate_local(175) 
  reac_source_local(27,175) = - reac_rate_local(175) 
  reac_source_local(01,176) = + reac_rate_local(176) 
  reac_source_local(02,176) = - reac_rate_local(176) 
  reac_source_local(21,176) = - reac_rate_local(176) 
  reac_source_local(27,176) = + reac_rate_local(176) 
  reac_source_local(02,177) = - reac_rate_local(177) 
  reac_source_local(03,177) = + reac_rate_local(177) 
  reac_source_local(21,177) = + reac_rate_local(177) 
  reac_source_local(27,177) = - reac_rate_local(177) 
  reac_source_local(02,178) = + reac_rate_local(178) 
  reac_source_local(03,178) = - reac_rate_local(178) 
  reac_source_local(21,178) = - reac_rate_local(178) 
  reac_source_local(27,178) = + reac_rate_local(178) 
  reac_source_local(03,179) = - reac_rate_local(179) 
  reac_source_local(04,179) = + reac_rate_local(179) 
  reac_source_local(21,179) = + reac_rate_local(179) 
  reac_source_local(27,179) = - reac_rate_local(179) 
  reac_source_local(03,180) = + reac_rate_local(180) 
  reac_source_local(04,180) = - reac_rate_local(180) 
  reac_source_local(21,180) = - reac_rate_local(180) 
  reac_source_local(27,180) = + reac_rate_local(180) 
  reac_source_local(04,181) = - reac_rate_local(181) 
  reac_source_local(05,181) = + reac_rate_local(181) 
  reac_source_local(21,181) = + reac_rate_local(181) 
  reac_source_local(27,181) = - reac_rate_local(181) 
  reac_source_local(04,182) = + reac_rate_local(182) 
  reac_source_local(05,182) = - reac_rate_local(182) 
  reac_source_local(21,182) = - reac_rate_local(182) 
  reac_source_local(27,182) = + reac_rate_local(182) 
  reac_source_local(05,183) = - reac_rate_local(183) 
  reac_source_local(06,183) = + reac_rate_local(183) 
  reac_source_local(21,183) = + reac_rate_local(183) 
  reac_source_local(27,183) = - reac_rate_local(183) 
  reac_source_local(05,184) = + reac_rate_local(184) 
  reac_source_local(06,184) = - reac_rate_local(184) 
  reac_source_local(21,184) = - reac_rate_local(184) 
  reac_source_local(27,184) = + reac_rate_local(184) 
  reac_source_local(06,185) = - reac_rate_local(185) 
  reac_source_local(07,185) = + reac_rate_local(185) 
  reac_source_local(21,185) = + reac_rate_local(185) 
  reac_source_local(27,185) = - reac_rate_local(185) 
  reac_source_local(06,186) = + reac_rate_local(186) 
  reac_source_local(07,186) = - reac_rate_local(186) 
  reac_source_local(21,186) = - reac_rate_local(186) 
  reac_source_local(27,186) = + reac_rate_local(186) 
  reac_source_local(07,187) = - reac_rate_local(187) 
  reac_source_local(08,187) = + reac_rate_local(187) 
  reac_source_local(21,187) = + reac_rate_local(187) 
  reac_source_local(27,187) = - reac_rate_local(187) 
  reac_source_local(07,188) = + reac_rate_local(188) 
  reac_source_local(08,188) = - reac_rate_local(188) 
  reac_source_local(21,188) = - reac_rate_local(188) 
  reac_source_local(27,188) = + reac_rate_local(188) 
  reac_source_local(08,189) = - reac_rate_local(189) 
  reac_source_local(09,189) = + reac_rate_local(189) 
  reac_source_local(21,189) = + reac_rate_local(189) 
  reac_source_local(27,189) = - reac_rate_local(189) 
  reac_source_local(08,190) = + reac_rate_local(190) 
  reac_source_local(09,190) = - reac_rate_local(190) 
  reac_source_local(21,190) = - reac_rate_local(190) 
  reac_source_local(27,190) = + reac_rate_local(190) 
  reac_source_local(01,191) = - reac_rate_local(191) 
  reac_source_local(02,191) = + reac_rate_local(191) 
  reac_source_local(27,191) = + reac_rate_local(191) 
  reac_source_local(28,191) = - reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(02,192) = - reac_rate_local(192) 
  reac_source_local(27,192) = - reac_rate_local(192) 
  reac_source_local(28,192) = + reac_rate_local(192) 
  reac_source_local(02,193) = - reac_rate_local(193) 
  reac_source_local(03,193) = + reac_rate_local(193) 
  reac_source_local(27,193) = + reac_rate_local(193) 
  reac_source_local(28,193) = - reac_rate_local(193) 
  reac_source_local(02,194) = + reac_rate_local(194) 
  reac_source_local(03,194) = - reac_rate_local(194) 
  reac_source_local(27,194) = - reac_rate_local(194) 
  reac_source_local(28,194) = + reac_rate_local(194) 
  reac_source_local(03,195) = - reac_rate_local(195) 
  reac_source_local(04,195) = + reac_rate_local(195) 
  reac_source_local(27,195) = + reac_rate_local(195) 
  reac_source_local(28,195) = - reac_rate_local(195) 
  reac_source_local(03,196) = + reac_rate_local(196) 
  reac_source_local(04,196) = - reac_rate_local(196) 
  reac_source_local(27,196) = - reac_rate_local(196) 
  reac_source_local(28,196) = + reac_rate_local(196) 
  reac_source_local(04,197) = - reac_rate_local(197) 
  reac_source_local(05,197) = + reac_rate_local(197) 
  reac_source_local(27,197) = + reac_rate_local(197) 
  reac_source_local(28,197) = - reac_rate_local(197) 
  reac_source_local(04,198) = + reac_rate_local(198) 
  reac_source_local(05,198) = - reac_rate_local(198) 
  reac_source_local(27,198) = - reac_rate_local(198) 
  reac_source_local(28,198) = + reac_rate_local(198) 
  reac_source_local(05,199) = - reac_rate_local(199) 
  reac_source_local(06,199) = + reac_rate_local(199) 
  reac_source_local(27,199) = + reac_rate_local(199) 
  reac_source_local(28,199) = - reac_rate_local(199) 
  reac_source_local(05,200) = + reac_rate_local(200) 
  reac_source_local(06,200) = - reac_rate_local(200) 
  reac_source_local(27,200) = - reac_rate_local(200) 
  reac_source_local(28,200) = + reac_rate_local(200) 
  reac_source_local(06,201) = - reac_rate_local(201) 
  reac_source_local(07,201) = + reac_rate_local(201) 
  reac_source_local(27,201) = + reac_rate_local(201) 
  reac_source_local(28,201) = - reac_rate_local(201) 
  reac_source_local(06,202) = + reac_rate_local(202) 
  reac_source_local(07,202) = - reac_rate_local(202) 
  reac_source_local(27,202) = - reac_rate_local(202) 
  reac_source_local(28,202) = + reac_rate_local(202) 
  reac_source_local(07,203) = - reac_rate_local(203) 
  reac_source_local(08,203) = + reac_rate_local(203) 
  reac_source_local(27,203) = + reac_rate_local(203) 
  reac_source_local(28,203) = - reac_rate_local(203) 
  reac_source_local(07,204) = + reac_rate_local(204) 
  reac_source_local(08,204) = - reac_rate_local(204) 
  reac_source_local(27,204) = - reac_rate_local(204) 
  reac_source_local(28,204) = + reac_rate_local(204) 
  reac_source_local(08,205) = - reac_rate_local(205) 
  reac_source_local(09,205) = + reac_rate_local(205) 
  reac_source_local(27,205) = + reac_rate_local(205) 
  reac_source_local(28,205) = - reac_rate_local(205) 
  reac_source_local(08,206) = + reac_rate_local(206) 
  reac_source_local(09,206) = - reac_rate_local(206) 
  reac_source_local(27,206) = - reac_rate_local(206) 
  reac_source_local(28,206) = + reac_rate_local(206) 
  reac_source_local(01,207) = - reac_rate_local(207) 
  reac_source_local(02,207) = + reac_rate_local(207) 
  reac_source_local(28,207) = + reac_rate_local(207) 
  reac_source_local(29,207) = - reac_rate_local(207) 
  reac_source_local(01,208) = + reac_rate_local(208) 
  reac_source_local(02,208) = - reac_rate_local(208) 
  reac_source_local(28,208) = - reac_rate_local(208) 
  reac_source_local(29,208) = + reac_rate_local(208) 
  reac_source_local(02,209) = - reac_rate_local(209) 
  reac_source_local(03,209) = + reac_rate_local(209) 
  reac_source_local(28,209) = + reac_rate_local(209) 
  reac_source_local(29,209) = - reac_rate_local(209) 
  reac_source_local(02,210) = + reac_rate_local(210) 
  reac_source_local(03,210) = - reac_rate_local(210) 
  reac_source_local(28,210) = - reac_rate_local(210) 
  reac_source_local(29,210) = + reac_rate_local(210) 
  reac_source_local(03,211) = - reac_rate_local(211) 
  reac_source_local(04,211) = + reac_rate_local(211) 
  reac_source_local(28,211) = + reac_rate_local(211) 
  reac_source_local(29,211) = - reac_rate_local(211) 
  reac_source_local(03,212) = + reac_rate_local(212) 
  reac_source_local(04,212) = - reac_rate_local(212) 
  reac_source_local(28,212) = - reac_rate_local(212) 
  reac_source_local(29,212) = + reac_rate_local(212) 
  reac_source_local(04,213) = - reac_rate_local(213) 
  reac_source_local(05,213) = + reac_rate_local(213) 
  reac_source_local(28,213) = + reac_rate_local(213) 
  reac_source_local(29,213) = - reac_rate_local(213) 
  reac_source_local(04,214) = + reac_rate_local(214) 
  reac_source_local(05,214) = - reac_rate_local(214) 
  reac_source_local(28,214) = - reac_rate_local(214) 
  reac_source_local(29,214) = + reac_rate_local(214) 
  reac_source_local(05,215) = - reac_rate_local(215) 
  reac_source_local(06,215) = + reac_rate_local(215) 
  reac_source_local(28,215) = + reac_rate_local(215) 
  reac_source_local(29,215) = - reac_rate_local(215) 
  reac_source_local(05,216) = + reac_rate_local(216) 
  reac_source_local(06,216) = - reac_rate_local(216) 
  reac_source_local(28,216) = - reac_rate_local(216) 
  reac_source_local(29,216) = + reac_rate_local(216) 
  reac_source_local(06,217) = - reac_rate_local(217) 
  reac_source_local(07,217) = + reac_rate_local(217) 
  reac_source_local(28,217) = + reac_rate_local(217) 
  reac_source_local(29,217) = - reac_rate_local(217) 
  reac_source_local(06,218) = + reac_rate_local(218) 
  reac_source_local(07,218) = - reac_rate_local(218) 
  reac_source_local(28,218) = - reac_rate_local(218) 
  reac_source_local(29,218) = + reac_rate_local(218) 
  reac_source_local(07,219) = - reac_rate_local(219) 
  reac_source_local(08,219) = + reac_rate_local(219) 
  reac_source_local(28,219) = + reac_rate_local(219) 
  reac_source_local(29,219) = - reac_rate_local(219) 
  reac_source_local(07,220) = + reac_rate_local(220) 
  reac_source_local(08,220) = - reac_rate_local(220) 
  reac_source_local(28,220) = - reac_rate_local(220) 
  reac_source_local(29,220) = + reac_rate_local(220) 
  reac_source_local(08,221) = - reac_rate_local(221) 
  reac_source_local(09,221) = + reac_rate_local(221) 
  reac_source_local(28,221) = + reac_rate_local(221) 
  reac_source_local(29,221) = - reac_rate_local(221) 
  reac_source_local(08,222) = + reac_rate_local(222) 
  reac_source_local(09,222) = - reac_rate_local(222) 
  reac_source_local(28,222) = - reac_rate_local(222) 
  reac_source_local(29,222) = + reac_rate_local(222) 
  reac_source_local(01,223) = + reac_rate_local(223) 
  reac_source_local(03,223) = - reac_rate_local(223) 
  reac_source_local(21,223) = - reac_rate_local(223) 
  reac_source_local(27,223) = + reac_rate_local(223) 
  reac_source_local(01,224) = - reac_rate_local(224) 
  reac_source_local(03,224) = + reac_rate_local(224) 
  reac_source_local(21,224) = + reac_rate_local(224) 
  reac_source_local(27,224) = - reac_rate_local(224) 
  reac_source_local(02,225) = + reac_rate_local(225) 
  reac_source_local(04,225) = - reac_rate_local(225) 
  reac_source_local(21,225) = - reac_rate_local(225) 
  reac_source_local(27,225) = + reac_rate_local(225) 
  reac_source_local(02,226) = - reac_rate_local(226) 
  reac_source_local(04,226) = + reac_rate_local(226) 
  reac_source_local(21,226) = + reac_rate_local(226) 
  reac_source_local(27,226) = - reac_rate_local(226) 
  reac_source_local(03,227) = + reac_rate_local(227) 
  reac_source_local(05,227) = - reac_rate_local(227) 
  reac_source_local(21,227) = - reac_rate_local(227) 
  reac_source_local(27,227) = + reac_rate_local(227) 
  reac_source_local(03,228) = - reac_rate_local(228) 
  reac_source_local(05,228) = + reac_rate_local(228) 
  reac_source_local(21,228) = + reac_rate_local(228) 
  reac_source_local(27,228) = - reac_rate_local(228) 
  reac_source_local(04,229) = + reac_rate_local(229) 
  reac_source_local(06,229) = - reac_rate_local(229) 
  reac_source_local(21,229) = - reac_rate_local(229) 
  reac_source_local(27,229) = + reac_rate_local(229) 
  reac_source_local(04,230) = - reac_rate_local(230) 
  reac_source_local(06,230) = + reac_rate_local(230) 
  reac_source_local(21,230) = + reac_rate_local(230) 
  reac_source_local(27,230) = - reac_rate_local(230) 
  reac_source_local(05,231) = + reac_rate_local(231) 
  reac_source_local(07,231) = - reac_rate_local(231) 
  reac_source_local(21,231) = - reac_rate_local(231) 
  reac_source_local(27,231) = + reac_rate_local(231) 
  reac_source_local(05,232) = - reac_rate_local(232) 
  reac_source_local(07,232) = + reac_rate_local(232) 
  reac_source_local(21,232) = + reac_rate_local(232) 
  reac_source_local(27,232) = - reac_rate_local(232) 
  reac_source_local(06,233) = + reac_rate_local(233) 
  reac_source_local(08,233) = - reac_rate_local(233) 
  reac_source_local(21,233) = - reac_rate_local(233) 
  reac_source_local(27,233) = + reac_rate_local(233) 
  reac_source_local(06,234) = - reac_rate_local(234) 
  reac_source_local(08,234) = + reac_rate_local(234) 
  reac_source_local(21,234) = + reac_rate_local(234) 
  reac_source_local(27,234) = - reac_rate_local(234) 
  reac_source_local(07,235) = + reac_rate_local(235) 
  reac_source_local(09,235) = - reac_rate_local(235) 
  reac_source_local(21,235) = - reac_rate_local(235) 
  reac_source_local(27,235) = + reac_rate_local(235) 
  reac_source_local(07,236) = - reac_rate_local(236) 
  reac_source_local(09,236) = + reac_rate_local(236) 
  reac_source_local(21,236) = + reac_rate_local(236) 
  reac_source_local(27,236) = - reac_rate_local(236) 
  reac_source_local(01,237) = + reac_rate_local(237) 
  reac_source_local(03,237) = - reac_rate_local(237) 
  reac_source_local(27,237) = - reac_rate_local(237) 
  reac_source_local(28,237) = + reac_rate_local(237) 
  reac_source_local(01,238) = - reac_rate_local(238) 
  reac_source_local(03,238) = + reac_rate_local(238) 
  reac_source_local(27,238) = + reac_rate_local(238) 
  reac_source_local(28,238) = - reac_rate_local(238) 
  reac_source_local(02,239) = + reac_rate_local(239) 
  reac_source_local(04,239) = - reac_rate_local(239) 
  reac_source_local(27,239) = - reac_rate_local(239) 
  reac_source_local(28,239) = + reac_rate_local(239) 
  reac_source_local(02,240) = - reac_rate_local(240) 
  reac_source_local(04,240) = + reac_rate_local(240) 
  reac_source_local(27,240) = + reac_rate_local(240) 
  reac_source_local(28,240) = - reac_rate_local(240) 
  reac_source_local(03,241) = + reac_rate_local(241) 
  reac_source_local(05,241) = - reac_rate_local(241) 
  reac_source_local(27,241) = - reac_rate_local(241) 
  reac_source_local(28,241) = + reac_rate_local(241) 
  reac_source_local(03,242) = - reac_rate_local(242) 
  reac_source_local(05,242) = + reac_rate_local(242) 
  reac_source_local(27,242) = + reac_rate_local(242) 
  reac_source_local(28,242) = - reac_rate_local(242) 
  reac_source_local(04,243) = + reac_rate_local(243) 
  reac_source_local(06,243) = - reac_rate_local(243) 
  reac_source_local(27,243) = - reac_rate_local(243) 
  reac_source_local(28,243) = + reac_rate_local(243) 
  reac_source_local(04,244) = - reac_rate_local(244) 
  reac_source_local(06,244) = + reac_rate_local(244) 
  reac_source_local(27,244) = + reac_rate_local(244) 
  reac_source_local(28,244) = - reac_rate_local(244) 
  reac_source_local(05,245) = + reac_rate_local(245) 
  reac_source_local(07,245) = - reac_rate_local(245) 
  reac_source_local(27,245) = - reac_rate_local(245) 
  reac_source_local(28,245) = + reac_rate_local(245) 
  reac_source_local(05,246) = - reac_rate_local(246) 
  reac_source_local(07,246) = + reac_rate_local(246) 
  reac_source_local(27,246) = + reac_rate_local(246) 
  reac_source_local(28,246) = - reac_rate_local(246) 
  reac_source_local(06,247) = + reac_rate_local(247) 
  reac_source_local(08,247) = - reac_rate_local(247) 
  reac_source_local(27,247) = - reac_rate_local(247) 
  reac_source_local(28,247) = + reac_rate_local(247) 
  reac_source_local(06,248) = - reac_rate_local(248) 
  reac_source_local(08,248) = + reac_rate_local(248) 
  reac_source_local(27,248) = + reac_rate_local(248) 
  reac_source_local(28,248) = - reac_rate_local(248) 
  reac_source_local(07,249) = + reac_rate_local(249) 
  reac_source_local(09,249) = - reac_rate_local(249) 
  reac_source_local(27,249) = - reac_rate_local(249) 
  reac_source_local(28,249) = + reac_rate_local(249) 
  reac_source_local(07,250) = - reac_rate_local(250) 
  reac_source_local(09,250) = + reac_rate_local(250) 
  reac_source_local(27,250) = + reac_rate_local(250) 
  reac_source_local(28,250) = - reac_rate_local(250) 
  reac_source_local(01,251) = + reac_rate_local(251) 
  reac_source_local(03,251) = - reac_rate_local(251) 
  reac_source_local(28,251) = - reac_rate_local(251) 
  reac_source_local(29,251) = + reac_rate_local(251) 
  reac_source_local(01,252) = - reac_rate_local(252) 
  reac_source_local(03,252) = + reac_rate_local(252) 
  reac_source_local(28,252) = + reac_rate_local(252) 
  reac_source_local(29,252) = - reac_rate_local(252) 
  reac_source_local(02,253) = + reac_rate_local(253) 
  reac_source_local(04,253) = - reac_rate_local(253) 
  reac_source_local(28,253) = - reac_rate_local(253) 
  reac_source_local(29,253) = + reac_rate_local(253) 
  reac_source_local(02,254) = - reac_rate_local(254) 
  reac_source_local(04,254) = + reac_rate_local(254) 
  reac_source_local(28,254) = + reac_rate_local(254) 
  reac_source_local(29,254) = - reac_rate_local(254) 
  reac_source_local(03,255) = + reac_rate_local(255) 
  reac_source_local(05,255) = - reac_rate_local(255) 
  reac_source_local(28,255) = - reac_rate_local(255) 
  reac_source_local(29,255) = + reac_rate_local(255) 
  reac_source_local(03,256) = - reac_rate_local(256) 
  reac_source_local(05,256) = + reac_rate_local(256) 
  reac_source_local(28,256) = + reac_rate_local(256) 
  reac_source_local(29,256) = - reac_rate_local(256) 
  reac_source_local(04,257) = + reac_rate_local(257) 
  reac_source_local(06,257) = - reac_rate_local(257) 
  reac_source_local(28,257) = - reac_rate_local(257) 
  reac_source_local(29,257) = + reac_rate_local(257) 
  reac_source_local(04,258) = - reac_rate_local(258) 
  reac_source_local(06,258) = + reac_rate_local(258) 
  reac_source_local(28,258) = + reac_rate_local(258) 
  reac_source_local(29,258) = - reac_rate_local(258) 
  reac_source_local(05,259) = + reac_rate_local(259) 
  reac_source_local(07,259) = - reac_rate_local(259) 
  reac_source_local(28,259) = - reac_rate_local(259) 
  reac_source_local(29,259) = + reac_rate_local(259) 
  reac_source_local(05,260) = - reac_rate_local(260) 
  reac_source_local(07,260) = + reac_rate_local(260) 
  reac_source_local(28,260) = + reac_rate_local(260) 
  reac_source_local(29,260) = - reac_rate_local(260) 
  reac_source_local(06,261) = + reac_rate_local(261) 
  reac_source_local(08,261) = - reac_rate_local(261) 
  reac_source_local(28,261) = - reac_rate_local(261) 
  reac_source_local(29,261) = + reac_rate_local(261) 
  reac_source_local(06,262) = - reac_rate_local(262) 
  reac_source_local(08,262) = + reac_rate_local(262) 
  reac_source_local(28,262) = + reac_rate_local(262) 
  reac_source_local(29,262) = - reac_rate_local(262) 
  reac_source_local(07,263) = + reac_rate_local(263) 
  reac_source_local(09,263) = - reac_rate_local(263) 
  reac_source_local(28,263) = - reac_rate_local(263) 
  reac_source_local(29,263) = + reac_rate_local(263) 
  reac_source_local(07,264) = - reac_rate_local(264) 
  reac_source_local(09,264) = + reac_rate_local(264) 
  reac_source_local(28,264) = + reac_rate_local(264) 
  reac_source_local(29,264) = - reac_rate_local(264) 
  reac_source_local(27,265) = + reac_rate_local(265) 
  reac_source_local(28,265) = - reac_rate_local(265) * 2.d0
  reac_source_local(29,265) = + reac_rate_local(265) 
  reac_source_local(27,266) = - reac_rate_local(266) 
  reac_source_local(28,266) = + reac_rate_local(266) * 2.d0
  reac_source_local(29,266) = - reac_rate_local(266) 
  reac_source_local(21,267) = + reac_rate_local(267) 
  reac_source_local(27,267) = - reac_rate_local(267) 
  reac_source_local(28,267) = - reac_rate_local(267) 
  reac_source_local(29,267) = + reac_rate_local(267) 
  reac_source_local(21,268) = - reac_rate_local(268) 
  reac_source_local(27,268) = + reac_rate_local(268) 
  reac_source_local(28,268) = + reac_rate_local(268) 
  reac_source_local(29,268) = - reac_rate_local(268) 
  reac_source_local(21,269) = + reac_rate_local(269) 
  reac_source_local(27,269) = - reac_rate_local(269) * 2.d0
  reac_source_local(28,269) = + reac_rate_local(269) 
  reac_source_local(21,270) = - reac_rate_local(270) 
  reac_source_local(27,270) = + reac_rate_local(270) * 2.d0
  reac_source_local(28,270) = - reac_rate_local(270) 
  reac_source_local(01,271) = + reac_rate_local(271) 
  reac_source_local(07,271) = - reac_rate_local(271) 
  reac_source_local(10,271) = - reac_rate_local(271) 
  reac_source_local(11,271) = + reac_rate_local(271) 
  reac_source_local(02,272) = + reac_rate_local(272) 
  reac_source_local(08,272) = - reac_rate_local(272) 
  reac_source_local(10,272) = - reac_rate_local(272) 
  reac_source_local(11,272) = + reac_rate_local(272) 
  reac_source_local(03,273) = + reac_rate_local(273) 
  reac_source_local(09,273) = - reac_rate_local(273) 
  reac_source_local(10,273) = - reac_rate_local(273) 
  reac_source_local(11,273) = + reac_rate_local(273) 
  reac_source_local(01,274) = - reac_rate_local(274) 
  reac_source_local(07,274) = + reac_rate_local(274) 
  reac_source_local(10,274) = + reac_rate_local(274) 
  reac_source_local(11,274) = - reac_rate_local(274) 
  reac_source_local(02,275) = - reac_rate_local(275) 
  reac_source_local(08,275) = + reac_rate_local(275) 
  reac_source_local(10,275) = + reac_rate_local(275) 
  reac_source_local(11,275) = - reac_rate_local(275) 
  reac_source_local(03,276) = - reac_rate_local(276) 
  reac_source_local(09,276) = + reac_rate_local(276) 
  reac_source_local(10,276) = + reac_rate_local(276) 
  reac_source_local(11,276) = - reac_rate_local(276) 
  reac_source_local(04,277) = - reac_rate_local(277) 
  reac_source_local(09,277) = + reac_rate_local(277) 
  reac_source_local(10,277) = + reac_rate_local(277) 
  reac_source_local(11,277) = - reac_rate_local(277) 
  reac_source_local(05,278) = - reac_rate_local(278) 
  reac_source_local(09,278) = + reac_rate_local(278) 
  reac_source_local(10,278) = + reac_rate_local(278) 
  reac_source_local(11,278) = - reac_rate_local(278) 
  reac_source_local(06,279) = - reac_rate_local(279) 
  reac_source_local(09,279) = + reac_rate_local(279) 
  reac_source_local(10,279) = + reac_rate_local(279) 
  reac_source_local(11,279) = - reac_rate_local(279) 
  reac_source_local(07,280) = - reac_rate_local(280) 
  reac_source_local(09,280) = + reac_rate_local(280) 
  reac_source_local(10,280) = + reac_rate_local(280) 
  reac_source_local(11,280) = - reac_rate_local(280) 
  reac_source_local(08,281) = - reac_rate_local(281) 
  reac_source_local(09,281) = + reac_rate_local(281) 
  reac_source_local(10,281) = + reac_rate_local(281) 
  reac_source_local(11,281) = - reac_rate_local(281) 
  reac_source_local(10,282) = + reac_rate_local(282) 
  reac_source_local(11,282) = - reac_rate_local(282) 
  reac_source_local(10,283) = + reac_rate_local(283) 
  reac_source_local(11,283) = - reac_rate_local(283) 
  reac_source_local(01,284) = + reac_rate_local(284) 
  reac_source_local(12,284) = - reac_rate_local(284) 
  reac_source_local(11,285) = + reac_rate_local(285) 
  reac_source_local(13,285) = - reac_rate_local(285) 
  reac_source_local(21,286) = - reac_rate_local(286) 
  reac_source_local(30,286) = + reac_rate_local(286) * 2.d0
  reac_source_local(01,287) = - reac_rate_local(287) 
  reac_source_local(14,287) = + reac_rate_local(287) * 2.d0
  reac_source_local(14,288) = + reac_rate_local(288) 
  reac_source_local(30,288) = + reac_rate_local(288) 
  reac_source_local(35,288) = - reac_rate_local(288) 
  reac_source_local(14,289) = + reac_rate_local(289) 
  reac_source_local(21,289) = + reac_rate_local(289) 
  reac_source_local(36,289) = - reac_rate_local(289) 
  reac_source_local(30,290) = + reac_rate_local(290) 
  reac_source_local(35,290) = + reac_rate_local(290) 
  reac_source_local(36,290) = - reac_rate_local(290) 
  reac_source_local(30,291) = + reac_rate_local(291) 
  reac_source_local(36,291) = + reac_rate_local(291) 
  reac_source_local(37,291) = - reac_rate_local(291) 
  reac_source_local(21,292) = + reac_rate_local(292) 
  reac_source_local(35,292) = + reac_rate_local(292) 
  reac_source_local(37,292) = - reac_rate_local(292) 
  reac_source_local(14,293) = + reac_rate_local(293) * 2.d0
  reac_source_local(18,293) = - reac_rate_local(293) 
  reac_source_local(43,293) = - reac_rate_local(293) 
  reac_source_local(14,294) = + reac_rate_local(294) 
  reac_source_local(15,294) = + reac_rate_local(294) 
  reac_source_local(18,294) = - reac_rate_local(294) 
  reac_source_local(43,294) = - reac_rate_local(294) 
  reac_source_local(14,295) = + reac_rate_local(295) 
  reac_source_local(16,295) = + reac_rate_local(295) 
  reac_source_local(18,295) = - reac_rate_local(295) 
  reac_source_local(43,295) = - reac_rate_local(295) 
  reac_source_local(01,296) = + reac_rate_local(296) 
  reac_source_local(14,296) = + reac_rate_local(296) 
  reac_source_local(19,296) = - reac_rate_local(296) 
  reac_source_local(43,296) = - reac_rate_local(296) 
  reac_source_local(01,297) = + reac_rate_local(297) * 2.d0
  reac_source_local(20,297) = - reac_rate_local(297) 
  reac_source_local(43,297) = - reac_rate_local(297) 
  reac_source_local(21,298) = - reac_rate_local(298) 
  reac_source_local(30,298) = + reac_rate_local(298) 
  reac_source_local(33,298) = + reac_rate_local(298) 
  reac_source_local(43,298) = + reac_rate_local(298) 
  reac_source_local(30,299) = + reac_rate_local(299) * 2.d0
  reac_source_local(31,299) = - reac_rate_local(299) 
  reac_source_local(43,299) = - reac_rate_local(299) 
  reac_source_local(30,300) = + reac_rate_local(300) * 3.d0
  reac_source_local(32,300) = - reac_rate_local(300) 
  reac_source_local(43,300) = - reac_rate_local(300) 
  reac_source_local(21,301) = + reac_rate_local(301) 
  reac_source_local(30,301) = + reac_rate_local(301) 
  reac_source_local(32,301) = - reac_rate_local(301) 
  reac_source_local(43,301) = - reac_rate_local(301) 
  reac_source_local(14,302) = + reac_rate_local(302) 
  reac_source_local(30,302) = + reac_rate_local(302) 
  reac_source_local(38,302) = - reac_rate_local(302) 
  reac_source_local(43,302) = - reac_rate_local(302) 
  reac_source_local(30,303) = + reac_rate_local(303) 
  reac_source_local(35,303) = + reac_rate_local(303) 
  reac_source_local(40,303) = - reac_rate_local(303) 
  reac_source_local(43,303) = - reac_rate_local(303) 
  reac_source_local(14,304) = + reac_rate_local(304) 
  reac_source_local(30,304) = + reac_rate_local(304) * 2.d0
  reac_source_local(40,304) = - reac_rate_local(304) 
  reac_source_local(43,304) = - reac_rate_local(304) 
  reac_source_local(30,305) = + reac_rate_local(305) * 2.d0
  reac_source_local(35,305) = + reac_rate_local(305) 
  reac_source_local(41,305) = - reac_rate_local(305) 
  reac_source_local(43,305) = - reac_rate_local(305) 
  reac_source_local(30,306) = + reac_rate_local(306) 
  reac_source_local(36,306) = + reac_rate_local(306) 
  reac_source_local(41,306) = - reac_rate_local(306) 
  reac_source_local(43,306) = - reac_rate_local(306) 
  reac_source_local(30,307) = + reac_rate_local(307) 
  reac_source_local(37,307) = + reac_rate_local(307) 
  reac_source_local(42,307) = - reac_rate_local(307) 
  reac_source_local(43,307) = - reac_rate_local(307) 
  reac_source_local(30,308) = + reac_rate_local(308) * 2.d0
  reac_source_local(36,308) = + reac_rate_local(308) 
  reac_source_local(42,308) = - reac_rate_local(308) 
  reac_source_local(43,308) = - reac_rate_local(308) 
  reac_source_local(01,309) = + reac_rate_local(309) 
  reac_source_local(30,309) = + reac_rate_local(309) 
  reac_source_local(39,309) = - reac_rate_local(309) 
  reac_source_local(43,309) = - reac_rate_local(309) 
  reac_source_local(18,310) = - reac_rate_local(310) 
  reac_source_local(21,310) = - reac_rate_local(310) 
  reac_source_local(30,310) = + reac_rate_local(310) 
  reac_source_local(39,310) = + reac_rate_local(310) 
  reac_source_local(10,311) = - reac_rate_local(311) 
  reac_source_local(14,311) = + reac_rate_local(311) 
  reac_source_local(18,311) = - reac_rate_local(311) 
  reac_source_local(19,311) = + reac_rate_local(311) 
  reac_source_local(01,312) = + reac_rate_local(312) 
  reac_source_local(14,312) = - reac_rate_local(312) 
  reac_source_local(17,312) = + reac_rate_local(312) 
  reac_source_local(18,312) = - reac_rate_local(312) 
  reac_source_local(14,313) = - reac_rate_local(313) 
  reac_source_local(18,313) = - reac_rate_local(313) 
  reac_source_local(19,313) = + reac_rate_local(313) 
  reac_source_local(01,314) = - reac_rate_local(314) 
  reac_source_local(18,314) = - reac_rate_local(314) 
  reac_source_local(20,314) = + reac_rate_local(314) 
  reac_source_local(01,315) = + reac_rate_local(315) 
  reac_source_local(18,315) = - reac_rate_local(315) 
  reac_source_local(37,315) = - reac_rate_local(315) 
  reac_source_local(41,315) = + reac_rate_local(315) 
  reac_source_local(01,316) = + reac_rate_local(316) 
  reac_source_local(14,316) = - reac_rate_local(316) 
  reac_source_local(18,316) = + reac_rate_local(316) 
  reac_source_local(19,316) = - reac_rate_local(316) 
  reac_source_local(01,317) = + reac_rate_local(317) 
  reac_source_local(18,317) = + reac_rate_local(317) 
  reac_source_local(20,317) = - reac_rate_local(317) 
  reac_source_local(01,318) = + reac_rate_local(318) * 2.d0
  reac_source_local(14,318) = - reac_rate_local(318) 
  reac_source_local(17,318) = + reac_rate_local(318) 
  reac_source_local(20,318) = - reac_rate_local(318) 
  reac_source_local(17,319) = - reac_rate_local(319) 
  reac_source_local(21,319) = - reac_rate_local(319) 
  reac_source_local(30,319) = + reac_rate_local(319) 
  reac_source_local(38,319) = + reac_rate_local(319) 
  reac_source_local(17,320) = - reac_rate_local(320) 
  reac_source_local(35,320) = + reac_rate_local(320) 
  reac_source_local(37,320) = - reac_rate_local(320) 
  reac_source_local(40,320) = + reac_rate_local(320) 
  reac_source_local(14,321) = + reac_rate_local(321) 
  reac_source_local(17,321) = - reac_rate_local(321) 
  reac_source_local(37,321) = - reac_rate_local(321) 
  reac_source_local(41,321) = + reac_rate_local(321) 
  reac_source_local(17,322) = - reac_rate_local(322) 
  reac_source_local(21,322) = + reac_rate_local(322) 
  reac_source_local(37,322) = - reac_rate_local(322) 
  reac_source_local(39,322) = + reac_rate_local(322) 
  reac_source_local(21,323) = + reac_rate_local(323) 
  reac_source_local(30,323) = - reac_rate_local(323) 
  reac_source_local(31,323) = - reac_rate_local(323) 
  reac_source_local(33,323) = + reac_rate_local(323) 
  reac_source_local(21,324) = - reac_rate_local(324) 
  reac_source_local(30,324) = + reac_rate_local(324) 
  reac_source_local(31,324) = - reac_rate_local(324) 
  reac_source_local(32,324) = + reac_rate_local(324) 
  reac_source_local(21,325) = + reac_rate_local(325) 
  reac_source_local(31,325) = - reac_rate_local(325) 
  reac_source_local(37,325) = - reac_rate_local(325) 
  reac_source_local(41,325) = + reac_rate_local(325) 
  reac_source_local(01,326) = - reac_rate_local(326) 
  reac_source_local(30,326) = + reac_rate_local(326) 
  reac_source_local(31,326) = - reac_rate_local(326) 
  reac_source_local(39,326) = + reac_rate_local(326) 
  reac_source_local(30,327) = + reac_rate_local(327) 
  reac_source_local(33,327) = - reac_rate_local(327) 
  reac_source_local(37,327) = - reac_rate_local(327) 
  reac_source_local(41,327) = + reac_rate_local(327) 
  reac_source_local(14,328) = + reac_rate_local(328) 
  reac_source_local(21,328) = - reac_rate_local(328) 
  reac_source_local(32,328) = + reac_rate_local(328) 
  reac_source_local(38,328) = - reac_rate_local(328) 
  reac_source_local(21,329) = - reac_rate_local(329) 
  reac_source_local(30,329) = + reac_rate_local(329) 
  reac_source_local(38,329) = - reac_rate_local(329) 
  reac_source_local(40,329) = + reac_rate_local(329) 
  reac_source_local(35,330) = + reac_rate_local(330) 
  reac_source_local(37,330) = - reac_rate_local(330) 
  reac_source_local(38,330) = - reac_rate_local(330) 
  reac_source_local(41,330) = + reac_rate_local(330) 
  reac_source_local(14,331) = + reac_rate_local(331) 
  reac_source_local(37,331) = - reac_rate_local(331) 
  reac_source_local(38,331) = - reac_rate_local(331) 
  reac_source_local(42,331) = + reac_rate_local(331) 
  reac_source_local(01,332) = - reac_rate_local(332) 
  reac_source_local(14,332) = + reac_rate_local(332) 
  reac_source_local(38,332) = - reac_rate_local(332) 
  reac_source_local(39,332) = + reac_rate_local(332) 
  reac_source_local(21,333) = - reac_rate_local(333) 
  reac_source_local(30,333) = + reac_rate_local(333) 
  reac_source_local(40,333) = - reac_rate_local(333) 
  reac_source_local(41,333) = + reac_rate_local(333) 
  reac_source_local(36,334) = + reac_rate_local(334) 
  reac_source_local(37,334) = - reac_rate_local(334) 
  reac_source_local(40,334) = - reac_rate_local(334) 
  reac_source_local(41,334) = + reac_rate_local(334) 
  reac_source_local(35,335) = + reac_rate_local(335) 
  reac_source_local(37,335) = - reac_rate_local(335) 
  reac_source_local(40,335) = - reac_rate_local(335) 
  reac_source_local(42,335) = + reac_rate_local(335) 
  reac_source_local(36,336) = + reac_rate_local(336) 
  reac_source_local(37,336) = - reac_rate_local(336) 
  reac_source_local(41,336) = - reac_rate_local(336) 
  reac_source_local(42,336) = + reac_rate_local(336) 
  reac_source_local(01,337) = + reac_rate_local(337) 
  reac_source_local(10,337) = - reac_rate_local(337) 
  reac_source_local(01,338) = + reac_rate_local(338) 
  reac_source_local(10,338) = - reac_rate_local(338) 
  reac_source_local(14,338) = - reac_rate_local(338) 
  reac_source_local(16,338) = + reac_rate_local(338) 
  reac_source_local(01,339) = + reac_rate_local(339) 
  reac_source_local(10,339) = - reac_rate_local(339) 
  reac_source_local(01,340) = + reac_rate_local(340) 
  reac_source_local(10,340) = - reac_rate_local(340) * 2.d0
  reac_source_local(11,340) = + reac_rate_local(340) 
  reac_source_local(01,341) = + reac_rate_local(341) 
  reac_source_local(10,341) = - reac_rate_local(341) * 2.d0
  reac_source_local(13,341) = + reac_rate_local(341) 
  reac_source_local(01,342) = + reac_rate_local(342) 
  reac_source_local(10,342) = - reac_rate_local(342) * 2.d0
  reac_source_local(14,342) = + reac_rate_local(342) * 2.d0
  reac_source_local(10,343) = + reac_rate_local(343) 
  reac_source_local(11,343) = - reac_rate_local(343) 
  reac_source_local(01,344) = + reac_rate_local(344) 
  reac_source_local(11,344) = - reac_rate_local(344) 
  reac_source_local(11,345) = + reac_rate_local(345) 
  reac_source_local(12,345) = - reac_rate_local(345) 
  reac_source_local(01,346) = + reac_rate_local(346) 
  reac_source_local(12,346) = - reac_rate_local(346) * 2.d0
  reac_source_local(18,346) = + reac_rate_local(346) 
  reac_source_local(43,346) = + reac_rate_local(346) 
  reac_source_local(12,347) = - reac_rate_local(347) * 2.d0
  reac_source_local(20,347) = + reac_rate_local(347) 
  reac_source_local(43,347) = + reac_rate_local(347) 
  reac_source_local(10,348) = - reac_rate_local(348) 
  reac_source_local(12,348) = - reac_rate_local(348) 
  reac_source_local(20,348) = + reac_rate_local(348) 
  reac_source_local(43,348) = + reac_rate_local(348) 
  reac_source_local(14,349) = + reac_rate_local(349) 
  reac_source_local(15,349) = - reac_rate_local(349) 
  reac_source_local(14,350) = + reac_rate_local(350) 
  reac_source_local(16,350) = - reac_rate_local(350) 
  reac_source_local(15,351) = + reac_rate_local(351) 
  reac_source_local(16,351) = - reac_rate_local(351) 
  reac_source_local(14,352) = + reac_rate_local(352) 
  reac_source_local(16,352) = - reac_rate_local(352) 
  reac_source_local(01,353) = + reac_rate_local(353) 
  reac_source_local(10,353) = - reac_rate_local(353) 
  reac_source_local(01,354) = + reac_rate_local(354) 
  reac_source_local(10,354) = - reac_rate_local(354) 
  reac_source_local(21,354) = - reac_rate_local(354) 
  reac_source_local(30,354) = + reac_rate_local(354) * 2.d0
  reac_source_local(01,355) = + reac_rate_local(355) 
  reac_source_local(10,355) = - reac_rate_local(355) 
  reac_source_local(10,356) = + reac_rate_local(356) 
  reac_source_local(11,356) = - reac_rate_local(356) 
  reac_source_local(01,357) = + reac_rate_local(357) 
  reac_source_local(12,357) = - reac_rate_local(357) 
  reac_source_local(01,358) = + reac_rate_local(358) 
  reac_source_local(12,358) = - reac_rate_local(358) 
  reac_source_local(21,358) = - reac_rate_local(358) 
  reac_source_local(30,358) = + reac_rate_local(358) * 2.d0
  reac_source_local(14,359) = - reac_rate_local(359) 
  reac_source_local(27,359) = - reac_rate_local(359) 
  reac_source_local(30,359) = + reac_rate_local(359) 
  reac_source_local(35,359) = + reac_rate_local(359) 
  reac_source_local(14,360) = - reac_rate_local(360) 
  reac_source_local(28,360) = - reac_rate_local(360) 
  reac_source_local(30,360) = + reac_rate_local(360) 
  reac_source_local(35,360) = + reac_rate_local(360) 
  reac_source_local(14,361) = - reac_rate_local(361) 
  reac_source_local(29,361) = - reac_rate_local(361) 
  reac_source_local(30,361) = + reac_rate_local(361) 
  reac_source_local(35,361) = + reac_rate_local(361) 
  reac_source_local(14,362) = - reac_rate_local(362) 
  reac_source_local(26,362) = - reac_rate_local(362) 
  reac_source_local(30,362) = + reac_rate_local(362) 
  reac_source_local(35,362) = + reac_rate_local(362) 
  reac_source_local(14,363) = - reac_rate_local(363) 
  reac_source_local(22,363) = - reac_rate_local(363) 
  reac_source_local(30,363) = + reac_rate_local(363) 
  reac_source_local(35,363) = + reac_rate_local(363) 
  reac_source_local(14,364) = - reac_rate_local(364) 
  reac_source_local(23,364) = - reac_rate_local(364) 
  reac_source_local(30,364) = + reac_rate_local(364) 
  reac_source_local(35,364) = + reac_rate_local(364) 
  reac_source_local(14,365) = - reac_rate_local(365) 
  reac_source_local(24,365) = - reac_rate_local(365) 
  reac_source_local(30,365) = + reac_rate_local(365) 
  reac_source_local(35,365) = + reac_rate_local(365) 
  reac_source_local(14,366) = - reac_rate_local(366) 
  reac_source_local(25,366) = - reac_rate_local(366) 
  reac_source_local(30,366) = + reac_rate_local(366) 
  reac_source_local(35,366) = + reac_rate_local(366) 
  reac_source_local(15,367) = - reac_rate_local(367) 
  reac_source_local(21,367) = - reac_rate_local(367) 
  reac_source_local(30,367) = + reac_rate_local(367) 
  reac_source_local(35,367) = + reac_rate_local(367) 
  reac_source_local(15,368) = - reac_rate_local(368) 
  reac_source_local(35,368) = + reac_rate_local(368) 
  reac_source_local(36,368) = + reac_rate_local(368) 
  reac_source_local(37,368) = - reac_rate_local(368) 
  reac_source_local(16,369) = - reac_rate_local(369) 
  reac_source_local(21,369) = - reac_rate_local(369) 
  reac_source_local(30,369) = + reac_rate_local(369) 
  reac_source_local(35,369) = + reac_rate_local(369) 
  reac_source_local(01,370) = + reac_rate_local(370) 
  reac_source_local(14,370) = - reac_rate_local(370) 
  reac_source_local(30,370) = + reac_rate_local(370) 
  reac_source_local(35,370) = - reac_rate_local(370) 
  reac_source_local(14,371) = + reac_rate_local(371) 
  reac_source_local(21,371) = + reac_rate_local(371) 
  reac_source_local(30,371) = - reac_rate_local(371) 
  reac_source_local(35,371) = - reac_rate_local(371) 
  reac_source_local(01,372) = + reac_rate_local(372) 
  reac_source_local(21,372) = + reac_rate_local(372) 
  reac_source_local(35,372) = - reac_rate_local(372) * 2.d0
  reac_source_local(14,373) = + reac_rate_local(373) 
  reac_source_local(35,373) = - reac_rate_local(373) * 2.d0
  reac_source_local(36,373) = + reac_rate_local(373) 
  reac_source_local(01,374) = + reac_rate_local(374) 
  reac_source_local(30,374) = + reac_rate_local(374) * 2.d0
  reac_source_local(35,374) = - reac_rate_local(374) * 2.d0
  reac_source_local(21,375) = + reac_rate_local(375) 
  reac_source_local(30,375) = - reac_rate_local(375) 
  reac_source_local(35,375) = + reac_rate_local(375) 
  reac_source_local(36,375) = - reac_rate_local(375) 
  reac_source_local(01,376) = + reac_rate_local(376) 
  reac_source_local(14,376) = - reac_rate_local(376) 
  reac_source_local(30,376) = + reac_rate_local(376) * 2.d0
  reac_source_local(36,376) = - reac_rate_local(376) 
  reac_source_local(01,377) = + reac_rate_local(377) 
  reac_source_local(14,377) = - reac_rate_local(377) 
  reac_source_local(21,377) = + reac_rate_local(377) 
  reac_source_local(36,377) = - reac_rate_local(377) 
  reac_source_local(14,378) = + reac_rate_local(378) 
  reac_source_local(35,378) = - reac_rate_local(378) 
  reac_source_local(36,378) = - reac_rate_local(378) 
  reac_source_local(37,378) = + reac_rate_local(378) 
  reac_source_local(21,379) = - reac_rate_local(379) 
  reac_source_local(30,379) = + reac_rate_local(379) 
  reac_source_local(36,379) = - reac_rate_local(379) 
  reac_source_local(37,379) = + reac_rate_local(379) 
  reac_source_local(21,380) = + reac_rate_local(380) 
  reac_source_local(30,380) = - reac_rate_local(380) 
  reac_source_local(36,380) = + reac_rate_local(380) 
  reac_source_local(37,380) = - reac_rate_local(380) 
  reac_source_local(10,381) = + reac_rate_local(381) 
  reac_source_local(14,381) = - reac_rate_local(381) * 2.d0
  reac_source_local(10,382) = + reac_rate_local(382) 
  reac_source_local(14,382) = - reac_rate_local(382) * 2.d0
  reac_source_local(10,383) = + reac_rate_local(383) 
  reac_source_local(14,383) = - reac_rate_local(383) * 2.d0
  reac_source_local(10,384) = + reac_rate_local(384) 
  reac_source_local(14,384) = - reac_rate_local(384) * 2.d0
  reac_source_local(10,385) = + reac_rate_local(385) 
  reac_source_local(14,385) = - reac_rate_local(385) * 2.d0
  reac_source_local(10,386) = + reac_rate_local(386) 
  reac_source_local(14,386) = - reac_rate_local(386) * 2.d0
  reac_source_local(10,387) = + reac_rate_local(387) 
  reac_source_local(14,387) = - reac_rate_local(387) * 2.d0
  reac_source_local(10,388) = + reac_rate_local(388) 
  reac_source_local(14,388) = - reac_rate_local(388) * 2.d0
  reac_source_local(10,389) = + reac_rate_local(389) 
  reac_source_local(14,389) = - reac_rate_local(389) * 2.d0
  reac_source_local(10,390) = + reac_rate_local(390) 
  reac_source_local(14,390) = - reac_rate_local(390) * 2.d0
  reac_source_local(10,391) = + reac_rate_local(391) 
  reac_source_local(14,391) = - reac_rate_local(391) * 2.d0
  reac_source_local(10,392) = + reac_rate_local(392) 
  reac_source_local(14,392) = - reac_rate_local(392) * 2.d0
  reac_source_local(10,393) = + reac_rate_local(393) 
  reac_source_local(14,393) = - reac_rate_local(393) * 2.d0
  reac_source_local(10,394) = + reac_rate_local(394) 
  reac_source_local(14,394) = - reac_rate_local(394) * 2.d0
  reac_source_local(10,395) = + reac_rate_local(395) 
  reac_source_local(14,395) = - reac_rate_local(395) * 2.d0
  reac_source_local(11,396) = + reac_rate_local(396) 
  reac_source_local(14,396) = - reac_rate_local(396) * 2.d0
  reac_source_local(11,397) = + reac_rate_local(397) 
  reac_source_local(14,397) = - reac_rate_local(397) * 2.d0
  reac_source_local(11,398) = + reac_rate_local(398) 
  reac_source_local(14,398) = - reac_rate_local(398) * 2.d0
  reac_source_local(11,399) = + reac_rate_local(399) 
  reac_source_local(14,399) = - reac_rate_local(399) * 2.d0
  reac_source_local(11,400) = + reac_rate_local(400) 
  reac_source_local(14,400) = - reac_rate_local(400) * 2.d0
  reac_source_local(11,401) = + reac_rate_local(401) 
  reac_source_local(14,401) = - reac_rate_local(401) * 2.d0
  reac_source_local(11,402) = + reac_rate_local(402) 
  reac_source_local(14,402) = - reac_rate_local(402) * 2.d0
  reac_source_local(11,403) = + reac_rate_local(403) 
  reac_source_local(14,403) = - reac_rate_local(403) * 2.d0
  reac_source_local(11,404) = + reac_rate_local(404) 
  reac_source_local(14,404) = - reac_rate_local(404) * 2.d0
  reac_source_local(11,405) = + reac_rate_local(405) 
  reac_source_local(14,405) = - reac_rate_local(405) * 2.d0
  reac_source_local(11,406) = + reac_rate_local(406) 
  reac_source_local(14,406) = - reac_rate_local(406) * 2.d0
  reac_source_local(11,407) = + reac_rate_local(407) 
  reac_source_local(14,407) = - reac_rate_local(407) * 2.d0
  reac_source_local(11,408) = + reac_rate_local(408) 
  reac_source_local(14,408) = - reac_rate_local(408) * 2.d0
  reac_source_local(11,409) = + reac_rate_local(409) 
  reac_source_local(14,409) = - reac_rate_local(409) * 2.d0
  reac_source_local(11,410) = + reac_rate_local(410) 
  reac_source_local(14,410) = - reac_rate_local(410) * 2.d0
  reac_source_local(01,411) = + reac_rate_local(411) 
  reac_source_local(14,411) = - reac_rate_local(411) * 2.d0
  reac_source_local(01,412) = + reac_rate_local(412) 
  reac_source_local(14,412) = - reac_rate_local(412) * 2.d0
  reac_source_local(21,413) = + reac_rate_local(413) 
  reac_source_local(30,413) = - reac_rate_local(413) * 2.d0
  reac_source_local(21,414) = + reac_rate_local(414) 
  reac_source_local(30,414) = - reac_rate_local(414) * 2.d0
  reac_source_local(14,415) = - reac_rate_local(415) 
  reac_source_local(30,415) = - reac_rate_local(415) 
  reac_source_local(35,415) = + reac_rate_local(415) 
  reac_source_local(14,416) = - reac_rate_local(416) 
  reac_source_local(30,416) = - reac_rate_local(416) 
  reac_source_local(35,416) = + reac_rate_local(416) 
  reac_source_local(14,417) = - reac_rate_local(417) 
  reac_source_local(21,417) = - reac_rate_local(417) 
  reac_source_local(36,417) = + reac_rate_local(417) 
  reac_source_local(14,418) = - reac_rate_local(418) 
  reac_source_local(21,418) = - reac_rate_local(418) 
  reac_source_local(36,418) = + reac_rate_local(418) 
  reac_source_local(30,419) = - reac_rate_local(419) 
  reac_source_local(35,419) = - reac_rate_local(419) 
  reac_source_local(36,419) = + reac_rate_local(419) 
  reac_source_local(30,420) = - reac_rate_local(420) 
  reac_source_local(35,420) = - reac_rate_local(420) 
  reac_source_local(36,420) = + reac_rate_local(420) 
  reac_source_local(30,421) = - reac_rate_local(421) 
  reac_source_local(36,421) = - reac_rate_local(421) 
  reac_source_local(37,421) = + reac_rate_local(421) 
  reac_source_local(30,422) = - reac_rate_local(422) 
  reac_source_local(36,422) = - reac_rate_local(422) 
  reac_source_local(37,422) = + reac_rate_local(422) 
  reac_source_local(21,423) = - reac_rate_local(423) 
  reac_source_local(35,423) = - reac_rate_local(423) 
  reac_source_local(37,423) = + reac_rate_local(423) 
  reac_source_local(21,424) = - reac_rate_local(424) 
  reac_source_local(35,424) = - reac_rate_local(424) 
  reac_source_local(37,424) = + reac_rate_local(424) 
  reac_source_local(01,425) = + reac_rate_local(425) 
  reac_source_local(10,425) = - reac_rate_local(425) 
  reac_source_local(11,426) = + reac_rate_local(426) 
  reac_source_local(12,426) = - reac_rate_local(426) 
  reac_source_local(01,427) = + reac_rate_local(427) 
  reac_source_local(02,427) = - reac_rate_local(427) 
  reac_source_local(02,428) = + reac_rate_local(428) 
  reac_source_local(03,428) = - reac_rate_local(428) 
  reac_source_local(03,429) = + reac_rate_local(429) 
  reac_source_local(04,429) = - reac_rate_local(429) 
  reac_source_local(04,430) = + reac_rate_local(430) 
  reac_source_local(05,430) = - reac_rate_local(430) 
  reac_source_local(05,431) = + reac_rate_local(431) 
  reac_source_local(06,431) = - reac_rate_local(431) 
  reac_source_local(06,432) = + reac_rate_local(432) 
  reac_source_local(07,432) = - reac_rate_local(432) 
  reac_source_local(07,433) = + reac_rate_local(433) 
  reac_source_local(08,433) = - reac_rate_local(433) 
  reac_source_local(08,434) = + reac_rate_local(434) 
  reac_source_local(09,434) = - reac_rate_local(434) 
  reac_source_local(21,435) = + reac_rate_local(435) 
  reac_source_local(22,435) = - reac_rate_local(435) 
  reac_source_local(21,436) = + reac_rate_local(436) 
  reac_source_local(23,436) = - reac_rate_local(436) 
  reac_source_local(21,437) = + reac_rate_local(437) 
  reac_source_local(24,437) = - reac_rate_local(437) 
  reac_source_local(21,438) = + reac_rate_local(438) 
  reac_source_local(25,438) = - reac_rate_local(438) 
  reac_source_local(21,439) = + reac_rate_local(439) 
  reac_source_local(27,439) = - reac_rate_local(439) 
  reac_source_local(27,440) = + reac_rate_local(440) 
  reac_source_local(28,440) = - reac_rate_local(440) 
  reac_source_local(28,441) = + reac_rate_local(441) 
  reac_source_local(29,441) = - reac_rate_local(441) 
  reac_source_local(26,442) = - reac_rate_local(442) 
  reac_source_local(30,442) = + reac_rate_local(442) 
  reac_source_local(34,442) = + reac_rate_local(442) 
  reac_source_local(43,442) = - reac_rate_local(442) 
  reac_source_local(30,443) = + reac_rate_local(443) * 3.d0
  reac_source_local(31,443) = - reac_rate_local(443) 
  reac_source_local(34,443) = - reac_rate_local(443) 
  reac_source_local(21,444) = + reac_rate_local(444) 
  reac_source_local(30,444) = + reac_rate_local(444) * 2.d0
  reac_source_local(32,444) = - reac_rate_local(444) 
  reac_source_local(34,444) = - reac_rate_local(444) 
  reac_source_local(01,445) = + reac_rate_local(445) 
  reac_source_local(18,445) = - reac_rate_local(445) 
  reac_source_local(30,445) = + reac_rate_local(445) 
  reac_source_local(34,445) = - reac_rate_local(445) 
  reac_source_local(01,446) = + reac_rate_local(446) * 2.d0
  reac_source_local(20,446) = - reac_rate_local(446) 
  reac_source_local(30,446) = + reac_rate_local(446) 
  reac_source_local(34,446) = - reac_rate_local(446) 
  reac_source_local(01,447) = + reac_rate_local(447) 
  reac_source_local(21,447) = + reac_rate_local(447) 
  reac_source_local(34,447) = - reac_rate_local(447) 
  reac_source_local(39,447) = - reac_rate_local(447) 
  reac_source_local(21,448) = + reac_rate_local(448) 
  reac_source_local(30,448) = + reac_rate_local(448) 
  reac_source_local(31,448) = - reac_rate_local(448) 
  reac_source_local(34,448) = - reac_rate_local(448) 
  reac_source_local(21,449) = + reac_rate_local(449) 
  reac_source_local(30,449) = + reac_rate_local(449) 
  reac_source_local(31,449) = - reac_rate_local(449) 
  reac_source_local(34,449) = - reac_rate_local(449) 
  reac_source_local(21,450) = + reac_rate_local(450) 
  reac_source_local(30,450) = + reac_rate_local(450) 
  reac_source_local(31,450) = - reac_rate_local(450) 
  reac_source_local(34,450) = - reac_rate_local(450) 
  reac_source_local(21,451) = + reac_rate_local(451) 
  reac_source_local(30,451) = + reac_rate_local(451) 
  reac_source_local(31,451) = - reac_rate_local(451) 
  reac_source_local(34,451) = - reac_rate_local(451) 
  reac_source_local(21,452) = + reac_rate_local(452) * 2.d0
  reac_source_local(32,452) = - reac_rate_local(452) 
  reac_source_local(34,452) = - reac_rate_local(452) 
  reac_source_local(21,453) = + reac_rate_local(453) * 2.d0
  reac_source_local(32,453) = - reac_rate_local(453) 
  reac_source_local(34,453) = - reac_rate_local(453) 
  reac_source_local(21,454) = + reac_rate_local(454) * 2.d0
  reac_source_local(32,454) = - reac_rate_local(454) 
  reac_source_local(34,454) = - reac_rate_local(454) 
  reac_source_local(21,455) = + reac_rate_local(455) * 2.d0
  reac_source_local(32,455) = - reac_rate_local(455) 
  reac_source_local(34,455) = - reac_rate_local(455) 
  reac_source_local(01,456) = + reac_rate_local(456) 
  reac_source_local(18,456) = - reac_rate_local(456) 
  reac_source_local(30,456) = + reac_rate_local(456) 
  reac_source_local(34,456) = - reac_rate_local(456) 
  reac_source_local(01,457) = + reac_rate_local(457) 
  reac_source_local(18,457) = - reac_rate_local(457) 
  reac_source_local(30,457) = + reac_rate_local(457) 
  reac_source_local(34,457) = - reac_rate_local(457) 
  reac_source_local(01,458) = + reac_rate_local(458) 
  reac_source_local(18,458) = - reac_rate_local(458) 
  reac_source_local(30,458) = + reac_rate_local(458) 
  reac_source_local(34,458) = - reac_rate_local(458) 
  reac_source_local(01,459) = + reac_rate_local(459) 
  reac_source_local(18,459) = - reac_rate_local(459) 
  reac_source_local(30,459) = + reac_rate_local(459) 
  reac_source_local(34,459) = - reac_rate_local(459) 
  reac_source_local(01,460) = + reac_rate_local(460) * 2.d0
  reac_source_local(20,460) = - reac_rate_local(460) 
  reac_source_local(30,460) = + reac_rate_local(460) 
  reac_source_local(34,460) = - reac_rate_local(460) 
  reac_source_local(01,461) = + reac_rate_local(461) * 2.d0
  reac_source_local(20,461) = - reac_rate_local(461) 
  reac_source_local(30,461) = + reac_rate_local(461) 
  reac_source_local(34,461) = - reac_rate_local(461) 
  reac_source_local(01,462) = + reac_rate_local(462) * 2.d0
  reac_source_local(20,462) = - reac_rate_local(462) 
  reac_source_local(30,462) = + reac_rate_local(462) 
  reac_source_local(34,462) = - reac_rate_local(462) 
  reac_source_local(01,463) = + reac_rate_local(463) * 2.d0
  reac_source_local(20,463) = - reac_rate_local(463) 
  reac_source_local(30,463) = + reac_rate_local(463) 
  reac_source_local(34,463) = - reac_rate_local(463) 
  reac_source_local(01,464) = + reac_rate_local(464) 
  reac_source_local(21,464) = + reac_rate_local(464) 
  reac_source_local(34,464) = - reac_rate_local(464) 
  reac_source_local(39,464) = - reac_rate_local(464) 
  reac_source_local(01,465) = + reac_rate_local(465) 
  reac_source_local(21,465) = + reac_rate_local(465) 
  reac_source_local(34,465) = - reac_rate_local(465) 
  reac_source_local(39,465) = - reac_rate_local(465) 
  reac_source_local(01,466) = + reac_rate_local(466) 
  reac_source_local(21,466) = + reac_rate_local(466) 
  reac_source_local(34,466) = - reac_rate_local(466) 
  reac_source_local(39,466) = - reac_rate_local(466) 
  reac_source_local(01,467) = + reac_rate_local(467) 
  reac_source_local(21,467) = + reac_rate_local(467) 
  reac_source_local(34,467) = - reac_rate_local(467) 
  reac_source_local(39,467) = - reac_rate_local(467) 
  reac_source_local(14,468) = - reac_rate_local(468) 
  reac_source_local(44,468) = - reac_rate_local(468) 
  reac_source_local(46,468) = + reac_rate_local(468) 
  reac_source_local(30,469) = - reac_rate_local(469) 
  reac_source_local(44,469) = - reac_rate_local(469) 
  reac_source_local(45,469) = + reac_rate_local(469) 
  reac_source_local(35,470) = - reac_rate_local(470) 
  reac_source_local(44,470) = - reac_rate_local(470) 
  reac_source_local(47,470) = + reac_rate_local(470) 
  reac_source_local(36,471) = - reac_rate_local(471) 
  reac_source_local(44,471) = - reac_rate_local(471) 
  reac_source_local(48,471) = + reac_rate_local(471) 
  reac_source_local(01,472) = + reac_rate_local(472) 
  reac_source_local(14,472) = - reac_rate_local(472) 
  reac_source_local(44,472) = + reac_rate_local(472) 
  reac_source_local(46,472) = - reac_rate_local(472) 
  reac_source_local(21,473) = + reac_rate_local(473) 
  reac_source_local(30,473) = - reac_rate_local(473) 
  reac_source_local(44,473) = + reac_rate_local(473) 
  reac_source_local(45,473) = - reac_rate_local(473) 
  reac_source_local(14,474) = - reac_rate_local(474) 
  reac_source_local(45,474) = - reac_rate_local(474) 
  reac_source_local(47,474) = + reac_rate_local(474) 
  reac_source_local(35,475) = - reac_rate_local(475) 
  reac_source_local(45,475) = - reac_rate_local(475) 
  reac_source_local(48,475) = + reac_rate_local(475) 
  reac_source_local(36,476) = - reac_rate_local(476) 
  reac_source_local(37,476) = + reac_rate_local(476) 
  reac_source_local(44,476) = + reac_rate_local(476) 
  reac_source_local(45,476) = - reac_rate_local(476) 
  reac_source_local(30,477) = - reac_rate_local(477) 
  reac_source_local(46,477) = - reac_rate_local(477) 
  reac_source_local(47,477) = + reac_rate_local(477) 
  reac_source_local(30,478) = - reac_rate_local(478) 
  reac_source_local(47,478) = - reac_rate_local(478) 
  reac_source_local(48,478) = + reac_rate_local(478) 
  reac_source_local(30,479) = - reac_rate_local(479) 
  reac_source_local(37,479) = + reac_rate_local(479) 
  reac_source_local(44,479) = + reac_rate_local(479) 
  reac_source_local(48,479) = - reac_rate_local(479) 
  reac_source_local(21,480) = - reac_rate_local(480) 
  reac_source_local(37,480) = + reac_rate_local(480) 
  reac_source_local(44,480) = + reac_rate_local(480) 
  reac_source_local(47,480) = - reac_rate_local(480) 
  reac_source_local(44,481) = + reac_rate_local(481) 
  reac_source_local(45,481) = - reac_rate_local(481) 
  reac_source_local(46,481) = - reac_rate_local(481) 
  reac_source_local(47,481) = + reac_rate_local(481) 
  reac_source_local(44,482) = + reac_rate_local(482) 
  reac_source_local(45,482) = - reac_rate_local(482) 
  reac_source_local(47,482) = - reac_rate_local(482) 
  reac_source_local(48,482) = + reac_rate_local(482) 
  reac_source_local(37,483) = + reac_rate_local(483) 
  reac_source_local(44,483) = + reac_rate_local(483) * 2.d0
  reac_source_local(45,483) = - reac_rate_local(483) 
  reac_source_local(48,483) = - reac_rate_local(483) 
  reac_source_local(01,484) = - reac_rate_local(484) 
  reac_source_local(44,484) = - reac_rate_local(484) * 2.d0
  reac_source_local(46,484) = + reac_rate_local(484) * 2.d0
  reac_source_local(02,485) = - reac_rate_local(485) 
  reac_source_local(44,485) = - reac_rate_local(485) * 2.d0
  reac_source_local(46,485) = + reac_rate_local(485) * 2.d0
  reac_source_local(03,486) = - reac_rate_local(486) 
  reac_source_local(44,486) = - reac_rate_local(486) * 2.d0
  reac_source_local(46,486) = + reac_rate_local(486) * 2.d0
  reac_source_local(04,487) = - reac_rate_local(487) 
  reac_source_local(44,487) = - reac_rate_local(487) * 2.d0
  reac_source_local(46,487) = + reac_rate_local(487) * 2.d0
  reac_source_local(05,488) = - reac_rate_local(488) 
  reac_source_local(44,488) = - reac_rate_local(488) * 2.d0
  reac_source_local(46,488) = + reac_rate_local(488) * 2.d0
  reac_source_local(06,489) = - reac_rate_local(489) 
  reac_source_local(44,489) = - reac_rate_local(489) * 2.d0
  reac_source_local(46,489) = + reac_rate_local(489) * 2.d0
  reac_source_local(07,490) = - reac_rate_local(490) 
  reac_source_local(44,490) = - reac_rate_local(490) * 2.d0
  reac_source_local(46,490) = + reac_rate_local(490) * 2.d0
  reac_source_local(08,491) = - reac_rate_local(491) 
  reac_source_local(44,491) = - reac_rate_local(491) * 2.d0
  reac_source_local(46,491) = + reac_rate_local(491) * 2.d0
  reac_source_local(09,492) = - reac_rate_local(492) 
  reac_source_local(44,492) = - reac_rate_local(492) * 2.d0
  reac_source_local(46,492) = + reac_rate_local(492) * 2.d0
  reac_source_local(10,493) = - reac_rate_local(493) 
  reac_source_local(44,493) = - reac_rate_local(493) * 2.d0
  reac_source_local(46,493) = + reac_rate_local(493) * 2.d0
  reac_source_local(11,494) = - reac_rate_local(494) 
  reac_source_local(44,494) = - reac_rate_local(494) * 2.d0
  reac_source_local(46,494) = + reac_rate_local(494) * 2.d0
  reac_source_local(12,495) = - reac_rate_local(495) 
  reac_source_local(44,495) = - reac_rate_local(495) * 2.d0
  reac_source_local(46,495) = + reac_rate_local(495) * 2.d0
  reac_source_local(13,496) = - reac_rate_local(496) 
  reac_source_local(44,496) = - reac_rate_local(496) * 2.d0
  reac_source_local(46,496) = + reac_rate_local(496) * 2.d0
  reac_source_local(21,497) = - reac_rate_local(497) 
  reac_source_local(44,497) = - reac_rate_local(497) * 2.d0
  reac_source_local(45,497) = + reac_rate_local(497) * 2.d0
  reac_source_local(27,498) = - reac_rate_local(498) 
  reac_source_local(44,498) = - reac_rate_local(498) * 2.d0
  reac_source_local(45,498) = + reac_rate_local(498) * 2.d0
  reac_source_local(28,499) = - reac_rate_local(499) 
  reac_source_local(44,499) = - reac_rate_local(499) * 2.d0
  reac_source_local(45,499) = + reac_rate_local(499) * 2.d0
  reac_source_local(29,500) = - reac_rate_local(500) 
  reac_source_local(44,500) = - reac_rate_local(500) * 2.d0
  reac_source_local(45,500) = + reac_rate_local(500) * 2.d0
  reac_source_local(22,501) = - reac_rate_local(501) 
  reac_source_local(44,501) = - reac_rate_local(501) * 2.d0
  reac_source_local(45,501) = + reac_rate_local(501) * 2.d0
  reac_source_local(23,502) = - reac_rate_local(502) 
  reac_source_local(44,502) = - reac_rate_local(502) * 2.d0
  reac_source_local(45,502) = + reac_rate_local(502) * 2.d0
  reac_source_local(24,503) = - reac_rate_local(503) 
  reac_source_local(44,503) = - reac_rate_local(503) * 2.d0
  reac_source_local(45,503) = + reac_rate_local(503) * 2.d0
  reac_source_local(25,504) = - reac_rate_local(504) 
  reac_source_local(44,504) = - reac_rate_local(504) * 2.d0
  reac_source_local(45,504) = + reac_rate_local(504) * 2.d0
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(49)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(43) 
  rrt(002) = rrt(002) * density(01) * density(43) 
  rrt(003) = rrt(003) * density(01) * density(43) 
  rrt(004) = rrt(004) * density(01) * density(43) 
  rrt(005) = rrt(005) * density(01) * density(43) 
  rrt(006) = rrt(006) * density(01) * density(43) 
  rrt(007) = rrt(007) * density(21) * density(43) 
  rrt(008) = rrt(008) * density(21) * density(43) 
  rrt(009) = rrt(009) * density(21) * density(43) 
  rrt(010) = rrt(010) * density(21) * density(43) 
  rrt(011) = rrt(011) * density(21) * density(43) 
  rrt(012) = rrt(012) * density(01) * density(43) 
  rrt(013) = rrt(013) * density(21) * density(43) 
  rrt(014) = rrt(014) * density(14) * density(43) 
  rrt(015) = rrt(015) * density(30) * density(43) 
  rrt(016) = rrt(016) * density(37) * density(43) 
  rrt(017) = rrt(017) * density(01) * density(43) 
  rrt(018) = rrt(018) * density(01) * density(43) 
  rrt(019) = rrt(019) * density(01) * density(43) 
  rrt(020) = rrt(020) * density(01) * density(43) 
  rrt(021) = rrt(021) * density(01) * density(43) 
  rrt(022) = rrt(022) * density(01) * density(43) 
  rrt(023) = rrt(023) * density(01) * density(43) 
  rrt(024) = rrt(024) * density(01) * density(43) 
  rrt(025) = rrt(025) * density(02) * density(43) 
  rrt(026) = rrt(026) * density(03) * density(43) 
  rrt(027) = rrt(027) * density(04) * density(43) 
  rrt(028) = rrt(028) * density(05) * density(43) 
  rrt(029) = rrt(029) * density(06) * density(43) 
  rrt(030) = rrt(030) * density(07) * density(43) 
  rrt(031) = rrt(031) * density(08) * density(43) 
  rrt(032) = rrt(032) * density(09) * density(43) 
  rrt(033) = rrt(033) * density(21) * density(43) 
  rrt(034) = rrt(034) * density(21) * density(43) 
  rrt(035) = rrt(035) * density(21) * density(43) 
  rrt(036) = rrt(036) * density(27) * density(43) 
  rrt(037) = rrt(037) * density(28) * density(43) 
  rrt(038) = rrt(038) * density(29) * density(43) 
  rrt(039) = rrt(039) * density(22) * density(43) 
  rrt(040) = rrt(040) * density(23) * density(43) 
  rrt(041) = rrt(041) * density(24) * density(43) 
  rrt(042) = rrt(042) * density(25) * density(43) 
  rrt(043) = rrt(043) * density(01) * density(02) 
  rrt(044) = rrt(044) * density(01)**2 
  rrt(045) = rrt(045) * density(01) * density(03) 
  rrt(046) = rrt(046) * density(01) * density(02) 
  rrt(047) = rrt(047) * density(01) * density(04) 
  rrt(048) = rrt(048) * density(01) * density(03) 
  rrt(049) = rrt(049) * density(01) * density(05) 
  rrt(050) = rrt(050) * density(01) * density(04) 
  rrt(051) = rrt(051) * density(01) * density(06) 
  rrt(052) = rrt(052) * density(01) * density(05) 
  rrt(053) = rrt(053) * density(01) * density(07) 
  rrt(054) = rrt(054) * density(01) * density(06) 
  rrt(055) = rrt(055) * density(01) * density(08) 
  rrt(056) = rrt(056) * density(01) * density(07) 
  rrt(057) = rrt(057) * density(01) * density(09) 
  rrt(058) = rrt(058) * density(01) * density(08) 
  rrt(059) = rrt(059) * density(02) * density(21) 
  rrt(060) = rrt(060) * density(01) * density(21) 
  rrt(061) = rrt(061) * density(03) * density(21) 
  rrt(062) = rrt(062) * density(02) * density(21) 
  rrt(063) = rrt(063) * density(04) * density(21) 
  rrt(064) = rrt(064) * density(03) * density(21) 
  rrt(065) = rrt(065) * density(05) * density(21) 
  rrt(066) = rrt(066) * density(04) * density(21) 
  rrt(067) = rrt(067) * density(06) * density(21) 
  rrt(068) = rrt(068) * density(05) * density(21) 
  rrt(069) = rrt(069) * density(07) * density(21) 
  rrt(070) = rrt(070) * density(06) * density(21) 
  rrt(071) = rrt(071) * density(08) * density(21) 
  rrt(072) = rrt(072) * density(07) * density(21) 
  rrt(073) = rrt(073) * density(09) * density(21) 
  rrt(074) = rrt(074) * density(08) * density(21) 
  rrt(075) = rrt(075) * density(02) * density(14) 
  rrt(076) = rrt(076) * density(01) * density(14) 
  rrt(077) = rrt(077) * density(03) * density(14) 
  rrt(078) = rrt(078) * density(02) * density(14) 
  rrt(079) = rrt(079) * density(04) * density(14) 
  rrt(080) = rrt(080) * density(03) * density(14) 
  rrt(081) = rrt(081) * density(05) * density(14) 
  rrt(082) = rrt(082) * density(04) * density(14) 
  rrt(083) = rrt(083) * density(06) * density(14) 
  rrt(084) = rrt(084) * density(05) * density(14) 
  rrt(085) = rrt(085) * density(07) * density(14) 
  rrt(086) = rrt(086) * density(06) * density(14) 
  rrt(087) = rrt(087) * density(08) * density(14) 
  rrt(088) = rrt(088) * density(07) * density(14) 
  rrt(089) = rrt(089) * density(09) * density(14) 
  rrt(090) = rrt(090) * density(08) * density(14) 
  rrt(091) = rrt(091) * density(02) * density(30) 
  rrt(092) = rrt(092) * density(01) * density(30) 
  rrt(093) = rrt(093) * density(03) * density(30) 
  rrt(094) = rrt(094) * density(02) * density(30) 
  rrt(095) = rrt(095) * density(04) * density(30) 
  rrt(096) = rrt(096) * density(03) * density(30) 
  rrt(097) = rrt(097) * density(05) * density(30) 
  rrt(098) = rrt(098) * density(04) * density(30) 
  rrt(099) = rrt(099) * density(06) * density(30) 
  rrt(100) = rrt(100) * density(05) * density(30) 
  rrt(101) = rrt(101) * density(07) * density(30) 
  rrt(102) = rrt(102) * density(06) * density(30) 
  rrt(103) = rrt(103) * density(08) * density(30) 
  rrt(104) = rrt(104) * density(07) * density(30) 
  rrt(105) = rrt(105) * density(09) * density(30) 
  rrt(106) = rrt(106) * density(08) * density(30) 
  rrt(107) = rrt(107) * density(21) * density(27) 
  rrt(108) = rrt(108) * density(21)**2 
  rrt(109) = rrt(109) * density(21) * density(28) 
  rrt(110) = rrt(110) * density(21) * density(27) 
  rrt(111) = rrt(111) * density(21) * density(29) 
  rrt(112) = rrt(112) * density(21) * density(28) 
  rrt(113) = rrt(113) * density(27) * density(30) 
  rrt(114) = rrt(114) * density(21) * density(30) 
  rrt(115) = rrt(115) * density(28) * density(30) 
  rrt(116) = rrt(116) * density(27) * density(30) 
  rrt(117) = rrt(117) * density(29) * density(30) 
  rrt(118) = rrt(118) * density(28) * density(30) 
  rrt(119) = rrt(119) * density(02)**2 
  rrt(120) = rrt(120) * density(01) * density(03) 
  rrt(121) = rrt(121) * density(02) * density(03) 
  rrt(122) = rrt(122) * density(01) * density(04) 
  rrt(123) = rrt(123) * density(02) * density(04) 
  rrt(124) = rrt(124) * density(01) * density(05) 
  rrt(125) = rrt(125) * density(02) * density(05) 
  rrt(126) = rrt(126) * density(01) * density(06) 
  rrt(127) = rrt(127) * density(02) * density(06) 
  rrt(128) = rrt(128) * density(01) * density(07) 
  rrt(129) = rrt(129) * density(02) * density(07) 
  rrt(130) = rrt(130) * density(01) * density(08) 
  rrt(131) = rrt(131) * density(02) * density(08) 
  rrt(132) = rrt(132) * density(01) * density(09) 
  rrt(133) = rrt(133) * density(03)**2 
  rrt(134) = rrt(134) * density(02) * density(04) 
  rrt(135) = rrt(135) * density(03) * density(04) 
  rrt(136) = rrt(136) * density(02) * density(05) 
  rrt(137) = rrt(137) * density(03) * density(05) 
  rrt(138) = rrt(138) * density(02) * density(06) 
  rrt(139) = rrt(139) * density(03) * density(06) 
  rrt(140) = rrt(140) * density(02) * density(07) 
  rrt(141) = rrt(141) * density(03) * density(07) 
  rrt(142) = rrt(142) * density(02) * density(08) 
  rrt(143) = rrt(143) * density(03) * density(08) 
  rrt(144) = rrt(144) * density(02) * density(09) 
  rrt(145) = rrt(145) * density(04)**2 
  rrt(146) = rrt(146) * density(03) * density(05) 
  rrt(147) = rrt(147) * density(04) * density(05) 
  rrt(148) = rrt(148) * density(03) * density(06) 
  rrt(149) = rrt(149) * density(04) * density(06) 
  rrt(150) = rrt(150) * density(03) * density(07) 
  rrt(151) = rrt(151) * density(04) * density(07) 
  rrt(152) = rrt(152) * density(03) * density(08) 
  rrt(153) = rrt(153) * density(04) * density(08) 
  rrt(154) = rrt(154) * density(03) * density(09) 
  rrt(155) = rrt(155) * density(05)**2 
  rrt(156) = rrt(156) * density(04) * density(06) 
  rrt(157) = rrt(157) * density(05) * density(06) 
  rrt(158) = rrt(158) * density(04) * density(07) 
  rrt(159) = rrt(159) * density(05) * density(07) 
  rrt(160) = rrt(160) * density(04) * density(08) 
  rrt(161) = rrt(161) * density(05) * density(08) 
  rrt(162) = rrt(162) * density(04) * density(09) 
  rrt(163) = rrt(163) * density(06)**2 
  rrt(164) = rrt(164) * density(05) * density(07) 
  rrt(165) = rrt(165) * density(06) * density(07) 
  rrt(166) = rrt(166) * density(05) * density(08) 
  rrt(167) = rrt(167) * density(06) * density(08) 
  rrt(168) = rrt(168) * density(05) * density(09) 
  rrt(169) = rrt(169) * density(07)**2 
  rrt(170) = rrt(170) * density(06) * density(08) 
  rrt(171) = rrt(171) * density(07) * density(08) 
  rrt(172) = rrt(172) * density(06) * density(09) 
  rrt(173) = rrt(173) * density(08)**2 
  rrt(174) = rrt(174) * density(07) * density(09) 
  rrt(175) = rrt(175) * density(01) * density(27) 
  rrt(176) = rrt(176) * density(02) * density(21) 
  rrt(177) = rrt(177) * density(02) * density(27) 
  rrt(178) = rrt(178) * density(03) * density(21) 
  rrt(179) = rrt(179) * density(03) * density(27) 
  rrt(180) = rrt(180) * density(04) * density(21) 
  rrt(181) = rrt(181) * density(04) * density(27) 
  rrt(182) = rrt(182) * density(05) * density(21) 
  rrt(183) = rrt(183) * density(05) * density(27) 
  rrt(184) = rrt(184) * density(06) * density(21) 
  rrt(185) = rrt(185) * density(06) * density(27) 
  rrt(186) = rrt(186) * density(07) * density(21) 
  rrt(187) = rrt(187) * density(07) * density(27) 
  rrt(188) = rrt(188) * density(08) * density(21) 
  rrt(189) = rrt(189) * density(08) * density(27) 
  rrt(190) = rrt(190) * density(09) * density(21) 
  rrt(191) = rrt(191) * density(01) * density(28) 
  rrt(192) = rrt(192) * density(02) * density(27) 
  rrt(193) = rrt(193) * density(02) * density(28) 
  rrt(194) = rrt(194) * density(03) * density(27) 
  rrt(195) = rrt(195) * density(03) * density(28) 
  rrt(196) = rrt(196) * density(04) * density(27) 
  rrt(197) = rrt(197) * density(04) * density(28) 
  rrt(198) = rrt(198) * density(05) * density(27) 
  rrt(199) = rrt(199) * density(05) * density(28) 
  rrt(200) = rrt(200) * density(06) * density(27) 
  rrt(201) = rrt(201) * density(06) * density(28) 
  rrt(202) = rrt(202) * density(07) * density(27) 
  rrt(203) = rrt(203) * density(07) * density(28) 
  rrt(204) = rrt(204) * density(08) * density(27) 
  rrt(205) = rrt(205) * density(08) * density(28) 
  rrt(206) = rrt(206) * density(09) * density(27) 
  rrt(207) = rrt(207) * density(01) * density(29) 
  rrt(208) = rrt(208) * density(02) * density(28) 
  rrt(209) = rrt(209) * density(02) * density(29) 
  rrt(210) = rrt(210) * density(03) * density(28) 
  rrt(211) = rrt(211) * density(03) * density(29) 
  rrt(212) = rrt(212) * density(04) * density(28) 
  rrt(213) = rrt(213) * density(04) * density(29) 
  rrt(214) = rrt(214) * density(05) * density(28) 
  rrt(215) = rrt(215) * density(05) * density(29) 
  rrt(216) = rrt(216) * density(06) * density(28) 
  rrt(217) = rrt(217) * density(06) * density(29) 
  rrt(218) = rrt(218) * density(07) * density(28) 
  rrt(219) = rrt(219) * density(07) * density(29) 
  rrt(220) = rrt(220) * density(08) * density(28) 
  rrt(221) = rrt(221) * density(08) * density(29) 
  rrt(222) = rrt(222) * density(09) * density(28) 
  rrt(223) = rrt(223) * density(03) * density(21) 
  rrt(224) = rrt(224) * density(01) * density(27) 
  rrt(225) = rrt(225) * density(04) * density(21) 
  rrt(226) = rrt(226) * density(02) * density(27) 
  rrt(227) = rrt(227) * density(05) * density(21) 
  rrt(228) = rrt(228) * density(03) * density(27) 
  rrt(229) = rrt(229) * density(06) * density(21) 
  rrt(230) = rrt(230) * density(04) * density(27) 
  rrt(231) = rrt(231) * density(07) * density(21) 
  rrt(232) = rrt(232) * density(05) * density(27) 
  rrt(233) = rrt(233) * density(08) * density(21) 
  rrt(234) = rrt(234) * density(06) * density(27) 
  rrt(235) = rrt(235) * density(09) * density(21) 
  rrt(236) = rrt(236) * density(07) * density(27) 
  rrt(237) = rrt(237) * density(03) * density(27) 
  rrt(238) = rrt(238) * density(01) * density(28) 
  rrt(239) = rrt(239) * density(04) * density(27) 
  rrt(240) = rrt(240) * density(02) * density(28) 
  rrt(241) = rrt(241) * density(05) * density(27) 
  rrt(242) = rrt(242) * density(03) * density(28) 
  rrt(243) = rrt(243) * density(06) * density(27) 
  rrt(244) = rrt(244) * density(04) * density(28) 
  rrt(245) = rrt(245) * density(07) * density(27) 
  rrt(246) = rrt(246) * density(05) * density(28) 
  rrt(247) = rrt(247) * density(08) * density(27) 
  rrt(248) = rrt(248) * density(06) * density(28) 
  rrt(249) = rrt(249) * density(09) * density(27) 
  rrt(250) = rrt(250) * density(07) * density(28) 
  rrt(251) = rrt(251) * density(03) * density(28) 
  rrt(252) = rrt(252) * density(01) * density(29) 
  rrt(253) = rrt(253) * density(04) * density(28) 
  rrt(254) = rrt(254) * density(02) * density(29) 
  rrt(255) = rrt(255) * density(05) * density(28) 
  rrt(256) = rrt(256) * density(03) * density(29) 
  rrt(257) = rrt(257) * density(06) * density(28) 
  rrt(258) = rrt(258) * density(04) * density(29) 
  rrt(259) = rrt(259) * density(07) * density(28) 
  rrt(260) = rrt(260) * density(05) * density(29) 
  rrt(261) = rrt(261) * density(08) * density(28) 
  rrt(262) = rrt(262) * density(06) * density(29) 
  rrt(263) = rrt(263) * density(09) * density(28) 
  rrt(264) = rrt(264) * density(07) * density(29) 
  rrt(265) = rrt(265) * density(28)**2 
  rrt(266) = rrt(266) * density(27) * density(29) 
  rrt(267) = rrt(267) * density(27) * density(28) 
  rrt(268) = rrt(268) * density(21) * density(29) 
  rrt(269) = rrt(269) * density(27)**2 
  rrt(270) = rrt(270) * density(21) * density(28) 
  rrt(271) = rrt(271) * density(07) * density(10) 
  rrt(272) = rrt(272) * density(08) * density(10) 
  rrt(273) = rrt(273) * density(09) * density(10) 
  rrt(274) = rrt(274) * density(01) * density(11) 
  rrt(275) = rrt(275) * density(02) * density(11) 
  rrt(276) = rrt(276) * density(03) * density(11) 
  rrt(277) = rrt(277) * density(04) * density(11) 
  rrt(278) = rrt(278) * density(05) * density(11) 
  rrt(279) = rrt(279) * density(06) * density(11) 
  rrt(280) = rrt(280) * density(07) * density(11) 
  rrt(281) = rrt(281) * density(08) * density(11) 
  rrt(282) = rrt(282) * density(09) * density(11) 
  rrt(283) = rrt(283) * density(11) 
  rrt(284) = rrt(284) * density(12) 
  rrt(285) = rrt(285) * density(13) 
  rrt(286) = rrt(286) * density(21) * density(43) 
  rrt(287) = rrt(287) * density(01) * density(43) 
  rrt(288) = rrt(288) * density(35) * density(43) 
  rrt(289) = rrt(289) * density(36) * density(43) 
  rrt(290) = rrt(290) * density(36) * density(43) 
  rrt(291) = rrt(291) * density(37) * density(43) 
  rrt(292) = rrt(292) * density(37) * density(43) 
  rrt(293) = rrt(293) * density(18) * density(43) 
  rrt(294) = rrt(294) * density(18) * density(43) 
  rrt(295) = rrt(295) * density(18) * density(43) 
  rrt(296) = rrt(296) * density(19) * density(43) 
  rrt(297) = rrt(297) * density(20) * density(43) 
  rrt(298) = rrt(298) * density(21) * density(43) 
  rrt(299) = rrt(299) * density(31) * density(43) 
  rrt(300) = rrt(300) * density(32) * density(43) 
  rrt(301) = rrt(301) * density(32) * density(43) 
  rrt(302) = rrt(302) * density(38) * density(43) 
  rrt(303) = rrt(303) * density(40) * density(43) 
  rrt(304) = rrt(304) * density(40) * density(43) 
  rrt(305) = rrt(305) * density(41) * density(43) 
  rrt(306) = rrt(306) * density(41) * density(43) 
  rrt(307) = rrt(307) * density(42) * density(43) 
  rrt(308) = rrt(308) * density(42) * density(43) 
  rrt(309) = rrt(309) * density(39) * density(43) 
  rrt(310) = rrt(310) * density(18) * density(21) 
  rrt(311) = rrt(311) * density(10) * density(18) 
  rrt(312) = rrt(312) * density(14) * density(18) 
  rrt(313) = rrt(313) * density(01) * density(14) * density(18) 
  rrt(314) = rrt(314) * density(01)**2 * density(18) 
  rrt(315) = rrt(315) * density(18) * density(37) 
  rrt(316) = rrt(316) * density(14) * density(19) 
  rrt(317) = rrt(317) * density(01) * density(20) 
  rrt(318) = rrt(318) * density(14) * density(20) 
  rrt(319) = rrt(319) * density(17) * density(21) 
  rrt(320) = rrt(320) * density(17) * density(37) 
  rrt(321) = rrt(321) * density(17) * density(37) 
  rrt(322) = rrt(322) * density(17) * density(37) 
  rrt(323) = rrt(323) * density(30) * density(31) 
  rrt(324) = rrt(324) * density(21) * density(31) 
  rrt(325) = rrt(325) * density(31) * density(37) 
  rrt(326) = rrt(326) * density(01) * density(31) 
  rrt(327) = rrt(327) * density(33) * density(37) 
  rrt(328) = rrt(328) * density(21) * density(38) 
  rrt(329) = rrt(329) * density(21) * density(38) 
  rrt(330) = rrt(330) * density(37) * density(38) 
  rrt(331) = rrt(331) * density(37) * density(38) 
  rrt(332) = rrt(332) * density(01) * density(38) 
  rrt(333) = rrt(333) * density(21) * density(40) 
  rrt(334) = rrt(334) * density(37) * density(40) 
  rrt(335) = rrt(335) * density(37) * density(40) 
  rrt(336) = rrt(336) * density(37) * density(41) 
  rrt(337) = rrt(337) * density(10) * density(14) 
  rrt(338) = rrt(338) * density(10) * density(14) 
  rrt(339) = rrt(339) * density(01) * density(10) 
  rrt(340) = rrt(340) * density(10)**2 
  rrt(341) = rrt(341) * density(10)**2 
  rrt(342) = rrt(342) * density(10)**2 
  rrt(343) = rrt(343) * density(01) * density(11) 
  rrt(344) = rrt(344) * density(01) * density(11) 
  rrt(345) = rrt(345) * density(01) * density(12) 
  rrt(346) = rrt(346) * density(12)**2 
  rrt(347) = rrt(347) * density(12)**2 
  rrt(348) = rrt(348) * density(10) * density(12) 
  rrt(349) = rrt(349) * density(01) * density(15) 
  rrt(350) = rrt(350) * density(14) * density(16) 
  rrt(351) = rrt(351) * density(14) * density(16) 
  rrt(352) = rrt(352) * density(01) * density(16) 
  rrt(353) = rrt(353) * density(10) * density(30) 
  rrt(354) = rrt(354) * density(10) * density(21) 
  rrt(355) = rrt(355) * density(10) * density(37) 
  rrt(356) = rrt(356) * density(11) * density(21) 
  rrt(357) = rrt(357) * density(12) * density(30) 
  rrt(358) = rrt(358) * density(12) * density(21) 
  rrt(359) = rrt(359) * density(14) * density(27) 
  rrt(360) = rrt(360) * density(14) * density(28) 
  rrt(361) = rrt(361) * density(14) * density(29) 
  rrt(362) = rrt(362) * density(14) * density(26) 
  rrt(363) = rrt(363) * density(14) * density(22) 
  rrt(364) = rrt(364) * density(14) * density(23) 
  rrt(365) = rrt(365) * density(14) * density(24) 
  rrt(366) = rrt(366) * density(14) * density(25) 
  rrt(367) = rrt(367) * density(15) * density(21) 
  rrt(368) = rrt(368) * density(15) * density(37) 
  rrt(369) = rrt(369) * density(16) * density(21) 
  rrt(370) = rrt(370) * density(14) * density(35) 
  rrt(371) = rrt(371) * density(30) * density(35) 
  rrt(372) = rrt(372) * density(35)**2 
  rrt(373) = rrt(373) * density(35)**2 
  rrt(374) = rrt(374) * density(35)**2 
  rrt(375) = rrt(375) * density(30) * density(36) 
  rrt(376) = rrt(376) * density(14) * density(36) 
  rrt(377) = rrt(377) * density(14) * density(36) 
  rrt(378) = rrt(378) * density(35) * density(36) 
  rrt(379) = rrt(379) * density(21) * density(36) 
  rrt(380) = rrt(380) * density(30) * density(37) 
  rrt(381) = rrt(381) * density(01) * density(14)**2 
  rrt(382) = rrt(382) * density(14)**2 * density(21) 
  rrt(383) = rrt(383) * density(02) * density(14)**2 
  rrt(384) = rrt(384) * density(03) * density(14)**2 
  rrt(385) = rrt(385) * density(04) * density(14)**2 
  rrt(386) = rrt(386) * density(05) * density(14)**2 
  rrt(387) = rrt(387) * density(06) * density(14)**2 
  rrt(388) = rrt(388) * density(07) * density(14)**2 
  rrt(389) = rrt(389) * density(08) * density(14)**2 
  rrt(390) = rrt(390) * density(09) * density(14)**2 
  rrt(391) = rrt(391) * density(14)**2 * density(27) 
  rrt(392) = rrt(392) * density(14)**2 * density(28) 
  rrt(393) = rrt(393) * density(14)**2 * density(29) 
  rrt(394) = rrt(394) * density(14)**3 
  rrt(395) = rrt(395) * density(14)**2 * density(30) 
  rrt(396) = rrt(396) * density(01) * density(14)**2 
  rrt(397) = rrt(397) * density(14)**2 * density(21) 
  rrt(398) = rrt(398) * density(02) * density(14)**2 
  rrt(399) = rrt(399) * density(03) * density(14)**2 
  rrt(400) = rrt(400) * density(04) * density(14)**2 
  rrt(401) = rrt(401) * density(05) * density(14)**2 
  rrt(402) = rrt(402) * density(06) * density(14)**2 
  rrt(403) = rrt(403) * density(07) * density(14)**2 
  rrt(404) = rrt(404) * density(08) * density(14)**2 
  rrt(405) = rrt(405) * density(09) * density(14)**2 
  rrt(406) = rrt(406) * density(14)**2 * density(27) 
  rrt(407) = rrt(407) * density(14)**2 * density(28) 
  rrt(408) = rrt(408) * density(14)**2 * density(29) 
  rrt(409) = rrt(409) * density(14)**3 
  rrt(410) = rrt(410) * density(14)**2 * density(30) 
  rrt(411) = rrt(411) * density(01) * density(14)**2 
  rrt(412) = rrt(412) * density(14)**2 * density(21) 
  rrt(413) = rrt(413) * density(21) * density(30)**2 
  rrt(414) = rrt(414) * density(01) * density(30)**2 
  rrt(415) = rrt(415) * density(01) * density(14) * density(30) 
  rrt(416) = rrt(416) * density(14) * density(21) * density(30) 
  rrt(417) = rrt(417) * density(01) * density(14) * density(21) 
  rrt(418) = rrt(418) * density(14) * density(21)**2 
  rrt(419) = rrt(419) * density(01) * density(30) * density(35) 
  rrt(420) = rrt(420) * density(21) * density(30) * density(35) 
  rrt(421) = rrt(421) * density(01) * density(30) * density(36) 
  rrt(422) = rrt(422) * density(21) * density(30) * density(36) 
  rrt(423) = rrt(423) * density(01) * density(21) * density(35) 
  rrt(424) = rrt(424) * density(21)**2 * density(35) 
  rrt(425) = rrt(425) * density(10) 
  rrt(426) = rrt(426) * density(12) 
  rrt(427) = rrt(427) * density(02) 
  rrt(428) = rrt(428) * density(03) 
  rrt(429) = rrt(429) * density(04) 
  rrt(430) = rrt(430) * density(05) 
  rrt(431) = rrt(431) * density(06) 
  rrt(432) = rrt(432) * density(07) 
  rrt(433) = rrt(433) * density(08) 
  rrt(434) = rrt(434) * density(09) 
  rrt(435) = rrt(435) * density(22) 
  rrt(436) = rrt(436) * density(23) 
  rrt(437) = rrt(437) * density(24) 
  rrt(438) = rrt(438) * density(25) 
  rrt(439) = rrt(439) * density(27) 
  rrt(440) = rrt(440) * density(28) 
  rrt(441) = rrt(441) * density(29) 
  rrt(442) = rrt(442) * density(26) * density(43) 
  rrt(443) = rrt(443) * density(31) * density(34) 
  rrt(444) = rrt(444) * density(32) * density(34) 
  rrt(445) = rrt(445) * density(18) * density(34) 
  rrt(446) = rrt(446) * density(20) * density(34) 
  rrt(447) = rrt(447) * density(34) * density(39) 
  rrt(448) = rrt(448) * density(01) * density(31) * density(34) 
  rrt(449) = rrt(449) * density(21) * density(31) * density(34) 
  rrt(450) = rrt(450) * density(14) * density(31) * density(34) 
  rrt(451) = rrt(451) * density(30) * density(31) * density(34) 
  rrt(452) = rrt(452) * density(01) * density(32) * density(34) 
  rrt(453) = rrt(453) * density(21) * density(32) * density(34) 
  rrt(454) = rrt(454) * density(14) * density(32) * density(34) 
  rrt(455) = rrt(455) * density(30) * density(32) * density(34) 
  rrt(456) = rrt(456) * density(01) * density(18) * density(34) 
  rrt(457) = rrt(457) * density(18) * density(21) * density(34) 
  rrt(458) = rrt(458) * density(14) * density(18) * density(34) 
  rrt(459) = rrt(459) * density(18) * density(30) * density(34) 
  rrt(460) = rrt(460) * density(01) * density(20) * density(34) 
  rrt(461) = rrt(461) * density(20) * density(21) * density(34) 
  rrt(462) = rrt(462) * density(14) * density(20) * density(34) 
  rrt(463) = rrt(463) * density(20) * density(30) * density(34) 
  rrt(464) = rrt(464) * density(01) * density(34) * density(39) 
  rrt(465) = rrt(465) * density(21) * density(34) * density(39) 
  rrt(466) = rrt(466) * density(14) * density(34) * density(39) 
  rrt(467) = rrt(467) * density(30) * density(34) * density(39) 
  rrt(468) = rrt(468) * density(14) * density(44) 
  rrt(469) = rrt(469) * density(30) * density(44) 
  rrt(470) = rrt(470) * density(35) * density(44) 
  rrt(471) = rrt(471) * density(36) * density(44) 
  rrt(472) = rrt(472) * density(14) * density(46) 
  rrt(473) = rrt(473) * density(30) * density(45) 
  rrt(474) = rrt(474) * density(14) * density(45) 
  rrt(475) = rrt(475) * density(35) * density(45) 
  rrt(476) = rrt(476) * density(36) * density(45) 
  rrt(477) = rrt(477) * density(30) * density(46) 
  rrt(478) = rrt(478) * density(30) * density(47) 
  rrt(479) = rrt(479) * density(30) * density(48) 
  rrt(480) = rrt(480) * density(21) * density(47) 
  rrt(481) = rrt(481) * density(45) * density(46) 
  rrt(482) = rrt(482) * density(45) * density(47) 
  rrt(483) = rrt(483) * density(45) * density(48) 
  rrt(484) = rrt(484) * density(01) * density(44)**2 
  rrt(485) = rrt(485) * density(02) * density(44)**2 
  rrt(486) = rrt(486) * density(03) * density(44)**2 
  rrt(487) = rrt(487) * density(04) * density(44)**2 
  rrt(488) = rrt(488) * density(05) * density(44)**2 
  rrt(489) = rrt(489) * density(06) * density(44)**2 
  rrt(490) = rrt(490) * density(07) * density(44)**2 
  rrt(491) = rrt(491) * density(08) * density(44)**2 
  rrt(492) = rrt(492) * density(09) * density(44)**2 
  rrt(493) = rrt(493) * density(10) * density(44)**2 
  rrt(494) = rrt(494) * density(11) * density(44)**2 
  rrt(495) = rrt(495) * density(12) * density(44)**2 
  rrt(496) = rrt(496) * density(13) * density(44)**2 
  rrt(497) = rrt(497) * density(21) * density(44)**2 
  rrt(498) = rrt(498) * density(27) * density(44)**2 
  rrt(499) = rrt(499) * density(28) * density(44)**2 
  rrt(500) = rrt(500) * density(29) * density(44)**2 
  rrt(501) = rrt(501) * density(22) * density(44)**2 
  rrt(502) = rrt(502) * density(23) * density(44)**2 
  rrt(503) = rrt(503) * density(24) * density(44)**2 
  rrt(504) = rrt(504) * density(25) * density(44)**2 
  ydot(01) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(012)-rrt(017)-rrt(018)-rrt(019)-rrt(020)-rrt(021)-rrt(022)&
             -rrt(023)-rrt(024)+rrt(025)+rrt(026)+rrt(027)+rrt(028)+rrt(029)+rrt(030)+rrt(031)+rrt(032)+rrt(043)-rrt(044)+rrt(059)&
             -rrt(060)+rrt(075)-rrt(076)+rrt(091)-rrt(092)+rrt(119)-rrt(120)+rrt(121)-rrt(122)+rrt(123)-rrt(124)+rrt(125)-rrt(126)&
             +rrt(127)-rrt(128)+rrt(129)-rrt(130)+rrt(131)-rrt(132)-rrt(175)+rrt(176)-rrt(191)+rrt(192)-rrt(207)+rrt(208)+rrt(223)&
             -rrt(224)+rrt(237)-rrt(238)+rrt(251)-rrt(252)+rrt(271)-rrt(274)+rrt(284)-rrt(287)+rrt(296)+  2.d0 * rrt(297)+rrt(309)&
             +rrt(312)-rrt(314)+rrt(315)+rrt(316)+rrt(317)+  2.d0 * rrt(318)-rrt(326)-rrt(332)+rrt(337)+rrt(338)+rrt(339)+rrt(340)&
             +rrt(341)+rrt(342)+rrt(344)+rrt(346)+rrt(353)+rrt(354)+rrt(355)+rrt(357)+rrt(358)+rrt(370)+rrt(372)+rrt(374)+rrt(376)&
             +rrt(377)+rrt(411)+rrt(412)+rrt(425)+rrt(427)+rrt(445)+  2.d0 * rrt(446)+rrt(447)+rrt(456)+rrt(457)+rrt(458)+rrt(459)&
             +  2.d0 * rrt(460)+  2.d0 * rrt(461)+  2.d0 * rrt(462)+  2.d0 * rrt(463)+rrt(464)+rrt(465)+rrt(466)+rrt(467)+rrt(472)&
             -rrt(484) 
  ydot(02) = +rrt(017)-rrt(025)-rrt(043)+rrt(044)+rrt(045)-rrt(046)-rrt(059)+rrt(060)+rrt(061)-rrt(062)-rrt(075)+rrt(076)+rrt(077)&
             -rrt(078)-rrt(091)+rrt(092)+rrt(093)-rrt(094)-  2.d0 * rrt(119)+  2.d0 * rrt(120)-rrt(121)+rrt(122)-rrt(123)+rrt(124)&
             -rrt(125)+rrt(126)-rrt(127)+rrt(128)-rrt(129)+rrt(130)-rrt(131)+rrt(132)+rrt(133)-rrt(134)+rrt(135)-rrt(136)+rrt(137)&
             -rrt(138)+rrt(139)-rrt(140)+rrt(141)-rrt(142)+rrt(143)-rrt(144)+rrt(175)-rrt(176)-rrt(177)+rrt(178)+rrt(191)-rrt(192)&
             -rrt(193)+rrt(194)+rrt(207)-rrt(208)-rrt(209)+rrt(210)+rrt(225)-rrt(226)+rrt(239)-rrt(240)+rrt(253)-rrt(254)+rrt(272)&
             -rrt(275)-rrt(427)+rrt(428)-rrt(485) 
  ydot(03) = +rrt(018)-rrt(026)-rrt(045)+rrt(046)+rrt(047)-rrt(048)-rrt(061)+rrt(062)+rrt(063)-rrt(064)-rrt(077)+rrt(078)+rrt(079)&
             -rrt(080)-rrt(093)+rrt(094)+rrt(095)-rrt(096)+rrt(119)-rrt(120)-rrt(121)+rrt(122)-  2.d0 * rrt(133)+  2.d0 * rrt(134)&
             -rrt(135)+rrt(136)-rrt(137)+rrt(138)-rrt(139)+rrt(140)-rrt(141)+rrt(142)-rrt(143)+rrt(144)+rrt(145)-rrt(146)+rrt(147)&
             -rrt(148)+rrt(149)-rrt(150)+rrt(151)-rrt(152)+rrt(153)-rrt(154)+rrt(177)-rrt(178)-rrt(179)+rrt(180)+rrt(193)-rrt(194)&
             -rrt(195)+rrt(196)+rrt(209)-rrt(210)-rrt(211)+rrt(212)-rrt(223)+rrt(224)+rrt(227)-rrt(228)-rrt(237)+rrt(238)+rrt(241)&
             -rrt(242)-rrt(251)+rrt(252)+rrt(255)-rrt(256)+rrt(273)-rrt(276)-rrt(428)+rrt(429)-rrt(486) 
  ydot(04) = +rrt(019)-rrt(027)-rrt(047)+rrt(048)+rrt(049)-rrt(050)-rrt(063)+rrt(064)+rrt(065)-rrt(066)-rrt(079)+rrt(080)+rrt(081)&
             -rrt(082)-rrt(095)+rrt(096)+rrt(097)-rrt(098)+rrt(121)-rrt(122)-rrt(123)+rrt(124)+rrt(133)-rrt(134)-rrt(135)+rrt(136)&
             -  2.d0 * rrt(145)+  2.d0 * rrt(146)-rrt(147)+rrt(148)-rrt(149)+rrt(150)-rrt(151)+rrt(152)-rrt(153)+rrt(154)+rrt(155)&
             -rrt(156)+rrt(157)-rrt(158)+rrt(159)-rrt(160)+rrt(161)-rrt(162)+rrt(179)-rrt(180)-rrt(181)+rrt(182)+rrt(195)-rrt(196)&
             -rrt(197)+rrt(198)+rrt(211)-rrt(212)-rrt(213)+rrt(214)-rrt(225)+rrt(226)+rrt(229)-rrt(230)-rrt(239)+rrt(240)+rrt(243)&
             -rrt(244)-rrt(253)+rrt(254)+rrt(257)-rrt(258)-rrt(277)-rrt(429)+rrt(430)-rrt(487) 
  ydot(05) = +rrt(020)-rrt(028)-rrt(049)+rrt(050)+rrt(051)-rrt(052)-rrt(065)+rrt(066)+rrt(067)-rrt(068)-rrt(081)+rrt(082)+rrt(083)&
             -rrt(084)-rrt(097)+rrt(098)+rrt(099)-rrt(100)+rrt(123)-rrt(124)-rrt(125)+rrt(126)+rrt(135)-rrt(136)-rrt(137)+rrt(138)&
             +rrt(145)-rrt(146)-rrt(147)+rrt(148)-  2.d0 * rrt(155)+  2.d0 * rrt(156)-rrt(157)+rrt(158)-rrt(159)+rrt(160)-rrt(161)&
             +rrt(162)+rrt(163)-rrt(164)+rrt(165)-rrt(166)+rrt(167)-rrt(168)+rrt(181)-rrt(182)-rrt(183)+rrt(184)+rrt(197)-rrt(198)&
             -rrt(199)+rrt(200)+rrt(213)-rrt(214)-rrt(215)+rrt(216)-rrt(227)+rrt(228)+rrt(231)-rrt(232)-rrt(241)+rrt(242)+rrt(245)&
             -rrt(246)-rrt(255)+rrt(256)+rrt(259)-rrt(260)-rrt(278)-rrt(430)+rrt(431)-rrt(488) 
  ydot(06) = +rrt(021)-rrt(029)-rrt(051)+rrt(052)+rrt(053)-rrt(054)-rrt(067)+rrt(068)+rrt(069)-rrt(070)-rrt(083)+rrt(084)+rrt(085)&
             -rrt(086)-rrt(099)+rrt(100)+rrt(101)-rrt(102)+rrt(125)-rrt(126)-rrt(127)+rrt(128)+rrt(137)-rrt(138)-rrt(139)+rrt(140)&
             +rrt(147)-rrt(148)-rrt(149)+rrt(150)+rrt(155)-rrt(156)-rrt(157)+rrt(158)-  2.d0 * rrt(163)+  2.d0 * rrt(164)-rrt(165)&
             +rrt(166)-rrt(167)+rrt(168)+rrt(169)-rrt(170)+rrt(171)-rrt(172)+rrt(183)-rrt(184)-rrt(185)+rrt(186)+rrt(199)-rrt(200)&
             -rrt(201)+rrt(202)+rrt(215)-rrt(216)-rrt(217)+rrt(218)-rrt(229)+rrt(230)+rrt(233)-rrt(234)-rrt(243)+rrt(244)+rrt(247)&
             -rrt(248)-rrt(257)+rrt(258)+rrt(261)-rrt(262)-rrt(279)-rrt(431)+rrt(432)-rrt(489) 
  ydot(07) = +rrt(022)-rrt(030)-rrt(053)+rrt(054)+rrt(055)-rrt(056)-rrt(069)+rrt(070)+rrt(071)-rrt(072)-rrt(085)+rrt(086)+rrt(087)&
             -rrt(088)-rrt(101)+rrt(102)+rrt(103)-rrt(104)+rrt(127)-rrt(128)-rrt(129)+rrt(130)+rrt(139)-rrt(140)-rrt(141)+rrt(142)&
             +rrt(149)-rrt(150)-rrt(151)+rrt(152)+rrt(157)-rrt(158)-rrt(159)+rrt(160)+rrt(163)-rrt(164)-rrt(165)+rrt(166)&
             -  2.d0 * rrt(169)+  2.d0 * rrt(170)-rrt(171)+rrt(172)+rrt(173)-rrt(174)+rrt(185)-rrt(186)-rrt(187)+rrt(188)+rrt(201)&
             -rrt(202)-rrt(203)+rrt(204)+rrt(217)-rrt(218)-rrt(219)+rrt(220)-rrt(231)+rrt(232)+rrt(235)-rrt(236)-rrt(245)+rrt(246)&
             +rrt(249)-rrt(250)-rrt(259)+rrt(260)+rrt(263)-rrt(264)-rrt(271)+rrt(274)-rrt(280)-rrt(432)+rrt(433)-rrt(490) 
  ydot(08) = +rrt(023)-rrt(031)-rrt(055)+rrt(056)+rrt(057)-rrt(058)-rrt(071)+rrt(072)+rrt(073)-rrt(074)-rrt(087)+rrt(088)+rrt(089)&
             -rrt(090)-rrt(103)+rrt(104)+rrt(105)-rrt(106)+rrt(129)-rrt(130)-rrt(131)+rrt(132)+rrt(141)-rrt(142)-rrt(143)+rrt(144)&
             +rrt(151)-rrt(152)-rrt(153)+rrt(154)+rrt(159)-rrt(160)-rrt(161)+rrt(162)+rrt(165)-rrt(166)-rrt(167)+rrt(168)+rrt(169)&
             -rrt(170)-rrt(171)+rrt(172)-  2.d0 * rrt(173)+  2.d0 * rrt(174)+rrt(187)-rrt(188)-rrt(189)+rrt(190)+rrt(203)-rrt(204)&
             -rrt(205)+rrt(206)+rrt(219)-rrt(220)-rrt(221)+rrt(222)-rrt(233)+rrt(234)-rrt(247)+rrt(248)-rrt(261)+rrt(262)-rrt(272)&
             +rrt(275)-rrt(281)-rrt(433)+rrt(434)-rrt(491) 
  ydot(09) = +rrt(024)-rrt(032)-rrt(057)+rrt(058)-rrt(073)+rrt(074)-rrt(089)+rrt(090)-rrt(105)+rrt(106)+rrt(131)-rrt(132)+rrt(143)&
             -rrt(144)+rrt(153)-rrt(154)+rrt(161)-rrt(162)+rrt(167)-rrt(168)+rrt(171)-rrt(172)+rrt(173)-rrt(174)+rrt(189)-rrt(190)&
             +rrt(205)-rrt(206)+rrt(221)-rrt(222)-rrt(235)+rrt(236)-rrt(249)+rrt(250)-rrt(263)+rrt(264)-rrt(273)+rrt(276)+rrt(277)&
             +rrt(278)+rrt(279)+rrt(280)+rrt(281)-rrt(434)-rrt(492) 
  ydot(10) = +rrt(001)+rrt(002)+rrt(003)-rrt(271)-rrt(272)-rrt(273)+rrt(274)+rrt(275)+rrt(276)+rrt(277)+rrt(278)+rrt(279)+rrt(280)&
             +rrt(281)+rrt(282)+rrt(283)-rrt(311)-rrt(337)-rrt(338)-rrt(339)-  2.d0 * rrt(340)-  2.d0 * rrt(341)-  2.d0 * rrt(342)&
             +rrt(343)-rrt(348)-rrt(353)-rrt(354)-rrt(355)+rrt(356)+rrt(381)+rrt(382)+rrt(383)+rrt(384)+rrt(385)+rrt(386)+rrt(387)&
             +rrt(388)+rrt(389)+rrt(390)+rrt(391)+rrt(392)+rrt(393)+rrt(394)+rrt(395)-rrt(425)-rrt(493) 
  ydot(11) = +rrt(004)+rrt(271)+rrt(272)+rrt(273)-rrt(274)-rrt(275)-rrt(276)-rrt(277)-rrt(278)-rrt(279)-rrt(280)-rrt(281)-rrt(282)&
             -rrt(283)+rrt(285)+rrt(340)-rrt(343)-rrt(344)+rrt(345)-rrt(356)+rrt(396)+rrt(397)+rrt(398)+rrt(399)+rrt(400)+rrt(401)&
             +rrt(402)+rrt(403)+rrt(404)+rrt(405)+rrt(406)+rrt(407)+rrt(408)+rrt(409)+rrt(410)+rrt(426)-rrt(494) 
  ydot(12) = +rrt(005)-rrt(284)-rrt(345)-  2.d0 * rrt(346)-  2.d0 * rrt(347)-rrt(348)-rrt(357)-rrt(358)-rrt(426)-rrt(495) 
  ydot(13) = +rrt(006)-rrt(285)+rrt(341)-rrt(496) 
  ydot(14) = -rrt(014)+  2.d0 * rrt(287)+rrt(288)+rrt(289)+  2.d0 * rrt(293)+rrt(294)+rrt(295)+rrt(296)+rrt(302)+rrt(304)+rrt(311)&
             -rrt(312)-rrt(313)-rrt(316)-rrt(318)+rrt(321)+rrt(328)+rrt(331)+rrt(332)-rrt(338)+  2.d0 * rrt(342)+rrt(349)+rrt(350)&
             +rrt(352)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)-rrt(370)+rrt(371)+rrt(373)-rrt(376)&
             -rrt(377)+rrt(378)-  2.d0 * rrt(381)-  2.d0 * rrt(382)-  2.d0 * rrt(383)-  2.d0 * rrt(384)-  2.d0 * rrt(385)&
             -  2.d0 * rrt(386)-  2.d0 * rrt(387)-  2.d0 * rrt(388)-  2.d0 * rrt(389)-  2.d0 * rrt(390)-  2.d0 * rrt(391)&
             -  2.d0 * rrt(392)-  2.d0 * rrt(393)-  2.d0 * rrt(394)-  2.d0 * rrt(395)-  2.d0 * rrt(396)-  2.d0 * rrt(397)&
             -  2.d0 * rrt(398)-  2.d0 * rrt(399)-  2.d0 * rrt(400)-  2.d0 * rrt(401)-  2.d0 * rrt(402)-  2.d0 * rrt(403)&
             -  2.d0 * rrt(404)-  2.d0 * rrt(405)-  2.d0 * rrt(406)-  2.d0 * rrt(407)-  2.d0 * rrt(408)-  2.d0 * rrt(409)&
             -  2.d0 * rrt(410)-  2.d0 * rrt(411)-  2.d0 * rrt(412)-rrt(415)-rrt(416)-rrt(417)-rrt(418)-rrt(468)-rrt(472)-rrt(474) 
  ydot(15) = +rrt(294)-rrt(349)+rrt(351)-rrt(367)-rrt(368) 
  ydot(16) = +rrt(295)+rrt(338)-rrt(350)-rrt(351)-rrt(352)-rrt(369) 
  ydot(17) = +rrt(014)+rrt(312)+rrt(318)-rrt(319)-rrt(320)-rrt(321)-rrt(322) 
  ydot(18) = +rrt(012)-rrt(293)-rrt(294)-rrt(295)-rrt(310)-rrt(311)-rrt(312)-rrt(313)-rrt(314)-rrt(315)+rrt(316)+rrt(317)+rrt(346)&
             -rrt(445)-rrt(456)-rrt(457)-rrt(458)-rrt(459) 
  ydot(19) = -rrt(296)+rrt(311)+rrt(313)-rrt(316) 
  ydot(20) = -rrt(297)+rrt(314)-rrt(317)-rrt(318)+rrt(347)+rrt(348)-rrt(446)-rrt(460)-rrt(461)-rrt(462)-rrt(463) 
  ydot(21) = -rrt(007)-rrt(008)-rrt(009)-rrt(010)-rrt(011)-rrt(013)-rrt(033)-rrt(034)-rrt(035)+rrt(036)+rrt(037)+rrt(038)+rrt(039)&
             +rrt(040)+rrt(041)+rrt(042)+rrt(107)-rrt(108)+rrt(113)-rrt(114)+rrt(175)-rrt(176)+rrt(177)-rrt(178)+rrt(179)-rrt(180)&
             +rrt(181)-rrt(182)+rrt(183)-rrt(184)+rrt(185)-rrt(186)+rrt(187)-rrt(188)+rrt(189)-rrt(190)-rrt(223)+rrt(224)-rrt(225)&
             +rrt(226)-rrt(227)+rrt(228)-rrt(229)+rrt(230)-rrt(231)+rrt(232)-rrt(233)+rrt(234)-rrt(235)+rrt(236)+rrt(267)-rrt(268)&
             +rrt(269)-rrt(270)-rrt(286)+rrt(289)+rrt(292)-rrt(298)+rrt(301)-rrt(310)-rrt(319)+rrt(322)+rrt(323)-rrt(324)+rrt(325)&
             -rrt(328)-rrt(329)-rrt(333)-rrt(354)-rrt(358)-rrt(367)-rrt(369)+rrt(371)+rrt(372)+rrt(375)+rrt(377)-rrt(379)+rrt(380)&
             +rrt(413)+rrt(414)-rrt(417)-rrt(418)-rrt(423)-rrt(424)+rrt(435)+rrt(436)+rrt(437)+rrt(438)+rrt(439)+rrt(444)+rrt(447)&
             +rrt(448)+rrt(449)+rrt(450)+rrt(451)+  2.d0 * rrt(452)+  2.d0 * rrt(453)+  2.d0 * rrt(454)+  2.d0 * rrt(455)+rrt(464)&
             +rrt(465)+rrt(466)+rrt(467)+rrt(473)-rrt(480)-rrt(497) 
  ydot(22) = +rrt(007)-rrt(039)-rrt(363)-rrt(435)-rrt(501) 
  ydot(23) = +rrt(008)-rrt(040)-rrt(364)-rrt(436)-rrt(502) 
  ydot(24) = +rrt(009)-rrt(041)-rrt(365)-rrt(437)-rrt(503) 
  ydot(25) = +rrt(010)-rrt(042)-rrt(366)-rrt(438)-rrt(504) 
  ydot(26) = +rrt(011)-rrt(362)-rrt(442) 
  ydot(27) = +rrt(033)-rrt(036)-rrt(107)+rrt(108)+rrt(109)-rrt(110)-rrt(113)+rrt(114)+rrt(115)-rrt(116)-rrt(175)+rrt(176)-rrt(177)&
             +rrt(178)-rrt(179)+rrt(180)-rrt(181)+rrt(182)-rrt(183)+rrt(184)-rrt(185)+rrt(186)-rrt(187)+rrt(188)-rrt(189)+rrt(190)&
             +rrt(191)-rrt(192)+rrt(193)-rrt(194)+rrt(195)-rrt(196)+rrt(197)-rrt(198)+rrt(199)-rrt(200)+rrt(201)-rrt(202)+rrt(203)&
             -rrt(204)+rrt(205)-rrt(206)+rrt(223)-rrt(224)+rrt(225)-rrt(226)+rrt(227)-rrt(228)+rrt(229)-rrt(230)+rrt(231)-rrt(232)&
             +rrt(233)-rrt(234)+rrt(235)-rrt(236)-rrt(237)+rrt(238)-rrt(239)+rrt(240)-rrt(241)+rrt(242)-rrt(243)+rrt(244)-rrt(245)&
             +rrt(246)-rrt(247)+rrt(248)-rrt(249)+rrt(250)+rrt(265)-rrt(266)-rrt(267)+rrt(268)-  2.d0 * rrt(269)+  2.d0 * rrt(270)&
             -rrt(359)-rrt(439)+rrt(440)-rrt(498) 
  ydot(28) = +rrt(034)-rrt(037)-rrt(109)+rrt(110)+rrt(111)-rrt(112)-rrt(115)+rrt(116)+rrt(117)-rrt(118)-rrt(191)+rrt(192)-rrt(193)&
             +rrt(194)-rrt(195)+rrt(196)-rrt(197)+rrt(198)-rrt(199)+rrt(200)-rrt(201)+rrt(202)-rrt(203)+rrt(204)-rrt(205)+rrt(206)&
             +rrt(207)-rrt(208)+rrt(209)-rrt(210)+rrt(211)-rrt(212)+rrt(213)-rrt(214)+rrt(215)-rrt(216)+rrt(217)-rrt(218)+rrt(219)&
             -rrt(220)+rrt(221)-rrt(222)+rrt(237)-rrt(238)+rrt(239)-rrt(240)+rrt(241)-rrt(242)+rrt(243)-rrt(244)+rrt(245)-rrt(246)&
             +rrt(247)-rrt(248)+rrt(249)-rrt(250)-rrt(251)+rrt(252)-rrt(253)+rrt(254)-rrt(255)+rrt(256)-rrt(257)+rrt(258)-rrt(259)&
             +rrt(260)-rrt(261)+rrt(262)-rrt(263)+rrt(264)-  2.d0 * rrt(265)+  2.d0 * rrt(266)-rrt(267)+rrt(268)+rrt(269)-rrt(270)&
             -rrt(360)-rrt(440)+rrt(441)-rrt(499) 
  ydot(29) = +rrt(035)-rrt(038)-rrt(111)+rrt(112)-rrt(117)+rrt(118)-rrt(207)+rrt(208)-rrt(209)+rrt(210)-rrt(211)+rrt(212)-rrt(213)&
             +rrt(214)-rrt(215)+rrt(216)-rrt(217)+rrt(218)-rrt(219)+rrt(220)-rrt(221)+rrt(222)+rrt(251)-rrt(252)+rrt(253)-rrt(254)&
             +rrt(255)-rrt(256)+rrt(257)-rrt(258)+rrt(259)-rrt(260)+rrt(261)-rrt(262)+rrt(263)-rrt(264)+rrt(265)-rrt(266)+rrt(267)&
             -rrt(268)-rrt(361)-rrt(441)-rrt(500) 
  ydot(30) = -rrt(015)+  2.d0 * rrt(286)+rrt(288)+rrt(290)+rrt(291)+rrt(298)+  2.d0 * rrt(299)+  3.d0 * rrt(300)+rrt(301)+rrt(302)&
             +rrt(303)+  2.d0 * rrt(304)+  2.d0 * rrt(305)+rrt(306)+rrt(307)+  2.d0 * rrt(308)+rrt(309)+rrt(310)+rrt(319)-rrt(323)&
             +rrt(324)+rrt(326)+rrt(327)+rrt(329)+rrt(333)+  2.d0 * rrt(354)+  2.d0 * rrt(358)+rrt(359)+rrt(360)+rrt(361)+rrt(362)&
             +rrt(363)+rrt(364)+rrt(365)+rrt(366)+rrt(367)+rrt(369)+rrt(370)-rrt(371)+  2.d0 * rrt(374)-rrt(375)+  2.d0 * rrt(376)&
             +rrt(379)-rrt(380)-  2.d0 * rrt(413)-  2.d0 * rrt(414)-rrt(415)-rrt(416)-rrt(419)-rrt(420)-rrt(421)-rrt(422)+rrt(442)&
             +  3.d0 * rrt(443)+  2.d0 * rrt(444)+rrt(445)+rrt(446)+rrt(448)+rrt(449)+rrt(450)+rrt(451)+rrt(456)+rrt(457)+rrt(458)&
             +rrt(459)+rrt(460)+rrt(461)+rrt(462)+rrt(463)-rrt(469)-rrt(473)-rrt(477)-rrt(478)-rrt(479) 
  ydot(31) = +rrt(013)-rrt(299)-rrt(323)-rrt(324)-rrt(325)-rrt(326)-rrt(443)-rrt(448)-rrt(449)-rrt(450)-rrt(451) 
  ydot(32) = -rrt(300)-rrt(301)+rrt(324)+rrt(328)-rrt(444)-rrt(452)-rrt(453)-rrt(454)-rrt(455) 
  ydot(33) = +rrt(015)+rrt(298)+rrt(323)-rrt(327) 
  ydot(34) = +rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)-rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(451)-rrt(452)-rrt(453)-rrt(454)&
             -rrt(455)-rrt(456)-rrt(457)-rrt(458)-rrt(459)-rrt(460)-rrt(461)-rrt(462)-rrt(463)-rrt(464)-rrt(465)-rrt(466)-rrt(467) 
  ydot(35) = -rrt(288)+rrt(290)+rrt(292)+rrt(303)+rrt(305)+rrt(320)+rrt(330)+rrt(335)+rrt(359)+rrt(360)+rrt(361)+rrt(362)+rrt(363)&
             +rrt(364)+rrt(365)+rrt(366)+rrt(367)+rrt(368)+rrt(369)-rrt(370)-rrt(371)-  2.d0 * rrt(372)-  2.d0 * rrt(373)&
             -  2.d0 * rrt(374)+rrt(375)-rrt(378)+rrt(415)+rrt(416)-rrt(419)-rrt(420)-rrt(423)-rrt(424)-rrt(470)-rrt(475) 
  ydot(36) = -rrt(289)-rrt(290)+rrt(291)+rrt(306)+rrt(308)+rrt(334)+rrt(336)+rrt(368)+rrt(373)-rrt(375)-rrt(376)-rrt(377)-rrt(378)&
             -rrt(379)+rrt(380)+rrt(417)+rrt(418)+rrt(419)+rrt(420)-rrt(421)-rrt(422)-rrt(471)-rrt(476) 
  ydot(37) = -rrt(016)-rrt(291)-rrt(292)+rrt(307)-rrt(315)-rrt(320)-rrt(321)-rrt(322)-rrt(325)-rrt(327)-rrt(330)-rrt(331)-rrt(334)&
             -rrt(335)-rrt(336)-rrt(368)+rrt(378)+rrt(379)-rrt(380)+rrt(421)+rrt(422)+rrt(423)+rrt(424)+rrt(476)+rrt(479)+rrt(480)&
             +rrt(483) 
  ydot(38) = -rrt(302)+rrt(319)-rrt(328)-rrt(329)-rrt(330)-rrt(331)-rrt(332) 
  ydot(39) = -rrt(309)+rrt(310)+rrt(322)+rrt(326)+rrt(332)-rrt(447)-rrt(464)-rrt(465)-rrt(466)-rrt(467) 
  ydot(40) = -rrt(303)-rrt(304)+rrt(320)+rrt(329)-rrt(333)-rrt(334)-rrt(335) 
  ydot(41) = +rrt(016)-rrt(305)-rrt(306)+rrt(315)+rrt(321)+rrt(325)+rrt(327)+rrt(330)+rrt(333)+rrt(334)-rrt(336) 
  ydot(42) = -rrt(307)-rrt(308)+rrt(331)+rrt(335)+rrt(336) 
  ydot(43) = +rrt(012)+rrt(013)+rrt(014)+rrt(015)+rrt(016)-rrt(293)-rrt(294)-rrt(295)-rrt(296)-rrt(297)+rrt(298)-rrt(299)-rrt(300)&
             -rrt(301)-rrt(302)-rrt(303)-rrt(304)-rrt(305)-rrt(306)-rrt(307)-rrt(308)-rrt(309)+rrt(346)+rrt(347)+rrt(348)-rrt(442) 
  ydot(44) = -rrt(468)-rrt(469)-rrt(470)-rrt(471)+rrt(472)+rrt(473)+rrt(476)+rrt(479)+rrt(480)+rrt(481)+rrt(482)+  2.d0 * rrt(483)&
             -  2.d0 * rrt(484)-  2.d0 * rrt(485)-  2.d0 * rrt(486)-  2.d0 * rrt(487)-  2.d0 * rrt(488)-  2.d0 * rrt(489)&
             -  2.d0 * rrt(490)-  2.d0 * rrt(491)-  2.d0 * rrt(492)-  2.d0 * rrt(493)-  2.d0 * rrt(494)-  2.d0 * rrt(495)&
             -  2.d0 * rrt(496)-  2.d0 * rrt(497)-  2.d0 * rrt(498)-  2.d0 * rrt(499)-  2.d0 * rrt(500)-  2.d0 * rrt(501)&
             -  2.d0 * rrt(502)-  2.d0 * rrt(503)-  2.d0 * rrt(504) 
  ydot(45) = +rrt(469)-rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(481)-rrt(482)-rrt(483)+  2.d0 * rrt(497)+  2.d0 * rrt(498)&
             +  2.d0 * rrt(499)+  2.d0 * rrt(500)+  2.d0 * rrt(501)+  2.d0 * rrt(502)+  2.d0 * rrt(503)+  2.d0 * rrt(504) 
  ydot(46) = +rrt(468)-rrt(472)-rrt(477)-rrt(481)+  2.d0 * rrt(484)+  2.d0 * rrt(485)+  2.d0 * rrt(486)+  2.d0 * rrt(487)&
             +  2.d0 * rrt(488)+  2.d0 * rrt(489)+  2.d0 * rrt(490)+  2.d0 * rrt(491)+  2.d0 * rrt(492)+  2.d0 * rrt(493)&
             +  2.d0 * rrt(494)+  2.d0 * rrt(495)+  2.d0 * rrt(496) 
  ydot(47) = +rrt(470)+rrt(474)+rrt(477)-rrt(478)-rrt(480)+rrt(481)-rrt(482) 
  ydot(48) = +rrt(471)+rrt(475)+rrt(478)-rrt(479)+rrt(482)-rrt(483) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(49) = 0.0d0
  if( lgas_heating ) then
    ydot(49) = ( ZDPlasKin_cfg(14)/k_B + ydot(49) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(49) = ydot(49) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(49)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(001) * density(43) 
  pd(01,43) = pd(01,43) - rrt(001) * density(01) 
  pd(10,01) = pd(10,01) + rrt(001) * density(43) 
  pd(10,43) = pd(10,43) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) - rrt(002) * density(43) 
  pd(01,43) = pd(01,43) - rrt(002) * density(01) 
  pd(10,01) = pd(10,01) + rrt(002) * density(43) 
  pd(10,43) = pd(10,43) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(43) 
  pd(01,43) = pd(01,43) - rrt(003) * density(01) 
  pd(10,01) = pd(10,01) + rrt(003) * density(43) 
  pd(10,43) = pd(10,43) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(43) 
  pd(01,43) = pd(01,43) - rrt(004) * density(01) 
  pd(11,01) = pd(11,01) + rrt(004) * density(43) 
  pd(11,43) = pd(11,43) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) - rrt(005) * density(43) 
  pd(01,43) = pd(01,43) - rrt(005) * density(01) 
  pd(12,01) = pd(12,01) + rrt(005) * density(43) 
  pd(12,43) = pd(12,43) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(43) 
  pd(01,43) = pd(01,43) - rrt(006) * density(01) 
  pd(13,01) = pd(13,01) + rrt(006) * density(43) 
  pd(13,43) = pd(13,43) + rrt(006) * density(01) 
  pd(21,21) = pd(21,21) - rrt(007) * density(43) 
  pd(21,43) = pd(21,43) - rrt(007) * density(21) 
  pd(22,21) = pd(22,21) + rrt(007) * density(43) 
  pd(22,43) = pd(22,43) + rrt(007) * density(21) 
  pd(21,21) = pd(21,21) - rrt(008) * density(43) 
  pd(21,43) = pd(21,43) - rrt(008) * density(21) 
  pd(23,21) = pd(23,21) + rrt(008) * density(43) 
  pd(23,43) = pd(23,43) + rrt(008) * density(21) 
  pd(21,21) = pd(21,21) - rrt(009) * density(43) 
  pd(21,43) = pd(21,43) - rrt(009) * density(21) 
  pd(24,21) = pd(24,21) + rrt(009) * density(43) 
  pd(24,43) = pd(24,43) + rrt(009) * density(21) 
  pd(21,21) = pd(21,21) - rrt(010) * density(43) 
  pd(21,43) = pd(21,43) - rrt(010) * density(21) 
  pd(25,21) = pd(25,21) + rrt(010) * density(43) 
  pd(25,43) = pd(25,43) + rrt(010) * density(21) 
  pd(21,21) = pd(21,21) - rrt(011) * density(43) 
  pd(21,43) = pd(21,43) - rrt(011) * density(21) 
  pd(26,21) = pd(26,21) + rrt(011) * density(43) 
  pd(26,43) = pd(26,43) + rrt(011) * density(21) 
  pd(01,01) = pd(01,01) - rrt(012) * density(43) 
  pd(01,43) = pd(01,43) - rrt(012) * density(01) 
  pd(18,01) = pd(18,01) + rrt(012) * density(43) 
  pd(18,43) = pd(18,43) + rrt(012) * density(01) 
  pd(43,01) = pd(43,01) + rrt(012) * density(43) 
  pd(43,43) = pd(43,43) + rrt(012) * density(01) 
  pd(21,21) = pd(21,21) - rrt(013) * density(43) 
  pd(21,43) = pd(21,43) - rrt(013) * density(21) 
  pd(31,21) = pd(31,21) + rrt(013) * density(43) 
  pd(31,43) = pd(31,43) + rrt(013) * density(21) 
  pd(43,21) = pd(43,21) + rrt(013) * density(43) 
  pd(43,43) = pd(43,43) + rrt(013) * density(21) 
  pd(14,14) = pd(14,14) - rrt(014) * density(43) 
  pd(14,43) = pd(14,43) - rrt(014) * density(14) 
  pd(17,14) = pd(17,14) + rrt(014) * density(43) 
  pd(17,43) = pd(17,43) + rrt(014) * density(14) 
  pd(43,14) = pd(43,14) + rrt(014) * density(43) 
  pd(43,43) = pd(43,43) + rrt(014) * density(14) 
  pd(30,30) = pd(30,30) - rrt(015) * density(43) 
  pd(30,43) = pd(30,43) - rrt(015) * density(30) 
  pd(33,30) = pd(33,30) + rrt(015) * density(43) 
  pd(33,43) = pd(33,43) + rrt(015) * density(30) 
  pd(43,30) = pd(43,30) + rrt(015) * density(43) 
  pd(43,43) = pd(43,43) + rrt(015) * density(30) 
  pd(37,37) = pd(37,37) - rrt(016) * density(43) 
  pd(37,43) = pd(37,43) - rrt(016) * density(37) 
  pd(41,37) = pd(41,37) + rrt(016) * density(43) 
  pd(41,43) = pd(41,43) + rrt(016) * density(37) 
  pd(43,37) = pd(43,37) + rrt(016) * density(43) 
  pd(43,43) = pd(43,43) + rrt(016) * density(37) 
  pd(01,01) = pd(01,01) - rrt(017) * density(43) 
  pd(01,43) = pd(01,43) - rrt(017) * density(01) 
  pd(02,01) = pd(02,01) + rrt(017) * density(43) 
  pd(02,43) = pd(02,43) + rrt(017) * density(01) 
  pd(01,01) = pd(01,01) - rrt(018) * density(43) 
  pd(01,43) = pd(01,43) - rrt(018) * density(01) 
  pd(03,01) = pd(03,01) + rrt(018) * density(43) 
  pd(03,43) = pd(03,43) + rrt(018) * density(01) 
  pd(01,01) = pd(01,01) - rrt(019) * density(43) 
  pd(01,43) = pd(01,43) - rrt(019) * density(01) 
  pd(04,01) = pd(04,01) + rrt(019) * density(43) 
  pd(04,43) = pd(04,43) + rrt(019) * density(01) 
  pd(01,01) = pd(01,01) - rrt(020) * density(43) 
  pd(01,43) = pd(01,43) - rrt(020) * density(01) 
  pd(05,01) = pd(05,01) + rrt(020) * density(43) 
  pd(05,43) = pd(05,43) + rrt(020) * density(01) 
  pd(01,01) = pd(01,01) - rrt(021) * density(43) 
  pd(01,43) = pd(01,43) - rrt(021) * density(01) 
  pd(06,01) = pd(06,01) + rrt(021) * density(43) 
  pd(06,43) = pd(06,43) + rrt(021) * density(01) 
  pd(01,01) = pd(01,01) - rrt(022) * density(43) 
  pd(01,43) = pd(01,43) - rrt(022) * density(01) 
  pd(07,01) = pd(07,01) + rrt(022) * density(43) 
  pd(07,43) = pd(07,43) + rrt(022) * density(01) 
  pd(01,01) = pd(01,01) - rrt(023) * density(43) 
  pd(01,43) = pd(01,43) - rrt(023) * density(01) 
  pd(08,01) = pd(08,01) + rrt(023) * density(43) 
  pd(08,43) = pd(08,43) + rrt(023) * density(01) 
  pd(01,01) = pd(01,01) - rrt(024) * density(43) 
  pd(01,43) = pd(01,43) - rrt(024) * density(01) 
  pd(09,01) = pd(09,01) + rrt(024) * density(43) 
  pd(09,43) = pd(09,43) + rrt(024) * density(01) 
  pd(01,02) = pd(01,02) + rrt(025) * density(43) 
  pd(01,43) = pd(01,43) + rrt(025) * density(02) 
  pd(02,02) = pd(02,02) - rrt(025) * density(43) 
  pd(02,43) = pd(02,43) - rrt(025) * density(02) 
  pd(01,03) = pd(01,03) + rrt(026) * density(43) 
  pd(01,43) = pd(01,43) + rrt(026) * density(03) 
  pd(03,03) = pd(03,03) - rrt(026) * density(43) 
  pd(03,43) = pd(03,43) - rrt(026) * density(03) 
  pd(01,04) = pd(01,04) + rrt(027) * density(43) 
  pd(01,43) = pd(01,43) + rrt(027) * density(04) 
  pd(04,04) = pd(04,04) - rrt(027) * density(43) 
  pd(04,43) = pd(04,43) - rrt(027) * density(04) 
  pd(01,05) = pd(01,05) + rrt(028) * density(43) 
  pd(01,43) = pd(01,43) + rrt(028) * density(05) 
  pd(05,05) = pd(05,05) - rrt(028) * density(43) 
  pd(05,43) = pd(05,43) - rrt(028) * density(05) 
  pd(01,06) = pd(01,06) + rrt(029) * density(43) 
  pd(01,43) = pd(01,43) + rrt(029) * density(06) 
  pd(06,06) = pd(06,06) - rrt(029) * density(43) 
  pd(06,43) = pd(06,43) - rrt(029) * density(06) 
  pd(01,07) = pd(01,07) + rrt(030) * density(43) 
  pd(01,43) = pd(01,43) + rrt(030) * density(07) 
  pd(07,07) = pd(07,07) - rrt(030) * density(43) 
  pd(07,43) = pd(07,43) - rrt(030) * density(07) 
  pd(01,08) = pd(01,08) + rrt(031) * density(43) 
  pd(01,43) = pd(01,43) + rrt(031) * density(08) 
  pd(08,08) = pd(08,08) - rrt(031) * density(43) 
  pd(08,43) = pd(08,43) - rrt(031) * density(08) 
  pd(01,09) = pd(01,09) + rrt(032) * density(43) 
  pd(01,43) = pd(01,43) + rrt(032) * density(09) 
  pd(09,09) = pd(09,09) - rrt(032) * density(43) 
  pd(09,43) = pd(09,43) - rrt(032) * density(09) 
  pd(21,21) = pd(21,21) - rrt(033) * density(43) 
  pd(21,43) = pd(21,43) - rrt(033) * density(21) 
  pd(27,21) = pd(27,21) + rrt(033) * density(43) 
  pd(27,43) = pd(27,43) + rrt(033) * density(21) 
  pd(21,21) = pd(21,21) - rrt(034) * density(43) 
  pd(21,43) = pd(21,43) - rrt(034) * density(21) 
  pd(28,21) = pd(28,21) + rrt(034) * density(43) 
  pd(28,43) = pd(28,43) + rrt(034) * density(21) 
  pd(21,21) = pd(21,21) - rrt(035) * density(43) 
  pd(21,43) = pd(21,43) - rrt(035) * density(21) 
  pd(29,21) = pd(29,21) + rrt(035) * density(43) 
  pd(29,43) = pd(29,43) + rrt(035) * density(21) 
  pd(21,27) = pd(21,27) + rrt(036) * density(43) 
  pd(21,43) = pd(21,43) + rrt(036) * density(27) 
  pd(27,27) = pd(27,27) - rrt(036) * density(43) 
  pd(27,43) = pd(27,43) - rrt(036) * density(27) 
  pd(21,28) = pd(21,28) + rrt(037) * density(43) 
  pd(21,43) = pd(21,43) + rrt(037) * density(28) 
  pd(28,28) = pd(28,28) - rrt(037) * density(43) 
  pd(28,43) = pd(28,43) - rrt(037) * density(28) 
  pd(21,29) = pd(21,29) + rrt(038) * density(43) 
  pd(21,43) = pd(21,43) + rrt(038) * density(29) 
  pd(29,29) = pd(29,29) - rrt(038) * density(43) 
  pd(29,43) = pd(29,43) - rrt(038) * density(29) 
  pd(21,22) = pd(21,22) + rrt(039) * density(43) 
  pd(21,43) = pd(21,43) + rrt(039) * density(22) 
  pd(22,22) = pd(22,22) - rrt(039) * density(43) 
  pd(22,43) = pd(22,43) - rrt(039) * density(22) 
  pd(21,23) = pd(21,23) + rrt(040) * density(43) 
  pd(21,43) = pd(21,43) + rrt(040) * density(23) 
  pd(23,23) = pd(23,23) - rrt(040) * density(43) 
  pd(23,43) = pd(23,43) - rrt(040) * density(23) 
  pd(21,24) = pd(21,24) + rrt(041) * density(43) 
  pd(21,43) = pd(21,43) + rrt(041) * density(24) 
  pd(24,24) = pd(24,24) - rrt(041) * density(43) 
  pd(24,43) = pd(24,43) - rrt(041) * density(24) 
  pd(21,25) = pd(21,25) + rrt(042) * density(43) 
  pd(21,43) = pd(21,43) + rrt(042) * density(25) 
  pd(25,25) = pd(25,25) - rrt(042) * density(43) 
  pd(25,43) = pd(25,43) - rrt(042) * density(25) 
  pd(01,01) = pd(01,01) + rrt(043) * density(02) 
  pd(01,02) = pd(01,02) + rrt(043) * density(01) 
  pd(02,01) = pd(02,01) - rrt(043) * density(02) 
  pd(02,02) = pd(02,02) - rrt(043) * density(01) 
  pd(01,01) = pd(01,01) - rrt(044) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(044) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(045) * density(03) 
  pd(02,03) = pd(02,03) + rrt(045) * density(01) 
  pd(03,01) = pd(03,01) - rrt(045) * density(03) 
  pd(03,03) = pd(03,03) - rrt(045) * density(01) 
  pd(02,01) = pd(02,01) - rrt(046) * density(02) 
  pd(02,02) = pd(02,02) - rrt(046) * density(01) 
  pd(03,01) = pd(03,01) + rrt(046) * density(02) 
  pd(03,02) = pd(03,02) + rrt(046) * density(01) 
  pd(03,01) = pd(03,01) + rrt(047) * density(04) 
  pd(03,04) = pd(03,04) + rrt(047) * density(01) 
  pd(04,01) = pd(04,01) - rrt(047) * density(04) 
  pd(04,04) = pd(04,04) - rrt(047) * density(01) 
  pd(03,01) = pd(03,01) - rrt(048) * density(03) 
  pd(03,03) = pd(03,03) - rrt(048) * density(01) 
  pd(04,01) = pd(04,01) + rrt(048) * density(03) 
  pd(04,03) = pd(04,03) + rrt(048) * density(01) 
  pd(04,01) = pd(04,01) + rrt(049) * density(05) 
  pd(04,05) = pd(04,05) + rrt(049) * density(01) 
  pd(05,01) = pd(05,01) - rrt(049) * density(05) 
  pd(05,05) = pd(05,05) - rrt(049) * density(01) 
  pd(04,01) = pd(04,01) - rrt(050) * density(04) 
  pd(04,04) = pd(04,04) - rrt(050) * density(01) 
  pd(05,01) = pd(05,01) + rrt(050) * density(04) 
  pd(05,04) = pd(05,04) + rrt(050) * density(01) 
  pd(05,01) = pd(05,01) + rrt(051) * density(06) 
  pd(05,06) = pd(05,06) + rrt(051) * density(01) 
  pd(06,01) = pd(06,01) - rrt(051) * density(06) 
  pd(06,06) = pd(06,06) - rrt(051) * density(01) 
  pd(05,01) = pd(05,01) - rrt(052) * density(05) 
  pd(05,05) = pd(05,05) - rrt(052) * density(01) 
  pd(06,01) = pd(06,01) + rrt(052) * density(05) 
  pd(06,05) = pd(06,05) + rrt(052) * density(01) 
  pd(06,01) = pd(06,01) + rrt(053) * density(07) 
  pd(06,07) = pd(06,07) + rrt(053) * density(01) 
  pd(07,01) = pd(07,01) - rrt(053) * density(07) 
  pd(07,07) = pd(07,07) - rrt(053) * density(01) 
  pd(06,01) = pd(06,01) - rrt(054) * density(06) 
  pd(06,06) = pd(06,06) - rrt(054) * density(01) 
  pd(07,01) = pd(07,01) + rrt(054) * density(06) 
  pd(07,06) = pd(07,06) + rrt(054) * density(01) 
  pd(07,01) = pd(07,01) + rrt(055) * density(08) 
  pd(07,08) = pd(07,08) + rrt(055) * density(01) 
  pd(08,01) = pd(08,01) - rrt(055) * density(08) 
  pd(08,08) = pd(08,08) - rrt(055) * density(01) 
  pd(07,01) = pd(07,01) - rrt(056) * density(07) 
  pd(07,07) = pd(07,07) - rrt(056) * density(01) 
  pd(08,01) = pd(08,01) + rrt(056) * density(07) 
  pd(08,07) = pd(08,07) + rrt(056) * density(01) 
  pd(08,01) = pd(08,01) + rrt(057) * density(09) 
  pd(08,09) = pd(08,09) + rrt(057) * density(01) 
  pd(09,01) = pd(09,01) - rrt(057) * density(09) 
  pd(09,09) = pd(09,09) - rrt(057) * density(01) 
  pd(08,01) = pd(08,01) - rrt(058) * density(08) 
  pd(08,08) = pd(08,08) - rrt(058) * density(01) 
  pd(09,01) = pd(09,01) + rrt(058) * density(08) 
  pd(09,08) = pd(09,08) + rrt(058) * density(01) 
  pd(01,02) = pd(01,02) + rrt(059) * density(21) 
  pd(01,21) = pd(01,21) + rrt(059) * density(02) 
  pd(02,02) = pd(02,02) - rrt(059) * density(21) 
  pd(02,21) = pd(02,21) - rrt(059) * density(02) 
  pd(01,01) = pd(01,01) - rrt(060) * density(21) 
  pd(01,21) = pd(01,21) - rrt(060) * density(01) 
  pd(02,01) = pd(02,01) + rrt(060) * density(21) 
  pd(02,21) = pd(02,21) + rrt(060) * density(01) 
  pd(02,03) = pd(02,03) + rrt(061) * density(21) 
  pd(02,21) = pd(02,21) + rrt(061) * density(03) 
  pd(03,03) = pd(03,03) - rrt(061) * density(21) 
  pd(03,21) = pd(03,21) - rrt(061) * density(03) 
  pd(02,02) = pd(02,02) - rrt(062) * density(21) 
  pd(02,21) = pd(02,21) - rrt(062) * density(02) 
  pd(03,02) = pd(03,02) + rrt(062) * density(21) 
  pd(03,21) = pd(03,21) + rrt(062) * density(02) 
  pd(03,04) = pd(03,04) + rrt(063) * density(21) 
  pd(03,21) = pd(03,21) + rrt(063) * density(04) 
  pd(04,04) = pd(04,04) - rrt(063) * density(21) 
  pd(04,21) = pd(04,21) - rrt(063) * density(04) 
  pd(03,03) = pd(03,03) - rrt(064) * density(21) 
  pd(03,21) = pd(03,21) - rrt(064) * density(03) 
  pd(04,03) = pd(04,03) + rrt(064) * density(21) 
  pd(04,21) = pd(04,21) + rrt(064) * density(03) 
  pd(04,05) = pd(04,05) + rrt(065) * density(21) 
  pd(04,21) = pd(04,21) + rrt(065) * density(05) 
  pd(05,05) = pd(05,05) - rrt(065) * density(21) 
  pd(05,21) = pd(05,21) - rrt(065) * density(05) 
  pd(04,04) = pd(04,04) - rrt(066) * density(21) 
  pd(04,21) = pd(04,21) - rrt(066) * density(04) 
  pd(05,04) = pd(05,04) + rrt(066) * density(21) 
  pd(05,21) = pd(05,21) + rrt(066) * density(04) 
  pd(05,06) = pd(05,06) + rrt(067) * density(21) 
  pd(05,21) = pd(05,21) + rrt(067) * density(06) 
  pd(06,06) = pd(06,06) - rrt(067) * density(21) 
  pd(06,21) = pd(06,21) - rrt(067) * density(06) 
  pd(05,05) = pd(05,05) - rrt(068) * density(21) 
  pd(05,21) = pd(05,21) - rrt(068) * density(05) 
  pd(06,05) = pd(06,05) + rrt(068) * density(21) 
  pd(06,21) = pd(06,21) + rrt(068) * density(05) 
  pd(06,07) = pd(06,07) + rrt(069) * density(21) 
  pd(06,21) = pd(06,21) + rrt(069) * density(07) 
  pd(07,07) = pd(07,07) - rrt(069) * density(21) 
  pd(07,21) = pd(07,21) - rrt(069) * density(07) 
  pd(06,06) = pd(06,06) - rrt(070) * density(21) 
  pd(06,21) = pd(06,21) - rrt(070) * density(06) 
  pd(07,06) = pd(07,06) + rrt(070) * density(21) 
  pd(07,21) = pd(07,21) + rrt(070) * density(06) 
  pd(07,08) = pd(07,08) + rrt(071) * density(21) 
  pd(07,21) = pd(07,21) + rrt(071) * density(08) 
  pd(08,08) = pd(08,08) - rrt(071) * density(21) 
  pd(08,21) = pd(08,21) - rrt(071) * density(08) 
  pd(07,07) = pd(07,07) - rrt(072) * density(21) 
  pd(07,21) = pd(07,21) - rrt(072) * density(07) 
  pd(08,07) = pd(08,07) + rrt(072) * density(21) 
  pd(08,21) = pd(08,21) + rrt(072) * density(07) 
  pd(08,09) = pd(08,09) + rrt(073) * density(21) 
  pd(08,21) = pd(08,21) + rrt(073) * density(09) 
  pd(09,09) = pd(09,09) - rrt(073) * density(21) 
  pd(09,21) = pd(09,21) - rrt(073) * density(09) 
  pd(08,08) = pd(08,08) - rrt(074) * density(21) 
  pd(08,21) = pd(08,21) - rrt(074) * density(08) 
  pd(09,08) = pd(09,08) + rrt(074) * density(21) 
  pd(09,21) = pd(09,21) + rrt(074) * density(08) 
  pd(01,02) = pd(01,02) + rrt(075) * density(14) 
  pd(01,14) = pd(01,14) + rrt(075) * density(02) 
  pd(02,02) = pd(02,02) - rrt(075) * density(14) 
  pd(02,14) = pd(02,14) - rrt(075) * density(02) 
  pd(01,01) = pd(01,01) - rrt(076) * density(14) 
  pd(01,14) = pd(01,14) - rrt(076) * density(01) 
  pd(02,01) = pd(02,01) + rrt(076) * density(14) 
  pd(02,14) = pd(02,14) + rrt(076) * density(01) 
  pd(02,03) = pd(02,03) + rrt(077) * density(14) 
  pd(02,14) = pd(02,14) + rrt(077) * density(03) 
  pd(03,03) = pd(03,03) - rrt(077) * density(14) 
  pd(03,14) = pd(03,14) - rrt(077) * density(03) 
  pd(02,02) = pd(02,02) - rrt(078) * density(14) 
  pd(02,14) = pd(02,14) - rrt(078) * density(02) 
  pd(03,02) = pd(03,02) + rrt(078) * density(14) 
  pd(03,14) = pd(03,14) + rrt(078) * density(02) 
  pd(03,04) = pd(03,04) + rrt(079) * density(14) 
  pd(03,14) = pd(03,14) + rrt(079) * density(04) 
  pd(04,04) = pd(04,04) - rrt(079) * density(14) 
  pd(04,14) = pd(04,14) - rrt(079) * density(04) 
  pd(03,03) = pd(03,03) - rrt(080) * density(14) 
  pd(03,14) = pd(03,14) - rrt(080) * density(03) 
  pd(04,03) = pd(04,03) + rrt(080) * density(14) 
  pd(04,14) = pd(04,14) + rrt(080) * density(03) 
  pd(04,05) = pd(04,05) + rrt(081) * density(14) 
  pd(04,14) = pd(04,14) + rrt(081) * density(05) 
  pd(05,05) = pd(05,05) - rrt(081) * density(14) 
  pd(05,14) = pd(05,14) - rrt(081) * density(05) 
  pd(04,04) = pd(04,04) - rrt(082) * density(14) 
  pd(04,14) = pd(04,14) - rrt(082) * density(04) 
  pd(05,04) = pd(05,04) + rrt(082) * density(14) 
  pd(05,14) = pd(05,14) + rrt(082) * density(04) 
  pd(05,06) = pd(05,06) + rrt(083) * density(14) 
  pd(05,14) = pd(05,14) + rrt(083) * density(06) 
  pd(06,06) = pd(06,06) - rrt(083) * density(14) 
  pd(06,14) = pd(06,14) - rrt(083) * density(06) 
  pd(05,05) = pd(05,05) - rrt(084) * density(14) 
  pd(05,14) = pd(05,14) - rrt(084) * density(05) 
  pd(06,05) = pd(06,05) + rrt(084) * density(14) 
  pd(06,14) = pd(06,14) + rrt(084) * density(05) 
  pd(06,07) = pd(06,07) + rrt(085) * density(14) 
  pd(06,14) = pd(06,14) + rrt(085) * density(07) 
  pd(07,07) = pd(07,07) - rrt(085) * density(14) 
  pd(07,14) = pd(07,14) - rrt(085) * density(07) 
  pd(06,06) = pd(06,06) - rrt(086) * density(14) 
  pd(06,14) = pd(06,14) - rrt(086) * density(06) 
  pd(07,06) = pd(07,06) + rrt(086) * density(14) 
  pd(07,14) = pd(07,14) + rrt(086) * density(06) 
  pd(07,08) = pd(07,08) + rrt(087) * density(14) 
  pd(07,14) = pd(07,14) + rrt(087) * density(08) 
  pd(08,08) = pd(08,08) - rrt(087) * density(14) 
  pd(08,14) = pd(08,14) - rrt(087) * density(08) 
  pd(07,07) = pd(07,07) - rrt(088) * density(14) 
  pd(07,14) = pd(07,14) - rrt(088) * density(07) 
  pd(08,07) = pd(08,07) + rrt(088) * density(14) 
  pd(08,14) = pd(08,14) + rrt(088) * density(07) 
  pd(08,09) = pd(08,09) + rrt(089) * density(14) 
  pd(08,14) = pd(08,14) + rrt(089) * density(09) 
  pd(09,09) = pd(09,09) - rrt(089) * density(14) 
  pd(09,14) = pd(09,14) - rrt(089) * density(09) 
  pd(08,08) = pd(08,08) - rrt(090) * density(14) 
  pd(08,14) = pd(08,14) - rrt(090) * density(08) 
  pd(09,08) = pd(09,08) + rrt(090) * density(14) 
  pd(09,14) = pd(09,14) + rrt(090) * density(08) 
  pd(01,02) = pd(01,02) + rrt(091) * density(30) 
  pd(01,30) = pd(01,30) + rrt(091) * density(02) 
  pd(02,02) = pd(02,02) - rrt(091) * density(30) 
  pd(02,30) = pd(02,30) - rrt(091) * density(02) 
  pd(01,01) = pd(01,01) - rrt(092) * density(30) 
  pd(01,30) = pd(01,30) - rrt(092) * density(01) 
  pd(02,01) = pd(02,01) + rrt(092) * density(30) 
  pd(02,30) = pd(02,30) + rrt(092) * density(01) 
  pd(02,03) = pd(02,03) + rrt(093) * density(30) 
  pd(02,30) = pd(02,30) + rrt(093) * density(03) 
  pd(03,03) = pd(03,03) - rrt(093) * density(30) 
  pd(03,30) = pd(03,30) - rrt(093) * density(03) 
  pd(02,02) = pd(02,02) - rrt(094) * density(30) 
  pd(02,30) = pd(02,30) - rrt(094) * density(02) 
  pd(03,02) = pd(03,02) + rrt(094) * density(30) 
  pd(03,30) = pd(03,30) + rrt(094) * density(02) 
  pd(03,04) = pd(03,04) + rrt(095) * density(30) 
  pd(03,30) = pd(03,30) + rrt(095) * density(04) 
  pd(04,04) = pd(04,04) - rrt(095) * density(30) 
  pd(04,30) = pd(04,30) - rrt(095) * density(04) 
  pd(03,03) = pd(03,03) - rrt(096) * density(30) 
  pd(03,30) = pd(03,30) - rrt(096) * density(03) 
  pd(04,03) = pd(04,03) + rrt(096) * density(30) 
  pd(04,30) = pd(04,30) + rrt(096) * density(03) 
  pd(04,05) = pd(04,05) + rrt(097) * density(30) 
  pd(04,30) = pd(04,30) + rrt(097) * density(05) 
  pd(05,05) = pd(05,05) - rrt(097) * density(30) 
  pd(05,30) = pd(05,30) - rrt(097) * density(05) 
  pd(04,04) = pd(04,04) - rrt(098) * density(30) 
  pd(04,30) = pd(04,30) - rrt(098) * density(04) 
  pd(05,04) = pd(05,04) + rrt(098) * density(30) 
  pd(05,30) = pd(05,30) + rrt(098) * density(04) 
  pd(05,06) = pd(05,06) + rrt(099) * density(30) 
  pd(05,30) = pd(05,30) + rrt(099) * density(06) 
  pd(06,06) = pd(06,06) - rrt(099) * density(30) 
  pd(06,30) = pd(06,30) - rrt(099) * density(06) 
  pd(05,05) = pd(05,05) - rrt(100) * density(30) 
  pd(05,30) = pd(05,30) - rrt(100) * density(05) 
  pd(06,05) = pd(06,05) + rrt(100) * density(30) 
  pd(06,30) = pd(06,30) + rrt(100) * density(05) 
  pd(06,07) = pd(06,07) + rrt(101) * density(30) 
  pd(06,30) = pd(06,30) + rrt(101) * density(07) 
  pd(07,07) = pd(07,07) - rrt(101) * density(30) 
  pd(07,30) = pd(07,30) - rrt(101) * density(07) 
  pd(06,06) = pd(06,06) - rrt(102) * density(30) 
  pd(06,30) = pd(06,30) - rrt(102) * density(06) 
  pd(07,06) = pd(07,06) + rrt(102) * density(30) 
  pd(07,30) = pd(07,30) + rrt(102) * density(06) 
  pd(07,08) = pd(07,08) + rrt(103) * density(30) 
  pd(07,30) = pd(07,30) + rrt(103) * density(08) 
  pd(08,08) = pd(08,08) - rrt(103) * density(30) 
  pd(08,30) = pd(08,30) - rrt(103) * density(08) 
  pd(07,07) = pd(07,07) - rrt(104) * density(30) 
  pd(07,30) = pd(07,30) - rrt(104) * density(07) 
  pd(08,07) = pd(08,07) + rrt(104) * density(30) 
  pd(08,30) = pd(08,30) + rrt(104) * density(07) 
  pd(08,09) = pd(08,09) + rrt(105) * density(30) 
  pd(08,30) = pd(08,30) + rrt(105) * density(09) 
  pd(09,09) = pd(09,09) - rrt(105) * density(30) 
  pd(09,30) = pd(09,30) - rrt(105) * density(09) 
  pd(08,08) = pd(08,08) - rrt(106) * density(30) 
  pd(08,30) = pd(08,30) - rrt(106) * density(08) 
  pd(09,08) = pd(09,08) + rrt(106) * density(30) 
  pd(09,30) = pd(09,30) + rrt(106) * density(08) 
  pd(21,21) = pd(21,21) + rrt(107) * density(27) 
  pd(21,27) = pd(21,27) + rrt(107) * density(21) 
  pd(27,21) = pd(27,21) - rrt(107) * density(27) 
  pd(27,27) = pd(27,27) - rrt(107) * density(21) 
  pd(21,21) = pd(21,21) - rrt(108) * density(21) * 2.0d0
  pd(27,21) = pd(27,21) + rrt(108) * density(21) * 2.0d0
  pd(27,21) = pd(27,21) + rrt(109) * density(28) 
  pd(27,28) = pd(27,28) + rrt(109) * density(21) 
  pd(28,21) = pd(28,21) - rrt(109) * density(28) 
  pd(28,28) = pd(28,28) - rrt(109) * density(21) 
  pd(27,21) = pd(27,21) - rrt(110) * density(27) 
  pd(27,27) = pd(27,27) - rrt(110) * density(21) 
  pd(28,21) = pd(28,21) + rrt(110) * density(27) 
  pd(28,27) = pd(28,27) + rrt(110) * density(21) 
  pd(28,21) = pd(28,21) + rrt(111) * density(29) 
  pd(28,29) = pd(28,29) + rrt(111) * density(21) 
  pd(29,21) = pd(29,21) - rrt(111) * density(29) 
  pd(29,29) = pd(29,29) - rrt(111) * density(21) 
  pd(28,21) = pd(28,21) - rrt(112) * density(28) 
  pd(28,28) = pd(28,28) - rrt(112) * density(21) 
  pd(29,21) = pd(29,21) + rrt(112) * density(28) 
  pd(29,28) = pd(29,28) + rrt(112) * density(21) 
  pd(21,27) = pd(21,27) + rrt(113) * density(30) 
  pd(21,30) = pd(21,30) + rrt(113) * density(27) 
  pd(27,27) = pd(27,27) - rrt(113) * density(30) 
  pd(27,30) = pd(27,30) - rrt(113) * density(27) 
  pd(21,21) = pd(21,21) - rrt(114) * density(30) 
  pd(21,30) = pd(21,30) - rrt(114) * density(21) 
  pd(27,21) = pd(27,21) + rrt(114) * density(30) 
  pd(27,30) = pd(27,30) + rrt(114) * density(21) 
  pd(27,28) = pd(27,28) + rrt(115) * density(30) 
  pd(27,30) = pd(27,30) + rrt(115) * density(28) 
  pd(28,28) = pd(28,28) - rrt(115) * density(30) 
  pd(28,30) = pd(28,30) - rrt(115) * density(28) 
  pd(27,27) = pd(27,27) - rrt(116) * density(30) 
  pd(27,30) = pd(27,30) - rrt(116) * density(27) 
  pd(28,27) = pd(28,27) + rrt(116) * density(30) 
  pd(28,30) = pd(28,30) + rrt(116) * density(27) 
  pd(28,29) = pd(28,29) + rrt(117) * density(30) 
  pd(28,30) = pd(28,30) + rrt(117) * density(29) 
  pd(29,29) = pd(29,29) - rrt(117) * density(30) 
  pd(29,30) = pd(29,30) - rrt(117) * density(29) 
  pd(28,28) = pd(28,28) - rrt(118) * density(30) 
  pd(28,30) = pd(28,30) - rrt(118) * density(28) 
  pd(29,28) = pd(29,28) + rrt(118) * density(30) 
  pd(29,30) = pd(29,30) + rrt(118) * density(28) 
  pd(01,02) = pd(01,02) + rrt(119) * density(02) * 2.0d0
  pd(02,02) = pd(02,02) - rrt(119) * density(02) * 4.0d0
  pd(03,02) = pd(03,02) + rrt(119) * density(02) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(120) * density(03) 
  pd(01,03) = pd(01,03) - rrt(120) * density(01) 
  pd(02,01) = pd(02,01) + rrt(120) * density(03) * 2.0d0
  pd(02,03) = pd(02,03) + rrt(120) * density(01) * 2.0d0
  pd(03,01) = pd(03,01) - rrt(120) * density(03) 
  pd(03,03) = pd(03,03) - rrt(120) * density(01) 
  pd(01,02) = pd(01,02) + rrt(121) * density(03) 
  pd(01,03) = pd(01,03) + rrt(121) * density(02) 
  pd(02,02) = pd(02,02) - rrt(121) * density(03) 
  pd(02,03) = pd(02,03) - rrt(121) * density(02) 
  pd(03,02) = pd(03,02) - rrt(121) * density(03) 
  pd(03,03) = pd(03,03) - rrt(121) * density(02) 
  pd(04,02) = pd(04,02) + rrt(121) * density(03) 
  pd(04,03) = pd(04,03) + rrt(121) * density(02) 
  pd(01,01) = pd(01,01) - rrt(122) * density(04) 
  pd(01,04) = pd(01,04) - rrt(122) * density(01) 
  pd(02,01) = pd(02,01) + rrt(122) * density(04) 
  pd(02,04) = pd(02,04) + rrt(122) * density(01) 
  pd(03,01) = pd(03,01) + rrt(122) * density(04) 
  pd(03,04) = pd(03,04) + rrt(122) * density(01) 
  pd(04,01) = pd(04,01) - rrt(122) * density(04) 
  pd(04,04) = pd(04,04) - rrt(122) * density(01) 
  pd(01,02) = pd(01,02) + rrt(123) * density(04) 
  pd(01,04) = pd(01,04) + rrt(123) * density(02) 
  pd(02,02) = pd(02,02) - rrt(123) * density(04) 
  pd(02,04) = pd(02,04) - rrt(123) * density(02) 
  pd(04,02) = pd(04,02) - rrt(123) * density(04) 
  pd(04,04) = pd(04,04) - rrt(123) * density(02) 
  pd(05,02) = pd(05,02) + rrt(123) * density(04) 
  pd(05,04) = pd(05,04) + rrt(123) * density(02) 
  pd(01,01) = pd(01,01) - rrt(124) * density(05) 
  pd(01,05) = pd(01,05) - rrt(124) * density(01) 
  pd(02,01) = pd(02,01) + rrt(124) * density(05) 
  pd(02,05) = pd(02,05) + rrt(124) * density(01) 
  pd(04,01) = pd(04,01) + rrt(124) * density(05) 
  pd(04,05) = pd(04,05) + rrt(124) * density(01) 
  pd(05,01) = pd(05,01) - rrt(124) * density(05) 
  pd(05,05) = pd(05,05) - rrt(124) * density(01) 
  pd(01,02) = pd(01,02) + rrt(125) * density(05) 
  pd(01,05) = pd(01,05) + rrt(125) * density(02) 
  pd(02,02) = pd(02,02) - rrt(125) * density(05) 
  pd(02,05) = pd(02,05) - rrt(125) * density(02) 
  pd(05,02) = pd(05,02) - rrt(125) * density(05) 
  pd(05,05) = pd(05,05) - rrt(125) * density(02) 
  pd(06,02) = pd(06,02) + rrt(125) * density(05) 
  pd(06,05) = pd(06,05) + rrt(125) * density(02) 
  pd(01,01) = pd(01,01) - rrt(126) * density(06) 
  pd(01,06) = pd(01,06) - rrt(126) * density(01) 
  pd(02,01) = pd(02,01) + rrt(126) * density(06) 
  pd(02,06) = pd(02,06) + rrt(126) * density(01) 
  pd(05,01) = pd(05,01) + rrt(126) * density(06) 
  pd(05,06) = pd(05,06) + rrt(126) * density(01) 
  pd(06,01) = pd(06,01) - rrt(126) * density(06) 
  pd(06,06) = pd(06,06) - rrt(126) * density(01) 
  pd(01,02) = pd(01,02) + rrt(127) * density(06) 
  pd(01,06) = pd(01,06) + rrt(127) * density(02) 
  pd(02,02) = pd(02,02) - rrt(127) * density(06) 
  pd(02,06) = pd(02,06) - rrt(127) * density(02) 
  pd(06,02) = pd(06,02) - rrt(127) * density(06) 
  pd(06,06) = pd(06,06) - rrt(127) * density(02) 
  pd(07,02) = pd(07,02) + rrt(127) * density(06) 
  pd(07,06) = pd(07,06) + rrt(127) * density(02) 
  pd(01,01) = pd(01,01) - rrt(128) * density(07) 
  pd(01,07) = pd(01,07) - rrt(128) * density(01) 
  pd(02,01) = pd(02,01) + rrt(128) * density(07) 
  pd(02,07) = pd(02,07) + rrt(128) * density(01) 
  pd(06,01) = pd(06,01) + rrt(128) * density(07) 
  pd(06,07) = pd(06,07) + rrt(128) * density(01) 
  pd(07,01) = pd(07,01) - rrt(128) * density(07) 
  pd(07,07) = pd(07,07) - rrt(128) * density(01) 
  pd(01,02) = pd(01,02) + rrt(129) * density(07) 
  pd(01,07) = pd(01,07) + rrt(129) * density(02) 
  pd(02,02) = pd(02,02) - rrt(129) * density(07) 
  pd(02,07) = pd(02,07) - rrt(129) * density(02) 
  pd(07,02) = pd(07,02) - rrt(129) * density(07) 
  pd(07,07) = pd(07,07) - rrt(129) * density(02) 
  pd(08,02) = pd(08,02) + rrt(129) * density(07) 
  pd(08,07) = pd(08,07) + rrt(129) * density(02) 
  pd(01,01) = pd(01,01) - rrt(130) * density(08) 
  pd(01,08) = pd(01,08) - rrt(130) * density(01) 
  pd(02,01) = pd(02,01) + rrt(130) * density(08) 
  pd(02,08) = pd(02,08) + rrt(130) * density(01) 
  pd(07,01) = pd(07,01) + rrt(130) * density(08) 
  pd(07,08) = pd(07,08) + rrt(130) * density(01) 
  pd(08,01) = pd(08,01) - rrt(130) * density(08) 
  pd(08,08) = pd(08,08) - rrt(130) * density(01) 
  pd(01,02) = pd(01,02) + rrt(131) * density(08) 
  pd(01,08) = pd(01,08) + rrt(131) * density(02) 
  pd(02,02) = pd(02,02) - rrt(131) * density(08) 
  pd(02,08) = pd(02,08) - rrt(131) * density(02) 
  pd(08,02) = pd(08,02) - rrt(131) * density(08) 
  pd(08,08) = pd(08,08) - rrt(131) * density(02) 
  pd(09,02) = pd(09,02) + rrt(131) * density(08) 
  pd(09,08) = pd(09,08) + rrt(131) * density(02) 
  pd(01,01) = pd(01,01) - rrt(132) * density(09) 
  pd(01,09) = pd(01,09) - rrt(132) * density(01) 
  pd(02,01) = pd(02,01) + rrt(132) * density(09) 
  pd(02,09) = pd(02,09) + rrt(132) * density(01) 
  pd(08,01) = pd(08,01) + rrt(132) * density(09) 
  pd(08,09) = pd(08,09) + rrt(132) * density(01) 
  pd(09,01) = pd(09,01) - rrt(132) * density(09) 
  pd(09,09) = pd(09,09) - rrt(132) * density(01) 
  pd(02,03) = pd(02,03) + rrt(133) * density(03) * 2.0d0
  pd(03,03) = pd(03,03) - rrt(133) * density(03) * 4.0d0
  pd(04,03) = pd(04,03) + rrt(133) * density(03) * 2.0d0
  pd(02,02) = pd(02,02) - rrt(134) * density(04) 
  pd(02,04) = pd(02,04) - rrt(134) * density(02) 
  pd(03,02) = pd(03,02) + rrt(134) * density(04) * 2.0d0
  pd(03,04) = pd(03,04) + rrt(134) * density(02) * 2.0d0
  pd(04,02) = pd(04,02) - rrt(134) * density(04) 
  pd(04,04) = pd(04,04) - rrt(134) * density(02) 
  pd(02,03) = pd(02,03) + rrt(135) * density(04) 
  pd(02,04) = pd(02,04) + rrt(135) * density(03) 
  pd(03,03) = pd(03,03) - rrt(135) * density(04) 
  pd(03,04) = pd(03,04) - rrt(135) * density(03) 
  pd(04,03) = pd(04,03) - rrt(135) * density(04) 
  pd(04,04) = pd(04,04) - rrt(135) * density(03) 
  pd(05,03) = pd(05,03) + rrt(135) * density(04) 
  pd(05,04) = pd(05,04) + rrt(135) * density(03) 
  pd(02,02) = pd(02,02) - rrt(136) * density(05) 
  pd(02,05) = pd(02,05) - rrt(136) * density(02) 
  pd(03,02) = pd(03,02) + rrt(136) * density(05) 
  pd(03,05) = pd(03,05) + rrt(136) * density(02) 
  pd(04,02) = pd(04,02) + rrt(136) * density(05) 
  pd(04,05) = pd(04,05) + rrt(136) * density(02) 
  pd(05,02) = pd(05,02) - rrt(136) * density(05) 
  pd(05,05) = pd(05,05) - rrt(136) * density(02) 
  pd(02,03) = pd(02,03) + rrt(137) * density(05) 
  pd(02,05) = pd(02,05) + rrt(137) * density(03) 
  pd(03,03) = pd(03,03) - rrt(137) * density(05) 
  pd(03,05) = pd(03,05) - rrt(137) * density(03) 
  pd(05,03) = pd(05,03) - rrt(137) * density(05) 
  pd(05,05) = pd(05,05) - rrt(137) * density(03) 
  pd(06,03) = pd(06,03) + rrt(137) * density(05) 
  pd(06,05) = pd(06,05) + rrt(137) * density(03) 
  pd(02,02) = pd(02,02) - rrt(138) * density(06) 
  pd(02,06) = pd(02,06) - rrt(138) * density(02) 
  pd(03,02) = pd(03,02) + rrt(138) * density(06) 
  pd(03,06) = pd(03,06) + rrt(138) * density(02) 
  pd(05,02) = pd(05,02) + rrt(138) * density(06) 
  pd(05,06) = pd(05,06) + rrt(138) * density(02) 
  pd(06,02) = pd(06,02) - rrt(138) * density(06) 
  pd(06,06) = pd(06,06) - rrt(138) * density(02) 
  pd(02,03) = pd(02,03) + rrt(139) * density(06) 
  pd(02,06) = pd(02,06) + rrt(139) * density(03) 
  pd(03,03) = pd(03,03) - rrt(139) * density(06) 
  pd(03,06) = pd(03,06) - rrt(139) * density(03) 
  pd(06,03) = pd(06,03) - rrt(139) * density(06) 
  pd(06,06) = pd(06,06) - rrt(139) * density(03) 
  pd(07,03) = pd(07,03) + rrt(139) * density(06) 
  pd(07,06) = pd(07,06) + rrt(139) * density(03) 
  pd(02,02) = pd(02,02) - rrt(140) * density(07) 
  pd(02,07) = pd(02,07) - rrt(140) * density(02) 
  pd(03,02) = pd(03,02) + rrt(140) * density(07) 
  pd(03,07) = pd(03,07) + rrt(140) * density(02) 
  pd(06,02) = pd(06,02) + rrt(140) * density(07) 
  pd(06,07) = pd(06,07) + rrt(140) * density(02) 
  pd(07,02) = pd(07,02) - rrt(140) * density(07) 
  pd(07,07) = pd(07,07) - rrt(140) * density(02) 
  pd(02,03) = pd(02,03) + rrt(141) * density(07) 
  pd(02,07) = pd(02,07) + rrt(141) * density(03) 
  pd(03,03) = pd(03,03) - rrt(141) * density(07) 
  pd(03,07) = pd(03,07) - rrt(141) * density(03) 
  pd(07,03) = pd(07,03) - rrt(141) * density(07) 
  pd(07,07) = pd(07,07) - rrt(141) * density(03) 
  pd(08,03) = pd(08,03) + rrt(141) * density(07) 
  pd(08,07) = pd(08,07) + rrt(141) * density(03) 
  pd(02,02) = pd(02,02) - rrt(142) * density(08) 
  pd(02,08) = pd(02,08) - rrt(142) * density(02) 
  pd(03,02) = pd(03,02) + rrt(142) * density(08) 
  pd(03,08) = pd(03,08) + rrt(142) * density(02) 
  pd(07,02) = pd(07,02) + rrt(142) * density(08) 
  pd(07,08) = pd(07,08) + rrt(142) * density(02) 
  pd(08,02) = pd(08,02) - rrt(142) * density(08) 
  pd(08,08) = pd(08,08) - rrt(142) * density(02) 
  pd(02,03) = pd(02,03) + rrt(143) * density(08) 
  pd(02,08) = pd(02,08) + rrt(143) * density(03) 
  pd(03,03) = pd(03,03) - rrt(143) * density(08) 
  pd(03,08) = pd(03,08) - rrt(143) * density(03) 
  pd(08,03) = pd(08,03) - rrt(143) * density(08) 
  pd(08,08) = pd(08,08) - rrt(143) * density(03) 
  pd(09,03) = pd(09,03) + rrt(143) * density(08) 
  pd(09,08) = pd(09,08) + rrt(143) * density(03) 
  pd(02,02) = pd(02,02) - rrt(144) * density(09) 
  pd(02,09) = pd(02,09) - rrt(144) * density(02) 
  pd(03,02) = pd(03,02) + rrt(144) * density(09) 
  pd(03,09) = pd(03,09) + rrt(144) * density(02) 
  pd(08,02) = pd(08,02) + rrt(144) * density(09) 
  pd(08,09) = pd(08,09) + rrt(144) * density(02) 
  pd(09,02) = pd(09,02) - rrt(144) * density(09) 
  pd(09,09) = pd(09,09) - rrt(144) * density(02) 
  pd(03,04) = pd(03,04) + rrt(145) * density(04) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(145) * density(04) * 4.0d0
  pd(05,04) = pd(05,04) + rrt(145) * density(04) * 2.0d0
  pd(03,03) = pd(03,03) - rrt(146) * density(05) 
  pd(03,05) = pd(03,05) - rrt(146) * density(03) 
  pd(04,03) = pd(04,03) + rrt(146) * density(05) * 2.0d0
  pd(04,05) = pd(04,05) + rrt(146) * density(03) * 2.0d0
  pd(05,03) = pd(05,03) - rrt(146) * density(05) 
  pd(05,05) = pd(05,05) - rrt(146) * density(03) 
  pd(03,04) = pd(03,04) + rrt(147) * density(05) 
  pd(03,05) = pd(03,05) + rrt(147) * density(04) 
  pd(04,04) = pd(04,04) - rrt(147) * density(05) 
  pd(04,05) = pd(04,05) - rrt(147) * density(04) 
  pd(05,04) = pd(05,04) - rrt(147) * density(05) 
  pd(05,05) = pd(05,05) - rrt(147) * density(04) 
  pd(06,04) = pd(06,04) + rrt(147) * density(05) 
  pd(06,05) = pd(06,05) + rrt(147) * density(04) 
  pd(03,03) = pd(03,03) - rrt(148) * density(06) 
  pd(03,06) = pd(03,06) - rrt(148) * density(03) 
  pd(04,03) = pd(04,03) + rrt(148) * density(06) 
  pd(04,06) = pd(04,06) + rrt(148) * density(03) 
  pd(05,03) = pd(05,03) + rrt(148) * density(06) 
  pd(05,06) = pd(05,06) + rrt(148) * density(03) 
  pd(06,03) = pd(06,03) - rrt(148) * density(06) 
  pd(06,06) = pd(06,06) - rrt(148) * density(03) 
  pd(03,04) = pd(03,04) + rrt(149) * density(06) 
  pd(03,06) = pd(03,06) + rrt(149) * density(04) 
  pd(04,04) = pd(04,04) - rrt(149) * density(06) 
  pd(04,06) = pd(04,06) - rrt(149) * density(04) 
  pd(06,04) = pd(06,04) - rrt(149) * density(06) 
  pd(06,06) = pd(06,06) - rrt(149) * density(04) 
  pd(07,04) = pd(07,04) + rrt(149) * density(06) 
  pd(07,06) = pd(07,06) + rrt(149) * density(04) 
  pd(03,03) = pd(03,03) - rrt(150) * density(07) 
  pd(03,07) = pd(03,07) - rrt(150) * density(03) 
  pd(04,03) = pd(04,03) + rrt(150) * density(07) 
  pd(04,07) = pd(04,07) + rrt(150) * density(03) 
  pd(06,03) = pd(06,03) + rrt(150) * density(07) 
  pd(06,07) = pd(06,07) + rrt(150) * density(03) 
  pd(07,03) = pd(07,03) - rrt(150) * density(07) 
  pd(07,07) = pd(07,07) - rrt(150) * density(03) 
  pd(03,04) = pd(03,04) + rrt(151) * density(07) 
  pd(03,07) = pd(03,07) + rrt(151) * density(04) 
  pd(04,04) = pd(04,04) - rrt(151) * density(07) 
  pd(04,07) = pd(04,07) - rrt(151) * density(04) 
  pd(07,04) = pd(07,04) - rrt(151) * density(07) 
  pd(07,07) = pd(07,07) - rrt(151) * density(04) 
  pd(08,04) = pd(08,04) + rrt(151) * density(07) 
  pd(08,07) = pd(08,07) + rrt(151) * density(04) 
  pd(03,03) = pd(03,03) - rrt(152) * density(08) 
  pd(03,08) = pd(03,08) - rrt(152) * density(03) 
  pd(04,03) = pd(04,03) + rrt(152) * density(08) 
  pd(04,08) = pd(04,08) + rrt(152) * density(03) 
  pd(07,03) = pd(07,03) + rrt(152) * density(08) 
  pd(07,08) = pd(07,08) + rrt(152) * density(03) 
  pd(08,03) = pd(08,03) - rrt(152) * density(08) 
  pd(08,08) = pd(08,08) - rrt(152) * density(03) 
  pd(03,04) = pd(03,04) + rrt(153) * density(08) 
  pd(03,08) = pd(03,08) + rrt(153) * density(04) 
  pd(04,04) = pd(04,04) - rrt(153) * density(08) 
  pd(04,08) = pd(04,08) - rrt(153) * density(04) 
  pd(08,04) = pd(08,04) - rrt(153) * density(08) 
  pd(08,08) = pd(08,08) - rrt(153) * density(04) 
  pd(09,04) = pd(09,04) + rrt(153) * density(08) 
  pd(09,08) = pd(09,08) + rrt(153) * density(04) 
  pd(03,03) = pd(03,03) - rrt(154) * density(09) 
  pd(03,09) = pd(03,09) - rrt(154) * density(03) 
  pd(04,03) = pd(04,03) + rrt(154) * density(09) 
  pd(04,09) = pd(04,09) + rrt(154) * density(03) 
  pd(08,03) = pd(08,03) + rrt(154) * density(09) 
  pd(08,09) = pd(08,09) + rrt(154) * density(03) 
  pd(09,03) = pd(09,03) - rrt(154) * density(09) 
  pd(09,09) = pd(09,09) - rrt(154) * density(03) 
  pd(04,05) = pd(04,05) + rrt(155) * density(05) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(155) * density(05) * 4.0d0
  pd(06,05) = pd(06,05) + rrt(155) * density(05) * 2.0d0
  pd(04,04) = pd(04,04) - rrt(156) * density(06) 
  pd(04,06) = pd(04,06) - rrt(156) * density(04) 
  pd(05,04) = pd(05,04) + rrt(156) * density(06) * 2.0d0
  pd(05,06) = pd(05,06) + rrt(156) * density(04) * 2.0d0
  pd(06,04) = pd(06,04) - rrt(156) * density(06) 
  pd(06,06) = pd(06,06) - rrt(156) * density(04) 
  pd(04,05) = pd(04,05) + rrt(157) * density(06) 
  pd(04,06) = pd(04,06) + rrt(157) * density(05) 
  pd(05,05) = pd(05,05) - rrt(157) * density(06) 
  pd(05,06) = pd(05,06) - rrt(157) * density(05) 
  pd(06,05) = pd(06,05) - rrt(157) * density(06) 
  pd(06,06) = pd(06,06) - rrt(157) * density(05) 
  pd(07,05) = pd(07,05) + rrt(157) * density(06) 
  pd(07,06) = pd(07,06) + rrt(157) * density(05) 
  pd(04,04) = pd(04,04) - rrt(158) * density(07) 
  pd(04,07) = pd(04,07) - rrt(158) * density(04) 
  pd(05,04) = pd(05,04) + rrt(158) * density(07) 
  pd(05,07) = pd(05,07) + rrt(158) * density(04) 
  pd(06,04) = pd(06,04) + rrt(158) * density(07) 
  pd(06,07) = pd(06,07) + rrt(158) * density(04) 
  pd(07,04) = pd(07,04) - rrt(158) * density(07) 
  pd(07,07) = pd(07,07) - rrt(158) * density(04) 
  pd(04,05) = pd(04,05) + rrt(159) * density(07) 
  pd(04,07) = pd(04,07) + rrt(159) * density(05) 
  pd(05,05) = pd(05,05) - rrt(159) * density(07) 
  pd(05,07) = pd(05,07) - rrt(159) * density(05) 
  pd(07,05) = pd(07,05) - rrt(159) * density(07) 
  pd(07,07) = pd(07,07) - rrt(159) * density(05) 
  pd(08,05) = pd(08,05) + rrt(159) * density(07) 
  pd(08,07) = pd(08,07) + rrt(159) * density(05) 
  pd(04,04) = pd(04,04) - rrt(160) * density(08) 
  pd(04,08) = pd(04,08) - rrt(160) * density(04) 
  pd(05,04) = pd(05,04) + rrt(160) * density(08) 
  pd(05,08) = pd(05,08) + rrt(160) * density(04) 
  pd(07,04) = pd(07,04) + rrt(160) * density(08) 
  pd(07,08) = pd(07,08) + rrt(160) * density(04) 
  pd(08,04) = pd(08,04) - rrt(160) * density(08) 
  pd(08,08) = pd(08,08) - rrt(160) * density(04) 
  pd(04,05) = pd(04,05) + rrt(161) * density(08) 
  pd(04,08) = pd(04,08) + rrt(161) * density(05) 
  pd(05,05) = pd(05,05) - rrt(161) * density(08) 
  pd(05,08) = pd(05,08) - rrt(161) * density(05) 
  pd(08,05) = pd(08,05) - rrt(161) * density(08) 
  pd(08,08) = pd(08,08) - rrt(161) * density(05) 
  pd(09,05) = pd(09,05) + rrt(161) * density(08) 
  pd(09,08) = pd(09,08) + rrt(161) * density(05) 
  pd(04,04) = pd(04,04) - rrt(162) * density(09) 
  pd(04,09) = pd(04,09) - rrt(162) * density(04) 
  pd(05,04) = pd(05,04) + rrt(162) * density(09) 
  pd(05,09) = pd(05,09) + rrt(162) * density(04) 
  pd(08,04) = pd(08,04) + rrt(162) * density(09) 
  pd(08,09) = pd(08,09) + rrt(162) * density(04) 
  pd(09,04) = pd(09,04) - rrt(162) * density(09) 
  pd(09,09) = pd(09,09) - rrt(162) * density(04) 
  pd(05,06) = pd(05,06) + rrt(163) * density(06) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(163) * density(06) * 4.0d0
  pd(07,06) = pd(07,06) + rrt(163) * density(06) * 2.0d0
  pd(05,05) = pd(05,05) - rrt(164) * density(07) 
  pd(05,07) = pd(05,07) - rrt(164) * density(05) 
  pd(06,05) = pd(06,05) + rrt(164) * density(07) * 2.0d0
  pd(06,07) = pd(06,07) + rrt(164) * density(05) * 2.0d0
  pd(07,05) = pd(07,05) - rrt(164) * density(07) 
  pd(07,07) = pd(07,07) - rrt(164) * density(05) 
  pd(05,06) = pd(05,06) + rrt(165) * density(07) 
  pd(05,07) = pd(05,07) + rrt(165) * density(06) 
  pd(06,06) = pd(06,06) - rrt(165) * density(07) 
  pd(06,07) = pd(06,07) - rrt(165) * density(06) 
  pd(07,06) = pd(07,06) - rrt(165) * density(07) 
  pd(07,07) = pd(07,07) - rrt(165) * density(06) 
  pd(08,06) = pd(08,06) + rrt(165) * density(07) 
  pd(08,07) = pd(08,07) + rrt(165) * density(06) 
  pd(05,05) = pd(05,05) - rrt(166) * density(08) 
  pd(05,08) = pd(05,08) - rrt(166) * density(05) 
  pd(06,05) = pd(06,05) + rrt(166) * density(08) 
  pd(06,08) = pd(06,08) + rrt(166) * density(05) 
  pd(07,05) = pd(07,05) + rrt(166) * density(08) 
  pd(07,08) = pd(07,08) + rrt(166) * density(05) 
  pd(08,05) = pd(08,05) - rrt(166) * density(08) 
  pd(08,08) = pd(08,08) - rrt(166) * density(05) 
  pd(05,06) = pd(05,06) + rrt(167) * density(08) 
  pd(05,08) = pd(05,08) + rrt(167) * density(06) 
  pd(06,06) = pd(06,06) - rrt(167) * density(08) 
  pd(06,08) = pd(06,08) - rrt(167) * density(06) 
  pd(08,06) = pd(08,06) - rrt(167) * density(08) 
  pd(08,08) = pd(08,08) - rrt(167) * density(06) 
  pd(09,06) = pd(09,06) + rrt(167) * density(08) 
  pd(09,08) = pd(09,08) + rrt(167) * density(06) 
  pd(05,05) = pd(05,05) - rrt(168) * density(09) 
  pd(05,09) = pd(05,09) - rrt(168) * density(05) 
  pd(06,05) = pd(06,05) + rrt(168) * density(09) 
  pd(06,09) = pd(06,09) + rrt(168) * density(05) 
  pd(08,05) = pd(08,05) + rrt(168) * density(09) 
  pd(08,09) = pd(08,09) + rrt(168) * density(05) 
  pd(09,05) = pd(09,05) - rrt(168) * density(09) 
  pd(09,09) = pd(09,09) - rrt(168) * density(05) 
  pd(06,07) = pd(06,07) + rrt(169) * density(07) * 2.0d0
  pd(07,07) = pd(07,07) - rrt(169) * density(07) * 4.0d0
  pd(08,07) = pd(08,07) + rrt(169) * density(07) * 2.0d0
  pd(06,06) = pd(06,06) - rrt(170) * density(08) 
  pd(06,08) = pd(06,08) - rrt(170) * density(06) 
  pd(07,06) = pd(07,06) + rrt(170) * density(08) * 2.0d0
  pd(07,08) = pd(07,08) + rrt(170) * density(06) * 2.0d0
  pd(08,06) = pd(08,06) - rrt(170) * density(08) 
  pd(08,08) = pd(08,08) - rrt(170) * density(06) 
  pd(06,07) = pd(06,07) + rrt(171) * density(08) 
  pd(06,08) = pd(06,08) + rrt(171) * density(07) 
  pd(07,07) = pd(07,07) - rrt(171) * density(08) 
  pd(07,08) = pd(07,08) - rrt(171) * density(07) 
  pd(08,07) = pd(08,07) - rrt(171) * density(08) 
  pd(08,08) = pd(08,08) - rrt(171) * density(07) 
  pd(09,07) = pd(09,07) + rrt(171) * density(08) 
  pd(09,08) = pd(09,08) + rrt(171) * density(07) 
  pd(06,06) = pd(06,06) - rrt(172) * density(09) 
  pd(06,09) = pd(06,09) - rrt(172) * density(06) 
  pd(07,06) = pd(07,06) + rrt(172) * density(09) 
  pd(07,09) = pd(07,09) + rrt(172) * density(06) 
  pd(08,06) = pd(08,06) + rrt(172) * density(09) 
  pd(08,09) = pd(08,09) + rrt(172) * density(06) 
  pd(09,06) = pd(09,06) - rrt(172) * density(09) 
  pd(09,09) = pd(09,09) - rrt(172) * density(06) 
  pd(07,08) = pd(07,08) + rrt(173) * density(08) * 2.0d0
  pd(08,08) = pd(08,08) - rrt(173) * density(08) * 4.0d0
  pd(09,08) = pd(09,08) + rrt(173) * density(08) * 2.0d0
  pd(07,07) = pd(07,07) - rrt(174) * density(09) 
  pd(07,09) = pd(07,09) - rrt(174) * density(07) 
  pd(08,07) = pd(08,07) + rrt(174) * density(09) * 2.0d0
  pd(08,09) = pd(08,09) + rrt(174) * density(07) * 2.0d0
  pd(09,07) = pd(09,07) - rrt(174) * density(09) 
  pd(09,09) = pd(09,09) - rrt(174) * density(07) 
  pd(01,01) = pd(01,01) - rrt(175) * density(27) 
  pd(01,27) = pd(01,27) - rrt(175) * density(01) 
  pd(02,01) = pd(02,01) + rrt(175) * density(27) 
  pd(02,27) = pd(02,27) + rrt(175) * density(01) 
  pd(21,01) = pd(21,01) + rrt(175) * density(27) 
  pd(21,27) = pd(21,27) + rrt(175) * density(01) 
  pd(27,01) = pd(27,01) - rrt(175) * density(27) 
  pd(27,27) = pd(27,27) - rrt(175) * density(01) 
  pd(01,02) = pd(01,02) + rrt(176) * density(21) 
  pd(01,21) = pd(01,21) + rrt(176) * density(02) 
  pd(02,02) = pd(02,02) - rrt(176) * density(21) 
  pd(02,21) = pd(02,21) - rrt(176) * density(02) 
  pd(21,02) = pd(21,02) - rrt(176) * density(21) 
  pd(21,21) = pd(21,21) - rrt(176) * density(02) 
  pd(27,02) = pd(27,02) + rrt(176) * density(21) 
  pd(27,21) = pd(27,21) + rrt(176) * density(02) 
  pd(02,02) = pd(02,02) - rrt(177) * density(27) 
  pd(02,27) = pd(02,27) - rrt(177) * density(02) 
  pd(03,02) = pd(03,02) + rrt(177) * density(27) 
  pd(03,27) = pd(03,27) + rrt(177) * density(02) 
  pd(21,02) = pd(21,02) + rrt(177) * density(27) 
  pd(21,27) = pd(21,27) + rrt(177) * density(02) 
  pd(27,02) = pd(27,02) - rrt(177) * density(27) 
  pd(27,27) = pd(27,27) - rrt(177) * density(02) 
  pd(02,03) = pd(02,03) + rrt(178) * density(21) 
  pd(02,21) = pd(02,21) + rrt(178) * density(03) 
  pd(03,03) = pd(03,03) - rrt(178) * density(21) 
  pd(03,21) = pd(03,21) - rrt(178) * density(03) 
  pd(21,03) = pd(21,03) - rrt(178) * density(21) 
  pd(21,21) = pd(21,21) - rrt(178) * density(03) 
  pd(27,03) = pd(27,03) + rrt(178) * density(21) 
  pd(27,21) = pd(27,21) + rrt(178) * density(03) 
  pd(03,03) = pd(03,03) - rrt(179) * density(27) 
  pd(03,27) = pd(03,27) - rrt(179) * density(03) 
  pd(04,03) = pd(04,03) + rrt(179) * density(27) 
  pd(04,27) = pd(04,27) + rrt(179) * density(03) 
  pd(21,03) = pd(21,03) + rrt(179) * density(27) 
  pd(21,27) = pd(21,27) + rrt(179) * density(03) 
  pd(27,03) = pd(27,03) - rrt(179) * density(27) 
  pd(27,27) = pd(27,27) - rrt(179) * density(03) 
  pd(03,04) = pd(03,04) + rrt(180) * density(21) 
  pd(03,21) = pd(03,21) + rrt(180) * density(04) 
  pd(04,04) = pd(04,04) - rrt(180) * density(21) 
  pd(04,21) = pd(04,21) - rrt(180) * density(04) 
  pd(21,04) = pd(21,04) - rrt(180) * density(21) 
  pd(21,21) = pd(21,21) - rrt(180) * density(04) 
  pd(27,04) = pd(27,04) + rrt(180) * density(21) 
  pd(27,21) = pd(27,21) + rrt(180) * density(04) 
  pd(04,04) = pd(04,04) - rrt(181) * density(27) 
  pd(04,27) = pd(04,27) - rrt(181) * density(04) 
  pd(05,04) = pd(05,04) + rrt(181) * density(27) 
  pd(05,27) = pd(05,27) + rrt(181) * density(04) 
  pd(21,04) = pd(21,04) + rrt(181) * density(27) 
  pd(21,27) = pd(21,27) + rrt(181) * density(04) 
  pd(27,04) = pd(27,04) - rrt(181) * density(27) 
  pd(27,27) = pd(27,27) - rrt(181) * density(04) 
  pd(04,05) = pd(04,05) + rrt(182) * density(21) 
  pd(04,21) = pd(04,21) + rrt(182) * density(05) 
  pd(05,05) = pd(05,05) - rrt(182) * density(21) 
  pd(05,21) = pd(05,21) - rrt(182) * density(05) 
  pd(21,05) = pd(21,05) - rrt(182) * density(21) 
  pd(21,21) = pd(21,21) - rrt(182) * density(05) 
  pd(27,05) = pd(27,05) + rrt(182) * density(21) 
  pd(27,21) = pd(27,21) + rrt(182) * density(05) 
  pd(05,05) = pd(05,05) - rrt(183) * density(27) 
  pd(05,27) = pd(05,27) - rrt(183) * density(05) 
  pd(06,05) = pd(06,05) + rrt(183) * density(27) 
  pd(06,27) = pd(06,27) + rrt(183) * density(05) 
  pd(21,05) = pd(21,05) + rrt(183) * density(27) 
  pd(21,27) = pd(21,27) + rrt(183) * density(05) 
  pd(27,05) = pd(27,05) - rrt(183) * density(27) 
  pd(27,27) = pd(27,27) - rrt(183) * density(05) 
  pd(05,06) = pd(05,06) + rrt(184) * density(21) 
  pd(05,21) = pd(05,21) + rrt(184) * density(06) 
  pd(06,06) = pd(06,06) - rrt(184) * density(21) 
  pd(06,21) = pd(06,21) - rrt(184) * density(06) 
  pd(21,06) = pd(21,06) - rrt(184) * density(21) 
  pd(21,21) = pd(21,21) - rrt(184) * density(06) 
  pd(27,06) = pd(27,06) + rrt(184) * density(21) 
  pd(27,21) = pd(27,21) + rrt(184) * density(06) 
  pd(06,06) = pd(06,06) - rrt(185) * density(27) 
  pd(06,27) = pd(06,27) - rrt(185) * density(06) 
  pd(07,06) = pd(07,06) + rrt(185) * density(27) 
  pd(07,27) = pd(07,27) + rrt(185) * density(06) 
  pd(21,06) = pd(21,06) + rrt(185) * density(27) 
  pd(21,27) = pd(21,27) + rrt(185) * density(06) 
  pd(27,06) = pd(27,06) - rrt(185) * density(27) 
  pd(27,27) = pd(27,27) - rrt(185) * density(06) 
  pd(06,07) = pd(06,07) + rrt(186) * density(21) 
  pd(06,21) = pd(06,21) + rrt(186) * density(07) 
  pd(07,07) = pd(07,07) - rrt(186) * density(21) 
  pd(07,21) = pd(07,21) - rrt(186) * density(07) 
  pd(21,07) = pd(21,07) - rrt(186) * density(21) 
  pd(21,21) = pd(21,21) - rrt(186) * density(07) 
  pd(27,07) = pd(27,07) + rrt(186) * density(21) 
  pd(27,21) = pd(27,21) + rrt(186) * density(07) 
  pd(07,07) = pd(07,07) - rrt(187) * density(27) 
  pd(07,27) = pd(07,27) - rrt(187) * density(07) 
  pd(08,07) = pd(08,07) + rrt(187) * density(27) 
  pd(08,27) = pd(08,27) + rrt(187) * density(07) 
  pd(21,07) = pd(21,07) + rrt(187) * density(27) 
  pd(21,27) = pd(21,27) + rrt(187) * density(07) 
  pd(27,07) = pd(27,07) - rrt(187) * density(27) 
  pd(27,27) = pd(27,27) - rrt(187) * density(07) 
  pd(07,08) = pd(07,08) + rrt(188) * density(21) 
  pd(07,21) = pd(07,21) + rrt(188) * density(08) 
  pd(08,08) = pd(08,08) - rrt(188) * density(21) 
  pd(08,21) = pd(08,21) - rrt(188) * density(08) 
  pd(21,08) = pd(21,08) - rrt(188) * density(21) 
  pd(21,21) = pd(21,21) - rrt(188) * density(08) 
  pd(27,08) = pd(27,08) + rrt(188) * density(21) 
  pd(27,21) = pd(27,21) + rrt(188) * density(08) 
  pd(08,08) = pd(08,08) - rrt(189) * density(27) 
  pd(08,27) = pd(08,27) - rrt(189) * density(08) 
  pd(09,08) = pd(09,08) + rrt(189) * density(27) 
  pd(09,27) = pd(09,27) + rrt(189) * density(08) 
  pd(21,08) = pd(21,08) + rrt(189) * density(27) 
  pd(21,27) = pd(21,27) + rrt(189) * density(08) 
  pd(27,08) = pd(27,08) - rrt(189) * density(27) 
  pd(27,27) = pd(27,27) - rrt(189) * density(08) 
  pd(08,09) = pd(08,09) + rrt(190) * density(21) 
  pd(08,21) = pd(08,21) + rrt(190) * density(09) 
  pd(09,09) = pd(09,09) - rrt(190) * density(21) 
  pd(09,21) = pd(09,21) - rrt(190) * density(09) 
  pd(21,09) = pd(21,09) - rrt(190) * density(21) 
  pd(21,21) = pd(21,21) - rrt(190) * density(09) 
  pd(27,09) = pd(27,09) + rrt(190) * density(21) 
  pd(27,21) = pd(27,21) + rrt(190) * density(09) 
  pd(01,01) = pd(01,01) - rrt(191) * density(28) 
  pd(01,28) = pd(01,28) - rrt(191) * density(01) 
  pd(02,01) = pd(02,01) + rrt(191) * density(28) 
  pd(02,28) = pd(02,28) + rrt(191) * density(01) 
  pd(27,01) = pd(27,01) + rrt(191) * density(28) 
  pd(27,28) = pd(27,28) + rrt(191) * density(01) 
  pd(28,01) = pd(28,01) - rrt(191) * density(28) 
  pd(28,28) = pd(28,28) - rrt(191) * density(01) 
  pd(01,02) = pd(01,02) + rrt(192) * density(27) 
  pd(01,27) = pd(01,27) + rrt(192) * density(02) 
  pd(02,02) = pd(02,02) - rrt(192) * density(27) 
  pd(02,27) = pd(02,27) - rrt(192) * density(02) 
  pd(27,02) = pd(27,02) - rrt(192) * density(27) 
  pd(27,27) = pd(27,27) - rrt(192) * density(02) 
  pd(28,02) = pd(28,02) + rrt(192) * density(27) 
  pd(28,27) = pd(28,27) + rrt(192) * density(02) 
  pd(02,02) = pd(02,02) - rrt(193) * density(28) 
  pd(02,28) = pd(02,28) - rrt(193) * density(02) 
  pd(03,02) = pd(03,02) + rrt(193) * density(28) 
  pd(03,28) = pd(03,28) + rrt(193) * density(02) 
  pd(27,02) = pd(27,02) + rrt(193) * density(28) 
  pd(27,28) = pd(27,28) + rrt(193) * density(02) 
  pd(28,02) = pd(28,02) - rrt(193) * density(28) 
  pd(28,28) = pd(28,28) - rrt(193) * density(02) 
  pd(02,03) = pd(02,03) + rrt(194) * density(27) 
  pd(02,27) = pd(02,27) + rrt(194) * density(03) 
  pd(03,03) = pd(03,03) - rrt(194) * density(27) 
  pd(03,27) = pd(03,27) - rrt(194) * density(03) 
  pd(27,03) = pd(27,03) - rrt(194) * density(27) 
  pd(27,27) = pd(27,27) - rrt(194) * density(03) 
  pd(28,03) = pd(28,03) + rrt(194) * density(27) 
  pd(28,27) = pd(28,27) + rrt(194) * density(03) 
  pd(03,03) = pd(03,03) - rrt(195) * density(28) 
  pd(03,28) = pd(03,28) - rrt(195) * density(03) 
  pd(04,03) = pd(04,03) + rrt(195) * density(28) 
  pd(04,28) = pd(04,28) + rrt(195) * density(03) 
  pd(27,03) = pd(27,03) + rrt(195) * density(28) 
  pd(27,28) = pd(27,28) + rrt(195) * density(03) 
  pd(28,03) = pd(28,03) - rrt(195) * density(28) 
  pd(28,28) = pd(28,28) - rrt(195) * density(03) 
  pd(03,04) = pd(03,04) + rrt(196) * density(27) 
  pd(03,27) = pd(03,27) + rrt(196) * density(04) 
  pd(04,04) = pd(04,04) - rrt(196) * density(27) 
  pd(04,27) = pd(04,27) - rrt(196) * density(04) 
  pd(27,04) = pd(27,04) - rrt(196) * density(27) 
  pd(27,27) = pd(27,27) - rrt(196) * density(04) 
  pd(28,04) = pd(28,04) + rrt(196) * density(27) 
  pd(28,27) = pd(28,27) + rrt(196) * density(04) 
  pd(04,04) = pd(04,04) - rrt(197) * density(28) 
  pd(04,28) = pd(04,28) - rrt(197) * density(04) 
  pd(05,04) = pd(05,04) + rrt(197) * density(28) 
  pd(05,28) = pd(05,28) + rrt(197) * density(04) 
  pd(27,04) = pd(27,04) + rrt(197) * density(28) 
  pd(27,28) = pd(27,28) + rrt(197) * density(04) 
  pd(28,04) = pd(28,04) - rrt(197) * density(28) 
  pd(28,28) = pd(28,28) - rrt(197) * density(04) 
  pd(04,05) = pd(04,05) + rrt(198) * density(27) 
  pd(04,27) = pd(04,27) + rrt(198) * density(05) 
  pd(05,05) = pd(05,05) - rrt(198) * density(27) 
  pd(05,27) = pd(05,27) - rrt(198) * density(05) 
  pd(27,05) = pd(27,05) - rrt(198) * density(27) 
  pd(27,27) = pd(27,27) - rrt(198) * density(05) 
  pd(28,05) = pd(28,05) + rrt(198) * density(27) 
  pd(28,27) = pd(28,27) + rrt(198) * density(05) 
  pd(05,05) = pd(05,05) - rrt(199) * density(28) 
  pd(05,28) = pd(05,28) - rrt(199) * density(05) 
  pd(06,05) = pd(06,05) + rrt(199) * density(28) 
  pd(06,28) = pd(06,28) + rrt(199) * density(05) 
  pd(27,05) = pd(27,05) + rrt(199) * density(28) 
  pd(27,28) = pd(27,28) + rrt(199) * density(05) 
  pd(28,05) = pd(28,05) - rrt(199) * density(28) 
  pd(28,28) = pd(28,28) - rrt(199) * density(05) 
  pd(05,06) = pd(05,06) + rrt(200) * density(27) 
  pd(05,27) = pd(05,27) + rrt(200) * density(06) 
  pd(06,06) = pd(06,06) - rrt(200) * density(27) 
  pd(06,27) = pd(06,27) - rrt(200) * density(06) 
  pd(27,06) = pd(27,06) - rrt(200) * density(27) 
  pd(27,27) = pd(27,27) - rrt(200) * density(06) 
  pd(28,06) = pd(28,06) + rrt(200) * density(27) 
  pd(28,27) = pd(28,27) + rrt(200) * density(06) 
  pd(06,06) = pd(06,06) - rrt(201) * density(28) 
  pd(06,28) = pd(06,28) - rrt(201) * density(06) 
  pd(07,06) = pd(07,06) + rrt(201) * density(28) 
  pd(07,28) = pd(07,28) + rrt(201) * density(06) 
  pd(27,06) = pd(27,06) + rrt(201) * density(28) 
  pd(27,28) = pd(27,28) + rrt(201) * density(06) 
  pd(28,06) = pd(28,06) - rrt(201) * density(28) 
  pd(28,28) = pd(28,28) - rrt(201) * density(06) 
  pd(06,07) = pd(06,07) + rrt(202) * density(27) 
  pd(06,27) = pd(06,27) + rrt(202) * density(07) 
  pd(07,07) = pd(07,07) - rrt(202) * density(27) 
  pd(07,27) = pd(07,27) - rrt(202) * density(07) 
  pd(27,07) = pd(27,07) - rrt(202) * density(27) 
  pd(27,27) = pd(27,27) - rrt(202) * density(07) 
  pd(28,07) = pd(28,07) + rrt(202) * density(27) 
  pd(28,27) = pd(28,27) + rrt(202) * density(07) 
  pd(07,07) = pd(07,07) - rrt(203) * density(28) 
  pd(07,28) = pd(07,28) - rrt(203) * density(07) 
  pd(08,07) = pd(08,07) + rrt(203) * density(28) 
  pd(08,28) = pd(08,28) + rrt(203) * density(07) 
  pd(27,07) = pd(27,07) + rrt(203) * density(28) 
  pd(27,28) = pd(27,28) + rrt(203) * density(07) 
  pd(28,07) = pd(28,07) - rrt(203) * density(28) 
  pd(28,28) = pd(28,28) - rrt(203) * density(07) 
  pd(07,08) = pd(07,08) + rrt(204) * density(27) 
  pd(07,27) = pd(07,27) + rrt(204) * density(08) 
  pd(08,08) = pd(08,08) - rrt(204) * density(27) 
  pd(08,27) = pd(08,27) - rrt(204) * density(08) 
  pd(27,08) = pd(27,08) - rrt(204) * density(27) 
  pd(27,27) = pd(27,27) - rrt(204) * density(08) 
  pd(28,08) = pd(28,08) + rrt(204) * density(27) 
  pd(28,27) = pd(28,27) + rrt(204) * density(08) 
  pd(08,08) = pd(08,08) - rrt(205) * density(28) 
  pd(08,28) = pd(08,28) - rrt(205) * density(08) 
  pd(09,08) = pd(09,08) + rrt(205) * density(28) 
  pd(09,28) = pd(09,28) + rrt(205) * density(08) 
  pd(27,08) = pd(27,08) + rrt(205) * density(28) 
  pd(27,28) = pd(27,28) + rrt(205) * density(08) 
  pd(28,08) = pd(28,08) - rrt(205) * density(28) 
  pd(28,28) = pd(28,28) - rrt(205) * density(08) 
  pd(08,09) = pd(08,09) + rrt(206) * density(27) 
  pd(08,27) = pd(08,27) + rrt(206) * density(09) 
  pd(09,09) = pd(09,09) - rrt(206) * density(27) 
  pd(09,27) = pd(09,27) - rrt(206) * density(09) 
  pd(27,09) = pd(27,09) - rrt(206) * density(27) 
  pd(27,27) = pd(27,27) - rrt(206) * density(09) 
  pd(28,09) = pd(28,09) + rrt(206) * density(27) 
  pd(28,27) = pd(28,27) + rrt(206) * density(09) 
  pd(01,01) = pd(01,01) - rrt(207) * density(29) 
  pd(01,29) = pd(01,29) - rrt(207) * density(01) 
  pd(02,01) = pd(02,01) + rrt(207) * density(29) 
  pd(02,29) = pd(02,29) + rrt(207) * density(01) 
  pd(28,01) = pd(28,01) + rrt(207) * density(29) 
  pd(28,29) = pd(28,29) + rrt(207) * density(01) 
  pd(29,01) = pd(29,01) - rrt(207) * density(29) 
  pd(29,29) = pd(29,29) - rrt(207) * density(01) 
  pd(01,02) = pd(01,02) + rrt(208) * density(28) 
  pd(01,28) = pd(01,28) + rrt(208) * density(02) 
  pd(02,02) = pd(02,02) - rrt(208) * density(28) 
  pd(02,28) = pd(02,28) - rrt(208) * density(02) 
  pd(28,02) = pd(28,02) - rrt(208) * density(28) 
  pd(28,28) = pd(28,28) - rrt(208) * density(02) 
  pd(29,02) = pd(29,02) + rrt(208) * density(28) 
  pd(29,28) = pd(29,28) + rrt(208) * density(02) 
  pd(02,02) = pd(02,02) - rrt(209) * density(29) 
  pd(02,29) = pd(02,29) - rrt(209) * density(02) 
  pd(03,02) = pd(03,02) + rrt(209) * density(29) 
  pd(03,29) = pd(03,29) + rrt(209) * density(02) 
  pd(28,02) = pd(28,02) + rrt(209) * density(29) 
  pd(28,29) = pd(28,29) + rrt(209) * density(02) 
  pd(29,02) = pd(29,02) - rrt(209) * density(29) 
  pd(29,29) = pd(29,29) - rrt(209) * density(02) 
  pd(02,03) = pd(02,03) + rrt(210) * density(28) 
  pd(02,28) = pd(02,28) + rrt(210) * density(03) 
  pd(03,03) = pd(03,03) - rrt(210) * density(28) 
  pd(03,28) = pd(03,28) - rrt(210) * density(03) 
  pd(28,03) = pd(28,03) - rrt(210) * density(28) 
  pd(28,28) = pd(28,28) - rrt(210) * density(03) 
  pd(29,03) = pd(29,03) + rrt(210) * density(28) 
  pd(29,28) = pd(29,28) + rrt(210) * density(03) 
  pd(03,03) = pd(03,03) - rrt(211) * density(29) 
  pd(03,29) = pd(03,29) - rrt(211) * density(03) 
  pd(04,03) = pd(04,03) + rrt(211) * density(29) 
  pd(04,29) = pd(04,29) + rrt(211) * density(03) 
  pd(28,03) = pd(28,03) + rrt(211) * density(29) 
  pd(28,29) = pd(28,29) + rrt(211) * density(03) 
  pd(29,03) = pd(29,03) - rrt(211) * density(29) 
  pd(29,29) = pd(29,29) - rrt(211) * density(03) 
  pd(03,04) = pd(03,04) + rrt(212) * density(28) 
  pd(03,28) = pd(03,28) + rrt(212) * density(04) 
  pd(04,04) = pd(04,04) - rrt(212) * density(28) 
  pd(04,28) = pd(04,28) - rrt(212) * density(04) 
  pd(28,04) = pd(28,04) - rrt(212) * density(28) 
  pd(28,28) = pd(28,28) - rrt(212) * density(04) 
  pd(29,04) = pd(29,04) + rrt(212) * density(28) 
  pd(29,28) = pd(29,28) + rrt(212) * density(04) 
  pd(04,04) = pd(04,04) - rrt(213) * density(29) 
  pd(04,29) = pd(04,29) - rrt(213) * density(04) 
  pd(05,04) = pd(05,04) + rrt(213) * density(29) 
  pd(05,29) = pd(05,29) + rrt(213) * density(04) 
  pd(28,04) = pd(28,04) + rrt(213) * density(29) 
  pd(28,29) = pd(28,29) + rrt(213) * density(04) 
  pd(29,04) = pd(29,04) - rrt(213) * density(29) 
  pd(29,29) = pd(29,29) - rrt(213) * density(04) 
  pd(04,05) = pd(04,05) + rrt(214) * density(28) 
  pd(04,28) = pd(04,28) + rrt(214) * density(05) 
  pd(05,05) = pd(05,05) - rrt(214) * density(28) 
  pd(05,28) = pd(05,28) - rrt(214) * density(05) 
  pd(28,05) = pd(28,05) - rrt(214) * density(28) 
  pd(28,28) = pd(28,28) - rrt(214) * density(05) 
  pd(29,05) = pd(29,05) + rrt(214) * density(28) 
  pd(29,28) = pd(29,28) + rrt(214) * density(05) 
  pd(05,05) = pd(05,05) - rrt(215) * density(29) 
  pd(05,29) = pd(05,29) - rrt(215) * density(05) 
  pd(06,05) = pd(06,05) + rrt(215) * density(29) 
  pd(06,29) = pd(06,29) + rrt(215) * density(05) 
  pd(28,05) = pd(28,05) + rrt(215) * density(29) 
  pd(28,29) = pd(28,29) + rrt(215) * density(05) 
  pd(29,05) = pd(29,05) - rrt(215) * density(29) 
  pd(29,29) = pd(29,29) - rrt(215) * density(05) 
  pd(05,06) = pd(05,06) + rrt(216) * density(28) 
  pd(05,28) = pd(05,28) + rrt(216) * density(06) 
  pd(06,06) = pd(06,06) - rrt(216) * density(28) 
  pd(06,28) = pd(06,28) - rrt(216) * density(06) 
  pd(28,06) = pd(28,06) - rrt(216) * density(28) 
  pd(28,28) = pd(28,28) - rrt(216) * density(06) 
  pd(29,06) = pd(29,06) + rrt(216) * density(28) 
  pd(29,28) = pd(29,28) + rrt(216) * density(06) 
  pd(06,06) = pd(06,06) - rrt(217) * density(29) 
  pd(06,29) = pd(06,29) - rrt(217) * density(06) 
  pd(07,06) = pd(07,06) + rrt(217) * density(29) 
  pd(07,29) = pd(07,29) + rrt(217) * density(06) 
  pd(28,06) = pd(28,06) + rrt(217) * density(29) 
  pd(28,29) = pd(28,29) + rrt(217) * density(06) 
  pd(29,06) = pd(29,06) - rrt(217) * density(29) 
  pd(29,29) = pd(29,29) - rrt(217) * density(06) 
  pd(06,07) = pd(06,07) + rrt(218) * density(28) 
  pd(06,28) = pd(06,28) + rrt(218) * density(07) 
  pd(07,07) = pd(07,07) - rrt(218) * density(28) 
  pd(07,28) = pd(07,28) - rrt(218) * density(07) 
  pd(28,07) = pd(28,07) - rrt(218) * density(28) 
  pd(28,28) = pd(28,28) - rrt(218) * density(07) 
  pd(29,07) = pd(29,07) + rrt(218) * density(28) 
  pd(29,28) = pd(29,28) + rrt(218) * density(07) 
  pd(07,07) = pd(07,07) - rrt(219) * density(29) 
  pd(07,29) = pd(07,29) - rrt(219) * density(07) 
  pd(08,07) = pd(08,07) + rrt(219) * density(29) 
  pd(08,29) = pd(08,29) + rrt(219) * density(07) 
  pd(28,07) = pd(28,07) + rrt(219) * density(29) 
  pd(28,29) = pd(28,29) + rrt(219) * density(07) 
  pd(29,07) = pd(29,07) - rrt(219) * density(29) 
  pd(29,29) = pd(29,29) - rrt(219) * density(07) 
  pd(07,08) = pd(07,08) + rrt(220) * density(28) 
  pd(07,28) = pd(07,28) + rrt(220) * density(08) 
  pd(08,08) = pd(08,08) - rrt(220) * density(28) 
  pd(08,28) = pd(08,28) - rrt(220) * density(08) 
  pd(28,08) = pd(28,08) - rrt(220) * density(28) 
  pd(28,28) = pd(28,28) - rrt(220) * density(08) 
  pd(29,08) = pd(29,08) + rrt(220) * density(28) 
  pd(29,28) = pd(29,28) + rrt(220) * density(08) 
  pd(08,08) = pd(08,08) - rrt(221) * density(29) 
  pd(08,29) = pd(08,29) - rrt(221) * density(08) 
  pd(09,08) = pd(09,08) + rrt(221) * density(29) 
  pd(09,29) = pd(09,29) + rrt(221) * density(08) 
  pd(28,08) = pd(28,08) + rrt(221) * density(29) 
  pd(28,29) = pd(28,29) + rrt(221) * density(08) 
  pd(29,08) = pd(29,08) - rrt(221) * density(29) 
  pd(29,29) = pd(29,29) - rrt(221) * density(08) 
  pd(08,09) = pd(08,09) + rrt(222) * density(28) 
  pd(08,28) = pd(08,28) + rrt(222) * density(09) 
  pd(09,09) = pd(09,09) - rrt(222) * density(28) 
  pd(09,28) = pd(09,28) - rrt(222) * density(09) 
  pd(28,09) = pd(28,09) - rrt(222) * density(28) 
  pd(28,28) = pd(28,28) - rrt(222) * density(09) 
  pd(29,09) = pd(29,09) + rrt(222) * density(28) 
  pd(29,28) = pd(29,28) + rrt(222) * density(09) 
  pd(01,03) = pd(01,03) + rrt(223) * density(21) 
  pd(01,21) = pd(01,21) + rrt(223) * density(03) 
  pd(03,03) = pd(03,03) - rrt(223) * density(21) 
  pd(03,21) = pd(03,21) - rrt(223) * density(03) 
  pd(21,03) = pd(21,03) - rrt(223) * density(21) 
  pd(21,21) = pd(21,21) - rrt(223) * density(03) 
  pd(27,03) = pd(27,03) + rrt(223) * density(21) 
  pd(27,21) = pd(27,21) + rrt(223) * density(03) 
  pd(01,01) = pd(01,01) - rrt(224) * density(27) 
  pd(01,27) = pd(01,27) - rrt(224) * density(01) 
  pd(03,01) = pd(03,01) + rrt(224) * density(27) 
  pd(03,27) = pd(03,27) + rrt(224) * density(01) 
  pd(21,01) = pd(21,01) + rrt(224) * density(27) 
  pd(21,27) = pd(21,27) + rrt(224) * density(01) 
  pd(27,01) = pd(27,01) - rrt(224) * density(27) 
  pd(27,27) = pd(27,27) - rrt(224) * density(01) 
  pd(02,04) = pd(02,04) + rrt(225) * density(21) 
  pd(02,21) = pd(02,21) + rrt(225) * density(04) 
  pd(04,04) = pd(04,04) - rrt(225) * density(21) 
  pd(04,21) = pd(04,21) - rrt(225) * density(04) 
  pd(21,04) = pd(21,04) - rrt(225) * density(21) 
  pd(21,21) = pd(21,21) - rrt(225) * density(04) 
  pd(27,04) = pd(27,04) + rrt(225) * density(21) 
  pd(27,21) = pd(27,21) + rrt(225) * density(04) 
  pd(02,02) = pd(02,02) - rrt(226) * density(27) 
  pd(02,27) = pd(02,27) - rrt(226) * density(02) 
  pd(04,02) = pd(04,02) + rrt(226) * density(27) 
  pd(04,27) = pd(04,27) + rrt(226) * density(02) 
  pd(21,02) = pd(21,02) + rrt(226) * density(27) 
  pd(21,27) = pd(21,27) + rrt(226) * density(02) 
  pd(27,02) = pd(27,02) - rrt(226) * density(27) 
  pd(27,27) = pd(27,27) - rrt(226) * density(02) 
  pd(03,05) = pd(03,05) + rrt(227) * density(21) 
  pd(03,21) = pd(03,21) + rrt(227) * density(05) 
  pd(05,05) = pd(05,05) - rrt(227) * density(21) 
  pd(05,21) = pd(05,21) - rrt(227) * density(05) 
  pd(21,05) = pd(21,05) - rrt(227) * density(21) 
  pd(21,21) = pd(21,21) - rrt(227) * density(05) 
  pd(27,05) = pd(27,05) + rrt(227) * density(21) 
  pd(27,21) = pd(27,21) + rrt(227) * density(05) 
  pd(03,03) = pd(03,03) - rrt(228) * density(27) 
  pd(03,27) = pd(03,27) - rrt(228) * density(03) 
  pd(05,03) = pd(05,03) + rrt(228) * density(27) 
  pd(05,27) = pd(05,27) + rrt(228) * density(03) 
  pd(21,03) = pd(21,03) + rrt(228) * density(27) 
  pd(21,27) = pd(21,27) + rrt(228) * density(03) 
  pd(27,03) = pd(27,03) - rrt(228) * density(27) 
  pd(27,27) = pd(27,27) - rrt(228) * density(03) 
  pd(04,06) = pd(04,06) + rrt(229) * density(21) 
  pd(04,21) = pd(04,21) + rrt(229) * density(06) 
  pd(06,06) = pd(06,06) - rrt(229) * density(21) 
  pd(06,21) = pd(06,21) - rrt(229) * density(06) 
  pd(21,06) = pd(21,06) - rrt(229) * density(21) 
  pd(21,21) = pd(21,21) - rrt(229) * density(06) 
  pd(27,06) = pd(27,06) + rrt(229) * density(21) 
  pd(27,21) = pd(27,21) + rrt(229) * density(06) 
  pd(04,04) = pd(04,04) - rrt(230) * density(27) 
  pd(04,27) = pd(04,27) - rrt(230) * density(04) 
  pd(06,04) = pd(06,04) + rrt(230) * density(27) 
  pd(06,27) = pd(06,27) + rrt(230) * density(04) 
  pd(21,04) = pd(21,04) + rrt(230) * density(27) 
  pd(21,27) = pd(21,27) + rrt(230) * density(04) 
  pd(27,04) = pd(27,04) - rrt(230) * density(27) 
  pd(27,27) = pd(27,27) - rrt(230) * density(04) 
  pd(05,07) = pd(05,07) + rrt(231) * density(21) 
  pd(05,21) = pd(05,21) + rrt(231) * density(07) 
  pd(07,07) = pd(07,07) - rrt(231) * density(21) 
  pd(07,21) = pd(07,21) - rrt(231) * density(07) 
  pd(21,07) = pd(21,07) - rrt(231) * density(21) 
  pd(21,21) = pd(21,21) - rrt(231) * density(07) 
  pd(27,07) = pd(27,07) + rrt(231) * density(21) 
  pd(27,21) = pd(27,21) + rrt(231) * density(07) 
  pd(05,05) = pd(05,05) - rrt(232) * density(27) 
  pd(05,27) = pd(05,27) - rrt(232) * density(05) 
  pd(07,05) = pd(07,05) + rrt(232) * density(27) 
  pd(07,27) = pd(07,27) + rrt(232) * density(05) 
  pd(21,05) = pd(21,05) + rrt(232) * density(27) 
  pd(21,27) = pd(21,27) + rrt(232) * density(05) 
  pd(27,05) = pd(27,05) - rrt(232) * density(27) 
  pd(27,27) = pd(27,27) - rrt(232) * density(05) 
  pd(06,08) = pd(06,08) + rrt(233) * density(21) 
  pd(06,21) = pd(06,21) + rrt(233) * density(08) 
  pd(08,08) = pd(08,08) - rrt(233) * density(21) 
  pd(08,21) = pd(08,21) - rrt(233) * density(08) 
  pd(21,08) = pd(21,08) - rrt(233) * density(21) 
  pd(21,21) = pd(21,21) - rrt(233) * density(08) 
  pd(27,08) = pd(27,08) + rrt(233) * density(21) 
  pd(27,21) = pd(27,21) + rrt(233) * density(08) 
  pd(06,06) = pd(06,06) - rrt(234) * density(27) 
  pd(06,27) = pd(06,27) - rrt(234) * density(06) 
  pd(08,06) = pd(08,06) + rrt(234) * density(27) 
  pd(08,27) = pd(08,27) + rrt(234) * density(06) 
  pd(21,06) = pd(21,06) + rrt(234) * density(27) 
  pd(21,27) = pd(21,27) + rrt(234) * density(06) 
  pd(27,06) = pd(27,06) - rrt(234) * density(27) 
  pd(27,27) = pd(27,27) - rrt(234) * density(06) 
  pd(07,09) = pd(07,09) + rrt(235) * density(21) 
  pd(07,21) = pd(07,21) + rrt(235) * density(09) 
  pd(09,09) = pd(09,09) - rrt(235) * density(21) 
  pd(09,21) = pd(09,21) - rrt(235) * density(09) 
  pd(21,09) = pd(21,09) - rrt(235) * density(21) 
  pd(21,21) = pd(21,21) - rrt(235) * density(09) 
  pd(27,09) = pd(27,09) + rrt(235) * density(21) 
  pd(27,21) = pd(27,21) + rrt(235) * density(09) 
  pd(07,07) = pd(07,07) - rrt(236) * density(27) 
  pd(07,27) = pd(07,27) - rrt(236) * density(07) 
  pd(09,07) = pd(09,07) + rrt(236) * density(27) 
  pd(09,27) = pd(09,27) + rrt(236) * density(07) 
  pd(21,07) = pd(21,07) + rrt(236) * density(27) 
  pd(21,27) = pd(21,27) + rrt(236) * density(07) 
  pd(27,07) = pd(27,07) - rrt(236) * density(27) 
  pd(27,27) = pd(27,27) - rrt(236) * density(07) 
  pd(01,03) = pd(01,03) + rrt(237) * density(27) 
  pd(01,27) = pd(01,27) + rrt(237) * density(03) 
  pd(03,03) = pd(03,03) - rrt(237) * density(27) 
  pd(03,27) = pd(03,27) - rrt(237) * density(03) 
  pd(27,03) = pd(27,03) - rrt(237) * density(27) 
  pd(27,27) = pd(27,27) - rrt(237) * density(03) 
  pd(28,03) = pd(28,03) + rrt(237) * density(27) 
  pd(28,27) = pd(28,27) + rrt(237) * density(03) 
  pd(01,01) = pd(01,01) - rrt(238) * density(28) 
  pd(01,28) = pd(01,28) - rrt(238) * density(01) 
  pd(03,01) = pd(03,01) + rrt(238) * density(28) 
  pd(03,28) = pd(03,28) + rrt(238) * density(01) 
  pd(27,01) = pd(27,01) + rrt(238) * density(28) 
  pd(27,28) = pd(27,28) + rrt(238) * density(01) 
  pd(28,01) = pd(28,01) - rrt(238) * density(28) 
  pd(28,28) = pd(28,28) - rrt(238) * density(01) 
  pd(02,04) = pd(02,04) + rrt(239) * density(27) 
  pd(02,27) = pd(02,27) + rrt(239) * density(04) 
  pd(04,04) = pd(04,04) - rrt(239) * density(27) 
  pd(04,27) = pd(04,27) - rrt(239) * density(04) 
  pd(27,04) = pd(27,04) - rrt(239) * density(27) 
  pd(27,27) = pd(27,27) - rrt(239) * density(04) 
  pd(28,04) = pd(28,04) + rrt(239) * density(27) 
  pd(28,27) = pd(28,27) + rrt(239) * density(04) 
  pd(02,02) = pd(02,02) - rrt(240) * density(28) 
  pd(02,28) = pd(02,28) - rrt(240) * density(02) 
  pd(04,02) = pd(04,02) + rrt(240) * density(28) 
  pd(04,28) = pd(04,28) + rrt(240) * density(02) 
  pd(27,02) = pd(27,02) + rrt(240) * density(28) 
  pd(27,28) = pd(27,28) + rrt(240) * density(02) 
  pd(28,02) = pd(28,02) - rrt(240) * density(28) 
  pd(28,28) = pd(28,28) - rrt(240) * density(02) 
  pd(03,05) = pd(03,05) + rrt(241) * density(27) 
  pd(03,27) = pd(03,27) + rrt(241) * density(05) 
  pd(05,05) = pd(05,05) - rrt(241) * density(27) 
  pd(05,27) = pd(05,27) - rrt(241) * density(05) 
  pd(27,05) = pd(27,05) - rrt(241) * density(27) 
  pd(27,27) = pd(27,27) - rrt(241) * density(05) 
  pd(28,05) = pd(28,05) + rrt(241) * density(27) 
  pd(28,27) = pd(28,27) + rrt(241) * density(05) 
  pd(03,03) = pd(03,03) - rrt(242) * density(28) 
  pd(03,28) = pd(03,28) - rrt(242) * density(03) 
  pd(05,03) = pd(05,03) + rrt(242) * density(28) 
  pd(05,28) = pd(05,28) + rrt(242) * density(03) 
  pd(27,03) = pd(27,03) + rrt(242) * density(28) 
  pd(27,28) = pd(27,28) + rrt(242) * density(03) 
  pd(28,03) = pd(28,03) - rrt(242) * density(28) 
  pd(28,28) = pd(28,28) - rrt(242) * density(03) 
  pd(04,06) = pd(04,06) + rrt(243) * density(27) 
  pd(04,27) = pd(04,27) + rrt(243) * density(06) 
  pd(06,06) = pd(06,06) - rrt(243) * density(27) 
  pd(06,27) = pd(06,27) - rrt(243) * density(06) 
  pd(27,06) = pd(27,06) - rrt(243) * density(27) 
  pd(27,27) = pd(27,27) - rrt(243) * density(06) 
  pd(28,06) = pd(28,06) + rrt(243) * density(27) 
  pd(28,27) = pd(28,27) + rrt(243) * density(06) 
  pd(04,04) = pd(04,04) - rrt(244) * density(28) 
  pd(04,28) = pd(04,28) - rrt(244) * density(04) 
  pd(06,04) = pd(06,04) + rrt(244) * density(28) 
  pd(06,28) = pd(06,28) + rrt(244) * density(04) 
  pd(27,04) = pd(27,04) + rrt(244) * density(28) 
  pd(27,28) = pd(27,28) + rrt(244) * density(04) 
  pd(28,04) = pd(28,04) - rrt(244) * density(28) 
  pd(28,28) = pd(28,28) - rrt(244) * density(04) 
  pd(05,07) = pd(05,07) + rrt(245) * density(27) 
  pd(05,27) = pd(05,27) + rrt(245) * density(07) 
  pd(07,07) = pd(07,07) - rrt(245) * density(27) 
  pd(07,27) = pd(07,27) - rrt(245) * density(07) 
  pd(27,07) = pd(27,07) - rrt(245) * density(27) 
  pd(27,27) = pd(27,27) - rrt(245) * density(07) 
  pd(28,07) = pd(28,07) + rrt(245) * density(27) 
  pd(28,27) = pd(28,27) + rrt(245) * density(07) 
  pd(05,05) = pd(05,05) - rrt(246) * density(28) 
  pd(05,28) = pd(05,28) - rrt(246) * density(05) 
  pd(07,05) = pd(07,05) + rrt(246) * density(28) 
  pd(07,28) = pd(07,28) + rrt(246) * density(05) 
  pd(27,05) = pd(27,05) + rrt(246) * density(28) 
  pd(27,28) = pd(27,28) + rrt(246) * density(05) 
  pd(28,05) = pd(28,05) - rrt(246) * density(28) 
  pd(28,28) = pd(28,28) - rrt(246) * density(05) 
  pd(06,08) = pd(06,08) + rrt(247) * density(27) 
  pd(06,27) = pd(06,27) + rrt(247) * density(08) 
  pd(08,08) = pd(08,08) - rrt(247) * density(27) 
  pd(08,27) = pd(08,27) - rrt(247) * density(08) 
  pd(27,08) = pd(27,08) - rrt(247) * density(27) 
  pd(27,27) = pd(27,27) - rrt(247) * density(08) 
  pd(28,08) = pd(28,08) + rrt(247) * density(27) 
  pd(28,27) = pd(28,27) + rrt(247) * density(08) 
  pd(06,06) = pd(06,06) - rrt(248) * density(28) 
  pd(06,28) = pd(06,28) - rrt(248) * density(06) 
  pd(08,06) = pd(08,06) + rrt(248) * density(28) 
  pd(08,28) = pd(08,28) + rrt(248) * density(06) 
  pd(27,06) = pd(27,06) + rrt(248) * density(28) 
  pd(27,28) = pd(27,28) + rrt(248) * density(06) 
  pd(28,06) = pd(28,06) - rrt(248) * density(28) 
  pd(28,28) = pd(28,28) - rrt(248) * density(06) 
  pd(07,09) = pd(07,09) + rrt(249) * density(27) 
  pd(07,27) = pd(07,27) + rrt(249) * density(09) 
  pd(09,09) = pd(09,09) - rrt(249) * density(27) 
  pd(09,27) = pd(09,27) - rrt(249) * density(09) 
  pd(27,09) = pd(27,09) - rrt(249) * density(27) 
  pd(27,27) = pd(27,27) - rrt(249) * density(09) 
  pd(28,09) = pd(28,09) + rrt(249) * density(27) 
  pd(28,27) = pd(28,27) + rrt(249) * density(09) 
  pd(07,07) = pd(07,07) - rrt(250) * density(28) 
  pd(07,28) = pd(07,28) - rrt(250) * density(07) 
  pd(09,07) = pd(09,07) + rrt(250) * density(28) 
  pd(09,28) = pd(09,28) + rrt(250) * density(07) 
  pd(27,07) = pd(27,07) + rrt(250) * density(28) 
  pd(27,28) = pd(27,28) + rrt(250) * density(07) 
  pd(28,07) = pd(28,07) - rrt(250) * density(28) 
  pd(28,28) = pd(28,28) - rrt(250) * density(07) 
  pd(01,03) = pd(01,03) + rrt(251) * density(28) 
  pd(01,28) = pd(01,28) + rrt(251) * density(03) 
  pd(03,03) = pd(03,03) - rrt(251) * density(28) 
  pd(03,28) = pd(03,28) - rrt(251) * density(03) 
  pd(28,03) = pd(28,03) - rrt(251) * density(28) 
  pd(28,28) = pd(28,28) - rrt(251) * density(03) 
  pd(29,03) = pd(29,03) + rrt(251) * density(28) 
  pd(29,28) = pd(29,28) + rrt(251) * density(03) 
  pd(01,01) = pd(01,01) - rrt(252) * density(29) 
  pd(01,29) = pd(01,29) - rrt(252) * density(01) 
  pd(03,01) = pd(03,01) + rrt(252) * density(29) 
  pd(03,29) = pd(03,29) + rrt(252) * density(01) 
  pd(28,01) = pd(28,01) + rrt(252) * density(29) 
  pd(28,29) = pd(28,29) + rrt(252) * density(01) 
  pd(29,01) = pd(29,01) - rrt(252) * density(29) 
  pd(29,29) = pd(29,29) - rrt(252) * density(01) 
  pd(02,04) = pd(02,04) + rrt(253) * density(28) 
  pd(02,28) = pd(02,28) + rrt(253) * density(04) 
  pd(04,04) = pd(04,04) - rrt(253) * density(28) 
  pd(04,28) = pd(04,28) - rrt(253) * density(04) 
  pd(28,04) = pd(28,04) - rrt(253) * density(28) 
  pd(28,28) = pd(28,28) - rrt(253) * density(04) 
  pd(29,04) = pd(29,04) + rrt(253) * density(28) 
  pd(29,28) = pd(29,28) + rrt(253) * density(04) 
  pd(02,02) = pd(02,02) - rrt(254) * density(29) 
  pd(02,29) = pd(02,29) - rrt(254) * density(02) 
  pd(04,02) = pd(04,02) + rrt(254) * density(29) 
  pd(04,29) = pd(04,29) + rrt(254) * density(02) 
  pd(28,02) = pd(28,02) + rrt(254) * density(29) 
  pd(28,29) = pd(28,29) + rrt(254) * density(02) 
  pd(29,02) = pd(29,02) - rrt(254) * density(29) 
  pd(29,29) = pd(29,29) - rrt(254) * density(02) 
  pd(03,05) = pd(03,05) + rrt(255) * density(28) 
  pd(03,28) = pd(03,28) + rrt(255) * density(05) 
  pd(05,05) = pd(05,05) - rrt(255) * density(28) 
  pd(05,28) = pd(05,28) - rrt(255) * density(05) 
  pd(28,05) = pd(28,05) - rrt(255) * density(28) 
  pd(28,28) = pd(28,28) - rrt(255) * density(05) 
  pd(29,05) = pd(29,05) + rrt(255) * density(28) 
  pd(29,28) = pd(29,28) + rrt(255) * density(05) 
  pd(03,03) = pd(03,03) - rrt(256) * density(29) 
  pd(03,29) = pd(03,29) - rrt(256) * density(03) 
  pd(05,03) = pd(05,03) + rrt(256) * density(29) 
  pd(05,29) = pd(05,29) + rrt(256) * density(03) 
  pd(28,03) = pd(28,03) + rrt(256) * density(29) 
  pd(28,29) = pd(28,29) + rrt(256) * density(03) 
  pd(29,03) = pd(29,03) - rrt(256) * density(29) 
  pd(29,29) = pd(29,29) - rrt(256) * density(03) 
  pd(04,06) = pd(04,06) + rrt(257) * density(28) 
  pd(04,28) = pd(04,28) + rrt(257) * density(06) 
  pd(06,06) = pd(06,06) - rrt(257) * density(28) 
  pd(06,28) = pd(06,28) - rrt(257) * density(06) 
  pd(28,06) = pd(28,06) - rrt(257) * density(28) 
  pd(28,28) = pd(28,28) - rrt(257) * density(06) 
  pd(29,06) = pd(29,06) + rrt(257) * density(28) 
  pd(29,28) = pd(29,28) + rrt(257) * density(06) 
  pd(04,04) = pd(04,04) - rrt(258) * density(29) 
  pd(04,29) = pd(04,29) - rrt(258) * density(04) 
  pd(06,04) = pd(06,04) + rrt(258) * density(29) 
  pd(06,29) = pd(06,29) + rrt(258) * density(04) 
  pd(28,04) = pd(28,04) + rrt(258) * density(29) 
  pd(28,29) = pd(28,29) + rrt(258) * density(04) 
  pd(29,04) = pd(29,04) - rrt(258) * density(29) 
  pd(29,29) = pd(29,29) - rrt(258) * density(04) 
  pd(05,07) = pd(05,07) + rrt(259) * density(28) 
  pd(05,28) = pd(05,28) + rrt(259) * density(07) 
  pd(07,07) = pd(07,07) - rrt(259) * density(28) 
  pd(07,28) = pd(07,28) - rrt(259) * density(07) 
  pd(28,07) = pd(28,07) - rrt(259) * density(28) 
  pd(28,28) = pd(28,28) - rrt(259) * density(07) 
  pd(29,07) = pd(29,07) + rrt(259) * density(28) 
  pd(29,28) = pd(29,28) + rrt(259) * density(07) 
  pd(05,05) = pd(05,05) - rrt(260) * density(29) 
  pd(05,29) = pd(05,29) - rrt(260) * density(05) 
  pd(07,05) = pd(07,05) + rrt(260) * density(29) 
  pd(07,29) = pd(07,29) + rrt(260) * density(05) 
  pd(28,05) = pd(28,05) + rrt(260) * density(29) 
  pd(28,29) = pd(28,29) + rrt(260) * density(05) 
  pd(29,05) = pd(29,05) - rrt(260) * density(29) 
  pd(29,29) = pd(29,29) - rrt(260) * density(05) 
  pd(06,08) = pd(06,08) + rrt(261) * density(28) 
  pd(06,28) = pd(06,28) + rrt(261) * density(08) 
  pd(08,08) = pd(08,08) - rrt(261) * density(28) 
  pd(08,28) = pd(08,28) - rrt(261) * density(08) 
  pd(28,08) = pd(28,08) - rrt(261) * density(28) 
  pd(28,28) = pd(28,28) - rrt(261) * density(08) 
  pd(29,08) = pd(29,08) + rrt(261) * density(28) 
  pd(29,28) = pd(29,28) + rrt(261) * density(08) 
  pd(06,06) = pd(06,06) - rrt(262) * density(29) 
  pd(06,29) = pd(06,29) - rrt(262) * density(06) 
  pd(08,06) = pd(08,06) + rrt(262) * density(29) 
  pd(08,29) = pd(08,29) + rrt(262) * density(06) 
  pd(28,06) = pd(28,06) + rrt(262) * density(29) 
  pd(28,29) = pd(28,29) + rrt(262) * density(06) 
  pd(29,06) = pd(29,06) - rrt(262) * density(29) 
  pd(29,29) = pd(29,29) - rrt(262) * density(06) 
  pd(07,09) = pd(07,09) + rrt(263) * density(28) 
  pd(07,28) = pd(07,28) + rrt(263) * density(09) 
  pd(09,09) = pd(09,09) - rrt(263) * density(28) 
  pd(09,28) = pd(09,28) - rrt(263) * density(09) 
  pd(28,09) = pd(28,09) - rrt(263) * density(28) 
  pd(28,28) = pd(28,28) - rrt(263) * density(09) 
  pd(29,09) = pd(29,09) + rrt(263) * density(28) 
  pd(29,28) = pd(29,28) + rrt(263) * density(09) 
  pd(07,07) = pd(07,07) - rrt(264) * density(29) 
  pd(07,29) = pd(07,29) - rrt(264) * density(07) 
  pd(09,07) = pd(09,07) + rrt(264) * density(29) 
  pd(09,29) = pd(09,29) + rrt(264) * density(07) 
  pd(28,07) = pd(28,07) + rrt(264) * density(29) 
  pd(28,29) = pd(28,29) + rrt(264) * density(07) 
  pd(29,07) = pd(29,07) - rrt(264) * density(29) 
  pd(29,29) = pd(29,29) - rrt(264) * density(07) 
  pd(27,28) = pd(27,28) + rrt(265) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(265) * density(28) * 4.0d0
  pd(29,28) = pd(29,28) + rrt(265) * density(28) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(266) * density(29) 
  pd(27,29) = pd(27,29) - rrt(266) * density(27) 
  pd(28,27) = pd(28,27) + rrt(266) * density(29) * 2.0d0
  pd(28,29) = pd(28,29) + rrt(266) * density(27) * 2.0d0
  pd(29,27) = pd(29,27) - rrt(266) * density(29) 
  pd(29,29) = pd(29,29) - rrt(266) * density(27) 
  pd(21,27) = pd(21,27) + rrt(267) * density(28) 
  pd(21,28) = pd(21,28) + rrt(267) * density(27) 
  pd(27,27) = pd(27,27) - rrt(267) * density(28) 
  pd(27,28) = pd(27,28) - rrt(267) * density(27) 
  pd(28,27) = pd(28,27) - rrt(267) * density(28) 
  pd(28,28) = pd(28,28) - rrt(267) * density(27) 
  pd(29,27) = pd(29,27) + rrt(267) * density(28) 
  pd(29,28) = pd(29,28) + rrt(267) * density(27) 
  pd(21,21) = pd(21,21) - rrt(268) * density(29) 
  pd(21,29) = pd(21,29) - rrt(268) * density(21) 
  pd(27,21) = pd(27,21) + rrt(268) * density(29) 
  pd(27,29) = pd(27,29) + rrt(268) * density(21) 
  pd(28,21) = pd(28,21) + rrt(268) * density(29) 
  pd(28,29) = pd(28,29) + rrt(268) * density(21) 
  pd(29,21) = pd(29,21) - rrt(268) * density(29) 
  pd(29,29) = pd(29,29) - rrt(268) * density(21) 
  pd(21,27) = pd(21,27) + rrt(269) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(269) * density(27) * 4.0d0
  pd(28,27) = pd(28,27) + rrt(269) * density(27) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(270) * density(28) 
  pd(21,28) = pd(21,28) - rrt(270) * density(21) 
  pd(27,21) = pd(27,21) + rrt(270) * density(28) * 2.0d0
  pd(27,28) = pd(27,28) + rrt(270) * density(21) * 2.0d0
  pd(28,21) = pd(28,21) - rrt(270) * density(28) 
  pd(28,28) = pd(28,28) - rrt(270) * density(21) 
  pd(01,07) = pd(01,07) + rrt(271) * density(10) 
  pd(01,10) = pd(01,10) + rrt(271) * density(07) 
  pd(07,07) = pd(07,07) - rrt(271) * density(10) 
  pd(07,10) = pd(07,10) - rrt(271) * density(07) 
  pd(10,07) = pd(10,07) - rrt(271) * density(10) 
  pd(10,10) = pd(10,10) - rrt(271) * density(07) 
  pd(11,07) = pd(11,07) + rrt(271) * density(10) 
  pd(11,10) = pd(11,10) + rrt(271) * density(07) 
  pd(02,08) = pd(02,08) + rrt(272) * density(10) 
  pd(02,10) = pd(02,10) + rrt(272) * density(08) 
  pd(08,08) = pd(08,08) - rrt(272) * density(10) 
  pd(08,10) = pd(08,10) - rrt(272) * density(08) 
  pd(10,08) = pd(10,08) - rrt(272) * density(10) 
  pd(10,10) = pd(10,10) - rrt(272) * density(08) 
  pd(11,08) = pd(11,08) + rrt(272) * density(10) 
  pd(11,10) = pd(11,10) + rrt(272) * density(08) 
  pd(03,09) = pd(03,09) + rrt(273) * density(10) 
  pd(03,10) = pd(03,10) + rrt(273) * density(09) 
  pd(09,09) = pd(09,09) - rrt(273) * density(10) 
  pd(09,10) = pd(09,10) - rrt(273) * density(09) 
  pd(10,09) = pd(10,09) - rrt(273) * density(10) 
  pd(10,10) = pd(10,10) - rrt(273) * density(09) 
  pd(11,09) = pd(11,09) + rrt(273) * density(10) 
  pd(11,10) = pd(11,10) + rrt(273) * density(09) 
  pd(01,01) = pd(01,01) - rrt(274) * density(11) 
  pd(01,11) = pd(01,11) - rrt(274) * density(01) 
  pd(07,01) = pd(07,01) + rrt(274) * density(11) 
  pd(07,11) = pd(07,11) + rrt(274) * density(01) 
  pd(10,01) = pd(10,01) + rrt(274) * density(11) 
  pd(10,11) = pd(10,11) + rrt(274) * density(01) 
  pd(11,01) = pd(11,01) - rrt(274) * density(11) 
  pd(11,11) = pd(11,11) - rrt(274) * density(01) 
  pd(02,02) = pd(02,02) - rrt(275) * density(11) 
  pd(02,11) = pd(02,11) - rrt(275) * density(02) 
  pd(08,02) = pd(08,02) + rrt(275) * density(11) 
  pd(08,11) = pd(08,11) + rrt(275) * density(02) 
  pd(10,02) = pd(10,02) + rrt(275) * density(11) 
  pd(10,11) = pd(10,11) + rrt(275) * density(02) 
  pd(11,02) = pd(11,02) - rrt(275) * density(11) 
  pd(11,11) = pd(11,11) - rrt(275) * density(02) 
  pd(03,03) = pd(03,03) - rrt(276) * density(11) 
  pd(03,11) = pd(03,11) - rrt(276) * density(03) 
  pd(09,03) = pd(09,03) + rrt(276) * density(11) 
  pd(09,11) = pd(09,11) + rrt(276) * density(03) 
  pd(10,03) = pd(10,03) + rrt(276) * density(11) 
  pd(10,11) = pd(10,11) + rrt(276) * density(03) 
  pd(11,03) = pd(11,03) - rrt(276) * density(11) 
  pd(11,11) = pd(11,11) - rrt(276) * density(03) 
  pd(04,04) = pd(04,04) - rrt(277) * density(11) 
  pd(04,11) = pd(04,11) - rrt(277) * density(04) 
  pd(09,04) = pd(09,04) + rrt(277) * density(11) 
  pd(09,11) = pd(09,11) + rrt(277) * density(04) 
  pd(10,04) = pd(10,04) + rrt(277) * density(11) 
  pd(10,11) = pd(10,11) + rrt(277) * density(04) 
  pd(11,04) = pd(11,04) - rrt(277) * density(11) 
  pd(11,11) = pd(11,11) - rrt(277) * density(04) 
  pd(05,05) = pd(05,05) - rrt(278) * density(11) 
  pd(05,11) = pd(05,11) - rrt(278) * density(05) 
  pd(09,05) = pd(09,05) + rrt(278) * density(11) 
  pd(09,11) = pd(09,11) + rrt(278) * density(05) 
  pd(10,05) = pd(10,05) + rrt(278) * density(11) 
  pd(10,11) = pd(10,11) + rrt(278) * density(05) 
  pd(11,05) = pd(11,05) - rrt(278) * density(11) 
  pd(11,11) = pd(11,11) - rrt(278) * density(05) 
  pd(06,06) = pd(06,06) - rrt(279) * density(11) 
  pd(06,11) = pd(06,11) - rrt(279) * density(06) 
  pd(09,06) = pd(09,06) + rrt(279) * density(11) 
  pd(09,11) = pd(09,11) + rrt(279) * density(06) 
  pd(10,06) = pd(10,06) + rrt(279) * density(11) 
  pd(10,11) = pd(10,11) + rrt(279) * density(06) 
  pd(11,06) = pd(11,06) - rrt(279) * density(11) 
  pd(11,11) = pd(11,11) - rrt(279) * density(06) 
  pd(07,07) = pd(07,07) - rrt(280) * density(11) 
  pd(07,11) = pd(07,11) - rrt(280) * density(07) 
  pd(09,07) = pd(09,07) + rrt(280) * density(11) 
  pd(09,11) = pd(09,11) + rrt(280) * density(07) 
  pd(10,07) = pd(10,07) + rrt(280) * density(11) 
  pd(10,11) = pd(10,11) + rrt(280) * density(07) 
  pd(11,07) = pd(11,07) - rrt(280) * density(11) 
  pd(11,11) = pd(11,11) - rrt(280) * density(07) 
  pd(08,08) = pd(08,08) - rrt(281) * density(11) 
  pd(08,11) = pd(08,11) - rrt(281) * density(08) 
  pd(09,08) = pd(09,08) + rrt(281) * density(11) 
  pd(09,11) = pd(09,11) + rrt(281) * density(08) 
  pd(10,08) = pd(10,08) + rrt(281) * density(11) 
  pd(10,11) = pd(10,11) + rrt(281) * density(08) 
  pd(11,08) = pd(11,08) - rrt(281) * density(11) 
  pd(11,11) = pd(11,11) - rrt(281) * density(08) 
  pd(10,09) = pd(10,09) + rrt(282) * density(11) 
  pd(10,11) = pd(10,11) + rrt(282) * density(09) 
  pd(11,09) = pd(11,09) - rrt(282) * density(11) 
  pd(11,11) = pd(11,11) - rrt(282) * density(09) 
  pd(10,11) = pd(10,11) + rrt(283) 
  pd(11,11) = pd(11,11) - rrt(283) 
  pd(01,12) = pd(01,12) + rrt(284) 
  pd(12,12) = pd(12,12) - rrt(284) 
  pd(11,13) = pd(11,13) + rrt(285) 
  pd(13,13) = pd(13,13) - rrt(285) 
  pd(21,21) = pd(21,21) - rrt(286) * density(43) 
  pd(21,43) = pd(21,43) - rrt(286) * density(21) 
  pd(30,21) = pd(30,21) + rrt(286) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) + rrt(286) * density(21) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(287) * density(43) 
  pd(01,43) = pd(01,43) - rrt(287) * density(01) 
  pd(14,01) = pd(14,01) + rrt(287) * density(43) * 2.0d0
  pd(14,43) = pd(14,43) + rrt(287) * density(01) * 2.0d0
  pd(14,35) = pd(14,35) + rrt(288) * density(43) 
  pd(14,43) = pd(14,43) + rrt(288) * density(35) 
  pd(30,35) = pd(30,35) + rrt(288) * density(43) 
  pd(30,43) = pd(30,43) + rrt(288) * density(35) 
  pd(35,35) = pd(35,35) - rrt(288) * density(43) 
  pd(35,43) = pd(35,43) - rrt(288) * density(35) 
  pd(14,36) = pd(14,36) + rrt(289) * density(43) 
  pd(14,43) = pd(14,43) + rrt(289) * density(36) 
  pd(21,36) = pd(21,36) + rrt(289) * density(43) 
  pd(21,43) = pd(21,43) + rrt(289) * density(36) 
  pd(36,36) = pd(36,36) - rrt(289) * density(43) 
  pd(36,43) = pd(36,43) - rrt(289) * density(36) 
  pd(30,36) = pd(30,36) + rrt(290) * density(43) 
  pd(30,43) = pd(30,43) + rrt(290) * density(36) 
  pd(35,36) = pd(35,36) + rrt(290) * density(43) 
  pd(35,43) = pd(35,43) + rrt(290) * density(36) 
  pd(36,36) = pd(36,36) - rrt(290) * density(43) 
  pd(36,43) = pd(36,43) - rrt(290) * density(36) 
  pd(30,37) = pd(30,37) + rrt(291) * density(43) 
  pd(30,43) = pd(30,43) + rrt(291) * density(37) 
  pd(36,37) = pd(36,37) + rrt(291) * density(43) 
  pd(36,43) = pd(36,43) + rrt(291) * density(37) 
  pd(37,37) = pd(37,37) - rrt(291) * density(43) 
  pd(37,43) = pd(37,43) - rrt(291) * density(37) 
  pd(21,37) = pd(21,37) + rrt(292) * density(43) 
  pd(21,43) = pd(21,43) + rrt(292) * density(37) 
  pd(35,37) = pd(35,37) + rrt(292) * density(43) 
  pd(35,43) = pd(35,43) + rrt(292) * density(37) 
  pd(37,37) = pd(37,37) - rrt(292) * density(43) 
  pd(37,43) = pd(37,43) - rrt(292) * density(37) 
  pd(14,18) = pd(14,18) + rrt(293) * density(43) * 2.0d0
  pd(14,43) = pd(14,43) + rrt(293) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(293) * density(43) 
  pd(18,43) = pd(18,43) - rrt(293) * density(18) 
  pd(43,18) = pd(43,18) - rrt(293) * density(43) 
  pd(43,43) = pd(43,43) - rrt(293) * density(18) 
  pd(14,18) = pd(14,18) + rrt(294) * density(43) 
  pd(14,43) = pd(14,43) + rrt(294) * density(18) 
  pd(15,18) = pd(15,18) + rrt(294) * density(43) 
  pd(15,43) = pd(15,43) + rrt(294) * density(18) 
  pd(18,18) = pd(18,18) - rrt(294) * density(43) 
  pd(18,43) = pd(18,43) - rrt(294) * density(18) 
  pd(43,18) = pd(43,18) - rrt(294) * density(43) 
  pd(43,43) = pd(43,43) - rrt(294) * density(18) 
  pd(14,18) = pd(14,18) + rrt(295) * density(43) 
  pd(14,43) = pd(14,43) + rrt(295) * density(18) 
  pd(16,18) = pd(16,18) + rrt(295) * density(43) 
  pd(16,43) = pd(16,43) + rrt(295) * density(18) 
  pd(18,18) = pd(18,18) - rrt(295) * density(43) 
  pd(18,43) = pd(18,43) - rrt(295) * density(18) 
  pd(43,18) = pd(43,18) - rrt(295) * density(43) 
  pd(43,43) = pd(43,43) - rrt(295) * density(18) 
  pd(01,19) = pd(01,19) + rrt(296) * density(43) 
  pd(01,43) = pd(01,43) + rrt(296) * density(19) 
  pd(14,19) = pd(14,19) + rrt(296) * density(43) 
  pd(14,43) = pd(14,43) + rrt(296) * density(19) 
  pd(19,19) = pd(19,19) - rrt(296) * density(43) 
  pd(19,43) = pd(19,43) - rrt(296) * density(19) 
  pd(43,19) = pd(43,19) - rrt(296) * density(43) 
  pd(43,43) = pd(43,43) - rrt(296) * density(19) 
  pd(01,20) = pd(01,20) + rrt(297) * density(43) * 2.0d0
  pd(01,43) = pd(01,43) + rrt(297) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(297) * density(43) 
  pd(20,43) = pd(20,43) - rrt(297) * density(20) 
  pd(43,20) = pd(43,20) - rrt(297) * density(43) 
  pd(43,43) = pd(43,43) - rrt(297) * density(20) 
  pd(21,21) = pd(21,21) - rrt(298) * density(43) 
  pd(21,43) = pd(21,43) - rrt(298) * density(21) 
  pd(30,21) = pd(30,21) + rrt(298) * density(43) 
  pd(30,43) = pd(30,43) + rrt(298) * density(21) 
  pd(33,21) = pd(33,21) + rrt(298) * density(43) 
  pd(33,43) = pd(33,43) + rrt(298) * density(21) 
  pd(43,21) = pd(43,21) + rrt(298) * density(43) 
  pd(43,43) = pd(43,43) + rrt(298) * density(21) 
  pd(30,31) = pd(30,31) + rrt(299) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) + rrt(299) * density(31) * 2.0d0
  pd(31,31) = pd(31,31) - rrt(299) * density(43) 
  pd(31,43) = pd(31,43) - rrt(299) * density(31) 
  pd(43,31) = pd(43,31) - rrt(299) * density(43) 
  pd(43,43) = pd(43,43) - rrt(299) * density(31) 
  pd(30,32) = pd(30,32) + rrt(300) * density(43) * 3.0d0
  pd(30,43) = pd(30,43) + rrt(300) * density(32) * 3.0d0
  pd(32,32) = pd(32,32) - rrt(300) * density(43) 
  pd(32,43) = pd(32,43) - rrt(300) * density(32) 
  pd(43,32) = pd(43,32) - rrt(300) * density(43) 
  pd(43,43) = pd(43,43) - rrt(300) * density(32) 
  pd(21,32) = pd(21,32) + rrt(301) * density(43) 
  pd(21,43) = pd(21,43) + rrt(301) * density(32) 
  pd(30,32) = pd(30,32) + rrt(301) * density(43) 
  pd(30,43) = pd(30,43) + rrt(301) * density(32) 
  pd(32,32) = pd(32,32) - rrt(301) * density(43) 
  pd(32,43) = pd(32,43) - rrt(301) * density(32) 
  pd(43,32) = pd(43,32) - rrt(301) * density(43) 
  pd(43,43) = pd(43,43) - rrt(301) * density(32) 
  pd(14,38) = pd(14,38) + rrt(302) * density(43) 
  pd(14,43) = pd(14,43) + rrt(302) * density(38) 
  pd(30,38) = pd(30,38) + rrt(302) * density(43) 
  pd(30,43) = pd(30,43) + rrt(302) * density(38) 
  pd(38,38) = pd(38,38) - rrt(302) * density(43) 
  pd(38,43) = pd(38,43) - rrt(302) * density(38) 
  pd(43,38) = pd(43,38) - rrt(302) * density(43) 
  pd(43,43) = pd(43,43) - rrt(302) * density(38) 
  pd(30,40) = pd(30,40) + rrt(303) * density(43) 
  pd(30,43) = pd(30,43) + rrt(303) * density(40) 
  pd(35,40) = pd(35,40) + rrt(303) * density(43) 
  pd(35,43) = pd(35,43) + rrt(303) * density(40) 
  pd(40,40) = pd(40,40) - rrt(303) * density(43) 
  pd(40,43) = pd(40,43) - rrt(303) * density(40) 
  pd(43,40) = pd(43,40) - rrt(303) * density(43) 
  pd(43,43) = pd(43,43) - rrt(303) * density(40) 
  pd(14,40) = pd(14,40) + rrt(304) * density(43) 
  pd(14,43) = pd(14,43) + rrt(304) * density(40) 
  pd(30,40) = pd(30,40) + rrt(304) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) + rrt(304) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(304) * density(43) 
  pd(40,43) = pd(40,43) - rrt(304) * density(40) 
  pd(43,40) = pd(43,40) - rrt(304) * density(43) 
  pd(43,43) = pd(43,43) - rrt(304) * density(40) 
  pd(30,41) = pd(30,41) + rrt(305) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) + rrt(305) * density(41) * 2.0d0
  pd(35,41) = pd(35,41) + rrt(305) * density(43) 
  pd(35,43) = pd(35,43) + rrt(305) * density(41) 
  pd(41,41) = pd(41,41) - rrt(305) * density(43) 
  pd(41,43) = pd(41,43) - rrt(305) * density(41) 
  pd(43,41) = pd(43,41) - rrt(305) * density(43) 
  pd(43,43) = pd(43,43) - rrt(305) * density(41) 
  pd(30,41) = pd(30,41) + rrt(306) * density(43) 
  pd(30,43) = pd(30,43) + rrt(306) * density(41) 
  pd(36,41) = pd(36,41) + rrt(306) * density(43) 
  pd(36,43) = pd(36,43) + rrt(306) * density(41) 
  pd(41,41) = pd(41,41) - rrt(306) * density(43) 
  pd(41,43) = pd(41,43) - rrt(306) * density(41) 
  pd(43,41) = pd(43,41) - rrt(306) * density(43) 
  pd(43,43) = pd(43,43) - rrt(306) * density(41) 
  pd(30,42) = pd(30,42) + rrt(307) * density(43) 
  pd(30,43) = pd(30,43) + rrt(307) * density(42) 
  pd(37,42) = pd(37,42) + rrt(307) * density(43) 
  pd(37,43) = pd(37,43) + rrt(307) * density(42) 
  pd(42,42) = pd(42,42) - rrt(307) * density(43) 
  pd(42,43) = pd(42,43) - rrt(307) * density(42) 
  pd(43,42) = pd(43,42) - rrt(307) * density(43) 
  pd(43,43) = pd(43,43) - rrt(307) * density(42) 
  pd(30,42) = pd(30,42) + rrt(308) * density(43) * 2.0d0
  pd(30,43) = pd(30,43) + rrt(308) * density(42) * 2.0d0
  pd(36,42) = pd(36,42) + rrt(308) * density(43) 
  pd(36,43) = pd(36,43) + rrt(308) * density(42) 
  pd(42,42) = pd(42,42) - rrt(308) * density(43) 
  pd(42,43) = pd(42,43) - rrt(308) * density(42) 
  pd(43,42) = pd(43,42) - rrt(308) * density(43) 
  pd(43,43) = pd(43,43) - rrt(308) * density(42) 
  pd(01,39) = pd(01,39) + rrt(309) * density(43) 
  pd(01,43) = pd(01,43) + rrt(309) * density(39) 
  pd(30,39) = pd(30,39) + rrt(309) * density(43) 
  pd(30,43) = pd(30,43) + rrt(309) * density(39) 
  pd(39,39) = pd(39,39) - rrt(309) * density(43) 
  pd(39,43) = pd(39,43) - rrt(309) * density(39) 
  pd(43,39) = pd(43,39) - rrt(309) * density(43) 
  pd(43,43) = pd(43,43) - rrt(309) * density(39) 
  pd(18,18) = pd(18,18) - rrt(310) * density(21) 
  pd(18,21) = pd(18,21) - rrt(310) * density(18) 
  pd(21,18) = pd(21,18) - rrt(310) * density(21) 
  pd(21,21) = pd(21,21) - rrt(310) * density(18) 
  pd(30,18) = pd(30,18) + rrt(310) * density(21) 
  pd(30,21) = pd(30,21) + rrt(310) * density(18) 
  pd(39,18) = pd(39,18) + rrt(310) * density(21) 
  pd(39,21) = pd(39,21) + rrt(310) * density(18) 
  pd(10,10) = pd(10,10) - rrt(311) * density(18) 
  pd(10,18) = pd(10,18) - rrt(311) * density(10) 
  pd(14,10) = pd(14,10) + rrt(311) * density(18) 
  pd(14,18) = pd(14,18) + rrt(311) * density(10) 
  pd(18,10) = pd(18,10) - rrt(311) * density(18) 
  pd(18,18) = pd(18,18) - rrt(311) * density(10) 
  pd(19,10) = pd(19,10) + rrt(311) * density(18) 
  pd(19,18) = pd(19,18) + rrt(311) * density(10) 
  pd(01,14) = pd(01,14) + rrt(312) * density(18) 
  pd(01,18) = pd(01,18) + rrt(312) * density(14) 
  pd(14,14) = pd(14,14) - rrt(312) * density(18) 
  pd(14,18) = pd(14,18) - rrt(312) * density(14) 
  pd(17,14) = pd(17,14) + rrt(312) * density(18) 
  pd(17,18) = pd(17,18) + rrt(312) * density(14) 
  pd(18,14) = pd(18,14) - rrt(312) * density(18) 
  pd(18,18) = pd(18,18) - rrt(312) * density(14) 
  pd(14,01) = pd(14,01) - rrt(313) * density(14) * density(18) 
  pd(14,14) = pd(14,14) - rrt(313) * density(01) * density(18) 
  pd(14,18) = pd(14,18) - rrt(313) * density(01) * density(14) 
  pd(18,01) = pd(18,01) - rrt(313) * density(14) * density(18) 
  pd(18,14) = pd(18,14) - rrt(313) * density(01) * density(18) 
  pd(18,18) = pd(18,18) - rrt(313) * density(01) * density(14) 
  pd(19,01) = pd(19,01) + rrt(313) * density(14) * density(18) 
  pd(19,14) = pd(19,14) + rrt(313) * density(01) * density(18) 
  pd(19,18) = pd(19,18) + rrt(313) * density(01) * density(14) 
  pd(01,01) = pd(01,01) - rrt(314) * density(01) * density(18) * 2.0d0
  pd(01,18) = pd(01,18) - rrt(314) * density(01)**2 
  pd(18,01) = pd(18,01) - rrt(314) * density(01) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(314) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(314) * density(01) * density(18) * 2.0d0
  pd(20,18) = pd(20,18) + rrt(314) * density(01)**2 
  pd(01,18) = pd(01,18) + rrt(315) * density(37) 
  pd(01,37) = pd(01,37) + rrt(315) * density(18) 
  pd(18,18) = pd(18,18) - rrt(315) * density(37) 
  pd(18,37) = pd(18,37) - rrt(315) * density(18) 
  pd(37,18) = pd(37,18) - rrt(315) * density(37) 
  pd(37,37) = pd(37,37) - rrt(315) * density(18) 
  pd(41,18) = pd(41,18) + rrt(315) * density(37) 
  pd(41,37) = pd(41,37) + rrt(315) * density(18) 
  pd(01,14) = pd(01,14) + rrt(316) * density(19) 
  pd(01,19) = pd(01,19) + rrt(316) * density(14) 
  pd(14,14) = pd(14,14) - rrt(316) * density(19) 
  pd(14,19) = pd(14,19) - rrt(316) * density(14) 
  pd(18,14) = pd(18,14) + rrt(316) * density(19) 
  pd(18,19) = pd(18,19) + rrt(316) * density(14) 
  pd(19,14) = pd(19,14) - rrt(316) * density(19) 
  pd(19,19) = pd(19,19) - rrt(316) * density(14) 
  pd(01,01) = pd(01,01) + rrt(317) * density(20) 
  pd(01,20) = pd(01,20) + rrt(317) * density(01) 
  pd(18,01) = pd(18,01) + rrt(317) * density(20) 
  pd(18,20) = pd(18,20) + rrt(317) * density(01) 
  pd(20,01) = pd(20,01) - rrt(317) * density(20) 
  pd(20,20) = pd(20,20) - rrt(317) * density(01) 
  pd(01,14) = pd(01,14) + rrt(318) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(318) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(318) * density(20) 
  pd(14,20) = pd(14,20) - rrt(318) * density(14) 
  pd(17,14) = pd(17,14) + rrt(318) * density(20) 
  pd(17,20) = pd(17,20) + rrt(318) * density(14) 
  pd(20,14) = pd(20,14) - rrt(318) * density(20) 
  pd(20,20) = pd(20,20) - rrt(318) * density(14) 
  pd(17,17) = pd(17,17) - rrt(319) * density(21) 
  pd(17,21) = pd(17,21) - rrt(319) * density(17) 
  pd(21,17) = pd(21,17) - rrt(319) * density(21) 
  pd(21,21) = pd(21,21) - rrt(319) * density(17) 
  pd(30,17) = pd(30,17) + rrt(319) * density(21) 
  pd(30,21) = pd(30,21) + rrt(319) * density(17) 
  pd(38,17) = pd(38,17) + rrt(319) * density(21) 
  pd(38,21) = pd(38,21) + rrt(319) * density(17) 
  pd(17,17) = pd(17,17) - rrt(320) * density(37) 
  pd(17,37) = pd(17,37) - rrt(320) * density(17) 
  pd(35,17) = pd(35,17) + rrt(320) * density(37) 
  pd(35,37) = pd(35,37) + rrt(320) * density(17) 
  pd(37,17) = pd(37,17) - rrt(320) * density(37) 
  pd(37,37) = pd(37,37) - rrt(320) * density(17) 
  pd(40,17) = pd(40,17) + rrt(320) * density(37) 
  pd(40,37) = pd(40,37) + rrt(320) * density(17) 
  pd(14,17) = pd(14,17) + rrt(321) * density(37) 
  pd(14,37) = pd(14,37) + rrt(321) * density(17) 
  pd(17,17) = pd(17,17) - rrt(321) * density(37) 
  pd(17,37) = pd(17,37) - rrt(321) * density(17) 
  pd(37,17) = pd(37,17) - rrt(321) * density(37) 
  pd(37,37) = pd(37,37) - rrt(321) * density(17) 
  pd(41,17) = pd(41,17) + rrt(321) * density(37) 
  pd(41,37) = pd(41,37) + rrt(321) * density(17) 
  pd(17,17) = pd(17,17) - rrt(322) * density(37) 
  pd(17,37) = pd(17,37) - rrt(322) * density(17) 
  pd(21,17) = pd(21,17) + rrt(322) * density(37) 
  pd(21,37) = pd(21,37) + rrt(322) * density(17) 
  pd(37,17) = pd(37,17) - rrt(322) * density(37) 
  pd(37,37) = pd(37,37) - rrt(322) * density(17) 
  pd(39,17) = pd(39,17) + rrt(322) * density(37) 
  pd(39,37) = pd(39,37) + rrt(322) * density(17) 
  pd(21,30) = pd(21,30) + rrt(323) * density(31) 
  pd(21,31) = pd(21,31) + rrt(323) * density(30) 
  pd(30,30) = pd(30,30) - rrt(323) * density(31) 
  pd(30,31) = pd(30,31) - rrt(323) * density(30) 
  pd(31,30) = pd(31,30) - rrt(323) * density(31) 
  pd(31,31) = pd(31,31) - rrt(323) * density(30) 
  pd(33,30) = pd(33,30) + rrt(323) * density(31) 
  pd(33,31) = pd(33,31) + rrt(323) * density(30) 
  pd(21,21) = pd(21,21) - rrt(324) * density(31) 
  pd(21,31) = pd(21,31) - rrt(324) * density(21) 
  pd(30,21) = pd(30,21) + rrt(324) * density(31) 
  pd(30,31) = pd(30,31) + rrt(324) * density(21) 
  pd(31,21) = pd(31,21) - rrt(324) * density(31) 
  pd(31,31) = pd(31,31) - rrt(324) * density(21) 
  pd(32,21) = pd(32,21) + rrt(324) * density(31) 
  pd(32,31) = pd(32,31) + rrt(324) * density(21) 
  pd(21,31) = pd(21,31) + rrt(325) * density(37) 
  pd(21,37) = pd(21,37) + rrt(325) * density(31) 
  pd(31,31) = pd(31,31) - rrt(325) * density(37) 
  pd(31,37) = pd(31,37) - rrt(325) * density(31) 
  pd(37,31) = pd(37,31) - rrt(325) * density(37) 
  pd(37,37) = pd(37,37) - rrt(325) * density(31) 
  pd(41,31) = pd(41,31) + rrt(325) * density(37) 
  pd(41,37) = pd(41,37) + rrt(325) * density(31) 
  pd(01,01) = pd(01,01) - rrt(326) * density(31) 
  pd(01,31) = pd(01,31) - rrt(326) * density(01) 
  pd(30,01) = pd(30,01) + rrt(326) * density(31) 
  pd(30,31) = pd(30,31) + rrt(326) * density(01) 
  pd(31,01) = pd(31,01) - rrt(326) * density(31) 
  pd(31,31) = pd(31,31) - rrt(326) * density(01) 
  pd(39,01) = pd(39,01) + rrt(326) * density(31) 
  pd(39,31) = pd(39,31) + rrt(326) * density(01) 
  pd(30,33) = pd(30,33) + rrt(327) * density(37) 
  pd(30,37) = pd(30,37) + rrt(327) * density(33) 
  pd(33,33) = pd(33,33) - rrt(327) * density(37) 
  pd(33,37) = pd(33,37) - rrt(327) * density(33) 
  pd(37,33) = pd(37,33) - rrt(327) * density(37) 
  pd(37,37) = pd(37,37) - rrt(327) * density(33) 
  pd(41,33) = pd(41,33) + rrt(327) * density(37) 
  pd(41,37) = pd(41,37) + rrt(327) * density(33) 
  pd(14,21) = pd(14,21) + rrt(328) * density(38) 
  pd(14,38) = pd(14,38) + rrt(328) * density(21) 
  pd(21,21) = pd(21,21) - rrt(328) * density(38) 
  pd(21,38) = pd(21,38) - rrt(328) * density(21) 
  pd(32,21) = pd(32,21) + rrt(328) * density(38) 
  pd(32,38) = pd(32,38) + rrt(328) * density(21) 
  pd(38,21) = pd(38,21) - rrt(328) * density(38) 
  pd(38,38) = pd(38,38) - rrt(328) * density(21) 
  pd(21,21) = pd(21,21) - rrt(329) * density(38) 
  pd(21,38) = pd(21,38) - rrt(329) * density(21) 
  pd(30,21) = pd(30,21) + rrt(329) * density(38) 
  pd(30,38) = pd(30,38) + rrt(329) * density(21) 
  pd(38,21) = pd(38,21) - rrt(329) * density(38) 
  pd(38,38) = pd(38,38) - rrt(329) * density(21) 
  pd(40,21) = pd(40,21) + rrt(329) * density(38) 
  pd(40,38) = pd(40,38) + rrt(329) * density(21) 
  pd(35,37) = pd(35,37) + rrt(330) * density(38) 
  pd(35,38) = pd(35,38) + rrt(330) * density(37) 
  pd(37,37) = pd(37,37) - rrt(330) * density(38) 
  pd(37,38) = pd(37,38) - rrt(330) * density(37) 
  pd(38,37) = pd(38,37) - rrt(330) * density(38) 
  pd(38,38) = pd(38,38) - rrt(330) * density(37) 
  pd(41,37) = pd(41,37) + rrt(330) * density(38) 
  pd(41,38) = pd(41,38) + rrt(330) * density(37) 
  pd(14,37) = pd(14,37) + rrt(331) * density(38) 
  pd(14,38) = pd(14,38) + rrt(331) * density(37) 
  pd(37,37) = pd(37,37) - rrt(331) * density(38) 
  pd(37,38) = pd(37,38) - rrt(331) * density(37) 
  pd(38,37) = pd(38,37) - rrt(331) * density(38) 
  pd(38,38) = pd(38,38) - rrt(331) * density(37) 
  pd(42,37) = pd(42,37) + rrt(331) * density(38) 
  pd(42,38) = pd(42,38) + rrt(331) * density(37) 
  pd(01,01) = pd(01,01) - rrt(332) * density(38) 
  pd(01,38) = pd(01,38) - rrt(332) * density(01) 
  pd(14,01) = pd(14,01) + rrt(332) * density(38) 
  pd(14,38) = pd(14,38) + rrt(332) * density(01) 
  pd(38,01) = pd(38,01) - rrt(332) * density(38) 
  pd(38,38) = pd(38,38) - rrt(332) * density(01) 
  pd(39,01) = pd(39,01) + rrt(332) * density(38) 
  pd(39,38) = pd(39,38) + rrt(332) * density(01) 
  pd(21,21) = pd(21,21) - rrt(333) * density(40) 
  pd(21,40) = pd(21,40) - rrt(333) * density(21) 
  pd(30,21) = pd(30,21) + rrt(333) * density(40) 
  pd(30,40) = pd(30,40) + rrt(333) * density(21) 
  pd(40,21) = pd(40,21) - rrt(333) * density(40) 
  pd(40,40) = pd(40,40) - rrt(333) * density(21) 
  pd(41,21) = pd(41,21) + rrt(333) * density(40) 
  pd(41,40) = pd(41,40) + rrt(333) * density(21) 
  pd(36,37) = pd(36,37) + rrt(334) * density(40) 
  pd(36,40) = pd(36,40) + rrt(334) * density(37) 
  pd(37,37) = pd(37,37) - rrt(334) * density(40) 
  pd(37,40) = pd(37,40) - rrt(334) * density(37) 
  pd(40,37) = pd(40,37) - rrt(334) * density(40) 
  pd(40,40) = pd(40,40) - rrt(334) * density(37) 
  pd(41,37) = pd(41,37) + rrt(334) * density(40) 
  pd(41,40) = pd(41,40) + rrt(334) * density(37) 
  pd(35,37) = pd(35,37) + rrt(335) * density(40) 
  pd(35,40) = pd(35,40) + rrt(335) * density(37) 
  pd(37,37) = pd(37,37) - rrt(335) * density(40) 
  pd(37,40) = pd(37,40) - rrt(335) * density(37) 
  pd(40,37) = pd(40,37) - rrt(335) * density(40) 
  pd(40,40) = pd(40,40) - rrt(335) * density(37) 
  pd(42,37) = pd(42,37) + rrt(335) * density(40) 
  pd(42,40) = pd(42,40) + rrt(335) * density(37) 
  pd(36,37) = pd(36,37) + rrt(336) * density(41) 
  pd(36,41) = pd(36,41) + rrt(336) * density(37) 
  pd(37,37) = pd(37,37) - rrt(336) * density(41) 
  pd(37,41) = pd(37,41) - rrt(336) * density(37) 
  pd(41,37) = pd(41,37) - rrt(336) * density(41) 
  pd(41,41) = pd(41,41) - rrt(336) * density(37) 
  pd(42,37) = pd(42,37) + rrt(336) * density(41) 
  pd(42,41) = pd(42,41) + rrt(336) * density(37) 
  pd(01,10) = pd(01,10) + rrt(337) * density(14) 
  pd(01,14) = pd(01,14) + rrt(337) * density(10) 
  pd(10,10) = pd(10,10) - rrt(337) * density(14) 
  pd(10,14) = pd(10,14) - rrt(337) * density(10) 
  pd(01,10) = pd(01,10) + rrt(338) * density(14) 
  pd(01,14) = pd(01,14) + rrt(338) * density(10) 
  pd(10,10) = pd(10,10) - rrt(338) * density(14) 
  pd(10,14) = pd(10,14) - rrt(338) * density(10) 
  pd(14,10) = pd(14,10) - rrt(338) * density(14) 
  pd(14,14) = pd(14,14) - rrt(338) * density(10) 
  pd(16,10) = pd(16,10) + rrt(338) * density(14) 
  pd(16,14) = pd(16,14) + rrt(338) * density(10) 
  pd(01,01) = pd(01,01) + rrt(339) * density(10) 
  pd(01,10) = pd(01,10) + rrt(339) * density(01) 
  pd(10,01) = pd(10,01) - rrt(339) * density(10) 
  pd(10,10) = pd(10,10) - rrt(339) * density(01) 
  pd(01,10) = pd(01,10) + rrt(340) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(340) * density(10) * 4.0d0
  pd(11,10) = pd(11,10) + rrt(340) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(341) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(341) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(341) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(342) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(342) * density(10) * 4.0d0
  pd(14,10) = pd(14,10) + rrt(342) * density(10) * 4.0d0
  pd(10,01) = pd(10,01) + rrt(343) * density(11) 
  pd(10,11) = pd(10,11) + rrt(343) * density(01) 
  pd(11,01) = pd(11,01) - rrt(343) * density(11) 
  pd(11,11) = pd(11,11) - rrt(343) * density(01) 
  pd(01,01) = pd(01,01) + rrt(344) * density(11) 
  pd(01,11) = pd(01,11) + rrt(344) * density(01) 
  pd(11,01) = pd(11,01) - rrt(344) * density(11) 
  pd(11,11) = pd(11,11) - rrt(344) * density(01) 
  pd(11,01) = pd(11,01) + rrt(345) * density(12) 
  pd(11,12) = pd(11,12) + rrt(345) * density(01) 
  pd(12,01) = pd(12,01) - rrt(345) * density(12) 
  pd(12,12) = pd(12,12) - rrt(345) * density(01) 
  pd(01,12) = pd(01,12) + rrt(346) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(346) * density(12) * 4.0d0
  pd(18,12) = pd(18,12) + rrt(346) * density(12) * 2.0d0
  pd(43,12) = pd(43,12) + rrt(346) * density(12) * 2.0d0
  pd(12,12) = pd(12,12) - rrt(347) * density(12) * 4.0d0
  pd(20,12) = pd(20,12) + rrt(347) * density(12) * 2.0d0
  pd(43,12) = pd(43,12) + rrt(347) * density(12) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(348) * density(12) 
  pd(10,12) = pd(10,12) - rrt(348) * density(10) 
  pd(12,10) = pd(12,10) - rrt(348) * density(12) 
  pd(12,12) = pd(12,12) - rrt(348) * density(10) 
  pd(20,10) = pd(20,10) + rrt(348) * density(12) 
  pd(20,12) = pd(20,12) + rrt(348) * density(10) 
  pd(43,10) = pd(43,10) + rrt(348) * density(12) 
  pd(43,12) = pd(43,12) + rrt(348) * density(10) 
  pd(14,01) = pd(14,01) + rrt(349) * density(15) 
  pd(14,15) = pd(14,15) + rrt(349) * density(01) 
  pd(15,01) = pd(15,01) - rrt(349) * density(15) 
  pd(15,15) = pd(15,15) - rrt(349) * density(01) 
  pd(14,14) = pd(14,14) + rrt(350) * density(16) 
  pd(14,16) = pd(14,16) + rrt(350) * density(14) 
  pd(16,14) = pd(16,14) - rrt(350) * density(16) 
  pd(16,16) = pd(16,16) - rrt(350) * density(14) 
  pd(15,14) = pd(15,14) + rrt(351) * density(16) 
  pd(15,16) = pd(15,16) + rrt(351) * density(14) 
  pd(16,14) = pd(16,14) - rrt(351) * density(16) 
  pd(16,16) = pd(16,16) - rrt(351) * density(14) 
  pd(14,01) = pd(14,01) + rrt(352) * density(16) 
  pd(14,16) = pd(14,16) + rrt(352) * density(01) 
  pd(16,01) = pd(16,01) - rrt(352) * density(16) 
  pd(16,16) = pd(16,16) - rrt(352) * density(01) 
  pd(01,10) = pd(01,10) + rrt(353) * density(30) 
  pd(01,30) = pd(01,30) + rrt(353) * density(10) 
  pd(10,10) = pd(10,10) - rrt(353) * density(30) 
  pd(10,30) = pd(10,30) - rrt(353) * density(10) 
  pd(01,10) = pd(01,10) + rrt(354) * density(21) 
  pd(01,21) = pd(01,21) + rrt(354) * density(10) 
  pd(10,10) = pd(10,10) - rrt(354) * density(21) 
  pd(10,21) = pd(10,21) - rrt(354) * density(10) 
  pd(21,10) = pd(21,10) - rrt(354) * density(21) 
  pd(21,21) = pd(21,21) - rrt(354) * density(10) 
  pd(30,10) = pd(30,10) + rrt(354) * density(21) * 2.0d0
  pd(30,21) = pd(30,21) + rrt(354) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(355) * density(37) 
  pd(01,37) = pd(01,37) + rrt(355) * density(10) 
  pd(10,10) = pd(10,10) - rrt(355) * density(37) 
  pd(10,37) = pd(10,37) - rrt(355) * density(10) 
  pd(10,11) = pd(10,11) + rrt(356) * density(21) 
  pd(10,21) = pd(10,21) + rrt(356) * density(11) 
  pd(11,11) = pd(11,11) - rrt(356) * density(21) 
  pd(11,21) = pd(11,21) - rrt(356) * density(11) 
  pd(01,12) = pd(01,12) + rrt(357) * density(30) 
  pd(01,30) = pd(01,30) + rrt(357) * density(12) 
  pd(12,12) = pd(12,12) - rrt(357) * density(30) 
  pd(12,30) = pd(12,30) - rrt(357) * density(12) 
  pd(01,12) = pd(01,12) + rrt(358) * density(21) 
  pd(01,21) = pd(01,21) + rrt(358) * density(12) 
  pd(12,12) = pd(12,12) - rrt(358) * density(21) 
  pd(12,21) = pd(12,21) - rrt(358) * density(12) 
  pd(21,12) = pd(21,12) - rrt(358) * density(21) 
  pd(21,21) = pd(21,21) - rrt(358) * density(12) 
  pd(30,12) = pd(30,12) + rrt(358) * density(21) * 2.0d0
  pd(30,21) = pd(30,21) + rrt(358) * density(12) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(359) * density(27) 
  pd(14,27) = pd(14,27) - rrt(359) * density(14) 
  pd(27,14) = pd(27,14) - rrt(359) * density(27) 
  pd(27,27) = pd(27,27) - rrt(359) * density(14) 
  pd(30,14) = pd(30,14) + rrt(359) * density(27) 
  pd(30,27) = pd(30,27) + rrt(359) * density(14) 
  pd(35,14) = pd(35,14) + rrt(359) * density(27) 
  pd(35,27) = pd(35,27) + rrt(359) * density(14) 
  pd(14,14) = pd(14,14) - rrt(360) * density(28) 
  pd(14,28) = pd(14,28) - rrt(360) * density(14) 
  pd(28,14) = pd(28,14) - rrt(360) * density(28) 
  pd(28,28) = pd(28,28) - rrt(360) * density(14) 
  pd(30,14) = pd(30,14) + rrt(360) * density(28) 
  pd(30,28) = pd(30,28) + rrt(360) * density(14) 
  pd(35,14) = pd(35,14) + rrt(360) * density(28) 
  pd(35,28) = pd(35,28) + rrt(360) * density(14) 
  pd(14,14) = pd(14,14) - rrt(361) * density(29) 
  pd(14,29) = pd(14,29) - rrt(361) * density(14) 
  pd(29,14) = pd(29,14) - rrt(361) * density(29) 
  pd(29,29) = pd(29,29) - rrt(361) * density(14) 
  pd(30,14) = pd(30,14) + rrt(361) * density(29) 
  pd(30,29) = pd(30,29) + rrt(361) * density(14) 
  pd(35,14) = pd(35,14) + rrt(361) * density(29) 
  pd(35,29) = pd(35,29) + rrt(361) * density(14) 
  pd(14,14) = pd(14,14) - rrt(362) * density(26) 
  pd(14,26) = pd(14,26) - rrt(362) * density(14) 
  pd(26,14) = pd(26,14) - rrt(362) * density(26) 
  pd(26,26) = pd(26,26) - rrt(362) * density(14) 
  pd(30,14) = pd(30,14) + rrt(362) * density(26) 
  pd(30,26) = pd(30,26) + rrt(362) * density(14) 
  pd(35,14) = pd(35,14) + rrt(362) * density(26) 
  pd(35,26) = pd(35,26) + rrt(362) * density(14) 
  pd(14,14) = pd(14,14) - rrt(363) * density(22) 
  pd(14,22) = pd(14,22) - rrt(363) * density(14) 
  pd(22,14) = pd(22,14) - rrt(363) * density(22) 
  pd(22,22) = pd(22,22) - rrt(363) * density(14) 
  pd(30,14) = pd(30,14) + rrt(363) * density(22) 
  pd(30,22) = pd(30,22) + rrt(363) * density(14) 
  pd(35,14) = pd(35,14) + rrt(363) * density(22) 
  pd(35,22) = pd(35,22) + rrt(363) * density(14) 
  pd(14,14) = pd(14,14) - rrt(364) * density(23) 
  pd(14,23) = pd(14,23) - rrt(364) * density(14) 
  pd(23,14) = pd(23,14) - rrt(364) * density(23) 
  pd(23,23) = pd(23,23) - rrt(364) * density(14) 
  pd(30,14) = pd(30,14) + rrt(364) * density(23) 
  pd(30,23) = pd(30,23) + rrt(364) * density(14) 
  pd(35,14) = pd(35,14) + rrt(364) * density(23) 
  pd(35,23) = pd(35,23) + rrt(364) * density(14) 
  pd(14,14) = pd(14,14) - rrt(365) * density(24) 
  pd(14,24) = pd(14,24) - rrt(365) * density(14) 
  pd(24,14) = pd(24,14) - rrt(365) * density(24) 
  pd(24,24) = pd(24,24) - rrt(365) * density(14) 
  pd(30,14) = pd(30,14) + rrt(365) * density(24) 
  pd(30,24) = pd(30,24) + rrt(365) * density(14) 
  pd(35,14) = pd(35,14) + rrt(365) * density(24) 
  pd(35,24) = pd(35,24) + rrt(365) * density(14) 
  pd(14,14) = pd(14,14) - rrt(366) * density(25) 
  pd(14,25) = pd(14,25) - rrt(366) * density(14) 
  pd(25,14) = pd(25,14) - rrt(366) * density(25) 
  pd(25,25) = pd(25,25) - rrt(366) * density(14) 
  pd(30,14) = pd(30,14) + rrt(366) * density(25) 
  pd(30,25) = pd(30,25) + rrt(366) * density(14) 
  pd(35,14) = pd(35,14) + rrt(366) * density(25) 
  pd(35,25) = pd(35,25) + rrt(366) * density(14) 
  pd(15,15) = pd(15,15) - rrt(367) * density(21) 
  pd(15,21) = pd(15,21) - rrt(367) * density(15) 
  pd(21,15) = pd(21,15) - rrt(367) * density(21) 
  pd(21,21) = pd(21,21) - rrt(367) * density(15) 
  pd(30,15) = pd(30,15) + rrt(367) * density(21) 
  pd(30,21) = pd(30,21) + rrt(367) * density(15) 
  pd(35,15) = pd(35,15) + rrt(367) * density(21) 
  pd(35,21) = pd(35,21) + rrt(367) * density(15) 
  pd(15,15) = pd(15,15) - rrt(368) * density(37) 
  pd(15,37) = pd(15,37) - rrt(368) * density(15) 
  pd(35,15) = pd(35,15) + rrt(368) * density(37) 
  pd(35,37) = pd(35,37) + rrt(368) * density(15) 
  pd(36,15) = pd(36,15) + rrt(368) * density(37) 
  pd(36,37) = pd(36,37) + rrt(368) * density(15) 
  pd(37,15) = pd(37,15) - rrt(368) * density(37) 
  pd(37,37) = pd(37,37) - rrt(368) * density(15) 
  pd(16,16) = pd(16,16) - rrt(369) * density(21) 
  pd(16,21) = pd(16,21) - rrt(369) * density(16) 
  pd(21,16) = pd(21,16) - rrt(369) * density(21) 
  pd(21,21) = pd(21,21) - rrt(369) * density(16) 
  pd(30,16) = pd(30,16) + rrt(369) * density(21) 
  pd(30,21) = pd(30,21) + rrt(369) * density(16) 
  pd(35,16) = pd(35,16) + rrt(369) * density(21) 
  pd(35,21) = pd(35,21) + rrt(369) * density(16) 
  pd(01,14) = pd(01,14) + rrt(370) * density(35) 
  pd(01,35) = pd(01,35) + rrt(370) * density(14) 
  pd(14,14) = pd(14,14) - rrt(370) * density(35) 
  pd(14,35) = pd(14,35) - rrt(370) * density(14) 
  pd(30,14) = pd(30,14) + rrt(370) * density(35) 
  pd(30,35) = pd(30,35) + rrt(370) * density(14) 
  pd(35,14) = pd(35,14) - rrt(370) * density(35) 
  pd(35,35) = pd(35,35) - rrt(370) * density(14) 
  pd(14,30) = pd(14,30) + rrt(371) * density(35) 
  pd(14,35) = pd(14,35) + rrt(371) * density(30) 
  pd(21,30) = pd(21,30) + rrt(371) * density(35) 
  pd(21,35) = pd(21,35) + rrt(371) * density(30) 
  pd(30,30) = pd(30,30) - rrt(371) * density(35) 
  pd(30,35) = pd(30,35) - rrt(371) * density(30) 
  pd(35,30) = pd(35,30) - rrt(371) * density(35) 
  pd(35,35) = pd(35,35) - rrt(371) * density(30) 
  pd(01,35) = pd(01,35) + rrt(372) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(372) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(372) * density(35) * 4.0d0
  pd(14,35) = pd(14,35) + rrt(373) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(373) * density(35) * 4.0d0
  pd(36,35) = pd(36,35) + rrt(373) * density(35) * 2.0d0
  pd(01,35) = pd(01,35) + rrt(374) * density(35) * 2.0d0
  pd(30,35) = pd(30,35) + rrt(374) * density(35) * 4.0d0
  pd(35,35) = pd(35,35) - rrt(374) * density(35) * 4.0d0
  pd(21,30) = pd(21,30) + rrt(375) * density(36) 
  pd(21,36) = pd(21,36) + rrt(375) * density(30) 
  pd(30,30) = pd(30,30) - rrt(375) * density(36) 
  pd(30,36) = pd(30,36) - rrt(375) * density(30) 
  pd(35,30) = pd(35,30) + rrt(375) * density(36) 
  pd(35,36) = pd(35,36) + rrt(375) * density(30) 
  pd(36,30) = pd(36,30) - rrt(375) * density(36) 
  pd(36,36) = pd(36,36) - rrt(375) * density(30) 
  pd(01,14) = pd(01,14) + rrt(376) * density(36) 
  pd(01,36) = pd(01,36) + rrt(376) * density(14) 
  pd(14,14) = pd(14,14) - rrt(376) * density(36) 
  pd(14,36) = pd(14,36) - rrt(376) * density(14) 
  pd(30,14) = pd(30,14) + rrt(376) * density(36) * 2.0d0
  pd(30,36) = pd(30,36) + rrt(376) * density(14) * 2.0d0
  pd(36,14) = pd(36,14) - rrt(376) * density(36) 
  pd(36,36) = pd(36,36) - rrt(376) * density(14) 
  pd(01,14) = pd(01,14) + rrt(377) * density(36) 
  pd(01,36) = pd(01,36) + rrt(377) * density(14) 
  pd(14,14) = pd(14,14) - rrt(377) * density(36) 
  pd(14,36) = pd(14,36) - rrt(377) * density(14) 
  pd(21,14) = pd(21,14) + rrt(377) * density(36) 
  pd(21,36) = pd(21,36) + rrt(377) * density(14) 
  pd(36,14) = pd(36,14) - rrt(377) * density(36) 
  pd(36,36) = pd(36,36) - rrt(377) * density(14) 
  pd(14,35) = pd(14,35) + rrt(378) * density(36) 
  pd(14,36) = pd(14,36) + rrt(378) * density(35) 
  pd(35,35) = pd(35,35) - rrt(378) * density(36) 
  pd(35,36) = pd(35,36) - rrt(378) * density(35) 
  pd(36,35) = pd(36,35) - rrt(378) * density(36) 
  pd(36,36) = pd(36,36) - rrt(378) * density(35) 
  pd(37,35) = pd(37,35) + rrt(378) * density(36) 
  pd(37,36) = pd(37,36) + rrt(378) * density(35) 
  pd(21,21) = pd(21,21) - rrt(379) * density(36) 
  pd(21,36) = pd(21,36) - rrt(379) * density(21) 
  pd(30,21) = pd(30,21) + rrt(379) * density(36) 
  pd(30,36) = pd(30,36) + rrt(379) * density(21) 
  pd(36,21) = pd(36,21) - rrt(379) * density(36) 
  pd(36,36) = pd(36,36) - rrt(379) * density(21) 
  pd(37,21) = pd(37,21) + rrt(379) * density(36) 
  pd(37,36) = pd(37,36) + rrt(379) * density(21) 
  pd(21,30) = pd(21,30) + rrt(380) * density(37) 
  pd(21,37) = pd(21,37) + rrt(380) * density(30) 
  pd(30,30) = pd(30,30) - rrt(380) * density(37) 
  pd(30,37) = pd(30,37) - rrt(380) * density(30) 
  pd(36,30) = pd(36,30) + rrt(380) * density(37) 
  pd(36,37) = pd(36,37) + rrt(380) * density(30) 
  pd(37,30) = pd(37,30) - rrt(380) * density(37) 
  pd(37,37) = pd(37,37) - rrt(380) * density(30) 
  pd(10,01) = pd(10,01) + rrt(381) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(381) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(381) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(381) * density(01) * density(14) * 4.0d0
  pd(10,14) = pd(10,14) + rrt(382) * density(14) * density(21) * 2.0d0
  pd(10,21) = pd(10,21) + rrt(382) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(382) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(382) * density(14)**2 * 2.0d0
  pd(10,02) = pd(10,02) + rrt(383) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(383) * density(02) * density(14) * 2.0d0
  pd(14,02) = pd(14,02) - rrt(383) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(383) * density(02) * density(14) * 4.0d0
  pd(10,03) = pd(10,03) + rrt(384) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(384) * density(03) * density(14) * 2.0d0
  pd(14,03) = pd(14,03) - rrt(384) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(384) * density(03) * density(14) * 4.0d0
  pd(10,04) = pd(10,04) + rrt(385) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(385) * density(04) * density(14) * 2.0d0
  pd(14,04) = pd(14,04) - rrt(385) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(385) * density(04) * density(14) * 4.0d0
  pd(10,05) = pd(10,05) + rrt(386) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(386) * density(05) * density(14) * 2.0d0
  pd(14,05) = pd(14,05) - rrt(386) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(386) * density(05) * density(14) * 4.0d0
  pd(10,06) = pd(10,06) + rrt(387) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(387) * density(06) * density(14) * 2.0d0
  pd(14,06) = pd(14,06) - rrt(387) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(387) * density(06) * density(14) * 4.0d0
  pd(10,07) = pd(10,07) + rrt(388) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(388) * density(07) * density(14) * 2.0d0
  pd(14,07) = pd(14,07) - rrt(388) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(388) * density(07) * density(14) * 4.0d0
  pd(10,08) = pd(10,08) + rrt(389) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(389) * density(08) * density(14) * 2.0d0
  pd(14,08) = pd(14,08) - rrt(389) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(389) * density(08) * density(14) * 4.0d0
  pd(10,09) = pd(10,09) + rrt(390) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(390) * density(09) * density(14) * 2.0d0
  pd(14,09) = pd(14,09) - rrt(390) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(390) * density(09) * density(14) * 4.0d0
  pd(10,14) = pd(10,14) + rrt(391) * density(14) * density(27) * 2.0d0
  pd(10,27) = pd(10,27) + rrt(391) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(391) * density(14) * density(27) * 4.0d0
  pd(14,27) = pd(14,27) - rrt(391) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(392) * density(14) * density(28) * 2.0d0
  pd(10,28) = pd(10,28) + rrt(392) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(392) * density(14) * density(28) * 4.0d0
  pd(14,28) = pd(14,28) - rrt(392) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(393) * density(14) * density(29) * 2.0d0
  pd(10,29) = pd(10,29) + rrt(393) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(393) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(393) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(394) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(394) * density(14)**2 * 6.0d0
  pd(10,14) = pd(10,14) + rrt(395) * density(14) * density(30) * 2.0d0
  pd(10,30) = pd(10,30) + rrt(395) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(395) * density(14) * density(30) * 4.0d0
  pd(14,30) = pd(14,30) - rrt(395) * density(14)**2 * 2.0d0
  pd(11,01) = pd(11,01) + rrt(396) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(396) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(396) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(396) * density(01) * density(14) * 4.0d0
  pd(11,14) = pd(11,14) + rrt(397) * density(14) * density(21) * 2.0d0
  pd(11,21) = pd(11,21) + rrt(397) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(397) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(397) * density(14)**2 * 2.0d0
  pd(11,02) = pd(11,02) + rrt(398) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(398) * density(02) * density(14) * 2.0d0
  pd(14,02) = pd(14,02) - rrt(398) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(398) * density(02) * density(14) * 4.0d0
  pd(11,03) = pd(11,03) + rrt(399) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(399) * density(03) * density(14) * 2.0d0
  pd(14,03) = pd(14,03) - rrt(399) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(399) * density(03) * density(14) * 4.0d0
  pd(11,04) = pd(11,04) + rrt(400) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(400) * density(04) * density(14) * 2.0d0
  pd(14,04) = pd(14,04) - rrt(400) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(400) * density(04) * density(14) * 4.0d0
  pd(11,05) = pd(11,05) + rrt(401) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(401) * density(05) * density(14) * 2.0d0
  pd(14,05) = pd(14,05) - rrt(401) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(401) * density(05) * density(14) * 4.0d0
  pd(11,06) = pd(11,06) + rrt(402) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(402) * density(06) * density(14) * 2.0d0
  pd(14,06) = pd(14,06) - rrt(402) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(402) * density(06) * density(14) * 4.0d0
  pd(11,07) = pd(11,07) + rrt(403) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(403) * density(07) * density(14) * 2.0d0
  pd(14,07) = pd(14,07) - rrt(403) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(403) * density(07) * density(14) * 4.0d0
  pd(11,08) = pd(11,08) + rrt(404) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(404) * density(08) * density(14) * 2.0d0
  pd(14,08) = pd(14,08) - rrt(404) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(404) * density(08) * density(14) * 4.0d0
  pd(11,09) = pd(11,09) + rrt(405) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(405) * density(09) * density(14) * 2.0d0
  pd(14,09) = pd(14,09) - rrt(405) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(405) * density(09) * density(14) * 4.0d0
  pd(11,14) = pd(11,14) + rrt(406) * density(14) * density(27) * 2.0d0
  pd(11,27) = pd(11,27) + rrt(406) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(406) * density(14) * density(27) * 4.0d0
  pd(14,27) = pd(14,27) - rrt(406) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(407) * density(14) * density(28) * 2.0d0
  pd(11,28) = pd(11,28) + rrt(407) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(407) * density(14) * density(28) * 4.0d0
  pd(14,28) = pd(14,28) - rrt(407) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(408) * density(14) * density(29) * 2.0d0
  pd(11,29) = pd(11,29) + rrt(408) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(408) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(408) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(409) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(409) * density(14)**2 * 6.0d0
  pd(11,14) = pd(11,14) + rrt(410) * density(14) * density(30) * 2.0d0
  pd(11,30) = pd(11,30) + rrt(410) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(410) * density(14) * density(30) * 4.0d0
  pd(14,30) = pd(14,30) - rrt(410) * density(14)**2 * 2.0d0
  pd(01,01) = pd(01,01) + rrt(411) * density(14)**2 
  pd(01,14) = pd(01,14) + rrt(411) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(411) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(411) * density(01) * density(14) * 4.0d0
  pd(01,14) = pd(01,14) + rrt(412) * density(14) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(412) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(412) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(412) * density(14)**2 * 2.0d0
  pd(21,21) = pd(21,21) + rrt(413) * density(30)**2 
  pd(21,30) = pd(21,30) + rrt(413) * density(21) * density(30) * 2.0d0
  pd(30,21) = pd(30,21) - rrt(413) * density(30)**2 * 2.0d0
  pd(30,30) = pd(30,30) - rrt(413) * density(21) * density(30) * 4.0d0
  pd(21,01) = pd(21,01) + rrt(414) * density(30)**2 
  pd(21,30) = pd(21,30) + rrt(414) * density(01) * density(30) * 2.0d0
  pd(30,01) = pd(30,01) - rrt(414) * density(30)**2 * 2.0d0
  pd(30,30) = pd(30,30) - rrt(414) * density(01) * density(30) * 4.0d0
  pd(14,01) = pd(14,01) - rrt(415) * density(14) * density(30) 
  pd(14,14) = pd(14,14) - rrt(415) * density(01) * density(30) 
  pd(14,30) = pd(14,30) - rrt(415) * density(01) * density(14) 
  pd(30,01) = pd(30,01) - rrt(415) * density(14) * density(30) 
  pd(30,14) = pd(30,14) - rrt(415) * density(01) * density(30) 
  pd(30,30) = pd(30,30) - rrt(415) * density(01) * density(14) 
  pd(35,01) = pd(35,01) + rrt(415) * density(14) * density(30) 
  pd(35,14) = pd(35,14) + rrt(415) * density(01) * density(30) 
  pd(35,30) = pd(35,30) + rrt(415) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(416) * density(21) * density(30) 
  pd(14,21) = pd(14,21) - rrt(416) * density(14) * density(30) 
  pd(14,30) = pd(14,30) - rrt(416) * density(14) * density(21) 
  pd(30,14) = pd(30,14) - rrt(416) * density(21) * density(30) 
  pd(30,21) = pd(30,21) - rrt(416) * density(14) * density(30) 
  pd(30,30) = pd(30,30) - rrt(416) * density(14) * density(21) 
  pd(35,14) = pd(35,14) + rrt(416) * density(21) * density(30) 
  pd(35,21) = pd(35,21) + rrt(416) * density(14) * density(30) 
  pd(35,30) = pd(35,30) + rrt(416) * density(14) * density(21) 
  pd(14,01) = pd(14,01) - rrt(417) * density(14) * density(21) 
  pd(14,14) = pd(14,14) - rrt(417) * density(01) * density(21) 
  pd(14,21) = pd(14,21) - rrt(417) * density(01) * density(14) 
  pd(21,01) = pd(21,01) - rrt(417) * density(14) * density(21) 
  pd(21,14) = pd(21,14) - rrt(417) * density(01) * density(21) 
  pd(21,21) = pd(21,21) - rrt(417) * density(01) * density(14) 
  pd(36,01) = pd(36,01) + rrt(417) * density(14) * density(21) 
  pd(36,14) = pd(36,14) + rrt(417) * density(01) * density(21) 
  pd(36,21) = pd(36,21) + rrt(417) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(418) * density(21)**2 
  pd(14,21) = pd(14,21) - rrt(418) * density(14) * density(21) * 2.0d0
  pd(21,14) = pd(21,14) - rrt(418) * density(21)**2 
  pd(21,21) = pd(21,21) - rrt(418) * density(14) * density(21) * 2.0d0
  pd(36,14) = pd(36,14) + rrt(418) * density(21)**2 
  pd(36,21) = pd(36,21) + rrt(418) * density(14) * density(21) * 2.0d0
  pd(30,01) = pd(30,01) - rrt(419) * density(30) * density(35) 
  pd(30,30) = pd(30,30) - rrt(419) * density(01) * density(35) 
  pd(30,35) = pd(30,35) - rrt(419) * density(01) * density(30) 
  pd(35,01) = pd(35,01) - rrt(419) * density(30) * density(35) 
  pd(35,30) = pd(35,30) - rrt(419) * density(01) * density(35) 
  pd(35,35) = pd(35,35) - rrt(419) * density(01) * density(30) 
  pd(36,01) = pd(36,01) + rrt(419) * density(30) * density(35) 
  pd(36,30) = pd(36,30) + rrt(419) * density(01) * density(35) 
  pd(36,35) = pd(36,35) + rrt(419) * density(01) * density(30) 
  pd(30,21) = pd(30,21) - rrt(420) * density(30) * density(35) 
  pd(30,30) = pd(30,30) - rrt(420) * density(21) * density(35) 
  pd(30,35) = pd(30,35) - rrt(420) * density(21) * density(30) 
  pd(35,21) = pd(35,21) - rrt(420) * density(30) * density(35) 
  pd(35,30) = pd(35,30) - rrt(420) * density(21) * density(35) 
  pd(35,35) = pd(35,35) - rrt(420) * density(21) * density(30) 
  pd(36,21) = pd(36,21) + rrt(420) * density(30) * density(35) 
  pd(36,30) = pd(36,30) + rrt(420) * density(21) * density(35) 
  pd(36,35) = pd(36,35) + rrt(420) * density(21) * density(30) 
  pd(30,01) = pd(30,01) - rrt(421) * density(30) * density(36) 
  pd(30,30) = pd(30,30) - rrt(421) * density(01) * density(36) 
  pd(30,36) = pd(30,36) - rrt(421) * density(01) * density(30) 
  pd(36,01) = pd(36,01) - rrt(421) * density(30) * density(36) 
  pd(36,30) = pd(36,30) - rrt(421) * density(01) * density(36) 
  pd(36,36) = pd(36,36) - rrt(421) * density(01) * density(30) 
  pd(37,01) = pd(37,01) + rrt(421) * density(30) * density(36) 
  pd(37,30) = pd(37,30) + rrt(421) * density(01) * density(36) 
  pd(37,36) = pd(37,36) + rrt(421) * density(01) * density(30) 
  pd(30,21) = pd(30,21) - rrt(422) * density(30) * density(36) 
  pd(30,30) = pd(30,30) - rrt(422) * density(21) * density(36) 
  pd(30,36) = pd(30,36) - rrt(422) * density(21) * density(30) 
  pd(36,21) = pd(36,21) - rrt(422) * density(30) * density(36) 
  pd(36,30) = pd(36,30) - rrt(422) * density(21) * density(36) 
  pd(36,36) = pd(36,36) - rrt(422) * density(21) * density(30) 
  pd(37,21) = pd(37,21) + rrt(422) * density(30) * density(36) 
  pd(37,30) = pd(37,30) + rrt(422) * density(21) * density(36) 
  pd(37,36) = pd(37,36) + rrt(422) * density(21) * density(30) 
  pd(21,01) = pd(21,01) - rrt(423) * density(21) * density(35) 
  pd(21,21) = pd(21,21) - rrt(423) * density(01) * density(35) 
  pd(21,35) = pd(21,35) - rrt(423) * density(01) * density(21) 
  pd(35,01) = pd(35,01) - rrt(423) * density(21) * density(35) 
  pd(35,21) = pd(35,21) - rrt(423) * density(01) * density(35) 
  pd(35,35) = pd(35,35) - rrt(423) * density(01) * density(21) 
  pd(37,01) = pd(37,01) + rrt(423) * density(21) * density(35) 
  pd(37,21) = pd(37,21) + rrt(423) * density(01) * density(35) 
  pd(37,35) = pd(37,35) + rrt(423) * density(01) * density(21) 
  pd(21,21) = pd(21,21) - rrt(424) * density(21) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) - rrt(424) * density(21)**2 
  pd(35,21) = pd(35,21) - rrt(424) * density(21) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(424) * density(21)**2 
  pd(37,21) = pd(37,21) + rrt(424) * density(21) * density(35) * 2.0d0
  pd(37,35) = pd(37,35) + rrt(424) * density(21)**2 
  pd(01,10) = pd(01,10) + rrt(425) 
  pd(10,10) = pd(10,10) - rrt(425) 
  pd(11,12) = pd(11,12) + rrt(426) 
  pd(12,12) = pd(12,12) - rrt(426) 
  pd(01,02) = pd(01,02) + rrt(427) 
  pd(02,02) = pd(02,02) - rrt(427) 
  pd(02,03) = pd(02,03) + rrt(428) 
  pd(03,03) = pd(03,03) - rrt(428) 
  pd(03,04) = pd(03,04) + rrt(429) 
  pd(04,04) = pd(04,04) - rrt(429) 
  pd(04,05) = pd(04,05) + rrt(430) 
  pd(05,05) = pd(05,05) - rrt(430) 
  pd(05,06) = pd(05,06) + rrt(431) 
  pd(06,06) = pd(06,06) - rrt(431) 
  pd(06,07) = pd(06,07) + rrt(432) 
  pd(07,07) = pd(07,07) - rrt(432) 
  pd(07,08) = pd(07,08) + rrt(433) 
  pd(08,08) = pd(08,08) - rrt(433) 
  pd(08,09) = pd(08,09) + rrt(434) 
  pd(09,09) = pd(09,09) - rrt(434) 
  pd(21,22) = pd(21,22) + rrt(435) 
  pd(22,22) = pd(22,22) - rrt(435) 
  pd(21,23) = pd(21,23) + rrt(436) 
  pd(23,23) = pd(23,23) - rrt(436) 
  pd(21,24) = pd(21,24) + rrt(437) 
  pd(24,24) = pd(24,24) - rrt(437) 
  pd(21,25) = pd(21,25) + rrt(438) 
  pd(25,25) = pd(25,25) - rrt(438) 
  pd(21,27) = pd(21,27) + rrt(439) 
  pd(27,27) = pd(27,27) - rrt(439) 
  pd(27,28) = pd(27,28) + rrt(440) 
  pd(28,28) = pd(28,28) - rrt(440) 
  pd(28,29) = pd(28,29) + rrt(441) 
  pd(29,29) = pd(29,29) - rrt(441) 
  pd(26,26) = pd(26,26) - rrt(442) * density(43) 
  pd(26,43) = pd(26,43) - rrt(442) * density(26) 
  pd(30,26) = pd(30,26) + rrt(442) * density(43) 
  pd(30,43) = pd(30,43) + rrt(442) * density(26) 
  pd(34,26) = pd(34,26) + rrt(442) * density(43) 
  pd(34,43) = pd(34,43) + rrt(442) * density(26) 
  pd(43,26) = pd(43,26) - rrt(442) * density(43) 
  pd(43,43) = pd(43,43) - rrt(442) * density(26) 
  pd(30,31) = pd(30,31) + rrt(443) * density(34) * 3.0d0
  pd(30,34) = pd(30,34) + rrt(443) * density(31) * 3.0d0
  pd(31,31) = pd(31,31) - rrt(443) * density(34) 
  pd(31,34) = pd(31,34) - rrt(443) * density(31) 
  pd(34,31) = pd(34,31) - rrt(443) * density(34) 
  pd(34,34) = pd(34,34) - rrt(443) * density(31) 
  pd(21,32) = pd(21,32) + rrt(444) * density(34) 
  pd(21,34) = pd(21,34) + rrt(444) * density(32) 
  pd(30,32) = pd(30,32) + rrt(444) * density(34) * 2.0d0
  pd(30,34) = pd(30,34) + rrt(444) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(444) * density(34) 
  pd(32,34) = pd(32,34) - rrt(444) * density(32) 
  pd(34,32) = pd(34,32) - rrt(444) * density(34) 
  pd(34,34) = pd(34,34) - rrt(444) * density(32) 
  pd(01,18) = pd(01,18) + rrt(445) * density(34) 
  pd(01,34) = pd(01,34) + rrt(445) * density(18) 
  pd(18,18) = pd(18,18) - rrt(445) * density(34) 
  pd(18,34) = pd(18,34) - rrt(445) * density(18) 
  pd(30,18) = pd(30,18) + rrt(445) * density(34) 
  pd(30,34) = pd(30,34) + rrt(445) * density(18) 
  pd(34,18) = pd(34,18) - rrt(445) * density(34) 
  pd(34,34) = pd(34,34) - rrt(445) * density(18) 
  pd(01,20) = pd(01,20) + rrt(446) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) + rrt(446) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(446) * density(34) 
  pd(20,34) = pd(20,34) - rrt(446) * density(20) 
  pd(30,20) = pd(30,20) + rrt(446) * density(34) 
  pd(30,34) = pd(30,34) + rrt(446) * density(20) 
  pd(34,20) = pd(34,20) - rrt(446) * density(34) 
  pd(34,34) = pd(34,34) - rrt(446) * density(20) 
  pd(01,34) = pd(01,34) + rrt(447) * density(39) 
  pd(01,39) = pd(01,39) + rrt(447) * density(34) 
  pd(21,34) = pd(21,34) + rrt(447) * density(39) 
  pd(21,39) = pd(21,39) + rrt(447) * density(34) 
  pd(34,34) = pd(34,34) - rrt(447) * density(39) 
  pd(34,39) = pd(34,39) - rrt(447) * density(34) 
  pd(39,34) = pd(39,34) - rrt(447) * density(39) 
  pd(39,39) = pd(39,39) - rrt(447) * density(34) 
  pd(21,01) = pd(21,01) + rrt(448) * density(31) * density(34) 
  pd(21,31) = pd(21,31) + rrt(448) * density(01) * density(34) 
  pd(21,34) = pd(21,34) + rrt(448) * density(01) * density(31) 
  pd(30,01) = pd(30,01) + rrt(448) * density(31) * density(34) 
  pd(30,31) = pd(30,31) + rrt(448) * density(01) * density(34) 
  pd(30,34) = pd(30,34) + rrt(448) * density(01) * density(31) 
  pd(31,01) = pd(31,01) - rrt(448) * density(31) * density(34) 
  pd(31,31) = pd(31,31) - rrt(448) * density(01) * density(34) 
  pd(31,34) = pd(31,34) - rrt(448) * density(01) * density(31) 
  pd(34,01) = pd(34,01) - rrt(448) * density(31) * density(34) 
  pd(34,31) = pd(34,31) - rrt(448) * density(01) * density(34) 
  pd(34,34) = pd(34,34) - rrt(448) * density(01) * density(31) 
  pd(21,21) = pd(21,21) + rrt(449) * density(31) * density(34) 
  pd(21,31) = pd(21,31) + rrt(449) * density(21) * density(34) 
  pd(21,34) = pd(21,34) + rrt(449) * density(21) * density(31) 
  pd(30,21) = pd(30,21) + rrt(449) * density(31) * density(34) 
  pd(30,31) = pd(30,31) + rrt(449) * density(21) * density(34) 
  pd(30,34) = pd(30,34) + rrt(449) * density(21) * density(31) 
  pd(31,21) = pd(31,21) - rrt(449) * density(31) * density(34) 
  pd(31,31) = pd(31,31) - rrt(449) * density(21) * density(34) 
  pd(31,34) = pd(31,34) - rrt(449) * density(21) * density(31) 
  pd(34,21) = pd(34,21) - rrt(449) * density(31) * density(34) 
  pd(34,31) = pd(34,31) - rrt(449) * density(21) * density(34) 
  pd(34,34) = pd(34,34) - rrt(449) * density(21) * density(31) 
  pd(21,14) = pd(21,14) + rrt(450) * density(31) * density(34) 
  pd(21,31) = pd(21,31) + rrt(450) * density(14) * density(34) 
  pd(21,34) = pd(21,34) + rrt(450) * density(14) * density(31) 
  pd(30,14) = pd(30,14) + rrt(450) * density(31) * density(34) 
  pd(30,31) = pd(30,31) + rrt(450) * density(14) * density(34) 
  pd(30,34) = pd(30,34) + rrt(450) * density(14) * density(31) 
  pd(31,14) = pd(31,14) - rrt(450) * density(31) * density(34) 
  pd(31,31) = pd(31,31) - rrt(450) * density(14) * density(34) 
  pd(31,34) = pd(31,34) - rrt(450) * density(14) * density(31) 
  pd(34,14) = pd(34,14) - rrt(450) * density(31) * density(34) 
  pd(34,31) = pd(34,31) - rrt(450) * density(14) * density(34) 
  pd(34,34) = pd(34,34) - rrt(450) * density(14) * density(31) 
  pd(21,30) = pd(21,30) + rrt(451) * density(31) * density(34) 
  pd(21,31) = pd(21,31) + rrt(451) * density(30) * density(34) 
  pd(21,34) = pd(21,34) + rrt(451) * density(30) * density(31) 
  pd(30,30) = pd(30,30) + rrt(451) * density(31) * density(34) 
  pd(30,31) = pd(30,31) + rrt(451) * density(30) * density(34) 
  pd(30,34) = pd(30,34) + rrt(451) * density(30) * density(31) 
  pd(31,30) = pd(31,30) - rrt(451) * density(31) * density(34) 
  pd(31,31) = pd(31,31) - rrt(451) * density(30) * density(34) 
  pd(31,34) = pd(31,34) - rrt(451) * density(30) * density(31) 
  pd(34,30) = pd(34,30) - rrt(451) * density(31) * density(34) 
  pd(34,31) = pd(34,31) - rrt(451) * density(30) * density(34) 
  pd(34,34) = pd(34,34) - rrt(451) * density(30) * density(31) 
  pd(21,01) = pd(21,01) + rrt(452) * density(32) * density(34) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(452) * density(01) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) + rrt(452) * density(01) * density(32) * 2.0d0
  pd(32,01) = pd(32,01) - rrt(452) * density(32) * density(34) 
  pd(32,32) = pd(32,32) - rrt(452) * density(01) * density(34) 
  pd(32,34) = pd(32,34) - rrt(452) * density(01) * density(32) 
  pd(34,01) = pd(34,01) - rrt(452) * density(32) * density(34) 
  pd(34,32) = pd(34,32) - rrt(452) * density(01) * density(34) 
  pd(34,34) = pd(34,34) - rrt(452) * density(01) * density(32) 
  pd(21,21) = pd(21,21) + rrt(453) * density(32) * density(34) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(453) * density(21) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) + rrt(453) * density(21) * density(32) * 2.0d0
  pd(32,21) = pd(32,21) - rrt(453) * density(32) * density(34) 
  pd(32,32) = pd(32,32) - rrt(453) * density(21) * density(34) 
  pd(32,34) = pd(32,34) - rrt(453) * density(21) * density(32) 
  pd(34,21) = pd(34,21) - rrt(453) * density(32) * density(34) 
  pd(34,32) = pd(34,32) - rrt(453) * density(21) * density(34) 
  pd(34,34) = pd(34,34) - rrt(453) * density(21) * density(32) 
  pd(21,14) = pd(21,14) + rrt(454) * density(32) * density(34) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(454) * density(14) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) + rrt(454) * density(14) * density(32) * 2.0d0
  pd(32,14) = pd(32,14) - rrt(454) * density(32) * density(34) 
  pd(32,32) = pd(32,32) - rrt(454) * density(14) * density(34) 
  pd(32,34) = pd(32,34) - rrt(454) * density(14) * density(32) 
  pd(34,14) = pd(34,14) - rrt(454) * density(32) * density(34) 
  pd(34,32) = pd(34,32) - rrt(454) * density(14) * density(34) 
  pd(34,34) = pd(34,34) - rrt(454) * density(14) * density(32) 
  pd(21,30) = pd(21,30) + rrt(455) * density(32) * density(34) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(455) * density(30) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) + rrt(455) * density(30) * density(32) * 2.0d0
  pd(32,30) = pd(32,30) - rrt(455) * density(32) * density(34) 
  pd(32,32) = pd(32,32) - rrt(455) * density(30) * density(34) 
  pd(32,34) = pd(32,34) - rrt(455) * density(30) * density(32) 
  pd(34,30) = pd(34,30) - rrt(455) * density(32) * density(34) 
  pd(34,32) = pd(34,32) - rrt(455) * density(30) * density(34) 
  pd(34,34) = pd(34,34) - rrt(455) * density(30) * density(32) 
  pd(01,01) = pd(01,01) + rrt(456) * density(18) * density(34) 
  pd(01,18) = pd(01,18) + rrt(456) * density(01) * density(34) 
  pd(01,34) = pd(01,34) + rrt(456) * density(01) * density(18) 
  pd(18,01) = pd(18,01) - rrt(456) * density(18) * density(34) 
  pd(18,18) = pd(18,18) - rrt(456) * density(01) * density(34) 
  pd(18,34) = pd(18,34) - rrt(456) * density(01) * density(18) 
  pd(30,01) = pd(30,01) + rrt(456) * density(18) * density(34) 
  pd(30,18) = pd(30,18) + rrt(456) * density(01) * density(34) 
  pd(30,34) = pd(30,34) + rrt(456) * density(01) * density(18) 
  pd(34,01) = pd(34,01) - rrt(456) * density(18) * density(34) 
  pd(34,18) = pd(34,18) - rrt(456) * density(01) * density(34) 
  pd(34,34) = pd(34,34) - rrt(456) * density(01) * density(18) 
  pd(01,18) = pd(01,18) + rrt(457) * density(21) * density(34) 
  pd(01,21) = pd(01,21) + rrt(457) * density(18) * density(34) 
  pd(01,34) = pd(01,34) + rrt(457) * density(18) * density(21) 
  pd(18,18) = pd(18,18) - rrt(457) * density(21) * density(34) 
  pd(18,21) = pd(18,21) - rrt(457) * density(18) * density(34) 
  pd(18,34) = pd(18,34) - rrt(457) * density(18) * density(21) 
  pd(30,18) = pd(30,18) + rrt(457) * density(21) * density(34) 
  pd(30,21) = pd(30,21) + rrt(457) * density(18) * density(34) 
  pd(30,34) = pd(30,34) + rrt(457) * density(18) * density(21) 
  pd(34,18) = pd(34,18) - rrt(457) * density(21) * density(34) 
  pd(34,21) = pd(34,21) - rrt(457) * density(18) * density(34) 
  pd(34,34) = pd(34,34) - rrt(457) * density(18) * density(21) 
  pd(01,14) = pd(01,14) + rrt(458) * density(18) * density(34) 
  pd(01,18) = pd(01,18) + rrt(458) * density(14) * density(34) 
  pd(01,34) = pd(01,34) + rrt(458) * density(14) * density(18) 
  pd(18,14) = pd(18,14) - rrt(458) * density(18) * density(34) 
  pd(18,18) = pd(18,18) - rrt(458) * density(14) * density(34) 
  pd(18,34) = pd(18,34) - rrt(458) * density(14) * density(18) 
  pd(30,14) = pd(30,14) + rrt(458) * density(18) * density(34) 
  pd(30,18) = pd(30,18) + rrt(458) * density(14) * density(34) 
  pd(30,34) = pd(30,34) + rrt(458) * density(14) * density(18) 
  pd(34,14) = pd(34,14) - rrt(458) * density(18) * density(34) 
  pd(34,18) = pd(34,18) - rrt(458) * density(14) * density(34) 
  pd(34,34) = pd(34,34) - rrt(458) * density(14) * density(18) 
  pd(01,18) = pd(01,18) + rrt(459) * density(30) * density(34) 
  pd(01,30) = pd(01,30) + rrt(459) * density(18) * density(34) 
  pd(01,34) = pd(01,34) + rrt(459) * density(18) * density(30) 
  pd(18,18) = pd(18,18) - rrt(459) * density(30) * density(34) 
  pd(18,30) = pd(18,30) - rrt(459) * density(18) * density(34) 
  pd(18,34) = pd(18,34) - rrt(459) * density(18) * density(30) 
  pd(30,18) = pd(30,18) + rrt(459) * density(30) * density(34) 
  pd(30,30) = pd(30,30) + rrt(459) * density(18) * density(34) 
  pd(30,34) = pd(30,34) + rrt(459) * density(18) * density(30) 
  pd(34,18) = pd(34,18) - rrt(459) * density(30) * density(34) 
  pd(34,30) = pd(34,30) - rrt(459) * density(18) * density(34) 
  pd(34,34) = pd(34,34) - rrt(459) * density(18) * density(30) 
  pd(01,01) = pd(01,01) + rrt(460) * density(20) * density(34) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(460) * density(01) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) + rrt(460) * density(01) * density(20) * 2.0d0
  pd(20,01) = pd(20,01) - rrt(460) * density(20) * density(34) 
  pd(20,20) = pd(20,20) - rrt(460) * density(01) * density(34) 
  pd(20,34) = pd(20,34) - rrt(460) * density(01) * density(20) 
  pd(30,01) = pd(30,01) + rrt(460) * density(20) * density(34) 
  pd(30,20) = pd(30,20) + rrt(460) * density(01) * density(34) 
  pd(30,34) = pd(30,34) + rrt(460) * density(01) * density(20) 
  pd(34,01) = pd(34,01) - rrt(460) * density(20) * density(34) 
  pd(34,20) = pd(34,20) - rrt(460) * density(01) * density(34) 
  pd(34,34) = pd(34,34) - rrt(460) * density(01) * density(20) 
  pd(01,20) = pd(01,20) + rrt(461) * density(21) * density(34) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(461) * density(20) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) + rrt(461) * density(20) * density(21) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(461) * density(21) * density(34) 
  pd(20,21) = pd(20,21) - rrt(461) * density(20) * density(34) 
  pd(20,34) = pd(20,34) - rrt(461) * density(20) * density(21) 
  pd(30,20) = pd(30,20) + rrt(461) * density(21) * density(34) 
  pd(30,21) = pd(30,21) + rrt(461) * density(20) * density(34) 
  pd(30,34) = pd(30,34) + rrt(461) * density(20) * density(21) 
  pd(34,20) = pd(34,20) - rrt(461) * density(21) * density(34) 
  pd(34,21) = pd(34,21) - rrt(461) * density(20) * density(34) 
  pd(34,34) = pd(34,34) - rrt(461) * density(20) * density(21) 
  pd(01,14) = pd(01,14) + rrt(462) * density(20) * density(34) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(462) * density(14) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) + rrt(462) * density(14) * density(20) * 2.0d0
  pd(20,14) = pd(20,14) - rrt(462) * density(20) * density(34) 
  pd(20,20) = pd(20,20) - rrt(462) * density(14) * density(34) 
  pd(20,34) = pd(20,34) - rrt(462) * density(14) * density(20) 
  pd(30,14) = pd(30,14) + rrt(462) * density(20) * density(34) 
  pd(30,20) = pd(30,20) + rrt(462) * density(14) * density(34) 
  pd(30,34) = pd(30,34) + rrt(462) * density(14) * density(20) 
  pd(34,14) = pd(34,14) - rrt(462) * density(20) * density(34) 
  pd(34,20) = pd(34,20) - rrt(462) * density(14) * density(34) 
  pd(34,34) = pd(34,34) - rrt(462) * density(14) * density(20) 
  pd(01,20) = pd(01,20) + rrt(463) * density(30) * density(34) * 2.0d0
  pd(01,30) = pd(01,30) + rrt(463) * density(20) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) + rrt(463) * density(20) * density(30) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(463) * density(30) * density(34) 
  pd(20,30) = pd(20,30) - rrt(463) * density(20) * density(34) 
  pd(20,34) = pd(20,34) - rrt(463) * density(20) * density(30) 
  pd(30,20) = pd(30,20) + rrt(463) * density(30) * density(34) 
  pd(30,30) = pd(30,30) + rrt(463) * density(20) * density(34) 
  pd(30,34) = pd(30,34) + rrt(463) * density(20) * density(30) 
  pd(34,20) = pd(34,20) - rrt(463) * density(30) * density(34) 
  pd(34,30) = pd(34,30) - rrt(463) * density(20) * density(34) 
  pd(34,34) = pd(34,34) - rrt(463) * density(20) * density(30) 
  pd(01,01) = pd(01,01) + rrt(464) * density(34) * density(39) 
  pd(01,34) = pd(01,34) + rrt(464) * density(01) * density(39) 
  pd(01,39) = pd(01,39) + rrt(464) * density(01) * density(34) 
  pd(21,01) = pd(21,01) + rrt(464) * density(34) * density(39) 
  pd(21,34) = pd(21,34) + rrt(464) * density(01) * density(39) 
  pd(21,39) = pd(21,39) + rrt(464) * density(01) * density(34) 
  pd(34,01) = pd(34,01) - rrt(464) * density(34) * density(39) 
  pd(34,34) = pd(34,34) - rrt(464) * density(01) * density(39) 
  pd(34,39) = pd(34,39) - rrt(464) * density(01) * density(34) 
  pd(39,01) = pd(39,01) - rrt(464) * density(34) * density(39) 
  pd(39,34) = pd(39,34) - rrt(464) * density(01) * density(39) 
  pd(39,39) = pd(39,39) - rrt(464) * density(01) * density(34) 
  pd(01,21) = pd(01,21) + rrt(465) * density(34) * density(39) 
  pd(01,34) = pd(01,34) + rrt(465) * density(21) * density(39) 
  pd(01,39) = pd(01,39) + rrt(465) * density(21) * density(34) 
  pd(21,21) = pd(21,21) + rrt(465) * density(34) * density(39) 
  pd(21,34) = pd(21,34) + rrt(465) * density(21) * density(39) 
  pd(21,39) = pd(21,39) + rrt(465) * density(21) * density(34) 
  pd(34,21) = pd(34,21) - rrt(465) * density(34) * density(39) 
  pd(34,34) = pd(34,34) - rrt(465) * density(21) * density(39) 
  pd(34,39) = pd(34,39) - rrt(465) * density(21) * density(34) 
  pd(39,21) = pd(39,21) - rrt(465) * density(34) * density(39) 
  pd(39,34) = pd(39,34) - rrt(465) * density(21) * density(39) 
  pd(39,39) = pd(39,39) - rrt(465) * density(21) * density(34) 
  pd(01,14) = pd(01,14) + rrt(466) * density(34) * density(39) 
  pd(01,34) = pd(01,34) + rrt(466) * density(14) * density(39) 
  pd(01,39) = pd(01,39) + rrt(466) * density(14) * density(34) 
  pd(21,14) = pd(21,14) + rrt(466) * density(34) * density(39) 
  pd(21,34) = pd(21,34) + rrt(466) * density(14) * density(39) 
  pd(21,39) = pd(21,39) + rrt(466) * density(14) * density(34) 
  pd(34,14) = pd(34,14) - rrt(466) * density(34) * density(39) 
  pd(34,34) = pd(34,34) - rrt(466) * density(14) * density(39) 
  pd(34,39) = pd(34,39) - rrt(466) * density(14) * density(34) 
  pd(39,14) = pd(39,14) - rrt(466) * density(34) * density(39) 
  pd(39,34) = pd(39,34) - rrt(466) * density(14) * density(39) 
  pd(39,39) = pd(39,39) - rrt(466) * density(14) * density(34) 
  pd(01,30) = pd(01,30) + rrt(467) * density(34) * density(39) 
  pd(01,34) = pd(01,34) + rrt(467) * density(30) * density(39) 
  pd(01,39) = pd(01,39) + rrt(467) * density(30) * density(34) 
  pd(21,30) = pd(21,30) + rrt(467) * density(34) * density(39) 
  pd(21,34) = pd(21,34) + rrt(467) * density(30) * density(39) 
  pd(21,39) = pd(21,39) + rrt(467) * density(30) * density(34) 
  pd(34,30) = pd(34,30) - rrt(467) * density(34) * density(39) 
  pd(34,34) = pd(34,34) - rrt(467) * density(30) * density(39) 
  pd(34,39) = pd(34,39) - rrt(467) * density(30) * density(34) 
  pd(39,30) = pd(39,30) - rrt(467) * density(34) * density(39) 
  pd(39,34) = pd(39,34) - rrt(467) * density(30) * density(39) 
  pd(39,39) = pd(39,39) - rrt(467) * density(30) * density(34) 
  pd(14,14) = pd(14,14) - rrt(468) * density(44) 
  pd(14,44) = pd(14,44) - rrt(468) * density(14) 
  pd(44,14) = pd(44,14) - rrt(468) * density(44) 
  pd(44,44) = pd(44,44) - rrt(468) * density(14) 
  pd(46,14) = pd(46,14) + rrt(468) * density(44) 
  pd(46,44) = pd(46,44) + rrt(468) * density(14) 
  pd(30,30) = pd(30,30) - rrt(469) * density(44) 
  pd(30,44) = pd(30,44) - rrt(469) * density(30) 
  pd(44,30) = pd(44,30) - rrt(469) * density(44) 
  pd(44,44) = pd(44,44) - rrt(469) * density(30) 
  pd(45,30) = pd(45,30) + rrt(469) * density(44) 
  pd(45,44) = pd(45,44) + rrt(469) * density(30) 
  pd(35,35) = pd(35,35) - rrt(470) * density(44) 
  pd(35,44) = pd(35,44) - rrt(470) * density(35) 
  pd(44,35) = pd(44,35) - rrt(470) * density(44) 
  pd(44,44) = pd(44,44) - rrt(470) * density(35) 
  pd(47,35) = pd(47,35) + rrt(470) * density(44) 
  pd(47,44) = pd(47,44) + rrt(470) * density(35) 
  pd(36,36) = pd(36,36) - rrt(471) * density(44) 
  pd(36,44) = pd(36,44) - rrt(471) * density(36) 
  pd(44,36) = pd(44,36) - rrt(471) * density(44) 
  pd(44,44) = pd(44,44) - rrt(471) * density(36) 
  pd(48,36) = pd(48,36) + rrt(471) * density(44) 
  pd(48,44) = pd(48,44) + rrt(471) * density(36) 
  pd(01,14) = pd(01,14) + rrt(472) * density(46) 
  pd(01,46) = pd(01,46) + rrt(472) * density(14) 
  pd(14,14) = pd(14,14) - rrt(472) * density(46) 
  pd(14,46) = pd(14,46) - rrt(472) * density(14) 
  pd(44,14) = pd(44,14) + rrt(472) * density(46) 
  pd(44,46) = pd(44,46) + rrt(472) * density(14) 
  pd(46,14) = pd(46,14) - rrt(472) * density(46) 
  pd(46,46) = pd(46,46) - rrt(472) * density(14) 
  pd(21,30) = pd(21,30) + rrt(473) * density(45) 
  pd(21,45) = pd(21,45) + rrt(473) * density(30) 
  pd(30,30) = pd(30,30) - rrt(473) * density(45) 
  pd(30,45) = pd(30,45) - rrt(473) * density(30) 
  pd(44,30) = pd(44,30) + rrt(473) * density(45) 
  pd(44,45) = pd(44,45) + rrt(473) * density(30) 
  pd(45,30) = pd(45,30) - rrt(473) * density(45) 
  pd(45,45) = pd(45,45) - rrt(473) * density(30) 
  pd(14,14) = pd(14,14) - rrt(474) * density(45) 
  pd(14,45) = pd(14,45) - rrt(474) * density(14) 
  pd(45,14) = pd(45,14) - rrt(474) * density(45) 
  pd(45,45) = pd(45,45) - rrt(474) * density(14) 
  pd(47,14) = pd(47,14) + rrt(474) * density(45) 
  pd(47,45) = pd(47,45) + rrt(474) * density(14) 
  pd(35,35) = pd(35,35) - rrt(475) * density(45) 
  pd(35,45) = pd(35,45) - rrt(475) * density(35) 
  pd(45,35) = pd(45,35) - rrt(475) * density(45) 
  pd(45,45) = pd(45,45) - rrt(475) * density(35) 
  pd(48,35) = pd(48,35) + rrt(475) * density(45) 
  pd(48,45) = pd(48,45) + rrt(475) * density(35) 
  pd(36,36) = pd(36,36) - rrt(476) * density(45) 
  pd(36,45) = pd(36,45) - rrt(476) * density(36) 
  pd(37,36) = pd(37,36) + rrt(476) * density(45) 
  pd(37,45) = pd(37,45) + rrt(476) * density(36) 
  pd(44,36) = pd(44,36) + rrt(476) * density(45) 
  pd(44,45) = pd(44,45) + rrt(476) * density(36) 
  pd(45,36) = pd(45,36) - rrt(476) * density(45) 
  pd(45,45) = pd(45,45) - rrt(476) * density(36) 
  pd(30,30) = pd(30,30) - rrt(477) * density(46) 
  pd(30,46) = pd(30,46) - rrt(477) * density(30) 
  pd(46,30) = pd(46,30) - rrt(477) * density(46) 
  pd(46,46) = pd(46,46) - rrt(477) * density(30) 
  pd(47,30) = pd(47,30) + rrt(477) * density(46) 
  pd(47,46) = pd(47,46) + rrt(477) * density(30) 
  pd(30,30) = pd(30,30) - rrt(478) * density(47) 
  pd(30,47) = pd(30,47) - rrt(478) * density(30) 
  pd(47,30) = pd(47,30) - rrt(478) * density(47) 
  pd(47,47) = pd(47,47) - rrt(478) * density(30) 
  pd(48,30) = pd(48,30) + rrt(478) * density(47) 
  pd(48,47) = pd(48,47) + rrt(478) * density(30) 
  pd(30,30) = pd(30,30) - rrt(479) * density(48) 
  pd(30,48) = pd(30,48) - rrt(479) * density(30) 
  pd(37,30) = pd(37,30) + rrt(479) * density(48) 
  pd(37,48) = pd(37,48) + rrt(479) * density(30) 
  pd(44,30) = pd(44,30) + rrt(479) * density(48) 
  pd(44,48) = pd(44,48) + rrt(479) * density(30) 
  pd(48,30) = pd(48,30) - rrt(479) * density(48) 
  pd(48,48) = pd(48,48) - rrt(479) * density(30) 
  pd(21,21) = pd(21,21) - rrt(480) * density(47) 
  pd(21,47) = pd(21,47) - rrt(480) * density(21) 
  pd(37,21) = pd(37,21) + rrt(480) * density(47) 
  pd(37,47) = pd(37,47) + rrt(480) * density(21) 
  pd(44,21) = pd(44,21) + rrt(480) * density(47) 
  pd(44,47) = pd(44,47) + rrt(480) * density(21) 
  pd(47,21) = pd(47,21) - rrt(480) * density(47) 
  pd(47,47) = pd(47,47) - rrt(480) * density(21) 
  pd(44,45) = pd(44,45) + rrt(481) * density(46) 
  pd(44,46) = pd(44,46) + rrt(481) * density(45) 
  pd(45,45) = pd(45,45) - rrt(481) * density(46) 
  pd(45,46) = pd(45,46) - rrt(481) * density(45) 
  pd(46,45) = pd(46,45) - rrt(481) * density(46) 
  pd(46,46) = pd(46,46) - rrt(481) * density(45) 
  pd(47,45) = pd(47,45) + rrt(481) * density(46) 
  pd(47,46) = pd(47,46) + rrt(481) * density(45) 
  pd(44,45) = pd(44,45) + rrt(482) * density(47) 
  pd(44,47) = pd(44,47) + rrt(482) * density(45) 
  pd(45,45) = pd(45,45) - rrt(482) * density(47) 
  pd(45,47) = pd(45,47) - rrt(482) * density(45) 
  pd(47,45) = pd(47,45) - rrt(482) * density(47) 
  pd(47,47) = pd(47,47) - rrt(482) * density(45) 
  pd(48,45) = pd(48,45) + rrt(482) * density(47) 
  pd(48,47) = pd(48,47) + rrt(482) * density(45) 
  pd(37,45) = pd(37,45) + rrt(483) * density(48) 
  pd(37,48) = pd(37,48) + rrt(483) * density(45) 
  pd(44,45) = pd(44,45) + rrt(483) * density(48) * 2.0d0
  pd(44,48) = pd(44,48) + rrt(483) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(483) * density(48) 
  pd(45,48) = pd(45,48) - rrt(483) * density(45) 
  pd(48,45) = pd(48,45) - rrt(483) * density(48) 
  pd(48,48) = pd(48,48) - rrt(483) * density(45) 
  pd(01,01) = pd(01,01) - rrt(484) * density(44)**2 
  pd(01,44) = pd(01,44) - rrt(484) * density(01) * density(44) * 2.0d0
  pd(44,01) = pd(44,01) - rrt(484) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(484) * density(01) * density(44) * 4.0d0
  pd(46,01) = pd(46,01) + rrt(484) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(484) * density(01) * density(44) * 4.0d0
  pd(02,02) = pd(02,02) - rrt(485) * density(44)**2 
  pd(02,44) = pd(02,44) - rrt(485) * density(02) * density(44) * 2.0d0
  pd(44,02) = pd(44,02) - rrt(485) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(485) * density(02) * density(44) * 4.0d0
  pd(46,02) = pd(46,02) + rrt(485) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(485) * density(02) * density(44) * 4.0d0
  pd(03,03) = pd(03,03) - rrt(486) * density(44)**2 
  pd(03,44) = pd(03,44) - rrt(486) * density(03) * density(44) * 2.0d0
  pd(44,03) = pd(44,03) - rrt(486) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(486) * density(03) * density(44) * 4.0d0
  pd(46,03) = pd(46,03) + rrt(486) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(486) * density(03) * density(44) * 4.0d0
  pd(04,04) = pd(04,04) - rrt(487) * density(44)**2 
  pd(04,44) = pd(04,44) - rrt(487) * density(04) * density(44) * 2.0d0
  pd(44,04) = pd(44,04) - rrt(487) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(487) * density(04) * density(44) * 4.0d0
  pd(46,04) = pd(46,04) + rrt(487) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(487) * density(04) * density(44) * 4.0d0
  pd(05,05) = pd(05,05) - rrt(488) * density(44)**2 
  pd(05,44) = pd(05,44) - rrt(488) * density(05) * density(44) * 2.0d0
  pd(44,05) = pd(44,05) - rrt(488) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(488) * density(05) * density(44) * 4.0d0
  pd(46,05) = pd(46,05) + rrt(488) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(488) * density(05) * density(44) * 4.0d0
  pd(06,06) = pd(06,06) - rrt(489) * density(44)**2 
  pd(06,44) = pd(06,44) - rrt(489) * density(06) * density(44) * 2.0d0
  pd(44,06) = pd(44,06) - rrt(489) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(489) * density(06) * density(44) * 4.0d0
  pd(46,06) = pd(46,06) + rrt(489) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(489) * density(06) * density(44) * 4.0d0
  pd(07,07) = pd(07,07) - rrt(490) * density(44)**2 
  pd(07,44) = pd(07,44) - rrt(490) * density(07) * density(44) * 2.0d0
  pd(44,07) = pd(44,07) - rrt(490) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(490) * density(07) * density(44) * 4.0d0
  pd(46,07) = pd(46,07) + rrt(490) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(490) * density(07) * density(44) * 4.0d0
  pd(08,08) = pd(08,08) - rrt(491) * density(44)**2 
  pd(08,44) = pd(08,44) - rrt(491) * density(08) * density(44) * 2.0d0
  pd(44,08) = pd(44,08) - rrt(491) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(491) * density(08) * density(44) * 4.0d0
  pd(46,08) = pd(46,08) + rrt(491) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(491) * density(08) * density(44) * 4.0d0
  pd(09,09) = pd(09,09) - rrt(492) * density(44)**2 
  pd(09,44) = pd(09,44) - rrt(492) * density(09) * density(44) * 2.0d0
  pd(44,09) = pd(44,09) - rrt(492) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(492) * density(09) * density(44) * 4.0d0
  pd(46,09) = pd(46,09) + rrt(492) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(492) * density(09) * density(44) * 4.0d0
  pd(10,10) = pd(10,10) - rrt(493) * density(44)**2 
  pd(10,44) = pd(10,44) - rrt(493) * density(10) * density(44) * 2.0d0
  pd(44,10) = pd(44,10) - rrt(493) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(493) * density(10) * density(44) * 4.0d0
  pd(46,10) = pd(46,10) + rrt(493) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(493) * density(10) * density(44) * 4.0d0
  pd(11,11) = pd(11,11) - rrt(494) * density(44)**2 
  pd(11,44) = pd(11,44) - rrt(494) * density(11) * density(44) * 2.0d0
  pd(44,11) = pd(44,11) - rrt(494) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(494) * density(11) * density(44) * 4.0d0
  pd(46,11) = pd(46,11) + rrt(494) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(494) * density(11) * density(44) * 4.0d0
  pd(12,12) = pd(12,12) - rrt(495) * density(44)**2 
  pd(12,44) = pd(12,44) - rrt(495) * density(12) * density(44) * 2.0d0
  pd(44,12) = pd(44,12) - rrt(495) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(495) * density(12) * density(44) * 4.0d0
  pd(46,12) = pd(46,12) + rrt(495) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(495) * density(12) * density(44) * 4.0d0
  pd(13,13) = pd(13,13) - rrt(496) * density(44)**2 
  pd(13,44) = pd(13,44) - rrt(496) * density(13) * density(44) * 2.0d0
  pd(44,13) = pd(44,13) - rrt(496) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(496) * density(13) * density(44) * 4.0d0
  pd(46,13) = pd(46,13) + rrt(496) * density(44)**2 * 2.0d0
  pd(46,44) = pd(46,44) + rrt(496) * density(13) * density(44) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(497) * density(44)**2 
  pd(21,44) = pd(21,44) - rrt(497) * density(21) * density(44) * 2.0d0
  pd(44,21) = pd(44,21) - rrt(497) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(497) * density(21) * density(44) * 4.0d0
  pd(45,21) = pd(45,21) + rrt(497) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(497) * density(21) * density(44) * 4.0d0
  pd(27,27) = pd(27,27) - rrt(498) * density(44)**2 
  pd(27,44) = pd(27,44) - rrt(498) * density(27) * density(44) * 2.0d0
  pd(44,27) = pd(44,27) - rrt(498) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(498) * density(27) * density(44) * 4.0d0
  pd(45,27) = pd(45,27) + rrt(498) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(498) * density(27) * density(44) * 4.0d0
  pd(28,28) = pd(28,28) - rrt(499) * density(44)**2 
  pd(28,44) = pd(28,44) - rrt(499) * density(28) * density(44) * 2.0d0
  pd(44,28) = pd(44,28) - rrt(499) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(499) * density(28) * density(44) * 4.0d0
  pd(45,28) = pd(45,28) + rrt(499) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(499) * density(28) * density(44) * 4.0d0
  pd(29,29) = pd(29,29) - rrt(500) * density(44)**2 
  pd(29,44) = pd(29,44) - rrt(500) * density(29) * density(44) * 2.0d0
  pd(44,29) = pd(44,29) - rrt(500) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(500) * density(29) * density(44) * 4.0d0
  pd(45,29) = pd(45,29) + rrt(500) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(500) * density(29) * density(44) * 4.0d0
  pd(22,22) = pd(22,22) - rrt(501) * density(44)**2 
  pd(22,44) = pd(22,44) - rrt(501) * density(22) * density(44) * 2.0d0
  pd(44,22) = pd(44,22) - rrt(501) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(501) * density(22) * density(44) * 4.0d0
  pd(45,22) = pd(45,22) + rrt(501) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(501) * density(22) * density(44) * 4.0d0
  pd(23,23) = pd(23,23) - rrt(502) * density(44)**2 
  pd(23,44) = pd(23,44) - rrt(502) * density(23) * density(44) * 2.0d0
  pd(44,23) = pd(44,23) - rrt(502) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(502) * density(23) * density(44) * 4.0d0
  pd(45,23) = pd(45,23) + rrt(502) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(502) * density(23) * density(44) * 4.0d0
  pd(24,24) = pd(24,24) - rrt(503) * density(44)**2 
  pd(24,44) = pd(24,44) - rrt(503) * density(24) * density(44) * 2.0d0
  pd(44,24) = pd(44,24) - rrt(503) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(503) * density(24) * density(44) * 4.0d0
  pd(45,24) = pd(45,24) + rrt(503) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(503) * density(24) * density(44) * 4.0d0
  pd(25,25) = pd(25,25) - rrt(504) * density(44)**2 
  pd(25,44) = pd(25,44) - rrt(504) * density(25) * density(44) * 2.0d0
  pd(44,25) = pd(44,25) - rrt(504) * density(44)**2 * 2.0d0
  pd(44,44) = pd(44,44) - rrt(504) * density(25) * density(44) * 4.0d0
  pd(45,25) = pd(45,25) + rrt(504) * density(44)**2 * 2.0d0
  pd(45,44) = pd(45,44) + rrt(504) * density(25) * density(44) * 4.0d0
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(49,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(49,:) = pd(49,:) * ZDPlasKin_cfg(13)
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
  double precision :: EN
  double precision :: Te
  double precision :: De
  DOUBLE PRECISION :: DTION, TEFFN2N, TEFFN2N2, TEFFN4N2 ! K
  DOUBLE PRECISION, PARAMETER :: DELTA_E_H2 = 107.7, E_10_H2 = 0.516/(1.38064852D-23*6.242D18), CHI_E_H2 = 0.02758, OMEGA_E_H2 = 4401. ! PARAMETERS FOR H2
  DOUBLE PRECISION, PARAMETER :: LLL = 0.2, H_2_MASS = 1.00784*2, E_V1_H2 = 0.515556157D0, E_V2_H2 = 1.001014015D0, E_V3_H2 = 1.456373572D0
  DOUBLE PRECISION, PARAMETER :: DELTA_E_N2 = 20.2, E_10_N2 = 0.288/(1.38064852D-23*6.242D18), CHI_E_N2 = 0.00613, OMEGA_E_N2 = 2359. ! PARAMETERS FOR N2
  DOUBLE PRECISION, PARAMETER :: N_2_MASS = 14.0067*2, E_V1_N2 = 0.288, E_V2_N2 = 0.573, E_V3_N2 = 0.855, E_V4_N2 = 1.133
  DOUBLE PRECISION, PARAMETER :: E_V5_N2 = 1.408, E_V6_N2 = 1.679, E_V7_N2 = 1.947, E_V8_N2 = 2.211
  DOUBLE PRECISION, PARAMETER :: E_20_N2 = 0.573/(1.38064852D-23*6.242D18)
  DOUBLE PRECISION, PARAMETER :: N_N2A = 1., M_N2A = 1., A_N2A = 7.8D-12, B_N2A = 2.18D2, C_N2A = 6.9D2, D_N2A = 1 ! PARAMETERS FOR N2-N2 (V-T)
  DOUBLE PRECISION :: REDUCED_MASS_N2_N2, K_10_VT_N2_N2, DELTA_VT_N2_N2_HI, DELTA_VT_N2_N2_LO
  DOUBLE PRECISION :: GAMMA_V1_N2_N2, GAMMA_V2_N2_N2, GAMMA_V3_N2_N2, GAMMA_V4_N2_N2, GAMMA_V5_N2_N2, GAMMA_V6_N2_N2, GAMMA_V7_N2_N2, GAMMA_V8_N2_N2
  DOUBLE PRECISION :: G_10_N2_N2, G_21_N2_N2, G_32_N2_N2, G_43_N2_N2, G_54_N2_N2, G_65_N2_N2, G_76_N2_N2, G_87_N2_N2 ! SINCE WE HAVE EIGHT REACTIONS FOR N2_N2 (V-T)
  DOUBLE PRECISION, PARAMETER :: N_H2A = 1., M_H2A = 2./3., A_H2A = 4.9D-12, B_H2A = 1.671D2, C_H2A = 3.94D2, D_H2A = 1 ! PARAMETERS FOR N2-H2 (V-T)
  DOUBLE PRECISION :: REDUCED_MASS_N2_H2, K_10_VT_N2_H2, DELTA_VT_N2_H2_HI, DELTA_VT_N2_H2_LO
  DOUBLE PRECISION :: GAMMA_V1_N2_H2, GAMMA_V2_N2_H2, GAMMA_V3_N2_H2, GAMMA_V4_N2_H2, GAMMA_V5_N2_H2, GAMMA_V6_N2_H2, GAMMA_V7_N2_H2, GAMMA_V8_N2_H2
  DOUBLE PRECISION :: G_10_N2_H2, G_21_N2_H2, G_32_N2_H2, G_43_N2_H2, G_54_N2_H2, G_65_N2_H2, G_76_N2_H2, G_87_N2_H2 ! SINCE WE HAVE EIGHT REACTIONS FOR N2_H2 (V-T)
  DOUBLE PRECISION, PARAMETER :: K_N2_ATOM_VT = 4D-10, E_A_N2_N = 7280 ! PARAMETERS FOR N2-N (V-T)
  DOUBLE PRECISION :: K_10_N2_N, K_21_N2_N, K_32_N2_N, K_43_N2_N, K_54_N2_N, K_65_N2_N, K_76_N2_N, K_87_N2_N ! SINCE WE HAVE EIGHT REACTIONS FOR N2_H2 (V-T)
  DOUBLE PRECISION, PARAMETER :: E_A_N2_H = 7500 ! PARAMETERS FOR N2-H (V-T)
  DOUBLE PRECISION :: K_10_N2_H, K_21_N2_H, K_32_N2_H, K_43_N2_H, K_54_N2_H, K_65_N2_H, K_76_N2_H, K_87_N2_H ! SINCE WE HAVE EIGHT REACTIONS FOR N2_H (V-T)
  DOUBLE PRECISION, PARAMETER :: N_H2B = 5.E-1, M_H2B = 0., A_H2B = 7.47D-12, B_H2B = 9.387D1, C_H2B = 0, D_H2B = 1 ! PARAMETERS FOR H2-H2 (V-T)
  DOUBLE PRECISION :: REDUCED_MASS_H2_H2, K_10_VT_H2_H2, GAMMA_V1_H2_H2, GAMMA_V2_H2_H2, GAMMA_V3_H2_H2, DELTA_VT_H2_H2_HI, DELTA_VT_H2_H2_LO
  DOUBLE PRECISION :: G_10_H2_H2, G_21_H2_H2, G_32_H2_H2 ! SINCE WE HAVE THREE REACTIONS FOR H2-H2 (V-T)
  DOUBLE PRECISION :: K_10_H2_H, K_21_H2_H, K_32_H2_H ! SINCE WE HAVE THREE REACTIONS FOR H2-H (V-T)
  DOUBLE PRECISION :: DELTA_10_H2_H, DELTA_21_H2_H, DELTA_32_H2_H ! SINCE WE HAVE THREE REACTIONS FOR H2-H (V-T)
  DOUBLE PRECISION :: DELTA_VV_0, Q_01_VV
  DOUBLE PRECISION :: DELTA_0_N2_VV, DELTA_1_N2_VV, DELTA_2_N2_VV, DELTA_3_N2_VV
  DOUBLE PRECISION :: DELTA_4_N2_VV, DELTA_5_N2_VV, DELTA_6_N2_VV, DELTA_7_N2_VV ! SINCE WE HAVE EIGHT LEVELS FOR N2(V)
  DOUBLE PRECISION :: VV_1102_N2_N2, VV_1203_N2_N2, VV_1304_N2_N2, VV_1405_N2_N2
  DOUBLE PRECISION :: VV_1506_N2_N2, VV_1607_N2_N2, VV_1708_N2_N2 ! SINCE WE HAVE SEVEN REACTIONS FOR N2(V1)-N2 (V-V)
  DOUBLE PRECISION :: F_1102_N2_N2, F_1203_N2_N2, F_1304_N2_N2, F_1405_N2_N2
  DOUBLE PRECISION :: F_1506_N2_N2, F_1607_N2_N2, F_1708_N2_N2
  DOUBLE PRECISION :: VV_2213_N2_N2, VV_2314_N2_N2, VV_2415_N2_N2
  DOUBLE PRECISION :: VV_2516_N2_N2, VV_2617_N2_N2, VV_2718_N2_N2 ! SINCE WE HAVE SIX REACTIONS FOR N2(V2)-N2 (V-V)
  DOUBLE PRECISION :: F_2213_N2_N2, F_2314_N2_N2, F_2415_N2_N2
  DOUBLE PRECISION :: F_2516_N2_N2, F_2617_N2_N2, F_2718_N2_N2
  DOUBLE PRECISION :: VV_3324_N2_N2, VV_3425_N2_N2, VV_3526_N2_N2
  DOUBLE PRECISION :: VV_3627_N2_N2, VV_3728_N2_N2 ! SINCE WE HAVE FIVE REACTIONS FOR N2(V3)-N2 (V-V)
  DOUBLE PRECISION :: F_3324_N2_N2, F_3425_N2_N2, F_3526_N2_N2
  DOUBLE PRECISION :: F_3627_N2_N2, F_3728_N2_N2
  DOUBLE PRECISION :: VV_4435_N2_N2, VV_4536_N2_N2
  DOUBLE PRECISION :: VV_4637_N2_N2, VV_4738_N2_N2 ! SINCE WE HAVE FOUR REACTIONS FOR N2(V4)-N2 (V-V)
  DOUBLE PRECISION :: F_4435_N2_N2, F_4536_N2_N2
  DOUBLE PRECISION :: F_4637_N2_N2, F_4738_N2_N2
  DOUBLE PRECISION :: VV_5546_N2_N2, VV_5647_N2_N2
  DOUBLE PRECISION :: VV_5748_N2_N2 ! SINCE WE HAVE THREE REACTIONS FOR N2(V5)-N2 (V-V)
  DOUBLE PRECISION :: F_5546_N2_N2, F_5647_N2_N2
  DOUBLE PRECISION :: F_5748_N2_N2
  DOUBLE PRECISION :: VV_6657_N2_N2, VV_6758_N2_N2 ! SINCE WE HAVE TWO REACTIONS FOR N2(V6)-N2 (V-V)
  DOUBLE PRECISION :: F_6657_N2_N2, F_6758_N2_N2
  DOUBLE PRECISION :: VV_7768_N2_N2 ! SINCE WE HAVE ONE REACTION FOR N2(V7)-N2 (V-V)
  DOUBLE PRECISION :: F_7768_N2_N2
  DOUBLE PRECISION :: DELTA_VV_H2_N2_S, Q_VV_H2_N2_S
  DOUBLE PRECISION :: DELTA_0_H2_VV, DELTA_1_H2_VV, DELTA_2_H2_VV ! SINCE WE HAVE THREE LEVELS FOR H2(V)
  DOUBLE PRECISION :: VV_1001_H2_N2, VV_1102_H2_N2, VV_1203_H2_N2, VV_1304_H2_N2
  DOUBLE PRECISION :: VV_1405_H2_N2, VV_1506_H2_N2, VV_1607_H2_N2, VV_1708_H2_N2 ! SINCE WE HAVE EIGHT REACTIONS FOR H2(V1)-N2 (V-V)
  DOUBLE PRECISION :: F_1001_H2_N2, F_1102_H2_N2, F_1203_H2_N2, F_1304_H2_N2
  DOUBLE PRECISION :: F_1405_H2_N2, F_1506_H2_N2, F_1607_H2_N2, F_1708_H2_N2
  DOUBLE PRECISION :: VV_2011_H2_N2, VV_2112_H2_N2, VV_2213_H2_N2, VV_2314_H2_N2
  DOUBLE PRECISION :: VV_2415_H2_N2, VV_2516_H2_N2, VV_2617_H2_N2, VV_2718_H2_N2 ! SINCE WE HAVE EIGHT REACTIONS FOR H2(V2)-N2 (V-V)
  DOUBLE PRECISION :: F_2011_H2_N2, F_2112_H2_N2, F_2213_H2_N2, F_2314_H2_N2
  DOUBLE PRECISION :: F_2415_H2_N2, F_2516_H2_N2, F_2617_H2_N2, F_2718_H2_N2
  DOUBLE PRECISION :: VV_3021_H2_N2, VV_3122_H2_N2, VV_3223_H2_N2, VV_3324_H2_N2
  DOUBLE PRECISION :: VV_3425_H2_N2, VV_3526_H2_N2, VV_3627_H2_N2, VV_3728_H2_N2 ! SINCE WE HAVE EIGHT REACTIONS FOR H2(V3)-N2 (V-V)
  DOUBLE PRECISION :: F_3021_H2_N2, F_3122_H2_N2, F_3223_H2_N2, F_3324_H2_N2
  DOUBLE PRECISION :: F_3425_H2_N2, F_3526_H2_N2, F_3627_H2_N2, F_3728_H2_N2
  DOUBLE PRECISION :: DELTA_VV_H2_N2_P, Q_VV_H2_N2_P
  DOUBLE PRECISION :: DELTA_0_N2_VV_P, DELTA_1_N2_VV_P, DELTA_2_N2_VV_P, DELTA_3_N2_VV_P
  DOUBLE PRECISION :: DELTA_4_N2_VV_P, DELTA_5_N2_VV_P, DELTA_6_N2_VV_P, DELTA_7_N2_VV_P ! SINCE WE HAVE EIGHT LEVELS FOR N2(V)-H2 POLY-QUANTUM
  DOUBLE PRECISION :: VV_2001_H2_N2, VV_3011_H2_N2, VV_4021_H2_N2, VV_5031_H2_N2
  DOUBLE PRECISION :: VV_6041_H2_N2, VV_7051_H2_N2, VV_8061_H2_N2 ! SINCE WE HAVE SEVEN REACTIONS FOR N2-H2(V0) (V-V) POLY-QUANTUM
  DOUBLE PRECISION :: F_2001_H2_N2, F_3011_H2_N2, F_4021_H2_N2, F_5031_H2_N2
  DOUBLE PRECISION :: F_6041_H2_N2, F_7051_H2_N2, F_8061_H2_N2
  DOUBLE PRECISION :: VV_2102_H2_N2, VV_3112_H2_N2, VV_4122_H2_N2, VV_5132_H2_N2
  DOUBLE PRECISION :: VV_6142_H2_N2, VV_7152_H2_N2, VV_8162_H2_N2 ! SINCE WE HAVE SEVEN REACTIONS FOR N2-H2(V1) (V-V) POLY-QUANTUM
  DOUBLE PRECISION :: F_2102_H2_N2, F_3112_H2_N2, F_4122_H2_N2, F_5132_H2_N2
  DOUBLE PRECISION :: F_6142_H2_N2, F_7152_H2_N2, F_8162_H2_N2
  DOUBLE PRECISION :: VV_2203_H2_N2, VV_3213_H2_N2, VV_4223_H2_N2, VV_5233_H2_N2
  DOUBLE PRECISION :: VV_6243_H2_N2, VV_7253_H2_N2, VV_8263_H2_N2 ! SINCE WE HAVE SEVEN REACTIONS FOR N2-H2(V2) (V-V) POLY-QUANTUM
  DOUBLE PRECISION :: F_2203_H2_N2, F_3213_H2_N2, F_4223_H2_N2, F_5233_H2_N2
  DOUBLE PRECISION :: F_6243_H2_N2, F_7253_H2_N2, F_8263_H2_N2
  DOUBLE PRECISION :: VV_2231_H2_H2, VV_2130_H2_H2, VV_1120_H2_H2 ! SINCE WE HAVE THREE REACTIONS FOR H2-H2 (V-V)
  DOUBLE PRECISION :: F_2231_H2_H2, F_2130_H2_H2, F_1120_H2_H2 !
  DOUBLE PRECISION, PARAMETER :: DIFF_L = 0.01, DIFF_COEFF_400 = 7.9D-1, MYPI = 16*ATAN(1./5.) - 4*ATAN(1./239.)
  DOUBLE PRECISION, PARAMETER :: H2_MASS_KG = 2.*1.6735575D-27, N2_MASS_KG = 2.*2.3258671D-26
  DOUBLE PRECISION, PARAMETER :: THE_V = 48.6, THE_AREA = 3370, ROUGHNESS = 2.1
  DOUBLE PRECISION, PARAMETER :: N2_WALL_V_GAMMA = 4.5D-4, N2_WALL_E_GAMMA = 1D-3
  DOUBLE PRECISION, PARAMETER :: H2_WALL_V_GAMMA = 1.0D-4, H2_WALL_E_GAMMA = 1D-3
  DOUBLE PRECISION :: GAMMA_D, N2_THERMAL_VEL, N2_WALL_SECOND_PART_V, N2_WALL_SECOND_PART_E
  DOUBLE PRECISION :: H2_THERMAL_VEL, H2_WALL_SECOND_PART_V, H2_WALL_SECOND_PART_E
  DOUBLE PRECISION, PARAMETER :: V_OVER_A = 7.0D-3, LAMBDA = 1.0D-2, SITE_DENSITY = 1.0D15
  DOUBLE PRECISION, PARAMETER :: EV_PER_J = 6.242D18, KB = 1.38064852D-23
  DOUBLE PRECISION, PARAMETER :: METAL_ED = 0.2D0, JUMP_FREQUENCY = 1.0D13
  DOUBLE PRECISION :: STVOL, DSURF, EZ, KADS_N, KADS_H, KADS_NH, KADS_NH2
  DOUBLE PRECISION :: KER_N, KER_H, KER_N_HS, KER_NH_HS, KER_NH2_HS
  DOUBLE PRECISION :: KER_H_NS, KER_H_NHS, KER_H_NH2S, KER_H2_NHS
  DOUBLE PRECISION :: KLH_S14, KLH_S15, KLH_S16
  DOUBLE PRECISION :: KDADS_N2_0, KDADS_N2_1, KDADS_N2_2, KDADS_N2_3, KDADS_N2_4
  DOUBLE PRECISION :: KDADS_N2_5, KDADS_N2_6, KDADS_N2_7, KDADS_N2_8, KDADS_N2_E
  DOUBLE PRECISION :: KDADS_H2_0, KDADS_H2_1, KDADS_H2_2, KDADS_H2_3, KDADS_H2_E
  DOUBLE PRECISION :: GAMMA_N2_0, GAMMA_N2_1, GAMMA_N2_2, GAMMA_N2_3, GAMMA_N2_4
  DOUBLE PRECISION :: GAMMA_N2_5, GAMMA_N2_6, GAMMA_N2_7, GAMMA_N2_8
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  De  = ZDPlasKin_cfg(6)
  DTION = 3.14D0 / ( 6.0D0 * 1.3807D-16 ) * 1.6605D-24 * ( 1.0D-17 * EN )**2
  TEFFN2N = TGAS + DTION * (28.0D0*14.0D0)/(28.0D0+14.0D0) * (1.60*2.45D19*(300/TGAS))**2
  TEFFN2N2 = TGAS + DTION * (28.0D0*28.0D0)/(28.0D0+28.0D0) * (1.60*2.45D19*(300/TGAS))**2
  TEFFN4N2 = TGAS + DTION * (56.0D0*28.0D0)/(28.0D0+56.0D0) * (2.32*2.45D19*(300/TGAS))**2
  K_10_VT_N2_N2 = A_N2A*(TGAS**N_N2A)*EXP(-B_N2A/(TGAS**(1./3.))+C_N2A/(TGAS**M_N2A))/(1-D_N2A*EXP(-E_10_N2/TGAS))
  REDUCED_MASS_N2_N2 = (N_2_MASS*N_2_MASS)/(N_2_MASS+N_2_MASS)
  DELTA_VT_N2_N2_LO = 0.427*DELTA_E_N2*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  DELTA_VT_N2_N2_HI = 1.87*DELTA_E_N2*(LLL**(2./3.))*((REDUCED_MASS_N2_N2/(E_10_N2*TGAS))**(1./3.))
  GAMMA_V1_N2_N2 = 0.32*(E_V1_N2-0.)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V1_N2_N2 .LT. 20) THEN; G_10_N2_N2= (0.+1.)*(1.-CHI_E_N2)/(1.-(0.+1.)*CHI_E_N2)*EXP(0.*DELTA_VT_N2_N2_LO); ELSE; G_10_N2_N2 = (0.+1.)*(1.-CHI_E_N2)/(1.-(0.+1.)*CHI_E_N2)*EXP(0.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V2_N2_N2 = 0.32*(E_V2_N2-E_V1_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V2_N2_N2 .LT. 20) THEN; G_21_N2_N2= (1.+1.)*(1.-CHI_E_N2)/(1.-(1.+1.)*CHI_E_N2)*EXP(1.*DELTA_VT_N2_N2_LO); ELSE; G_21_N2_N2 = (1.+1.)*(1.-CHI_E_N2)/(1.-(1.+1.)*CHI_E_N2)*EXP(1.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V3_N2_N2 = 0.32*(E_V3_N2-E_V2_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V3_N2_N2 .LT. 20) THEN; G_32_N2_N2= (2.+1.)*(1.-CHI_E_N2)/(1.-(2.+1.)*CHI_E_N2)*EXP(2.*DELTA_VT_N2_N2_LO); ELSE; G_32_N2_N2 = (2.+1.)*(1.-CHI_E_N2)/(1.-(2.+1.)*CHI_E_N2)*EXP(2.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V4_N2_N2 = 0.32*(E_V4_N2-E_V3_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V4_N2_N2 .LT. 20) THEN; G_43_N2_N2= (3.+1.)*(1.-CHI_E_N2)/(1.-(3.+1.)*CHI_E_N2)*EXP(3.*DELTA_VT_N2_N2_LO); ELSE; G_43_N2_N2 = (3.+1.)*(1.-CHI_E_N2)/(1.-(3.+1.)*CHI_E_N2)*EXP(3.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V5_N2_N2 = 0.32*(E_V5_N2-E_V4_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V5_N2_N2 .LT. 20) THEN; G_54_N2_N2= (4.+1.)*(1.-CHI_E_N2)/(1.-(4.+1.)*CHI_E_N2)*EXP(4.*DELTA_VT_N2_N2_LO); ELSE; G_54_N2_N2 = (4.+1.)*(1.-CHI_E_N2)/(1.-(4.+1.)*CHI_E_N2)*EXP(4.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V6_N2_N2 = 0.32*(E_V6_N2-E_V5_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V6_N2_N2 .LT. 20) THEN; G_65_N2_N2= (5.+1.)*(1.-CHI_E_N2)/(1.-(5.+1.)*CHI_E_N2)*EXP(5.*DELTA_VT_N2_N2_LO); ELSE; G_65_N2_N2 = (5.+1.)*(1.-CHI_E_N2)/(1.-(5.+1.)*CHI_E_N2)*EXP(5.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V7_N2_N2 = 0.32*(E_V7_N2-E_V6_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V7_N2_N2 .LT. 20) THEN; G_76_N2_N2= (6.+1.)*(1.-CHI_E_N2)/(1.-(6.+1.)*CHI_E_N2)*EXP(6.*DELTA_VT_N2_N2_LO); ELSE; G_76_N2_N2 = (6.+1.)*(1.-CHI_E_N2)/(1.-(6.+1.)*CHI_E_N2)*EXP(6.*DELTA_VT_N2_N2_HI);END IF
  GAMMA_V8_N2_N2 = 0.32*(E_V8_N2-E_V6_N2)*LLL*SQRT(REDUCED_MASS_N2_N2/TGAS)
  IF (GAMMA_V8_N2_N2 .LT. 20) THEN; G_87_N2_N2= (7.+1.)*(1.-CHI_E_N2)/(1.-(7.+1.)*CHI_E_N2)*EXP(7.*DELTA_VT_N2_N2_LO); ELSE; G_87_N2_N2 = (7.+1.)*(1.-CHI_E_N2)/(1.-(7.+1.)*CHI_E_N2)*EXP(7.*DELTA_VT_N2_N2_HI);END IF
  K_10_VT_N2_H2 = A_H2A*(TGAS**N_H2A)*EXP(-B_H2A/(TGAS**(1./3.))+C_H2A/(TGAS**M_H2A))/(1-D_H2A*EXP(-E_10_N2/TGAS))
  REDUCED_MASS_N2_H2 = (N_2_MASS*H_2_MASS)/(N_2_MASS+H_2_MASS)
  DELTA_VT_N2_H2_LO = 0.427*DELTA_E_N2*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  DELTA_VT_N2_H2_HI = 1.87*DELTA_E_N2*(LLL**(2./3.))*((REDUCED_MASS_N2_H2/(E_10_N2*TGAS))**(1./3.))
  GAMMA_V1_N2_H2 = 0.32*(E_V1_N2-0.)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V1_N2_H2 .LT. 20) THEN; G_10_N2_H2= (0.+1.)*(1.-CHI_E_N2)/(1.-(0.+1.)*CHI_E_N2)*EXP(0.*DELTA_VT_N2_H2_LO); ELSE; G_10_N2_H2 = (0.+1.)*(1.-CHI_E_N2)/(1.-(0.+1.)*CHI_E_N2)*EXP(0.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V2_N2_H2 = 0.32*(E_V2_N2-E_V1_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V2_N2_H2 .LT. 20) THEN; G_21_N2_H2= (1.+1.)*(1.-CHI_E_N2)/(1.-(1.+1.)*CHI_E_N2)*EXP(1.*DELTA_VT_N2_H2_LO); ELSE; G_21_N2_H2 = (1.+1.)*(1.-CHI_E_N2)/(1.-(1.+1.)*CHI_E_N2)*EXP(1.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V3_N2_H2 = 0.32*(E_V3_N2-E_V2_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V3_N2_H2 .LT. 20) THEN; G_32_N2_H2= (2.+1.)*(1.-CHI_E_N2)/(1.-(2.+1.)*CHI_E_N2)*EXP(2.*DELTA_VT_N2_H2_LO); ELSE; G_32_N2_H2 = (2.+1.)*(1.-CHI_E_N2)/(1.-(2.+1.)*CHI_E_N2)*EXP(2.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V4_N2_H2 = 0.32*(E_V4_N2-E_V3_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V4_N2_H2 .LT. 20) THEN; G_43_N2_H2= (3.+1.)*(1.-CHI_E_N2)/(1.-(3.+1.)*CHI_E_N2)*EXP(3.*DELTA_VT_N2_H2_LO); ELSE; G_43_N2_H2 = (3.+1.)*(1.-CHI_E_N2)/(1.-(3.+1.)*CHI_E_N2)*EXP(3.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V5_N2_H2 = 0.32*(E_V5_N2-E_V4_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V5_N2_H2 .LT. 20) THEN; G_54_N2_H2= (4.+1.)*(1.-CHI_E_N2)/(1.-(4.+1.)*CHI_E_N2)*EXP(4.*DELTA_VT_N2_H2_LO); ELSE; G_54_N2_H2 = (4.+1.)*(1.-CHI_E_N2)/(1.-(4.+1.)*CHI_E_N2)*EXP(4.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V6_N2_H2 = 0.32*(E_V6_N2-E_V5_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V6_N2_H2 .LT. 20) THEN; G_65_N2_H2= (5.+1.)*(1.-CHI_E_N2)/(1.-(5.+1.)*CHI_E_N2)*EXP(5.*DELTA_VT_N2_H2_LO); ELSE; G_54_N2_H2 = (5.+1.)*(1.-CHI_E_N2)/(1.-(5.+1.)*CHI_E_N2)*EXP(5.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V7_N2_H2 = 0.32*(E_V7_N2-E_V6_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V7_N2_H2 .LT. 20) THEN; G_76_N2_H2= (6.+1.)*(1.-CHI_E_N2)/(1.-(6.+1.)*CHI_E_N2)*EXP(6.*DELTA_VT_N2_H2_LO); ELSE; G_76_N2_H2 = (6.+1.)*(1.-CHI_E_N2)/(1.-(6.+1.)*CHI_E_N2)*EXP(6.*DELTA_VT_N2_H2_HI);END IF
  GAMMA_V8_N2_H2 = 0.32*(E_V8_N2-E_V7_N2)*LLL*SQRT(REDUCED_MASS_N2_H2/TGAS)
  IF (GAMMA_V8_N2_H2 .LT. 20) THEN; G_87_N2_H2= (7.+1.)*(1.-CHI_E_N2)/(1.-(7.+1.)*CHI_E_N2)*EXP(7.*DELTA_VT_N2_H2_LO); ELSE; G_87_N2_H2 = (7.+1.)*(1.-CHI_E_N2)/(1.-(7.+1.)*CHI_E_N2)*EXP(7.*DELTA_VT_N2_H2_HI);END IF
  IF (0.065*E_V1_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_10_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V1_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_10_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V2_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_21_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V2_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_21_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V3_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_32_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V3_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_32_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V4_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_43_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V4_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_43_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V5_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_54_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V5_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_54_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V6_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_65_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V6_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_65_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V7_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_76_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V7_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_76_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.065*E_V8_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_N) THEN; K_87_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_N/TGAS+(0.065*E_V8_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_87_N2_N = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V1_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_10_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V1_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_10_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V2_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_21_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V2_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_21_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V3_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_32_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V3_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_32_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V4_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_43_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V4_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_43_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V5_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_54_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V5_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_54_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V6_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_65_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V6_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_65_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V7_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_76_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V7_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_76_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  IF (0.105*E_V8_N2/(1.38064852D-23*6.242D18) .LE. E_A_N2_H) THEN; K_87_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.)*EXP(-E_A_N2_H/TGAS+(0.105*E_V8_N2/(1.38064852D-23*6.242D18))/TGAS); ELSE; K_87_N2_H = K_N2_ATOM_VT*SQRT(TGAS/300.) ;END IF
  K_10_VT_H2_H2 = A_H2B*(TGAS**N_H2B)*EXP(-B_H2B/(TGAS**(1./3.))+C_H2B/(TGAS**M_H2B))/(1-D_H2B*EXP(-E_10_H2/TGAS))
  REDUCED_MASS_H2_H2 = (H_2_MASS*H_2_MASS)/(H_2_MASS+H_2_MASS)
  DELTA_VT_H2_H2_LO = 0.427*DELTA_E_H2*LLL*SQRT(REDUCED_MASS_H2_H2/TGAS)
  DELTA_VT_H2_H2_HI = 1.87*DELTA_E_H2*(LLL**(2./3.))*((REDUCED_MASS_H2_H2/(E_10_H2*TGAS))**(1./3.))
  GAMMA_V1_H2_H2 = 0.32*(E_V1_H2-0.)*LLL*SQRT(REDUCED_MASS_H2_H2/TGAS)
  IF (GAMMA_V1_H2_H2 .LT. 20) THEN; G_10_H2_H2 = (0.+1.)*(1.-CHI_E_H2)/(1.-(0.+1.)*CHI_E_H2)*EXP(0.*DELTA_VT_H2_H2_LO); ELSE; G_10_H2_H2 = (0.+1.)*(1.-CHI_E_H2)/(1.-(0.+1.)*CHI_E_H2)*EXP(0.*DELTA_VT_H2_H2_HI);END IF
  GAMMA_V2_H2_H2 = 0.32*(E_V2_H2-E_V1_H2)*LLL*SQRT(REDUCED_MASS_H2_H2/TGAS)
  IF (GAMMA_V2_H2_H2 .LT. 20) THEN; G_21_H2_H2 = (1.+1.)*(1.-CHI_E_H2)/(1.-(1.+1.)*CHI_E_H2)*EXP(1.*DELTA_VT_H2_H2_LO); ELSE; G_21_H2_H2 = (1.+1.)*(1.-CHI_E_H2)/(1.-(1.+1.)*CHI_E_H2)*EXP(1.*DELTA_VT_H2_H2_HI);END IF
  GAMMA_V3_H2_H2 = 0.32*(E_V3_H2-E_V2_H2)*LLL*SQRT(REDUCED_MASS_H2_H2/TGAS)
  IF (GAMMA_V3_H2_H2 .LT. 20) THEN; G_32_H2_H2 = (2.+1.)*(1.-CHI_E_H2)/(1.-(2.+1.)*CHI_E_H2)*EXP(2.*DELTA_VT_H2_H2_LO); ELSE; G_32_H2_H2 = (2.+1.)*(1.-CHI_E_H2)/(1.-(2.+1.)*CHI_E_H2)*EXP(1.*DELTA_VT_H2_H2_HI);END IF
  DELTA_10_H2_H = (0.+1.)*(1.+0.*DELTA_E_H2/E_10_H2)
  K_10_H2_H = DELTA_10_H2_H*((1.-0.*5.67D-2)**2.66)*2.4D-8*EXP(-(162.6/((TGAS)**(1./3.)))*(1-0.*5.67D-2)**0.681)
  DELTA_21_H2_H = (1.+1.)*(1.+1.*DELTA_E_H2/E_10_H2)
  K_21_H2_H = DELTA_21_H2_H*((1.-1.*5.67D-2)**2.66)*2.4D-8*EXP(-(162.6/((TGAS)**(1./3.)))*(1-1.*5.67D-2)**0.681)
  DELTA_32_H2_H = (2.+1.)*(1.+2.*DELTA_E_H2/E_10_H2)
  K_32_H2_H = DELTA_32_H2_H*((1.-2.*5.67D-2)**2.66)*2.4D-8*EXP(-(162.6/((TGAS)**(1./3.)))*(1-2.*5.67D-2)**0.681)
  DELTA_VV_0 = 6.8/SQRT(TGAS)
  Q_01_VV = 2.5D-14*(TGAS**(3./2.))
  DELTA_0_N2_VV = (0.+1.)*(1.+0.*DELTA_E_N2/E_10_N2)
  DELTA_1_N2_VV = (1.+1.)*(1.+1.*DELTA_E_N2/E_10_N2)
  DELTA_2_N2_VV = (2.+1.)*(1.+2.*DELTA_E_N2/E_10_N2)
  DELTA_3_N2_VV = (3.+1.)*(1.+3.*DELTA_E_N2/E_10_N2)
  DELTA_4_N2_VV = (4.+1.)*(1.+4.*DELTA_E_N2/E_10_N2)
  DELTA_5_N2_VV = (5.+1.)*(1.+5.*DELTA_E_N2/E_10_N2)
  DELTA_6_N2_VV = (6.+1.)*(1.+6.*DELTA_E_N2/E_10_N2)
  DELTA_7_N2_VV = (7.+1.)*(1.+7.*DELTA_E_N2/E_10_N2)
  F_1102_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V2_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1102_N2_N2 = Q_01_VV*DELTA_1_N2_VV*DELTA_0_N2_VV*F_1102_N2_N2*(1.5-0.5*F_1102_N2_N2)
  F_1203_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V3_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1203_N2_N2 = Q_01_VV*DELTA_2_N2_VV*DELTA_0_N2_VV*F_1203_N2_N2*(1.5-0.5*F_1203_N2_N2)
  F_1304_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1304_N2_N2 = Q_01_VV*DELTA_3_N2_VV*DELTA_0_N2_VV*F_1304_N2_N2*(1.5-0.5*F_1304_N2_N2)
  F_1405_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1405_N2_N2 = Q_01_VV*DELTA_4_N2_VV*DELTA_0_N2_VV*F_1405_N2_N2*(1.5-0.5*F_1405_N2_N2)
  F_1506_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1506_N2_N2 = Q_01_VV*DELTA_5_N2_VV*DELTA_0_N2_VV*F_1506_N2_N2*(1.5-0.5*F_1506_N2_N2)
  F_1607_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1607_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_0_N2_VV*F_1607_N2_N2*(1.5-0.5*F_1607_N2_N2)
  F_1708_N2_N2 = EXP((-DELTA_VV_0*((E_V1_N2-0.)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_1708_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_0_N2_VV*F_1708_N2_N2*(1.5-0.5*F_1708_N2_N2)
  F_2213_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V3_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2213_N2_N2 = Q_01_VV*DELTA_2_N2_VV*DELTA_1_N2_VV*F_2213_N2_N2*(1.5-0.5*F_2213_N2_N2)
  F_2314_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2314_N2_N2 = Q_01_VV*DELTA_3_N2_VV*DELTA_1_N2_VV*F_2314_N2_N2*(1.5-0.5*F_2314_N2_N2)
  F_2415_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2415_N2_N2 = Q_01_VV*DELTA_4_N2_VV*DELTA_1_N2_VV*F_2415_N2_N2*(1.5-0.5*F_2415_N2_N2)
  F_2516_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2516_N2_N2 = Q_01_VV*DELTA_5_N2_VV*DELTA_1_N2_VV*F_2516_N2_N2*(1.5-0.5*F_2516_N2_N2)
  F_2617_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2617_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_1_N2_VV*F_2617_N2_N2*(1.5-0.5*F_2617_N2_N2)
  F_2718_N2_N2 = EXP((-DELTA_VV_0*((E_V2_N2-E_V1_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_2718_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_1_N2_VV*F_2718_N2_N2*(1.5-0.5*F_2718_N2_N2)
  F_3324_N2_N2 = EXP((-DELTA_VV_0*((E_V3_N2-E_V2_N2)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_3324_N2_N2 = Q_01_VV*DELTA_3_N2_VV*DELTA_2_N2_VV*F_3324_N2_N2*(1.5-0.5*F_3324_N2_N2)
  F_3425_N2_N2 = EXP((-DELTA_VV_0*((E_V3_N2-E_V2_N2)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_3425_N2_N2 = Q_01_VV*DELTA_4_N2_VV*DELTA_2_N2_VV*F_3425_N2_N2*(1.5-0.5*F_3425_N2_N2)
  F_3526_N2_N2 = EXP((-DELTA_VV_0*((E_V3_N2-E_V2_N2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_3526_N2_N2 = Q_01_VV*DELTA_5_N2_VV*DELTA_2_N2_VV*F_3526_N2_N2*(1.5-0.5*F_3526_N2_N2)
  F_3627_N2_N2 = EXP((-DELTA_VV_0*((E_V3_N2-E_V2_N2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_3627_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_2_N2_VV*F_3627_N2_N2*(1.5-0.5*F_3627_N2_N2)
  F_3728_N2_N2 = EXP((-DELTA_VV_0*((E_V3_N2-E_V2_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_3728_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_2_N2_VV*F_3728_N2_N2*(1.5-0.5*F_3728_N2_N2)
  F_4435_N2_N2 = EXP((-DELTA_VV_0*((E_V4_N2-E_V3_N2)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_4435_N2_N2 = Q_01_VV*DELTA_4_N2_VV*DELTA_3_N2_VV*F_4435_N2_N2*(1.5-0.5*F_4435_N2_N2)
  F_4536_N2_N2 = EXP((-DELTA_VV_0*((E_V4_N2-E_V3_N2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_4536_N2_N2 = Q_01_VV*DELTA_5_N2_VV*DELTA_3_N2_VV*F_4536_N2_N2*(1.5-0.5*F_4536_N2_N2)
  F_4637_N2_N2 = EXP((-DELTA_VV_0*((E_V4_N2-E_V3_N2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_4637_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_3_N2_VV*F_4637_N2_N2*(1.5-0.5*F_4637_N2_N2)
  F_4738_N2_N2 = EXP((-DELTA_VV_0*((E_V4_N2-E_V3_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_4738_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_3_N2_VV*F_4738_N2_N2*(1.5-0.5*F_4738_N2_N2)
  F_5546_N2_N2 = EXP((-DELTA_VV_0*((E_V5_N2-E_V4_N2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_5546_N2_N2 = Q_01_VV*DELTA_5_N2_VV*DELTA_4_N2_VV*F_5546_N2_N2*(1.5-0.5*F_5546_N2_N2)
  F_5647_N2_N2 = EXP((-DELTA_VV_0*((E_V5_N2-E_V4_N2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_5647_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_4_N2_VV*F_5647_N2_N2*(1.5-0.5*F_5647_N2_N2)
  F_5748_N2_N2 = EXP((-DELTA_VV_0*((E_V5_N2-E_V4_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_5748_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_4_N2_VV*F_5748_N2_N2*(1.5-0.5*F_5748_N2_N2)
  F_6657_N2_N2 = EXP((-DELTA_VV_0*((E_V6_N2-E_V5_N2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_6657_N2_N2 = Q_01_VV*DELTA_6_N2_VV*DELTA_5_N2_VV*F_6657_N2_N2*(1.5-0.5*F_6657_N2_N2)
  F_6758_N2_N2 = EXP((-DELTA_VV_0*((E_V6_N2-E_V5_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_6758_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_5_N2_VV*F_6758_N2_N2*(1.5-0.5*F_6758_N2_N2)
  F_7768_N2_N2 = EXP((-DELTA_VV_0*((E_V7_N2-E_V6_N2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_N2))
  VV_7768_N2_N2 = Q_01_VV*DELTA_7_N2_VV*DELTA_6_N2_VV*F_7768_N2_N2*(1.5-0.5*F_7768_N2_N2)
  DELTA_VV_H2_N2_S = 144./SQRT(TGAS)
  Q_VV_H2_N2_S = 1.9D-13*(TGAS**(3./2.))
  DELTA_0_H2_VV = (0.+1.)*(1.+0.*DELTA_E_H2/E_10_H2)
  DELTA_1_H2_VV = (1.+1.)*(1.+1.*DELTA_E_H2/E_10_H2)
  DELTA_2_H2_VV = (2.+1.)*(1.+2.*DELTA_E_H2/E_10_H2)
  F_1001_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V1_N2-0.))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1001_H2_N2 = Q_VV_H2_N2_S*DELTA_0_N2_VV*DELTA_0_H2_VV*F_1001_H2_N2*(1.5-0.5*F_1001_H2_N2)
  F_1102_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V2_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1102_H2_N2 = Q_VV_H2_N2_S*DELTA_1_N2_VV*DELTA_0_H2_VV*F_1102_H2_N2*(1.5-0.5*F_1102_H2_N2)
  F_1203_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V3_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1203_H2_N2 = Q_VV_H2_N2_S*DELTA_2_N2_VV*DELTA_0_H2_VV*F_1203_H2_N2*(1.5-0.5*F_1203_H2_N2)
  F_1304_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1304_H2_N2 = Q_VV_H2_N2_S*DELTA_3_N2_VV*DELTA_0_H2_VV*F_1304_H2_N2*(1.5-0.5*F_1304_H2_N2)
  F_1405_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1405_H2_N2 = Q_VV_H2_N2_S*DELTA_4_N2_VV*DELTA_0_H2_VV*F_1405_H2_N2*(1.5-0.5*F_1405_H2_N2)
  F_1506_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1506_H2_N2 = Q_VV_H2_N2_S*DELTA_5_N2_VV*DELTA_0_H2_VV*F_1506_H2_N2*(1.5-0.5*F_1506_H2_N2)
  F_1607_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1607_H2_N2 = Q_VV_H2_N2_S*DELTA_6_N2_VV*DELTA_0_H2_VV*F_1607_H2_N2*(1.5-0.5*F_1607_H2_N2)
  F_1708_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V1_H2-0.)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_1708_H2_N2 = Q_VV_H2_N2_S*DELTA_7_N2_VV*DELTA_0_H2_VV*F_1708_H2_N2*(1.5-0.5*F_1708_H2_N2)
  F_2011_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V1_N2-0.))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2011_H2_N2 = Q_VV_H2_N2_S*DELTA_0_N2_VV*DELTA_1_H2_VV*F_2011_H2_N2*(1.5-0.5*F_2011_H2_N2)
  F_2112_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V2_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2112_H2_N2 = Q_VV_H2_N2_S*DELTA_1_N2_VV*DELTA_1_H2_VV*F_2112_H2_N2*(1.5-0.5*F_2112_H2_N2)
  F_2213_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V3_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2213_H2_N2 = Q_VV_H2_N2_S*DELTA_2_N2_VV*DELTA_1_H2_VV*F_2213_H2_N2*(1.5-0.5*F_2213_H2_N2)
  F_2314_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2314_H2_N2 = Q_VV_H2_N2_S*DELTA_3_N2_VV*DELTA_1_H2_VV*F_2314_H2_N2*(1.5-0.5*F_2314_H2_N2)
  F_2415_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2415_H2_N2 = Q_VV_H2_N2_S*DELTA_4_N2_VV*DELTA_1_H2_VV*F_2415_H2_N2*(1.5-0.5*F_2415_H2_N2)
  F_2516_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2516_H2_N2 = Q_VV_H2_N2_S*DELTA_5_N2_VV*DELTA_1_H2_VV*F_2516_H2_N2*(1.5-0.5*F_2516_H2_N2)
  F_2617_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2617_H2_N2 = Q_VV_H2_N2_S*DELTA_6_N2_VV*DELTA_1_H2_VV*F_2617_H2_N2*(1.5-0.5*F_2617_H2_N2)
  F_2718_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V2_H2-E_V1_H2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_2718_H2_N2 = Q_VV_H2_N2_S*DELTA_7_N2_VV*DELTA_1_H2_VV*F_2718_H2_N2*(1.5-0.5*F_2718_H2_N2)
  F_3021_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V1_N2-0.))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3021_H2_N2 = Q_VV_H2_N2_S*DELTA_0_N2_VV*DELTA_2_H2_VV*F_3021_H2_N2*(1.5-0.5*F_3021_H2_N2)
  F_3122_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V2_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3122_H2_N2 = Q_VV_H2_N2_S*DELTA_1_N2_VV*DELTA_2_H2_VV*F_3122_H2_N2*(1.5-0.5*F_3122_H2_N2)
  F_3223_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V3_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3223_H2_N2 = Q_VV_H2_N2_S*DELTA_2_N2_VV*DELTA_2_H2_VV*F_3223_H2_N2*(1.5-0.5*F_3223_H2_N2)
  F_3324_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V4_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3324_H2_N2 = Q_VV_H2_N2_S*DELTA_3_N2_VV*DELTA_2_H2_VV*F_3324_H2_N2*(1.5-0.5*F_3324_H2_N2)
  F_3425_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V5_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3425_H2_N2 = Q_VV_H2_N2_S*DELTA_4_N2_VV*DELTA_2_H2_VV*F_3425_H2_N2*(1.5-0.5*F_3425_H2_N2)
  F_3526_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V6_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3526_H2_N2 = Q_VV_H2_N2_S*DELTA_5_N2_VV*DELTA_2_H2_VV*F_3526_H2_N2*(1.5-0.5*F_3526_H2_N2)
  F_3627_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3627_H2_N2 = Q_VV_H2_N2_S*DELTA_6_N2_VV*DELTA_2_H2_VV*F_3627_H2_N2*(1.5-0.5*F_3627_H2_N2)
  F_3728_H2_N2 = EXP((-DELTA_VV_H2_N2_S*ABS((E_V3_H2-E_V2_H2)-(E_V8_N2-E_V7_N2))/(1.38064852D-23*6.242D18))/(E_10_H2-E_10_N2))
  VV_3728_H2_N2 = Q_VV_H2_N2_S*DELTA_7_N2_VV*DELTA_2_H2_VV*F_3728_H2_N2*(1.5-0.5*F_3728_H2_N2)
  DELTA_VV_H2_N2_P = 37./SQRT(TGAS)
  Q_VV_H2_N2_P = 2.4D-15*(TGAS**(3./2.))
  DELTA_0_N2_VV_P = 0.5*DELTA_0_N2_VV*(DELTA_0_N2_VV+1)
  DELTA_1_N2_VV_P = 0.5*DELTA_1_N2_VV*(DELTA_1_N2_VV+1)
  DELTA_2_N2_VV_P = 0.5*DELTA_2_N2_VV*(DELTA_2_N2_VV+1)
  DELTA_3_N2_VV_P = 0.5*DELTA_3_N2_VV*(DELTA_3_N2_VV+1)
  DELTA_4_N2_VV_P = 0.5*DELTA_4_N2_VV*(DELTA_4_N2_VV+1)
  DELTA_5_N2_VV_P = 0.5*DELTA_5_N2_VV*(DELTA_5_N2_VV+1)
  DELTA_6_N2_VV_P = 0.5*DELTA_6_N2_VV*(DELTA_6_N2_VV+1)
  DELTA_7_N2_VV_P = 0.5*DELTA_7_N2_VV*(DELTA_7_N2_VV+1)
  F_2001_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V2_N2-0.))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_2001_H2_N2 = Q_VV_H2_N2_P*DELTA_0_N2_VV_P*DELTA_0_H2_VV*F_2001_H2_N2*(1.5-0.5*F_2001_H2_N2)
  F_3011_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V3_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_3011_H2_N2 = Q_VV_H2_N2_P*DELTA_1_N2_VV_P*DELTA_0_H2_VV*F_3011_H2_N2*(1.5-0.5*F_3011_H2_N2)
  F_4021_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V4_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_4021_H2_N2 = Q_VV_H2_N2_P*DELTA_2_N2_VV_P*DELTA_0_H2_VV*F_4021_H2_N2*(1.5-0.5*F_4021_H2_N2)
  F_5031_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V5_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_5031_H2_N2 = Q_VV_H2_N2_P*DELTA_3_N2_VV_P*DELTA_0_H2_VV*F_5031_H2_N2*(1.5-0.5*F_5031_H2_N2)
  F_6041_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V6_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_6041_H2_N2 = Q_VV_H2_N2_P*DELTA_4_N2_VV_P*DELTA_0_H2_VV*F_6041_H2_N2*(1.5-0.5*F_6041_H2_N2)
  F_7051_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V7_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_7051_H2_N2 = Q_VV_H2_N2_P*DELTA_5_N2_VV_P*DELTA_0_H2_VV*F_7051_H2_N2*(1.5-0.5*F_7051_H2_N2)
  F_8061_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V1_H2-0.)-(E_V8_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_8061_H2_N2 = Q_VV_H2_N2_P*DELTA_6_N2_VV_P*DELTA_0_H2_VV*F_8061_H2_N2*(1.5-0.5*F_8061_H2_N2)
  F_2102_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V2_N2-0.))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_2102_H2_N2 = Q_VV_H2_N2_P*DELTA_0_N2_VV_P*DELTA_1_H2_VV*F_2102_H2_N2*(1.5-0.5*F_2102_H2_N2)
  F_3112_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V3_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_3112_H2_N2 = Q_VV_H2_N2_P*DELTA_1_N2_VV_P*DELTA_1_H2_VV*F_3112_H2_N2*(1.5-0.5*F_3112_H2_N2)
  F_4122_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V4_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_4122_H2_N2 = Q_VV_H2_N2_P*DELTA_2_N2_VV_P*DELTA_1_H2_VV*F_4122_H2_N2*(1.5-0.5*F_4122_H2_N2)
  F_5132_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V5_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_5132_H2_N2 = Q_VV_H2_N2_P*DELTA_3_N2_VV_P*DELTA_1_H2_VV*F_5132_H2_N2*(1.5-0.5*F_5132_H2_N2)
  F_6142_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V6_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_6142_H2_N2 = Q_VV_H2_N2_P*DELTA_4_N2_VV_P*DELTA_1_H2_VV*F_6142_H2_N2*(1.5-0.5*F_6142_H2_N2)
  F_7152_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V7_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_7152_H2_N2 = Q_VV_H2_N2_P*DELTA_5_N2_VV_P*DELTA_1_H2_VV*F_7152_H2_N2*(1.5-0.5*F_7152_H2_N2)
  F_8162_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V2_H2-E_V1_H2)-(E_V7_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_8162_H2_N2 = Q_VV_H2_N2_P*DELTA_6_N2_VV_P*DELTA_1_H2_VV*F_8162_H2_N2*(1.5-0.5*F_8162_H2_N2)
  F_2203_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V2_N2-0.))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_2203_H2_N2 = Q_VV_H2_N2_P*DELTA_0_N2_VV_P*DELTA_2_H2_VV*F_2203_H2_N2*(1.5-0.5*F_2203_H2_N2)
  F_3213_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V3_N2-E_V1_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_3213_H2_N2 = Q_VV_H2_N2_P*DELTA_1_N2_VV_P*DELTA_2_H2_VV*F_3213_H2_N2*(1.5-0.5*F_3213_H2_N2)
  F_4223_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V4_N2-E_V2_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_4223_H2_N2 = Q_VV_H2_N2_P*DELTA_2_N2_VV_P*DELTA_2_H2_VV*F_4223_H2_N2*(1.5-0.5*F_4223_H2_N2)
  F_5233_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V5_N2-E_V3_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_5233_H2_N2 = Q_VV_H2_N2_P*DELTA_3_N2_VV_P*DELTA_2_H2_VV*F_5233_H2_N2*(1.5-0.5*F_5233_H2_N2)
  F_6243_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V6_N2-E_V4_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_6243_H2_N2 = Q_VV_H2_N2_P*DELTA_4_N2_VV_P*DELTA_2_H2_VV*F_6243_H2_N2*(1.5-0.5*F_6243_H2_N2)
  F_7253_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V7_N2-E_V5_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_7253_H2_N2 = Q_VV_H2_N2_P*DELTA_5_N2_VV_P*DELTA_2_H2_VV*F_7253_H2_N2*(1.5-0.5*F_7253_H2_N2)
  F_8263_H2_N2 = EXP((-DELTA_VV_H2_N2_P*ABS((E_V3_H2-E_V2_H2)-(E_V8_N2-E_V6_N2))/(1.38064852D-23*6.242D18))/ABS(E_20_N2-E_10_H2))
  VV_8263_H2_N2 = Q_VV_H2_N2_P*DELTA_6_N2_VV_P*DELTA_2_H2_VV*F_8263_H2_N2*(1.5-0.5*F_8263_H2_N2)
  F_2231_H2_H2 = EXP((-DELTA_VV_0*((E_V2_H2-E_V1_H2)-(E_V3_H2-E_V2_H2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_H2))
  VV_2231_H2_H2 = Q_01_VV*DELTA_2_H2_VV*DELTA_1_H2_VV*F_2231_H2_H2*(1.5-0.5*F_2231_H2_H2)
  F_2130_H2_H2 = EXP((-DELTA_VV_0*((E_V1_H2-0.)-(E_V3_H2-E_V2_H2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_H2))
  VV_2130_H2_H2 = Q_01_VV*DELTA_2_H2_VV*DELTA_0_H2_VV*F_2130_H2_H2*(1.5-0.5*F_2130_H2_H2)
  F_1120_H2_H2 = EXP((-DELTA_VV_0*((E_V1_H2-0.)-(E_V2_H2-E_V1_H2))/(1.38064852D-23*6.242D18))/(2*DELTA_E_H2))
  VV_1120_H2_H2 = Q_01_VV*DELTA_1_H2_VV*DELTA_0_H2_VV*F_1120_H2_H2*(1.5-0.5*F_1120_H2_H2)
  GAMMA_D = DIFF_L**2/(DIFF_COEFF_400*((TGAS/300)**(3./2.)))
  N2_THERMAL_VEL = SQRT(8*1.38064852D-23*TGAS/(N2_MASS_KG*MYPI))*100
  N2_WALL_SECOND_PART_V = (THE_V/(THE_AREA*ROUGHNESS))*(2./N2_THERMAL_VEL)*(2.-N2_WALL_V_GAMMA)/N2_WALL_V_GAMMA
  N2_WALL_SECOND_PART_E = (THE_V/(THE_AREA*ROUGHNESS))*(2./N2_THERMAL_VEL)*(2.-N2_WALL_E_GAMMA)/N2_WALL_E_GAMMA
  H2_THERMAL_VEL = SQRT(8*1.38064852D-23*TGAS/(H2_MASS_KG*MYPI))*100
  H2_WALL_SECOND_PART_V = (THE_V/(THE_AREA*ROUGHNESS))*(2./H2_THERMAL_VEL)*(2.-H2_WALL_V_GAMMA)/H2_WALL_V_GAMMA
  H2_WALL_SECOND_PART_E = (THE_V/(THE_AREA*ROUGHNESS))*(2./H2_THERMAL_VEL)*(2.-H2_WALL_E_GAMMA)/H2_WALL_E_GAMMA
  STVOL = SITE_DENSITY / V_OVER_A
  DSURF = 7.9D-1*(TGAS/300.D0)**1.5D0
  EZ = KB*TGAS*EV_PER_J
  KADS_N = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N_THERMAL_VEL )*(2.D0-1.D0)/1.D0)/STVOL
  KADS_H = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H_THERMAL_VEL )*(2.D0-1.D0)/1.D0)/STVOL
  KADS_NH = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/NH_THERMAL_VEL)*(2.D0-1.D0)/1.D0)/STVOL
  KADS_NH2 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/NH2_THERMAL_VEL)*(2.D0-1.D0)/1.D0)/STVOL
  KER_N = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N_THERMAL_VEL )*(2.D0-6.D-3)/6.D-3)/STVOL
  KER_H = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H_THERMAL_VEL )*(2.D0-1.5D-3)/1.5D-3)/STVOL
  KER_N_HS = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N_THERMAL_VEL )*(2.D0-1.D-2)/1.D-2)/STVOL
  KER_NH_HS = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/NH_THERMAL_VEL)*(2.D0-1.D-2)/1.D-2)/STVOL
  KER_NH2_HS = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/NH2_THERMAL_VEL)*(2.D0-1.D-2)/1.D-2)/STVOL
  KER_H_NS = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H_THERMAL_VEL )*(2.D0-8.D-3)/8.D-3)/STVOL
  KER_H_NHS = KER_H_NS
  KER_H_NH2S = KER_H_NS
  KER_H2_NHS = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-8.D-4)/8.D-4)/STVOL
  KLH_S14 = JUMP_FREQUENCY/(4.D0*STVOL)*EXP(-(METAL_ED+1.099D0)/(KB*TGAS*EV_PER_J))
  KLH_S15 = JUMP_FREQUENCY/(4.D0*STVOL)*EXP(-(METAL_ED+0.3D0)/(KB*TGAS*EV_PER_J))
  KLH_S16 = JUMP_FREQUENCY/(4.D0*STVOL)*EXP(-(METAL_ED+0.2D0)/(KB*TGAS*EV_PER_J))
  GAMMA_N2_0 = 10.D0**(-6.616D0+6.258D0*(1.D0-EXP(-1.821D0*EZ))*(1.D0-EXP( 0.606D0*EZ)))
  GAMMA_N2_1 = 10.D0**(-5.768D0+4.581D0*(1.D0-EXP(-1.237D0*EZ))*(1.D0-EXP(-4.462D0*EZ)))
  GAMMA_N2_2 = 10.D0**(-4.894D0+9.073D0*(1.D0-EXP(-0.411D0*EZ))*(1.D0-EXP(-6.115D0*EZ)))
  GAMMA_N2_3 = 10.D0**(-3.939D0+2.923D0*(1.D0-EXP(-2.120D0*EZ))*(1.D0-EXP(-0.627D0*EZ)))
  GAMMA_N2_4 = 10.D0**(-3.387D0+1.511D0*(1.D0-EXP(-6.768D0*EZ))*(1.D0-EXP(-1.444D0*EZ)))
  GAMMA_N2_5 = 10.D0**(-2.736D0+1.751D0*(1.D0-EXP(-9.403D0*EZ))*(1.D0-EXP( 0.248D0*EZ)))
  GAMMA_N2_6 = 10.D0**(-2.503D0+1.755D0*(1.D0-EXP(-7.913D0*EZ))*(1.D0-EXP( 0.526D0*EZ)))
  GAMMA_N2_7 = 10.D0**(-2.265D0+1.505D0*(1.D0-EXP(-8.096D0*EZ))*(1.D0-EXP( 0.430D0*EZ)))
  GAMMA_N2_8 = 10.D0**(-2.040D0+1.376D0*(1.D0-EXP(-7.239D0*EZ))*(1.D0-EXP( 0.489D0*EZ)))
  KDADS_N2_0 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_0)/GAMMA_N2_0)/(STVOL*STVOL)
  KDADS_N2_1 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_1)/GAMMA_N2_1)/(STVOL*STVOL)
  KDADS_N2_2 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_2)/GAMMA_N2_2)/(STVOL*STVOL)
  KDADS_N2_3 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_3)/GAMMA_N2_3)/(STVOL*STVOL)
  KDADS_N2_4 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_4)/GAMMA_N2_4)/(STVOL*STVOL)
  KDADS_N2_5 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_5)/GAMMA_N2_5)/(STVOL*STVOL)
  KDADS_N2_6 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_6)/GAMMA_N2_6)/(STVOL*STVOL)
  KDADS_N2_7 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_7)/GAMMA_N2_7)/(STVOL*STVOL)
  KDADS_N2_8 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-GAMMA_N2_8)/GAMMA_N2_8)/(STVOL*STVOL)
  KDADS_N2_E = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/N2_THERMAL_VEL)*(2.D0-1.D-1)/1.D-1)/(STVOL*STVOL)
  KDADS_H2_0 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-1.D-3)/1.D-3)/(STVOL*STVOL)
  KDADS_H2_1 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-1.D-2)/1.D-2)/(STVOL*STVOL)
  KDADS_H2_2 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-5.D-2)/5.D-2)/(STVOL*STVOL)
  KDADS_H2_3 = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-1.D-1)/1.D-1)/(STVOL*STVOL)
  KDADS_H2_E = 1.D0/(LAMBDA**2/DSURF + V_OVER_A*(2.D0/H2_THERMAL_VEL)*(2.D0-1.D-2)/1.D-2)/(STVOL*STVOL)
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = 1.26D-10*(TE/11600)-1.72D-10*(TE/11600)**2+6.51D-11*(TE/11600)**3-5.75D-12*(TE/11600)**4+1.71D-13*(TE/11600)**5
  rrt(015) = 6.50D-9*((TE/11600)**0.49)*EXP(-12.89/(TE/11600))
  rrt(016) = 1.53D-10*(TE/11600)-2.24D-10*(TE/11600)**2+9.37D-11*(TE/11600)**3-9.97D-12*(TE/11600)**4+3.33D-13*(TE/11600)**5
  rrt(017) = bolsig_rates(bolsig_pointer(14))
  rrt(018) = bolsig_rates(bolsig_pointer(15))
  rrt(019) = bolsig_rates(bolsig_pointer(16))
  rrt(020) = bolsig_rates(bolsig_pointer(17))
  rrt(021) = bolsig_rates(bolsig_pointer(18))
  rrt(022) = bolsig_rates(bolsig_pointer(19))
  rrt(023) = bolsig_rates(bolsig_pointer(20))
  rrt(024) = bolsig_rates(bolsig_pointer(21))
  rrt(025) = bolsig_rates(bolsig_pointer(22))
  rrt(026) = bolsig_rates(bolsig_pointer(23))
  rrt(027) = bolsig_rates(bolsig_pointer(24))
  rrt(028) = bolsig_rates(bolsig_pointer(25))
  rrt(029) = bolsig_rates(bolsig_pointer(26))
  rrt(030) = bolsig_rates(bolsig_pointer(27))
  rrt(031) = bolsig_rates(bolsig_pointer(28))
  rrt(032) = bolsig_rates(bolsig_pointer(29))
  rrt(033) = bolsig_rates(bolsig_pointer(30))
  rrt(034) = bolsig_rates(bolsig_pointer(31))
  rrt(035) = bolsig_rates(bolsig_pointer(32))
  rrt(036) = bolsig_rates(bolsig_pointer(33))
  rrt(037) = bolsig_rates(bolsig_pointer(34))
  rrt(038) = bolsig_rates(bolsig_pointer(35))
  rrt(039) = bolsig_rates(bolsig_pointer(36))
  rrt(040) = bolsig_rates(bolsig_pointer(37))
  rrt(041) = bolsig_rates(bolsig_pointer(38))
  rrt(042) = bolsig_rates(bolsig_pointer(39))
  rrt(043) = K_10_VT_N2_N2*G_10_N2_N2
  rrt(044) = K_10_VT_N2_N2*G_10_N2_N2*EXP(-(E_V1_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(045) = K_10_VT_N2_N2*G_21_N2_N2
  rrt(046) = K_10_VT_N2_N2*G_21_N2_N2*EXP(-(E_V2_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(047) = K_10_VT_N2_N2*G_32_N2_N2
  rrt(048) = K_10_VT_N2_N2*G_32_N2_N2*EXP(-(E_V3_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(049) = K_10_VT_N2_N2*G_43_N2_N2
  rrt(050) = K_10_VT_N2_N2*G_43_N2_N2*EXP(-(E_V4_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(051) = K_10_VT_N2_N2*G_54_N2_N2
  rrt(052) = K_10_VT_N2_N2*G_54_N2_N2*EXP(-(E_V5_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(053) = K_10_VT_N2_N2*G_65_N2_N2
  rrt(054) = K_10_VT_N2_N2*G_65_N2_N2*EXP(-(E_V6_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(055) = K_10_VT_N2_N2*G_76_N2_N2
  rrt(056) = K_10_VT_N2_N2*G_76_N2_N2*EXP(-(E_V7_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(057) = K_10_VT_N2_N2*G_87_N2_N2
  rrt(058) = K_10_VT_N2_N2*G_87_N2_N2*EXP(-(E_V8_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(059) = K_10_VT_N2_H2*G_10_N2_H2
  rrt(060) = K_10_VT_N2_H2*G_10_N2_H2*EXP(-(E_V1_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(061) = K_10_VT_N2_H2*G_21_N2_H2
  rrt(062) = K_10_VT_N2_H2*G_21_N2_H2*EXP(-(E_V2_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(063) = K_10_VT_N2_H2*G_32_N2_H2
  rrt(064) = K_10_VT_N2_H2*G_32_N2_H2*EXP(-(E_V3_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(065) = K_10_VT_N2_H2*G_43_N2_H2
  rrt(066) = K_10_VT_N2_H2*G_43_N2_H2*EXP(-(E_V4_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(067) = K_10_VT_N2_H2*G_54_N2_H2
  rrt(068) = K_10_VT_N2_H2*G_54_N2_H2*EXP(-(E_V5_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(069) = K_10_VT_N2_H2*G_65_N2_H2
  rrt(070) = K_10_VT_N2_H2*G_65_N2_H2*EXP(-(E_V6_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(071) = K_10_VT_N2_H2*G_76_N2_H2
  rrt(072) = K_10_VT_N2_H2*G_76_N2_H2*EXP(-(E_V7_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(073) = K_10_VT_N2_H2*G_87_N2_H2
  rrt(074) = K_10_VT_N2_H2*G_87_N2_H2*EXP(-(E_V8_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(075) = K_10_N2_N
  rrt(076) = K_10_N2_N*EXP(-(E_V1_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(077) = K_21_N2_N
  rrt(078) = K_21_N2_N*EXP(-(E_V2_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(079) = K_32_N2_N
  rrt(080) = K_32_N2_N*EXP(-(E_V3_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(081) = K_43_N2_N
  rrt(082) = K_43_N2_N*EXP(-(E_V4_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(083) = K_54_N2_N
  rrt(084) = K_54_N2_N*EXP(-(E_V5_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(085) = K_65_N2_N
  rrt(086) = K_65_N2_N*EXP(-(E_V6_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(087) = K_76_N2_N
  rrt(088) = K_76_N2_N*EXP(-(E_V7_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(089) = K_87_N2_N
  rrt(090) = K_87_N2_N*EXP(-(E_V8_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(091) = K_10_N2_H
  rrt(092) = 1.9D-21*EXP(-(E_V1_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(093) = K_21_N2_H
  rrt(094) = 4.9D-21*EXP(-(E_V2_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(095) = K_32_N2_H
  rrt(096) = 1.3D-20*EXP(-(E_V3_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(097) = K_43_N2_H
  rrt(098) = 3.3D-20*EXP(-(E_V4_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(099) = K_54_N2_H
  rrt(100) = 8.3D-20*EXP(-(E_V5_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(101) = K_65_N2_H
  rrt(102) = 2.1D-19*EXP(-(E_V6_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(103) = K_76_N2_H
  rrt(104) = 5.3D-19*EXP(-(E_V7_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(105) = K_87_N2_H
  rrt(106) = 1.3D-18*EXP(-(E_V8_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(107) = K_10_VT_H2_H2*G_10_H2_H2
  rrt(108) = K_10_VT_H2_H2*G_10_H2_H2*EXP(-(E_V1_H2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(109) = K_10_VT_H2_H2*G_21_H2_H2
  rrt(110) = K_10_VT_H2_H2*G_21_H2_H2*EXP(-(E_V2_H2-E_V1_H2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(111) = K_10_VT_H2_H2*G_32_H2_H2
  rrt(112) = K_10_VT_H2_H2*G_32_H2_H2*EXP(-(E_V3_H2-E_V2_H2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(113) = K_10_H2_H
  rrt(114) = K_10_H2_H*EXP(-(E_V1_H2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(115) = K_21_H2_H
  rrt(116) = K_21_H2_H*EXP(-(E_V2_H2-E_V1_H2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(117) = K_32_H2_H
  rrt(118) = K_32_H2_H*EXP(-(E_V3_H2-E_V2_H2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(119) = VV_1102_N2_N2
  rrt(120) = VV_1102_N2_N2*EXP(-(E_V1_N2-0.+E_V1_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(121) = VV_1203_N2_N2
  rrt(122) = VV_1203_N2_N2*EXP(-(E_V1_N2-0.+E_V2_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(123) = VV_1304_N2_N2
  rrt(124) = VV_1304_N2_N2*EXP(-(E_V1_N2-0.+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(125) = VV_1405_N2_N2
  rrt(126) = VV_1405_N2_N2*EXP(-(E_V1_N2-0.+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(127) = VV_1506_N2_N2
  rrt(128) = VV_1506_N2_N2*EXP(-(E_V1_N2-0.+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(129) = VV_1607_N2_N2
  rrt(130) = VV_1607_N2_N2*EXP(-(E_V1_N2-0.+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(131) = VV_1708_N2_N2
  rrt(132) = VV_1708_N2_N2*EXP(-(E_V1_N2-0.+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(133) = VV_2213_N2_N2
  rrt(134) = VV_2213_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V2_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(135) = VV_2314_N2_N2
  rrt(136) = VV_2314_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(137) = VV_2415_N2_N2
  rrt(138) = VV_2415_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(139) = VV_2516_N2_N2
  rrt(140) = VV_2516_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(141) = VV_2617_N2_N2
  rrt(142) = VV_2617_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(143) = VV_2718_N2_N2
  rrt(144) = VV_2718_N2_N2*EXP(-(E_V2_N2-E_V1_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(145) = VV_3324_N2_N2
  rrt(146) = VV_3324_N2_N2*EXP(-(E_V3_N2-E_V2_N2+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(147) = VV_3425_N2_N2
  rrt(148) = VV_3425_N2_N2*EXP(-(E_V3_N2-E_V2_N2+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(149) = VV_3526_N2_N2
  rrt(150) = VV_3526_N2_N2*EXP(-(E_V3_N2-E_V2_N2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(151) = VV_3627_N2_N2
  rrt(152) = VV_3627_N2_N2*EXP(-(E_V3_N2-E_V2_N2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(153) = VV_3728_N2_N2
  rrt(154) = VV_3728_N2_N2*EXP(-(E_V3_N2-E_V2_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(155) = VV_4435_N2_N2
  rrt(156) = VV_4435_N2_N2*EXP(-(E_V4_N2-E_V3_N2+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(157) = VV_4536_N2_N2
  rrt(158) = VV_4536_N2_N2*EXP(-(E_V4_N2-E_V3_N2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(159) = VV_4637_N2_N2
  rrt(160) = VV_4637_N2_N2*EXP(-(E_V4_N2-E_V3_N2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(161) = VV_4738_N2_N2
  rrt(162) = VV_4738_N2_N2*EXP(-(E_V4_N2-E_V3_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(163) = VV_5546_N2_N2
  rrt(164) = VV_5546_N2_N2*EXP(-(E_V5_N2-E_V4_N2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(165) = VV_5647_N2_N2
  rrt(166) = VV_5647_N2_N2*EXP(-(E_V5_N2-E_V4_N2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(167) = VV_5748_N2_N2
  rrt(168) = VV_5748_N2_N2*EXP(-(E_V5_N2-E_V4_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(169) = VV_6657_N2_N2
  rrt(170) = VV_6657_N2_N2*EXP(-(E_V6_N2-E_V5_N2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(171) = VV_6758_N2_N2
  rrt(172) = VV_6758_N2_N2*EXP(-(E_V6_N2-E_V5_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(173) = VV_7768_N2_N2
  rrt(174) = VV_7768_N2_N2*EXP(-(E_V7_N2-E_V6_N2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(175) = VV_1001_H2_N2
  rrt(176) = VV_1001_H2_N2*EXP(-(E_V1_H2-0.+0.-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(177) = VV_1102_H2_N2
  rrt(178) = VV_1102_H2_N2*EXP(-(E_V1_H2-0.+E_V1_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(179) = VV_1203_H2_N2
  rrt(180) = VV_1203_H2_N2*EXP(-(E_V1_H2-0.+E_V2_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(181) = VV_1304_H2_N2
  rrt(182) = VV_1304_H2_N2*EXP(-(E_V1_H2-0.+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(183) = VV_1405_H2_N2
  rrt(184) = VV_1405_H2_N2*EXP(-(E_V1_H2-0.+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(185) = VV_1506_H2_N2
  rrt(186) = VV_1506_H2_N2*EXP(-(E_V1_H2-0.+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(187) = VV_1607_H2_N2
  rrt(188) = VV_1607_H2_N2*EXP(-(E_V1_H2-0.+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(189) = VV_1708_H2_N2
  rrt(190) = VV_1708_H2_N2*EXP(-(E_V1_H2-0.+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(191) = VV_2011_H2_N2
  rrt(192) = VV_2011_H2_N2*EXP(-(E_V2_H2-E_V1_H2+0.-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(193) = VV_2112_H2_N2
  rrt(194) = VV_2112_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V1_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(195) = VV_2213_H2_N2
  rrt(196) = VV_2213_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V2_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(197) = VV_2314_H2_N2
  rrt(198) = VV_2314_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(199) = VV_2415_H2_N2
  rrt(200) = VV_2415_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(201) = VV_2516_H2_N2
  rrt(202) = VV_2516_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(203) = VV_2617_H2_N2
  rrt(204) = VV_2617_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(205) = VV_2718_H2_N2
  rrt(206) = VV_2718_H2_N2*EXP(-(E_V2_H2-E_V1_H2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(207) = VV_3021_H2_N2
  rrt(208) = VV_3021_H2_N2*EXP(-(E_V3_H2-E_V2_H2+0.-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(209) = VV_3122_H2_N2
  rrt(210) = VV_3122_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V1_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(211) = VV_3223_H2_N2
  rrt(212) = VV_3223_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V2_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(213) = VV_3324_H2_N2
  rrt(214) = VV_3324_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V3_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(215) = VV_3425_H2_N2
  rrt(216) = VV_3425_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V4_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(217) = VV_3526_H2_N2
  rrt(218) = VV_3526_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V5_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(219) = VV_3627_H2_N2
  rrt(220) = VV_3627_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V6_N2-E_V7_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(221) = VV_3728_H2_N2
  rrt(222) = VV_3728_H2_N2*EXP(-(E_V3_H2-E_V2_H2+E_V7_N2-E_V8_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(223) = VV_2001_H2_N2
  rrt(224) = VV_2001_H2_N2*EXP(-(0.-E_V1_H2+E_V2_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(225) = VV_3011_H2_N2
  rrt(226) = VV_3011_H2_N2*EXP(-(0.-E_V1_H2+E_V3_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(227) = VV_4021_H2_N2
  rrt(228) = VV_4021_H2_N2*EXP(-(0.-E_V1_H2+E_V4_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(229) = VV_5031_H2_N2
  rrt(230) = VV_5031_H2_N2*EXP(-(0.-E_V1_H2+E_V5_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(231) = VV_6041_H2_N2
  rrt(232) = VV_6041_H2_N2*EXP(-(0.-E_V1_H2+E_V6_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(233) = VV_7051_H2_N2
  rrt(234) = VV_7051_H2_N2*EXP(-(0.-E_V1_H2+E_V7_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(235) = VV_8061_H2_N2
  rrt(236) = VV_8061_H2_N2*EXP(-(0.-E_V1_H2+E_V8_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(237) = VV_2102_H2_N2
  rrt(238) = VV_2102_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V2_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(239) = VV_3112_H2_N2
  rrt(240) = VV_3112_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V3_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(241) = VV_4122_H2_N2
  rrt(242) = VV_4122_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V4_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(243) = VV_5132_H2_N2
  rrt(244) = VV_5132_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V5_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(245) = VV_6142_H2_N2
  rrt(246) = VV_6142_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V6_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(247) = VV_7152_H2_N2
  rrt(248) = VV_7152_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V7_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(249) = VV_8162_H2_N2
  rrt(250) = VV_8162_H2_N2*EXP(-(E_V1_H2-E_V2_H2+E_V8_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(251) = VV_2203_H2_N2
  rrt(252) = VV_2203_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V2_N2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(253) = VV_3213_H2_N2
  rrt(254) = VV_3213_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V3_N2-E_V1_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(255) = VV_4223_H2_N2
  rrt(256) = VV_4223_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V4_N2-E_V2_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(257) = VV_5233_H2_N2
  rrt(258) = VV_5233_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V5_N2-E_V3_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(259) = VV_6243_H2_N2
  rrt(260) = VV_6243_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V6_N2-E_V4_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(261) = VV_7253_H2_N2
  rrt(262) = VV_7253_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V7_N2-E_V5_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(263) = VV_8263_H2_N2
  rrt(264) = VV_8263_H2_N2*EXP(-(E_V2_H2-E_V3_H2+E_V8_N2-E_V6_N2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(265) = VV_2231_H2_H2
  rrt(266) = VV_2231_H2_H2*EXP(-(E_V2_H2-E_V3_H2+E_V2_H2-E_V1_H2)/(1.38064852D-23*TGAS*6.242D18))
  rrt(267) = F_2130_H2_H2
  rrt(268) = F_2130_H2_H2*EXP(-(E_V2_H2-E_V3_H2+E_V1_H2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(269) = VV_1120_H2_H2
  rrt(270) = VV_1120_H2_H2*EXP(-(E_V1_H2-E_V2_H2+E_V1_H2-0.)/(1.38064852D-23*TGAS*6.242D18))
  rrt(271) = 3.0D-11
  rrt(272) = 3.0D-11
  rrt(273) = 3.0D-11
  rrt(274) = 3.0D-11
  rrt(275) = 3.0D-11
  rrt(276) = 3.0D-11
  rrt(277) = 3.0D-11
  rrt(278) = 3.0D-11
  rrt(279) = 3.0D-11
  rrt(280) = 3.0D-11
  rrt(281) = 3.0D-11
  rrt(282) = 3.0D-11
  rrt(283) = 1.34D5
  rrt(284) = 1.0D2
  rrt(285) = 2.45D7
  rrt(286) = 1.75D-7*((TE/11600)**(-1.24))*(EXP(-12.59/(TE/11600)))
  rrt(287) = bolsig_rates(bolsig_pointer(40))
  rrt(288) = 5.0D-8*((TE/11600)**0.5)*(EXP(-8.6/(TE/11600)))
  rrt(289) = 5.0D-8*((TE/11600)**0.5)*(EXP(-7.6/(TE/11600)))
  rrt(290) = rrt(289)
  rrt(291) = 5.0D-8*((TE/11600)**0.5)*(EXP(-4.4/(TE/11600)))
  rrt(292) = 5.0D-8*((TE/11600)**0.5)*(EXP(-5.5/(TE/11600)))
  rrt(293) = 9.0D-8*((300.0D0/TE)**0.39)
  rrt(294) = 8.1D-8*((300.0D0/TE)**0.39)
  rrt(295) = 9.0D-9*((300.0D0/TE)**0.39)
  rrt(296) = 2.0D-7*((300.0D0/TE)**0.5)
  rrt(297) = 2.3D-6*((300.0D0/TE)**0.53)
  rrt(298) = 3.00D-8*((TE/11600)**0.44)*EXP(-37.72/(TE/11600))
  rrt(299) = 7.51D-9-1.12D-9*(TE/11600)+1.03D-10*((TE/11600)**2)-4.15D-12*((TE/11600)**3)+5.86D-14*((TE/11600)**4)
  rrt(300) = &
  0.5*8.39D-9+3.02D-9*(TE/11600)-3.80D-10*((TE/11600)**2)+1.31D-11*((TE/11600)**3)+2.42D-13*((TE/11600)**4)-2.30D-14*((TE/11600)**5)+3.55D-16*((TE/11600)**6)
  rrt(301) = rrt(300)
  rrt(302) = 4.30D-8*((0.026/(TE/11600))**0.5)
  rrt(303) = 1.02D-7*((0.026/(TE/11600))**0.4)
  rrt(304) = 1.98D-7*((0.026/(TE/11600))**0.4)
  rrt(305) = 1.55D-7*((0.026/(TE/11600))**0.5)
  rrt(306) = rrt(305)
  rrt(307) = 8.01D-7*((0.026/(TE/11600))**0.605)
  rrt(308) = 1.23D-7*((0.026/(TE/11600))**0.605)
  rrt(309) = 7.1D-7*((0.026/(TE/11600))**0.72)
  rrt(310) = 2.00D-9
  rrt(311) = 3.0D-10
  rrt(312) = 7.2D-13*(TEFFN2N/300)**2.2
  rrt(313) = 9.0D-30*EXP(400.0D0/TEFFN2N)
  rrt(314) = 5.2D-29*(300.0D0/TEFFN2N2)**2.2
  rrt(315) = 1.95D-9
  rrt(316) = 6.6D-11
  rrt(317) = MIN(2.1D-16*EXP(TEFFN4N2/121.0D0),5.0D-10)
  rrt(318) = 1.0D-11
  rrt(319) = 5.00D-10
  rrt(320) = 0.20*2.35D-9
  rrt(321) = 0.71*2.35D-9
  rrt(322) = 0.09*2.35D-9
  rrt(323) = 6.4D-10
  rrt(324) = 2.00D-9
  rrt(325) = 5.70D-9
  rrt(326) = 2.00D-9
  rrt(327) = 5.20D-9
  rrt(328) = 0.15*1.23D-9
  rrt(329) = 0.85*1.23D-9
  rrt(330) = 0.75*2.40D-9
  rrt(331) = 0.25*2.40D-9
  rrt(332) = 6.50D-10
  rrt(333) = 1.95D-10
  rrt(334) = 0.5*2.30D-9
  rrt(335) = rrt(334)
  rrt(336) = 2.10D-9
  rrt(337) = 2.0D-12
  rrt(338) = 4.0D-11*(300.0D0/TGAS)**0.667
  rrt(339) = 3.0D-16
  rrt(340) = 3.0D-10
  rrt(341) = 1.5D-10
  rrt(342) = 3.0D-11
  rrt(343) = 3.0D-11
  rrt(344) = 2.0D-12
  rrt(345) = 1.9D-13
  rrt(346) = 1.0D-11
  rrt(347) = 1.0D-11
  rrt(348) = 4.0D-12
  rrt(349) = 2.3D-14*EXP(-510.0D0/TGAS)
  rrt(350) = 1.8D-12
  rrt(351) = 6.0D-13
  rrt(352) = 6.0D-14
  rrt(353) = 5D-11
  rrt(354) = 2D-10*EXP(-3500/TGAS)
  rrt(355) = 1.6D-10
  rrt(356) = 2.5D-11
  rrt(357) = 1.5D-11
  rrt(358) = 2.6D-11
  rrt(359) = 4.0D-10*(TGAS/300)**0.5*EXP(-16600/TGAS+0.3*E_V1_H2/(1.38064852D-23*6.242D18*TGAS))
  rrt(360) = 4.0D-10*(TGAS/300)**0.5*EXP(-16600/TGAS+0.3*E_V2_H2/(1.38064852D-23*6.242D18*TGAS))
  rrt(361) = 4.0D-10*(TGAS/300)**0.5*EXP(-16600/TGAS+0.3*E_V3_H2/(1.38064852D-23*6.242D18*TGAS))
  rrt(362) = 4.0D-10*(TGAS/300)**0.5
  rrt(363) = rrt(362)
  rrt(364) = rrt(362)
  rrt(365) = rrt(362)
  rrt(366) = rrt(362)
  rrt(367) = 2.3D-12
  rrt(368) = 1.1D-10
  rrt(369) = 2.5D-14
  rrt(370) = 5D-11
  rrt(371) = 5.4D-11*EXP(-165/TGAS)
  rrt(372) = 5D-14*(TGAS/300)
  rrt(373) = 1.7D-12*(TGAS/300)**1.5
  rrt(374) = 8.5D-11
  rrt(375) = 6.6D-11*EXP(-1840/TGAS)
  rrt(376) = 1.2D-10
  rrt(377) = 1.2D-10
  rrt(378) = 1.66D-12
  rrt(379) = 5.4D-11*EXP(-6492/TGAS)
  rrt(380) = 8.4D-14*(TGAS/300)**4.1*EXP(-4760/TGAS)
  rrt(381) = 4.4D-36
  rrt(382) = 4.4D-36
  rrt(383) = 4.4D-36
  rrt(384) = 4.4D-36
  rrt(385) = 4.4D-36
  rrt(386) = 4.4D-36
  rrt(387) = 4.4D-36
  rrt(388) = 4.4D-36
  rrt(389) = 4.4D-36
  rrt(390) = 4.4D-36
  rrt(391) = 4.4D-36
  rrt(392) = 4.4D-36
  rrt(393) = 4.4D-36
  rrt(394) = 2.6D-35
  rrt(395) = 2.6D-35
  rrt(396) = 6.2D-36
  rrt(397) = 6.2D-36
  rrt(398) = 6.2D-36
  rrt(399) = 6.2D-36
  rrt(400) = 6.2D-36
  rrt(401) = 6.2D-36
  rrt(402) = 6.2D-36
  rrt(403) = 6.2D-36
  rrt(404) = 6.2D-36
  rrt(405) = 6.2D-36
  rrt(406) = 6.2D-36
  rrt(407) = 6.2D-36
  rrt(408) = 6.2D-36
  rrt(409) = 3.6D-35
  rrt(410) = 3.6D-35
  rrt(411) = 2.4D-36*EXP(500/TGAS)
  rrt(412) = rrt(411)
  rrt(413) = 2.3D-35*(300/TGAS)**0.6
  rrt(414) = 2.2D-35*(300/TGAS)
  rrt(415) = 2.6D-36
  rrt(416) = 2.6D-36
  rrt(417) = 2.6D-37
  rrt(418) = 2.6D-37
  rrt(419) = 2.6D-35
  rrt(420) = 2.6D-35
  rrt(421) = 1.43D-32
  rrt(422) = 1.43D-32
  rrt(423) = 6.5D-38*(TGAS/300)*EXP(1700/TGAS)
  rrt(424) = rrt(423)
  rrt(425) = 1./(GAMMA_D+N2_WALL_SECOND_PART_E)
  rrt(426) = rrt(425)
  rrt(427) = 1./(GAMMA_D+N2_WALL_SECOND_PART_V)
  rrt(428) = rrt(427)
  rrt(429) = rrt(427)
  rrt(430) = rrt(427)
  rrt(431) = rrt(427)
  rrt(432) = rrt(427)
  rrt(433) = rrt(427)
  rrt(434) = rrt(427)
  rrt(435) = 1./(GAMMA_D+H2_WALL_SECOND_PART_E)
  rrt(436) = rrt(435)
  rrt(437) = rrt(435)
  rrt(438) = rrt(435)
  rrt(439) = 1./(GAMMA_D+H2_WALL_SECOND_PART_V)
  rrt(440) = rrt(439)
  rrt(441) = rrt(439)
  rrt(442) = bolsig_rates(bolsig_pointer(41))
  rrt(443) = 2D-7*(300/TGAS)
  rrt(444) = rrt(443)
  rrt(445) = rrt(443)
  rrt(446) = rrt(443)
  rrt(447) = rrt(443)
  rrt(448) = 2D-25*(300/TGAS)**2.5
  rrt(449) = rrt(448)
  rrt(450) = rrt(448)
  rrt(451) = rrt(448)
  rrt(452) = rrt(448)
  rrt(453) = rrt(448)
  rrt(454) = rrt(448)
  rrt(455) = rrt(448)
  rrt(456) = rrt(448)
  rrt(457) = rrt(448)
  rrt(458) = rrt(448)
  rrt(459) = rrt(448)
  rrt(460) = rrt(448)
  rrt(461) = rrt(448)
  rrt(462) = rrt(448)
  rrt(463) = rrt(448)
  rrt(464) = rrt(448)
  rrt(465) = rrt(448)
  rrt(466) = rrt(448)
  rrt(467) = rrt(448)
  rrt(468) = KADS_N
  rrt(469) = KADS_H
  rrt(470) = KADS_NH
  rrt(471) = KADS_NH2
  rrt(472) = KER_N
  rrt(473) = KER_H
  rrt(474) = KER_N_HS
  rrt(475) = KER_NH_HS
  rrt(476) = KER_NH2_HS
  rrt(477) = KER_H_NS
  rrt(478) = KER_H_NHS
  rrt(479) = KER_H_NH2S
  rrt(480) = KER_H2_NHS
  rrt(481) = KLH_S14
  rrt(482) = KLH_S15
  rrt(483) = KLH_S16
  rrt(484) = KDADS_N2_0
  rrt(485) = KDADS_N2_1
  rrt(486) = KDADS_N2_2
  rrt(487) = KDADS_N2_3
  rrt(488) = KDADS_N2_4
  rrt(489) = KDADS_N2_5
  rrt(490) = KDADS_N2_6
  rrt(491) = KDADS_N2_7
  rrt(492) = KDADS_N2_8
  rrt(493) = KDADS_N2_E
  rrt(494) = KDADS_N2_E
  rrt(495) = KDADS_N2_E
  rrt(496) = KDADS_N2_E
  rrt(497) = KDADS_H2_0
  rrt(498) = KDADS_H2_1
  rrt(499) = KDADS_H2_2
  rrt(500) = KDADS_H2_3
  rrt(501) = KDADS_H2_E
  rrt(502) = KDADS_H2_E
  rrt(503) = KDADS_H2_E
  rrt(504) = KDADS_H2_E
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------

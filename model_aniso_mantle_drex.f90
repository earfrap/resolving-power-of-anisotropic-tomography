!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!
!
!  input is (r, theta, phi), output is the matrix cij(6x6)
!  0 <= r <= 1, 0 <= phi <= 80 , 70 <= theta <= 110
!
!  read cij from text file and interpolate them to the grid
!
!
!--------------------------------------------------------------------------------------------------

  module model_aniso_mantle_drex_par

  ! model_aniso_mantle_drex_variables
  double precision,dimension(:,:,:,:),allocatable :: AMM_V_Cij
  double precision,dimension(:,:),allocatable :: AMM_V_Cijp  
  double precision,dimension(:),allocatable :: AMM_V_pro,AMM_V_prop
  double precision, dimension(:), allocatable :: AMM_V_lon
  double precision, dimension(:), allocatable :: AMM_V_colat
  integer :: nx,ny,nz,nzp
 end module model_aniso_mantle_drex_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_aniso_mantle_drex_broadcast()

! standard routine to setup model

  use constants, only: myrank,IMAIN
  use model_aniso_mantle_drex_par

  implicit none

  ! local parameters
  integer :: ier


  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: drex_model (aniso_mantle from drex) model'
    call flush_IMAIN()
   endif


!for drex nx801 ny401 nz71
!for prem nx25 ny25 nz101

  ! allocates model arrays
  ! modify these values according to your needs
    allocate(AMM_V_Cij(22,801,401,71), & 
             AMM_V_pro(71),&
             AMM_V_Cijp(3,20),&
             AMM_V_prop(20),&
             AMM_V_lon(801),&
             AMM_V_colat(401), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating AMM_V arrays')

  ! the variables read are declared and stored in structure AMM_V
  if (myrank == 0) call read_aniso_mantle_model_drex()


  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(nx)
  call bcast_all_singlei(ny)
  call bcast_all_singlei(nz)
  call bcast_all_singlei(nzp)
  call bcast_all_dp(AMM_V_Cij,22*801*401*71) !22*nx*ny*nz
  call bcast_all_dp(AMM_V_pro,71)  !nz
  call bcast_all_dp(AMM_V_lon,801)  !nx
  call bcast_all_dp(AMM_V_colat,401)  !ny
  call bcast_all_dp(AMM_V_Cijp,3*20)  !3*nzp
  call bcast_all_dp(AMM_V_prop,20)  !nzp


   end subroutine model_aniso_mantle_drex_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_aniso_mantle_drex(r,theta,phi,&
                                rho,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use constants, only: PI,GRAV,RHOAV,DEGREES_TO_RADIANS,R_EARTH,R_EARTH_KM,R_UNIT_SPHERE,ZERO
  use shared_parameters, only: R670

  use model_aniso_mantle_drex_par

  implicit none

  double precision,intent(in) :: r,theta,phi
  double precision,intent(out) :: rho  !(inout)
  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  ! local parameters
  double precision :: depth,lon,colat,A,C,F,L,N,colat_const,lon_const
  double precision :: vp,vs
  double precision :: anispara(3,2),elpar(3)
  double precision :: dprof1,scale_Pa,scaleval
  integer :: idep,ipar,icz0,icz1


!inizialize values
    A = ZERO
    F = ZERO
    C = ZERO
    L = ZERO
    N = ZERO

 ! (change to CMB if desired)
 if (r <= R670/R_EARTH) then  ! if r < 0.8


! dimensionalize
 depth = R_EARTH_KM*(R_UNIT_SPHERE - r)
 if (depth < AMM_V_prop(nzp) .or. depth > AMM_V_prop(1)) print *, 'depth',depth
 if (depth < AMM_V_prop(nzp) .or. depth > AMM_V_prop(1)) print *, 'AMM_V_prop(nzp)',AMM_V_prop(nzp)
 if (depth < AMM_V_prop(nzp) .or. depth > AMM_V_prop(1)) print *, 'AMM_V_prop(1)',AMM_V_prop(1)
 if (depth < AMM_V_prop(nzp) .or. depth > AMM_V_prop(1)) call exit_MPI_without_rank('r out of range in lower mantle')

    icz0 = 0
    do idep = 1,nzp   
      if (AMM_V_prop(idep) > depth) icz0 = icz0 + 1
    enddo
    icz1 = icz0 + 1


    if (icz0 < 1 .or. icz0 > nzp) call exit_MPI_without_rank('icz0 out of range in lower mantle')
    if (icz1 < 1 .or. icz1 > nzp) call exit_MPI_without_rank('icz1 out of range in lower mantle')
    
    do ipar = 1,3 
      anispara(ipar,1) = AMM_V_Cijp(ipar,icz0)
      anispara(ipar,2) = AMM_V_Cijp(ipar,icz1)
    enddo

 
    dprof1 = (depth - AMM_V_prop(icz1))/(AMM_V_prop(icz0) - AMM_V_prop(icz1)) 

    do ipar = 1,3  
       elpar(ipar) = anispara(ipar,1)*dprof1 + anispara(ipar,2)*(1.0-dprof1)
    enddo

  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_Pa =(RHOAV)*((R_EARTH*scaleval)**2)

!inizialize values

    rho=ZERO
    vp=ZERO
    vs=ZERO

!!    c11 = ZERO
!!    c12 = ZERO
!!    c13 = ZERO
!!    c14 = ZERO
!!    c15 = ZERO
!!    c16 = ZERO
!!    c22 = ZERO
!!    c23 = ZERO
!!    c24 = ZERO
!!    c25 = ZERO
!!    c26 = ZERO
!!    c33 = ZERO
!!    c34 = ZERO
!!    c35 = ZERO
!!    c36 = ZERO
!!    c44 = ZERO
!!    c45 = ZERO
!!    c46 = ZERO
!!    c55 = ZERO
!!    c56 = ZERO
!!    c66 = ZERO


    rho = elpar(1)
    vp = elpar(2)
    vs = elpar(3)
    !print *,'rho kg/m3',rho,'vp m/s',vp,'vs m/s',vs
    ! c11 from vp,vs and rho in m/s and kg/m3 is in Pa 
    c11 = (rho*vp*vp)/scale_Pa
    c12 = (rho*(vp*vp-2.d0*vs*vs))/scale_Pa
    c13 = c12
    c14 = 0.d0
    c15 = 0.d0
    c16 = 0.d0
    c22 = c11
    c23 = c12
    c24 = 0.d0
    c25 = 0.d0
    c26 = 0.d0
    c33 = c11
    c34 = 0.d0
    c35 = 0.d0
    c36 = 0.d0
    c44 = (rho*vs*vs)/scale_Pa
    c45 = 0.d0
    c46 = 0.d0
    c55 = c44
    c56 = 0.d0
    c66 = c44

! non dimensionalized parameters
  rho = rho/RHOAV  ! rho from text file is kg/m3, we need it non-dimens
! vp from text file is in m/s, we need it non-dimens 
  vp = vp/(scaleval * R_EARTH) 
  vs = vs/(scaleval * R_EARTH)


 else
    ! above 670 discontinuity
    ! from ~24.4 to 670 km, 0.8<r<1
    ! converts colat/lon to degrees if necessary
    lon = phi / DEGREES_TO_RADIANS
    colat = theta / DEGREES_TO_RADIANS

    ! it reads the model parameters from text file and interpolate them to the grid
    call build_cij_drex(AMM_V_pro,rho,r,colat,lon, &
                   c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66)

        ! DONT USE IT IF YOU USE ROTATION (RAD to GLOB) IN MESHFEM3D_MODELS.f90
        ! 24/08/2020  in the case of 1D model from matlab or if drex produces tensor in
        ! local (radial) coordinates system use rad_to_glob
        !! call rotate_tensor_radial_to_global(theta,phi,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
        !!                                   c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
        !!                                 c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
        !!                               c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)


! ----------------------------------------------------------
! 02/09/2020 modification to create a 1D or 3D radially anisotropic model
 !activate it only if you want to calculate approximated tensor
 
        IF(1==0) THEN  !calculation radial anisotropic model

        ! select constant values only for 1D cases and put them into build_cij_drex
        colat_const = 90.d0
        lon_const = 40.d0

        call build_cij_drex(AMM_V_pro,rho,r,colat,lon, &
                   c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36, &
                   c44,c45,c46,c55,c56,c66) 
        

        ! calculate Love Parameters (Love 1927)
        A = (3.d0/8.d0)*(c11+c22)+0.25d0*c12+0.5d0*c66
        C = c33
        F = 0.5d0*(c13+c23)
        L = 0.5d0*(c44+c55)
        N = (1.d0/8.d0)*(c11+c22)-0.25d0*c12+0.5d0*c66

        ! radially anisotropic elastic tensor:
        !      A       A-2N   F    0        0        0
        !      A-2N     A    F     0        0        0
        !      F        F    C     0        0        0
        !      0        0    0     L        0        0
        !      0        0    0     0        L        0
        !      0        0    0     0        0        N

! calculate back the 21 components of which 5 are independent
c11=A
c12=A-(2.d0*N)
c13=F
c14=0.d0
c15=0.d0
c16=0.d0
c22=A
c23=F
c24=0.d0
c25=0.d0
c26=0.d0
c33=C
c34=0.d0
c35=0.d0
c36=0.d0
c44=L
c45=0.d0
c46=0.d0
c55=L
c56=0.d0
c66=N

        
        END IF  !end of calculation radial anisotropic model
! ---------------------------------------------------------

endif

  end subroutine model_aniso_mantle_drex

!--------------------------------------------------------------------

  subroutine build_cij_drex(pro,rho,r,colat,lon, &
                       d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,d33,d34,d35,d36, &
                       d44,d45,d46,d55,d56,d66)

  use constants, only: ZERO,R_EARTH,R_EARTH_KM,R_UNIT_SPHERE,DEGREES_TO_RADIANS,PI,GRAV,RHOAV
  use model_aniso_mantle_drex_par

  implicit none

  double precision,intent(in) :: pro(nz)  

  double precision,intent(out) :: rho

   double precision,intent(in) :: r,colat,lon
  double precision,intent(out) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66  

  ! local parameters
 
 double precision :: tet,ph,dtet,dph,depth,PwaveMod,SwaveMod
  double precision :: d1,d2,d3,d4,pc1,pc2,pc3,pc4, &
                   dpr1,dpr2,scale_GPa,scaleval
  double precision :: anispara(22,2,4),elpar(22)
  integer :: idep,ipar,icolat,ilon,icz0, &
          ict0,ict1,icp0,icp1,icz1 


  tet = colat !in degrees 
  ph = lon   ! ""


! avoid edge effects
! the grid has to be inside the model chunck
if (tet < AMM_V_colat(1)) tet = AMM_V_colat(1)+0.00001
if (ph < AMM_V_lon(1)) ph = AMM_V_lon(1)+0.00001
if (tet > AMM_V_colat(ny)) tet = AMM_V_colat(ny)-0.00001
if (ph > AMM_V_lon(nx)) ph = AMM_V_lon(nx)-0.00001



! dimensionalize
  depth = R_EARTH_KM*(R_UNIT_SPHERE - r)

  if (depth <= pro(nz) .or. depth >= pro(1)) print *, 'depth',depth,'pro(nz)',pro(nz),'pro(1)',pro(1)
  if (depth <= pro(nz) .or. depth >= pro(1)) call exit_MPI_without_rank('r out of range in build_cij')


  ict0 = 0
  do icolat = 1,ny-1  ! I add -1 to avoid icolat=ny and ict0>ny
    if (AMM_V_colat(icolat) < tet) ict0 = ict0 + 1 
  enddo

  icp0 = 0
  do ilon = 1,nx-1
    if (AMM_V_lon(ilon) < ph) icp0 = icp0 + 1  
  enddo

  icz0 = 0
  do idep = 1,nz     
    if (pro(idep) > depth) icz0 = icz0 + 1
  enddo



  ict1 = ict0 + 1 

  icp1 = icp0 + 1

  icz1 = icz0 + 1


  if (icp0 < 1 .or. icp0 > nx) print *,ph,icp0,nx,AMM_V_lon(1),AMM_V_lon(nx)
  if (ict0 < 1 .or. ict0 > ny) print *,tet,ict0,ny,AMM_V_colat(1),AMM_V_colat(ny)
  if (ict1 < 1 .or. ict1 > ny) print *,tet,ict1,ny,AMM_V_colat(1),AMM_V_colat(ny)
! check that parameters make sense
  if (ict0 < 1 .or. ict0 > ny) call exit_MPI_without_rank('ict0 out of range')
  if (ict1 < 1 .or. ict1 > ny) call exit_MPI_without_rank('ict1 out of range')
  if (icp0 < 1 .or. icp0 > nx) call exit_MPI_without_rank('icp0 out of range')
  if (icp1 < 1 .or. icp1 > nx) call exit_MPI_without_rank('icp1 out of range')
  if (icz0 < 1 .or. icz0 > nz) call exit_MPI_without_rank('icz0 out of range')
  if (icz1 < 1 .or. icz1 > nz) call exit_MPI_without_rank('icz1 out of range')


! intepolate the 22 parameters from AMM_V_Cij to 8 nodes of the integration cell
! I do that 22 times 

  do ipar = 1,22
    anispara(ipar,1,1) = AMM_V_Cij(ipar,icp0,ict0,icz0)
    anispara(ipar,2,1) = AMM_V_Cij(ipar,icp1,ict0,icz0)
    anispara(ipar,1,2) = AMM_V_Cij(ipar,icp0,ict0,icz1)
    anispara(ipar,2,2) = AMM_V_Cij(ipar,icp1,ict0,icz1)
    anispara(ipar,1,3) = AMM_V_Cij(ipar,icp0,ict1,icz0)
    anispara(ipar,2,3) = AMM_V_Cij(ipar,icp1,ict1,icz0)
    anispara(ipar,1,4) = AMM_V_Cij(ipar,icp0,ict1,icz1)
    anispara(ipar,2,4) = AMM_V_Cij(ipar,icp1,ict1,icz1)
  enddo

!
! calculation of distances between the selected point and grid points
!
    dtet = (tet - AMM_V_colat(ict0))/(AMM_V_colat(ict1) - AMM_V_colat(ict0))
    dph  = (ph  - AMM_V_lon(icp0))/(AMM_V_lon(icp1) - AMM_V_lon(icp0)) 
 ! check that parameters make sense
    if (dtet<0.d0)  call exit_MPI_without_rank('dtet < 0')
    if (dph<0.d0)  call exit_MPI_without_rank('dph < 0')
    if (dtet>1.d0)  call exit_MPI_without_rank('dtet > 1')
    if (dph>1.d0)  call exit_MPI_without_rank('dph > 1')


  d1 = (1.0 - dtet)*(1.0 - dph)

  d2 = dtet*(1.0 - dph)

  d3 = (1.0 - dtet)*dph

  d4 = dtet*dph

 dpr1 = (depth - pro(icz1))/(pro(icz0) - pro(icz1))
 dpr2 = 1.0 - dpr1

!Manuele 25/05/2020
  do ipar = 1,22
     pc1 = anispara(ipar,1,1)*dpr1+anispara(ipar,1,2)*dpr2
     pc2 = anispara(ipar,1,3)*dpr1+anispara(ipar,1,4)*dpr2
     pc3 = anispara(ipar,2,1)*dpr1+anispara(ipar,2,2)*dpr2
     pc4 = anispara(ipar,2,3)*dpr1+anispara(ipar,2,4)*dpr2
     elpar(ipar) = pc1*d1 + pc2*d2 + pc3*d3 + pc4*d4
  enddo

  d11 = ZERO
  d12 = ZERO
  d13 = ZERO
  d14 = ZERO
  d15 = ZERO
  d16 = ZERO
  d22 = ZERO
  d23 = ZERO
  d24 = ZERO
  d25 = ZERO
  d26 = ZERO
  d33 = ZERO
  d34 = ZERO
  d35 = ZERO
  d36 = ZERO
  d44 = ZERO
  d45 = ZERO
  d46 = ZERO
  d55 = ZERO
  d56 = ZERO
  d66 = ZERO
!
!   create dij
!
  rho = elpar(1)
  d11 = elpar(2)
  d12 = elpar(3)
  d13 = elpar(4)
  d14 = elpar(5)
  d15 = elpar(6)
  d16 = elpar(7)
  d22 = elpar(8)
  d23 = elpar(9)
  d24 = elpar(10)
  d25 = elpar(11)
  d26 = elpar(12)
  d33 = elpar(13)
  d34 = elpar(14)
  d35 = elpar(15)
  d36 = elpar(16)
  d44 = elpar(17)
  d45 = elpar(18)
  d46 = elpar(19)
  d55 = elpar(20)
  d56 = elpar(21)
  d66 = elpar(22)


! -----------------------------------------------------------------------
! 04/06/2020 calculate isotropic elastic tensor if necessary

  IF(1==0) THEN  
        ! elastic tensor for hexagonal symmetry in reduced notation:
        !      c11 c12 c13  0   0        0
        !      c12 c11 c13  0   0        0
        !      c13 c13 c33  0   0        0
        !       0   0   0  c44  0        0
        !       0   0   0   0  c44       0
        !       0   0   0   0   0  c66=(c11-c12)/2
        ! where 
        ! c11=k+4/3G
        ! c22=k+4/3G
        ! c33=k+4/3G
        ! c44=G
        ! c55=G
        ! c66=G
        ! c12=k-2/3G
        ! c13=k-2/3G
        ! c23=k-2/3G

! first - calculate PwaveMod and SwaveMod
          PwaveMod = (3.d0/15.d0)*(d11+d22+d33)+(2.d0/15.d0)*(d23+d13+d12)+(4.d0/15.d0)*(d44+d55+d66)
          SwaveMod = (1.d0/15.d0)*(d11+d22+d33)-(1.d0/15.d0)*(d23+d13+d12)+(3.d0/15.d0)*(d44+d55+d66)
! second - calculate isotropic component

  d11 = ZERO
  d12 = ZERO
  d13 = ZERO
  d14 = ZERO
  d15 = ZERO
  d16 = ZERO
  d22 = ZERO
  d23 = ZERO
  d24 = ZERO
  d25 = ZERO
  d26 = ZERO
  d33 = ZERO
  d34 = ZERO
  d35 = ZERO
  d36 = ZERO
  d44 = ZERO
  d45 = ZERO
  d46 = ZERO
  d55 = ZERO
  d56 = ZERO
  d66 = ZERO


  d11 = PwaveMod
  d12 = PwaveMod-2.d0*SwaveMod
  d13 = PwaveMod-2.d0*SwaveMod
  d14 = 0.d0
  d15 = 0.d0
  d16 = 0.d0
  d22 = PwaveMod
  d23 = PwaveMod-2.d0*SwaveMod
  d24 = 0.d0
  d25 = 0.d0
  d26 = 0.d0
  d33 = PwaveMod
  d34 = 0.d0
  d35 = 0.d0
  d36 = 0.d0
  d44 = SwaveMod
  d45 = 0.d0
  d46 = 0.d0
  d55 = SwaveMod
  d56 = 0.d0
  d66 = SwaveMod

  END IF
! -------------------------------------------------------------

! non-dimensionalize the elastic coefficients using
! the scale of GPa--[g/cm^3][(km/s)^2]
  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_GPa =(RHOAV/1000.d0)*((R_EARTH*scaleval/1000.d0)**2)


  d11 = d11/scale_GPa
  d12 = d12/scale_GPa
  d13 = d13/scale_GPa
  d14 = d14/scale_GPa
  d15 = d15/scale_GPa
  d16 = d16/scale_GPa
  d22 = d22/scale_GPa
  d23 = d23/scale_GPa
  d24 = d24/scale_GPa
  d25 = d25/scale_GPa
  d26 = d26/scale_GPa
  d33 = d33/scale_GPa
  d34 = d34/scale_GPa
  d35 = d35/scale_GPa
  d36 = d36/scale_GPa
  d44 = d44/scale_GPa
  d45 = d45/scale_GPa
  d46 = d46/scale_GPa
  d55 = d55/scale_GPa
  d56 = d56/scale_GPa
  d66 = d66/scale_GPa

! ---------------------------------------------------------------------------------
! non-dimensionalize
  rho = rho/RHOAV ! rho*1000.d0/RHOAV  rho [Kg/m3] RHOAV 5514.3 [kg/m3]

  end subroutine build_cij_drex

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_aniso_mantle_model_drex()

  use constants, only: IIN,DEGREES_TO_RADIANS 
  use shared_parameters, only:RMOHO
  use model_aniso_mantle_drex_par

  implicit none

  ! local parameters
  integer :: i,k,l,param,paramp
  integer :: ier
  character(len=*), parameter :: drex_model = 'DATA/DREX/drex_model.xyz'
  character(len=*), parameter :: iso_prem = 'DATA/DREX/iso_prem.dat'
  
! read the model iso_prem model

  open(IIN,file = iso_prem,status='old',action='read',iostat=ier)

  if (ier /= 0 ) stop 'Error opening file iso_prem'

! read the number of nodes in the three dimensions
  read(IIN, '(a)',end = 88)
  read(IIN, *,end = 88) nzp
  read(IIN, '(a)',end = 88)
  read(IIN, *,end = 88) AMM_V_prop


  AMM_V_prop = AMM_V_prop/1000.00 ! convert from meters to km
 

  read(IIN, '(a)',end = 88)
  do l = 1,nzp

     read(IIN,*,end = 88) (AMM_V_Cijp(paramp,l),paramp=1,3)

  enddo


88 close(IIN)

  !Read DREX model
  open(IIN,file=drex_model,status='old',action='read',iostat=ier)

  if (ier /= 0 ) stop 'Error opening file drex_model.xyz'

! read the number of nodes in the three dimensions
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) nx, ny, nz
  !print *, 'nx',nx,'ny',ny,'nz',nz

  !nz  Number of layers in the DREX model
 

  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_lon
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_colat
  read(IIN, '(a)',end = 888)
  read(IIN, *,end = 888) AMM_V_pro

  ! REMEMBER to modify this part according to your needs
  AMM_V_pro = (AMM_V_pro + (6371000.00 - RMOHO) - 500)/1000.00  ! convert meters to Km and add 24.4km of moho
  !AMM_V_pro = (AMM_V_pro) / 1000.00  ! we want it in km  
  AMM_V_lon = (AMM_V_lon / DEGREES_TO_RADIANS) !+ 80.00  ! +40.00 convert radians to degree 
  AMM_V_colat = (AMM_V_colat / DEGREES_TO_RADIANS)  


  read(IIN, '(a)',end = 888)

  do l = 1,nz
    do i = 1,ny
      do k = 1,nx
           read(IIN, *,end = 888) (AMM_V_Cij(param,k,i,l), param=1,22)
           if (AMM_V_Cij(1,k,i,l) == 0 .and. AMM_V_Cij(2,k,i,l) == 0 .and. AMM_V_Cij(3,k,i,l) == 0 & 
           .and. AMM_V_Cij(4,k,i,l) == 0) then 
                  do param =1,22 
                     IF(k>1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k-1,i,l)
                     IF(k==1 .AND. i> 1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k,i-1,l)
                     IF(k==1 .AND. i==1 .AND. l > 1) AMM_V_Cij(param,k,i,l)=AMM_V_Cij(param,k,i,l-1)
                  enddo
           endif
           ! if an entire row of parameters is zero make it equal to the
           ! previous row
       enddo
    enddo
  enddo

888 close(IIN)


  end subroutine read_aniso_mantle_model_drex

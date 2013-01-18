!****************************************************************************
!****************************************************************************
!*                                      *
!*  ModuleFile : Global_variables                         *
!*  Purpose    : This module includes global variables, which are used to  *
!*        peridynamics computations and data storages             *
!*                                      *
!*  PROGRAMMER : Dr. Youn Doh Ha                      *
!*                                      *
!****************************************************************************
module Global_variables

  use Objects

  !_____________________ Refinement Global variables _____________________
  type(node), allocatable,target :: nodes(:), original_nodes(:) ! global vector of nodes
  integer(4) ::nnodes ! number of nodes
  integer(4) ::nnodes_original ! number of nodes in the initial uniform grid
  integer(4) ::nbc, nforce ! number of nodes for prescribing boundary conditions and tractions  
    
  type(element), allocatable, target :: elements(:) ! global vector of  leath elements
  type(element), allocatable :: old_elements(:) ! old global vector of leath elements
  type(element), allocatable, target :: quadtree_elements(:), original_quadtree_elements(:)! initial vector of elements( input elements)
  integer(4) :: nelements_leath ! global number of elements
  integer(4) :: nroot_elements
  integer(4) :: global_depth
  type(element) :: aux_elements(8)

  logical :: grid_ok =.true. ! In the case of two neighbor elements with difference in the level of refinement larger than 2,
                           ! this variable recevies .false. 
  
  integer(4), allocatable :: temp(:)
  integer(4) :: max_depth  ! maximum depth of an element on the grid
  integer(4) :: min_depth  ! minimum depth of an element on the grid
  integer(4) :: max_size  ! maximum size of an element on the grid
  integer(4) :: min_size  ! minimum size of an element on the grid
  integer(4) :: dim     ! dimension of the problem: if dim=2, then is a 2D problem
  integer(4) :: iteration=0
  integer(4) :: iterations
  integer(4) :: mnode
  real(8) :: x_min_lim, y_min_lim, z_min_lim, x_max_lim, y_max_lim, z_max_lim,cont

  integer(4) :: nbox_x,nbox_y,nboxes
  real(8) :: max_strain_energy
  real(8) :: damage_index_cri
  
  logical :: modified_mesh=.false. ! This variable tells us if the mesh have been modified( refined or coarsed) 

  !_____________________ Solver Global variables _____________________
  integer(4), dimension(:,:), allocatable :: bond
  real(8), dimension(:), allocatable :: critical_strain
  real(8) :: critical_stretch
  real(8) :: fracture_energy
  real(8) :: dc1, dc2, dc3
  integer(4) :: myn,mbd,nbd1,nbd2
  real(8) :: delta, dt
    logical :: brokij, brokij_old, broke_first

  logical, dimension(:), allocatable :: mine(:)

  integer(4), dimension(:,:), allocatable :: broke
  integer(4) :: nbroke, ncyc
  logical, dimension(:), allocatable :: omitt

    real(8), dimension(:), allocatable ::  edt, edt_add
!    , ut1_old, ut2_old
    real(8), dimension(:), allocatable ::  ut1new, ut2new, ut1, ut2
    real(8), dimension(:), allocatable ::  ut1old, ut2old
    real(8), dimension(:), allocatable ::  f1, f2
  real(8), dimension(:), allocatable ::  wt

!  real(8), dimension(:), allocatable :: k1_x,k1_y,k2_x,k2_y,k3_x,k3_y,k4_x,k4_y ! Variables used in the Runge-Kutta time integration

    real(8), dimension(:), allocatable :: acc1, acc2  ! Acceleration of the nodes: acc1 = acceleration in the x-direction, acc2 = acceleration in the y-direction 

  real(8), dimension(:), allocatable :: lc1bd,lc2bd, u1bd,u2bd,v1bd,v2bd, tendbd

    real(8), dimension(:), allocatable :: vis
  integer(4), dimension(:), allocatable :: node_type, nofail

  real(8), dimension(:), allocatable :: fcoefm
! coefficents of bond
  integer(4), dimension(:), allocatable :: fnorm
! normalization of bond force
  real(8), dimension(:), allocatable :: crit_exten
! critical extension of bonds
  real(8), dimension(:), allocatable :: damage_index   ! Each node has a value for damage_index, which is the ration of broken bonds for a node to the total number of bonds for the same node in the reference configuration
! critical extension of bonds
  real(8), dimension(:), allocatable :: volnod
! volume of node
!  real(8), dimension(:,:), allocatable :: ffint
! coefficients
  integer(4), dimension(:), allocatable :: nodbd, mbdtyp
  real(8), dimension(:), allocatable :: spsum
  real(8):: young
  real(8):: denst
    integer(4),parameter :: ncell1=10
    integer(4),parameter :: ncell2=10
    integer(4),dimension(ncell1,ncell2) :: cell_location_in_bin, number_of_nodes_in_cell
    integer(4), parameter :: maxfam=2000 
  integer(4) :: mfam, m_def_fam
  integer(4), parameter :: max_global_nodes=2000000   ! Maximum number of nodes 
    integer(4), dimension(:), allocatable :: family
  integer(4), dimension(:), allocatable :: nodes_in_bin
    integer(4),parameter :: ncell_def1=10
    integer(4),parameter :: ncell_def2=10
    integer(4),dimension(ncell_def1,ncell_def2) :: cell_location_in_bin_def, number_of_nodes_in_cell_def
    integer(4), dimension(:), allocatable :: def_family
  integer(4), dimension(:), allocatable :: nodes_in_bin_def
  real(8), dimension(:), allocatable :: radnod
  real(8) :: rad_search, radnod_max

  ! parameters
    real(8) :: x1_lo, x1_hi
    real(8) :: x2_lo, x2_hi
    real(8) :: y1_lo, y1_hi
    real(8) :: y2_lo, y2_hi
    real(8), dimension(:), allocatable :: xt1, xt2, yt1, yt2
 
    real(8), dimension(:), allocatable :: a,b
    integer(4), dimension(:), allocatable :: jsnr,irnr
  integer(4) :: mctr
!    integer ::  n_elements,info_element(21000,4)
    integer(4) ::  n_elements
    integer(4), dimension(:,:), allocatable :: info_element
! Number of deleted regions
  integer(4),parameter  :: ndel=1
  integer(4) :: mdtype(ndel)
    real(8) :: x1dlo(ndel), x1dhi(ndel), x2dlo(ndel), x2dhi(ndel)
  real(8) ::  x1dcen(ndel), x2dcen(ndel), raddl(ndel)
  type(line),allocatable :: Crack_Line(:)
  integer(4) :: ncracks=0
  logical :: newiteration

! Time steps
    integer(4) :: iter
  integer(4) :: nt

  integer(4):: snapshots_frequence
! external forces
 !   real(8) :: ext1(max_global_nodes), ext2(max_global_nodes), interval(max_global_nodes,2), RadHorizon(max_global_nodes)
 
    real(8), dimension(:), allocatable :: xtempbotton(:), ytempbotton(:), xtemptop(:), ytemptop(:)
    integer(4), dimension(:), allocatable :: ncorrespond_botton(:), ncorrespond_top(:)

    real(8), dimension(:), allocatable :: ext1, ext2, RadHorizon
    real(8), dimension(:,:), allocatable :: interval
  real(8) :: shortrange_dist_fac_nom(10),  shortrange_dist_fac_init(10),  shortrange_force_fac(10) 
  real(8) :: speed_iter=0
  real(8),dimension(:), allocatable :: top_crack_speed,botton_crack_speed  
  real(8) :: bc_velocity
  real(8) :: force_mag    !>> 03062009_YounDoh; external force magnitude
  real(8) :: visco, visdk  !>> 03072009_YounDoh; damage viscosity parameters
  real(8) :: horizon_factor
! Decision variables
  integer(4) :: do_adap_refinement
  integer(4) :: introduce_precrack
  integer(4) :: dynamic_or_static
  integer(4) :: sw_micro
! Initial Crack
  real(8) :: initial_crack_x0,initial_crack_y0,initial_crack_x1,initial_crack_y1
   
end module Global_variables

module Input_subroutines

 use Objects
 use Global_variables
 ! use Refine_subroutines

contains


subroutine read_input_file(filename2)

  character(19), intent(in) :: filename2
  ! integer :: number_of_iterations
  ! integer :: nlomitt, i, j
  ! integer nomitt(2)
  integer   :: w
  character(100):: char, char2, char3

  !  Check whether the file exists
  logical :: file_exists
  inquire(file=filename2, exist=file_exists)
  if (.not. file_exists) then
    print *, "Input file", filename2, "does not exist"
  endif
  
  w = 11 
  open (w, FILE= filename2, STATUS='OLD')
    
  print*,"......Input Main File....."  
  read(w,*) char
  read(w,*) char
  read(w,*) char
  read(w,*) char ! nodes input file
  print *, "Input node file = ", char
    
  read(w,*) char2
  read(w,*) char2 ! elements input file
  print *, "Input element file = ", char

  !read(w,*) char3
  !read(w,*) do_adap_refinement ! Flag to decide if adaptive refinement procedure will be applied: 1 = yes, 0 = no

  read(w,*) char3
  read(w,*) dynamic_or_static ! Flag to decide if the analysis will be dynamic or static: 1 = dynamic, 0 = static  

  read(w,*) char3
  read(w,*) sw_micro ! Flag to decide what kind of micromodulus is used : 0 = constant, 1 = conical micromodulus

  if (do_adap_refinement == 1) then
    mnode = 4
  else
    mnode = 1
  endif

  call assign_nodes(char)
  call assign_elements(char2)
  allocate(omitt(nnodes*mnode))
  omitt=.FALSE.

  read(w,*) char3
  read(w,*) nt  ! number of iterations
  read(w,*) char3
  read(w,*) dt  ! step size
  read(w,*) char3
  read(w,*) snapshots_frequence  ! Number of iterations between snapshots
  read(w,*) char3
  read(w,*) young  ! Young's modulus
  read(w,*) char3
  read(w,*) denst  ! density
  read(w,*) char3
! read(w,*) critical_stretch  ! critical stretch
  read(w,*) fracture_energy  ! fracture energy
  read(w,*) char3
  read(w,*) horizon_factor !Horizon's scalar factor: horizon is calculated as the delta_x(distance between nodes) times a scalar
  read(w,*) char3
  read(w,*) introduce_precrack ! Flag to decide if the following precrack will be applied: 1 = yes, 0 = no
  read(w,*) char3
  read(w,*) initial_crack_x0,initial_crack_y0,initial_crack_x1,initial_crack_y1 ! Initial crack (x0,y0) (x1,y1)
! read(w,*) char3
! read(w,*) nlomitt ! Number of layers for omitted nodes
! read(w,*) char3
! if (nlomitt.gt.0) then
!   do i=1, nlomitt
!     read(w,*) nomitt(1),nomitt(2)
!     if (nomitt(1).ge.nomitt(2)) then
!       do j=nomitt(2), nomitt(1)
!         omitt(j) = .TRUE.
!       enddo
!     else
!       do j=nomitt(1), nomitt(2)
!         omitt(j) = .TRUE.
!       enddo
!     endif
!   enddo
! endif
  read(w,*) char3
  read(w,*) bc_velocity ! Value of the applied velocity for boundary conditions
!<<
  read(w,*) char3
  read(w,*) force_mag ! Value of the external force magnitude
!>> 03062009_YounDoh
!<<
  read(w,*) char3
  read(w,*) visco, visdk ! Damage viscosity parameters
!>> 03072009_YounDoh
  read(w,*) char3
  read(w,*) damage_index_cri ! damage index criterion for adaptive refinement
  read(w,*) char3
  read(w,*) dc1, dc2, dc3 ! damage stretch coefficients (dc1, dc2, dc3)
  
    
! call assign_nodes(char)
!! if (do_adap_refinement == 1)then ! If  is  =1 then the refinement and coarsening will be performed
!   call assign_elements(char2)
!! endif
  
end subroutine read_input_file


!***** method that reads the input node file and allocate these nodes to the global array of nodes : nodes(i)
subroutine assign_nodes (filename)
    
  character(25), intent(in) :: filename
  integer:: n_nodes, i, j, idnode, u
  real(8):: strain_energy = 0
  character(90):: char
    
  u = 63
  open (u, FILE= filename, STATUS='OLD')
    
  print*,"......Input Nodes....."
    
  read(u,*) char,dim

  read(u,*) char,n_nodes
  read(u,*) char,nbc
  read(u,*) char,nforce
  read(u,*) char
  nnodes = n_nodes
  nnodes_original = nnodes

  print *, "n_nodes = ", n_nodes, " mnode = ", mnode
  allocate(nodes(n_nodes*mnode))
  if (allocated(nodes)) then
    print *, "nodes allocated"
  endif
  allocate(original_nodes(nnodes_original*1))
  if (allocated(original_nodes)) then
    print *, "original_nodes allocated"
  endif

  do i = 1,n_nodes
    allocate(nodes(i)%pos(dim))
    read(u,*) idnode,(nodes(i)%pos(j), j=1, dim)
    nodes(i)%id = idnode
    allocate(nodes(i)%disp(dim))
    nodes(i)%strain_energy = strain_energy
    
    nodes(i)%iflag = .false.
    allocate(nodes(i)%old_disp(dim))
    allocate(nodes(i)%new_disp(dim))
    allocate(nodes(i)%veloc(dim))
    allocate(nodes(i)%new_veloc(dim))
    allocate(nodes(i)%accel(dim))


    !**original_nodes array:
    allocate(original_nodes(i)%pos(dim))
    allocate(original_nodes(i)%disp(dim))
    allocate(original_nodes(i)%new_disp(dim))
    allocate(original_nodes(i)%veloc(dim))
    allocate(original_nodes(i)%new_veloc(dim))
    allocate(original_nodes(i)%accel(dim))
    
      
    original_nodes(i)%id = idnode
    do j=1, dim
      original_nodes(i)%pos(j) = nodes(i)%pos(j)
      original_nodes(i)%disp(j) = nodes(i)%disp(j)
    enddo
    original_nodes(i)%strain_energy = nodes(i)%strain_energy
    original_nodes(i)%iflag = nodes(i)%iflag

    allocate(nodes(i)%bc(dim))
    allocate(nodes(i)%force(dim))
  end do

  do i = 1,nbc
    read(u,*) idnode,(nodes(idnode)%bc(j), j=1, dim)
  end do

  do i = 1,nforce
    read(u,*) idnode,(nodes(idnode)%force(j), j=1, dim)
  end do

  close(u)

end subroutine assign_nodes



!*** method that read the input element file and allocate these elements to the global array of elements: elements(i)
subroutine assign_elements(filename)
    
  character(len=*), intent(in) :: filename
  integer:: i,j, idelement, idnode1, idnode2, idnode3, idnode4, u, status
  character(90):: char
  real(8)::size,x1,x2,x3,x4,y1,y2,y3,y4
    
  u=64
  open (u, FILE= filename, STATUS='OLD')
    
  print*,"......Input Elements....."
    
  read(u,*) char
  read(u,*) n_elements
  read(u,*) char
  nroot_elements = n_elements
  nelements_leath=n_elements
  max_depth=1
  min_depth=1
  
  print *, "n_elements = ", n_elements
  allocate(info_element(n_elements*2,4), stat=status)
  if (allocated(info_element)) then
    print *, "info_element allocated %d", n_elements
  endif
  allocate(elements(n_elements*2), stat=status)
  allocate(quadtree_elements(n_elements*2), stat=status)
  allocate(original_quadtree_elements(n_elements*2), stat=status)
  allocate(temp(0), stat=status)
    
  min_size = 1000000000
  max_size = -1
  global_depth=0

    
  do i = 1,n_elements  ! Find, based on the node id number, correct nodes for the element  
        
    read(u,*) idelement,idnode1, idnode2, idnode3, idnode4      
        
    allocate(elements(i)%elementnodes(4))
    allocate(elements(i)%node1)
    allocate(elements(i)%node2)
    allocate(elements(i)%node3)
    allocate(elements(i)%node4)

    allocate(original_quadtree_elements(i)%elementnodes(4))
    allocate(original_quadtree_elements(i)%node1)
    allocate(original_quadtree_elements(i)%node2)
    allocate(original_quadtree_elements(i)%node3)
    allocate(original_quadtree_elements(i)%node4)

    info_element(i,1)=idnode1
    info_element(i,2)=idnode2
    info_element(i,3)=idnode3
    info_element(i,4)=idnode4

        
    elements(i)%id = idelement
    elements(i)%leath = .true.
    elements(i)%depth = 1
    elements(i)%root = .true.

    original_quadtree_elements(i)%id = idelement
    original_quadtree_elements(i)%leath = .true.
    original_quadtree_elements(i)%depth = 1
    original_quadtree_elements(i)%root = .true.

    do j = 1,nnodes
      if(nodes(j)%id == idnode1)then      
        elements(i)%node1 => nodes(j)
        elements(i)%elementnodes(1) = nodes(j)
            
        original_quadtree_elements(i)%node1 =>original_nodes(j)
        original_quadtree_elements(i)%elementnodes(1) =original_nodes(j)

      end if   
      if(nodes(j)%id == idnode2)then          
        elements(i)%node2 => nodes(j)
        elements(i)%elementnodes(2) = nodes(j)

        original_quadtree_elements(i)%node2 =>original_nodes(j)
        original_quadtree_elements(i)%elementnodes(2) =original_nodes(j)
          
      end if
      if(nodes(j)%id == idnode3)then
        elements(i)%node3 => nodes(j)
        elements(i)%elementnodes(3) = nodes(j)

        original_quadtree_elements(i)%node3 =>original_nodes(j)
        original_quadtree_elements(i)%elementnodes(3) =original_nodes(j)
  
      end if
      if(nodes(j)%id == idnode4)then
        elements(i)%node4 => nodes(j)
        elements(i)%elementnodes(4) = nodes(j)

        original_quadtree_elements(i)%node4 =>original_nodes(j)
        original_quadtree_elements(i)%elementnodes(4) =original_nodes(j)

      end if
          
    end do

    !  quadtree_elements(i) = elements(i)
    call Element_atribution(quadtree_elements(i),elements(i))
    x1 = elements(i)%elementnodes(1)%pos(1)
    x2 = elements(i)%elementnodes(2)%pos(1)
    x3 = elements(i)%elementnodes(3)%pos(1)
    x4 = elements(i)%elementnodes(4)%pos(1)
    y1 = elements(i)%elementnodes(1)%pos(2)
    y2 = elements(i)%elementnodes(2)%pos(2)
    y3 = elements(i)%elementnodes(3)%pos(2)
    y4 = elements(i)%elementnodes(4)%pos(2)
        

    size = dabs(0.5d0*((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y4-x4*y3)+(x4*y1-x1*y4)))
                
    if(int(size)<min_size)then
      min_size = int(size)
    endif

    if(int(size)>max_size)then
      max_size = int(size)
    endif
  
              
  end do
  close(u)

end subroutine assign_elements


end module Input_subroutines

module Objects
    

  !********* this structure defines the node type
    type node
      !type(mat_model), pointer  :: mat
      integer(4) :: id
      real(8) :: density
        integer(4) :: mat_type
      real(8) :: young
      real(8) :: strain_energy =0
      real(8) :: damage_index
      real(8), allocatable :: pos(:)
        real(8), allocatable :: disp(:)
      real(8), allocatable :: veloc(:)
      real(8), allocatable :: accel(:)
      real(8), allocatable :: new_veloc(:)
      real(8), allocatable :: new_disp(:)
      real(8), allocatable :: old_disp(:)
      real(8), allocatable :: force(:)
      integer(4), allocatable :: bc(:)

      
    
      real(8) :: horizon_size
      logical :: iflag       ! iflag = 1: node is a hanging node; iflag = 0: node is not a hanging node
      real(8) :: volume 

    
      integer(4) :: nnodeelements
      integer(4) :: nodeelements(10) ! The id of the elements adjacent to this node, can be 1,2,3 or 4 elements
      real(8) :: nodeelements_size(10) ! The size(area) of the elements adjacent to this node, can be 1,2,3 or 4 elements
      integer(4) :: nodeelements_depth(10) ! The depth of each element adjacent to this node
      integer(4) :: nodeelements_nhanging_nodes(10) ! The number of hanging nodes for each element connected to this node
    end type node

    type element
      
      !   Node position on an element
      !    _______
      !  1|       |2 
      !  |       |
      !  |    |
      !  4|_______|3
      !      
      
      integer(4) :: id
      type(node), allocatable :: elementnodes(:)
    
      type(node),pointer :: node1
      type(node),pointer :: node2 
      type(node),pointer :: node3
      type(node),pointer :: node4 
      
      integer(4) :: depth
      logical :: leath  = .true.    ! tell us if this element is a leaf of the quadtree stucture
      logical :: root = .false.
      logical :: dif_level_refine = .false. ! flag that tell us if the element should be refined because a dif. level of refinment >2 situation happened
          logical :: strain_energy_refine 
      type(element),pointer :: child1    ! we can improve this by changing this pointer to point to a id number instead of a data structure  
      type(element),pointer :: child2    ! this id number is from the global array of elements
      type(element),pointer :: child3 
      type(element),pointer :: child4

      type(element),pointer :: father

      type(element),pointer :: quad_element
      

      !type(element),pointer :: next  ! next is a pointer to the begin of the array of elements_neighbors ids 
      integer(4),allocatable :: neighborhood(:)
      integer(4) :: n_neighbors 
    

    end type element

      
    type line
    
      real(8) :: x1,y1,x2,y2 
    
    end type line

  
  contains
    ! THis function is created to assign one element to another. For some Fortran compilers, you need to assign field by field
    ! of the structure when this type has a field that is a pointer
    subroutine Element_atribution(element1,element2)
      type(element) element1,element2
      integer i
      
      element1%id = element2%id
      allocate(element1%elementnodes(4))
      allocate(element1%neighborhood(element2%n_neighbors))
      allocate(element1%node1)
      allocate(element1%node2)
      allocate(element1%node3)
      allocate(element1%node4)
      element1%elementnodes(1) = element2%elementnodes(1)
      element1%elementnodes(2) = element2%elementnodes(2)
      element1%elementnodes(3) = element2%elementnodes(3)
      element1%elementnodes(4) = element2%elementnodes(4)
      element1%depth = element2%depth
      element1%leath = element2%leath
      element1%root = element2%root
      element1%child1 => element2%child1
      element1%child2 => element2%child2
      element1%child3 => element2%child3
      element1%child4 => element2%child4
      element1%father => element2%father

      element1%node1 => element2%node1 
      element1%node2 => element2%node2 
      element1%node3 => element2%node3 
      element1%node4 => element2%node4 

      element1%quad_element => element2%quad_element
      do i =1,element2%n_neighbors
        element1%neighborhood(i) = element2%neighborhood(i)
      end do
      element1%n_neighbors = element2%n_neighbors



    end subroutine Element_atribution

  
  
  end module Objects
  
  
  

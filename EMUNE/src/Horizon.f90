module Horizon

  use Global_variables
  use Input_subroutines
  use precision

  contains

  ! ****************************************************************************************************
  ! This subroutine will calculate the horizon)size for each node of any uniform and non-uniform grid. 
  ! This value is the largest edge connceted to this node times an scalar defined by the user
  !*****************************************************************************************************
  subroutine Calculate_Horizon()
    real(8) :: max_length,x1,x2,x3,x4,y1,y2,y3,y4, xlength,ylength
    ! real(8) :: xl1,xl2
    real(8) :: max_interval_x, max_interval_y
    ! real(8) ::  radnod_min
    integer(4)::cont, i, j, k
    type(node)::node1
    type(node)::node2
    do i=1, nnodes
      max_length=zero
      max_interval_x=zero
      max_interval_y=zero
!     xl1=0.0
!     xl2=0.0
      cont = 0
      nodes(i)%nodeelements_nhanging_nodes = 0
      nodes(i)%horizon_size = 0
      do j=1, nelements_leath  
        
        do k=1,4
          
!         if(elements(j)%elementnodes(k)%pos(1)==nodes(i)%pos(1) .and. elements(j)%elementnodes(k)%pos(2)==nodes(i)%pos(2)) then
          if(elements(j)%elementnodes(k)%id==nodes(i)%id) then
            cont = cont +1
            nodes(i)%nodeelements(cont) = elements(j)%id
            !For the calculation of the area of the elements connected to this node, we need to be carefull because
            ! these elements can be irregular polygons, so we are going to use the general area calculation formula
            x1 = elements(j)%node1%pos(1)
            x2 = elements(j)%node2%pos(1)
            x3 = elements(j)%node3%pos(1)
            x4 = elements(j)%node4%pos(1)
            y1 = elements(j)%node1%pos(2)
            y2 = elements(j)%node2%pos(2)
            y3 = elements(j)%node3%pos(2)
            y4 = elements(j)%node4%pos(2)
            
            
            ! If this node is a hanging node
            if(elements(j)%node1%iflag .eqv. .TRUE.)then
              nodes(i)%nodeelements_nhanging_nodes(cont) = nodes(i)%nodeelements_nhanging_nodes(cont) +1
            endif
            if(elements(j)%node2%iflag .eqv. .TRUE.)then
              nodes(i)%nodeelements_nhanging_nodes(cont) = nodes(i)%nodeelements_nhanging_nodes(cont) +1
            endif
            if(elements(j)%node3%iflag .eqv. .TRUE.)then
              nodes(i)%nodeelements_nhanging_nodes(cont) = nodes(i)%nodeelements_nhanging_nodes(cont) +1
            endif
            if(elements(j)%node4%iflag .eqv. .TRUE.)then
              nodes(i)%nodeelements_nhanging_nodes(cont) = nodes(i)%nodeelements_nhanging_nodes(cont) +1
            endif

            nodes(i)%nodeelements_depth(cont) = elements(j)%depth  
            nodes(i)%nodeelements_size(cont) = dabs(0.5*((x1*y2-x2*y1)+(x2*y3-x3*y2)+(x3*y4-x4*y3)+(x4*y1-x1*y4) ))
           
            if(k.gt.1) then
              node1= elements(j)%elementnodes(k-1)
            else
              node1=elements(j)%elementnodes(4)
            endif
            if(k.lt.4) then
              node2=elements(j)%elementnodes(k+1)
            else
              node2=elements(j)%elementnodes(1)
            endif
          
!<<
!           xlength=sqrt((node1%pos(1)- node2%pos(1))*(node1%pos(1)-node2%pos(1))+(node1%pos(2)-node2%pos(2))*(node1%pos(2)-node2%pos(2)))
            xlength = dabs(node1%pos(1)-node2%pos(1))
            ylength = dabs(node1%pos(2)-node2%pos(2))
!>>02272009_YounDoh
!           max_length=max(max_length,xlength)
            max_interval_x=dmax1(max_interval_x, xlength)            
            max_interval_y=dmax1(max_interval_y, ylength)            
          endif
        enddo
      enddo
      nodes(i)%nnodeelements = cont

      interval(i,1) = max_interval_x
      interval(i,2) = max_interval_y
      max_length = dmax1(max_interval_x,max_interval_y)

      nodes(i)%horizon_size = horizon_factor*max_length
!<<
      radnod(i) = 0.5d0*max_length
!>>02272009_YounDoh
      
    enddo

    radnod_max = -1.0d-6
    do i = 1, nnodes
      if(radnod(i) .gt. radnod_max) then
        radnod_max = radnod(i)
      endif
    enddo

  end subroutine Calculate_Horizon


  !**** method to calculate the horizon of each node in a uniform grid
  subroutine cal_interval  

    real(8) :: max_interval_x, max_interval_y
    ! real(8) :: radnod_min
    integer(4) :: i,j,k,kn,kp
      

    do i = 1, nnodes
      max_interval_x=1.0d-6
      max_interval_y=1.0d-6
      do j = 1, n_elements
        do k = 1, 4
          kn = info_element(j,k)
          if(i == kn) then
            if(k.lt.4) then
              kp=info_element(j,k+1)
            else
              kp=info_element(j,1)
            endif
            max_interval_x=dmax1(dabs(nodes(kp)%pos(1)-nodes(kn)%pos(1)), max_interval_x)
            max_interval_y=dmax1(dabs(nodes(kp)%pos(2)-nodes(kn)%pos(2)), max_interval_y)
          endif
        enddo
      enddo
      interval(i,1)=max_interval_x
      interval(i,2)=max_interval_y
      nodes(i)%horizon_size= horizon_factor*dmax1(max_interval_x, max_interval_y)
!<<
      radnod(i) = 0.5d0*dmax1(max_interval_x, max_interval_y)
!>>02272009_YounDoh
    enddo
!<<
    radnod_max = -1.0d-6
    do i = 1, nnodes
      if(radnod(i) .gt. radnod_max) then
        radnod_max = radnod(i)
      endif
    enddo
!>>02272009_YounDoh

  end subroutine cal_interval

end module Horizon

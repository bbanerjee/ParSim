module volume_partition
  
  use Global_variables
  use precision
  contains

  !**** general method to calculate the volume of each node in any type of grid( uniform and non-uniform)
  subroutine Calculate_Volume() ! O(n)
        
    real(8) :: vol(10),depth(10),depth_max,depth_min
    ! real(8) :: vol_max,vol_min
    integer(4):: num,num_depth_min
    ! logical::t

    do i=1,nnodes
      num=nodes(i)%nnodeelements
      num_depth_min =0
      num_depth_max =0
      depth_min = 1
      nodes(i)%volume = zero

      do j =1,num
        vol(j) = nodes(i)%nodeelements_size(j)
        depth(j) = nodes(i)%nodeelements_depth(j)
           
      enddo
        
      depth_max = -100000
      depth_min = 10000000
      do j=1,num
        if(depth(j) >= depth_max)then
          depth_max= depth(j)
        endif
        if(depth_min >= depth(j))then
          depth_min = depth(j)
        endif
      enddo
        
      if (nodes(i)%iflag .eqv. .false.) then
          
        if(depth_max  == depth_min)then
            
          do j =1,num
            nodes(i)%volume = nodes(i)%volume + vol(j)*0.25d0
          enddo

        endif

        if(depth_max /= depth_min)then
            
          do j=1,num
            
            if(nodes(i)%nodeelements_nhanging_nodes(j) == 0)then
              nodes(i)%volume = nodes(i)%volume + vol(j)*0.25d0 
            endif

            if(nodes(i)%nodeelements_nhanging_nodes(j) == 1)then
              nodes(i)%volume = nodes(i)%volume + vol(j)*0.5d0
            endif  

            if(nodes(i)%nodeelements_nhanging_nodes(j) == 2)then
              nodes(i)%volume = nodes(i)%volume + vol(j)*0.75d0
            endif

          enddo
          
        endif

      else
        nodes(i)%volume =zero
      endif

    enddo

  end subroutine Calculate_Volume

  !**** method to calculate the volume of each node in a uniform grid
  subroutine cal_volume
    integer(4) :: i,j,k,kp
    real(8) :: area
   
    do i = 1, nnodes
      nodes(i)%volume = zero
    enddo

    do i = 1, n_elements
      area=zero
      do j = 1, 4
        k = info_element(i,j)
        if(j.lt.4) then
          kp=info_element(i,j+1)
          area=area+nodes(k)%pos(1)*nodes(kp)%pos(2)-nodes(kp)%pos(1)*nodes(k)%pos(2)
        else
          kp=info_element(i,1)
          area=area+nodes(k)%pos(1)*nodes(kp)%pos(2)-nodes(kp)%pos(1)*nodes(k)%pos(2)
        endif
      enddo
      if (i .eq. 1) then
        print *, 'Element = ', i, ' Area = ', area
      endif
      do j=1,4
        k = info_element(i,j)
        nodes(k)%volume=nodes(k)%volume+area*0.125d0
        !if (i .eq. 1) then
        !  print *, '   Node = ', k, ' Volume = ', nodes(k)%volume
        !endif
      enddo
  
    enddo
    do ii = 1, nnodes
      if (ii == 1 .or. ii == 2 .or. ii == 104 .or. ii == 105) then
        print *, '   Node = ', ii, ' Volume = ', nodes(k)%volume
      endif
    enddo

  end subroutine cal_volume

end module volume_partition

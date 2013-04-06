#include <BondFamilyComputer.h> 

using namespace Emu2DC;

BondFamilyComputer::BondFamilyComputer(const Domain& domain,
                                       const NodePArray& nodeList)
{
}

BondFamilyComputer::~BondFamilyComputer()
{
}

// Finds the family of node m: Find all the nodes inside the horizon of node m
 
// input:
//   m...                          node number (global)
//   horiz_p...                    horizon of node m
//   xt1_global,...                node positions
//   total_global_nodes...         number of nodes in the node position arrays
//   ncell1,ncell2,ncell3...       number of cells in each direction for sorting
//   maxfam...                     dimension of the family array
//   cell_location_in_bin...       pointer to the element in the bin corresponding to a cell
//   number_of_nodes_in_cell...    number of nodes whose position in each cell
//   nodes_in_bin...               node numbers sorted by cell compressed into the bin
//   x1_lo,x1_hi,...               upper and lower coordinates of search region
 
// output:
//   family...                     array containing the family of node m
//   mfam...                       number of nodes in the family of node m
void 
BondFamilyComputer::getFamilyReference(const Domain& domain,
                               const Node* node,
                               NodeArray& family)
{
 
    ! parameters
    int m
 
    ! local variables
    int cell1, cell2
    int bin_loc, bin_loc_start, bin_loc_end
    double dist_sq
    double dx1_cell, dx2_cell
    int search_lo_1, search_lo_2
    int search_hi_1, search_hi_2
    int search_cell_1, search_cell_2
    int search_node
    int i, j, n1, n2, n3
    ! static arrays used by sort_ref
 
    ! initialize
    family = 0

    ! set the cell limits and widths.
    dx1_cell = (x1_hi-x1_lo)/ncell1
    dx2_cell = (x2_hi-x2_lo)/ncell2
 
    ! find the cell number for node m.
    cell1 = 1+(nodes(m)%pos(1)-x1_lo)/dx1_cell
    cell2 = 1+(nodes(m)%pos(2)-x2_lo)/dx2_cell
 
    ! find the search distance for cells that might contain family members.
    search_lo_1 = 1+(nodes(m)%pos(1)-delta-x1_lo)/dx1_cell
    search_hi_1 = 1+(nodes(m)%pos(1)+delta-x1_lo)/dx1_cell
    search_lo_2 = 1+(nodes(m)%pos(2)-delta-x2_lo)/dx2_cell
    search_hi_2 = 1+(nodes(m)%pos(2)+delta-x2_lo)/dx2_cell
 
    search_lo_1 = max(1, search_lo_1)
    search_lo_2 = max(1, search_lo_2)
    search_hi_1 = min(ncell1, search_hi_1)
    search_hi_2 = min(ncell2, search_hi_2)
    ! do the search.
    ! write(*,*) m, delta, xt1(m),xt2(m), dx1_cell, dx2_cell, search_lo_1, search_hi_1
    ! pause 2222

    mfam = 0
    do search_cell_1 = search_lo_1, search_hi_1
      do search_cell_2 = search_lo_2, search_hi_2
        if(number_of_nodes_in_cell(search_cell_1,search_cell_2)>0) then
          bin_loc_start = cell_location_in_bin(search_cell_1,search_cell_2)
          bin_loc_end = bin_loc_start + number_of_nodes_in_cell(search_cell_1,search_cell_2)-1
          do bin_loc = bin_loc_start, bin_loc_end
            search_node = nodes_in_bin(bin_loc)
            dist_sq = (nodes(search_node)%pos(1)-nodes(m)%pos(1))**2 + (nodes(search_node)%pos(2)-nodes(m)%pos(2))**2 
            if(m/=search_node .and. dist_sq < delta**2) then
              mfam = mfam+1
              if(mfam <= maxfam) then
                family(mfam) = search_node
              endif
            endif
          enddo
        endif
      enddo
    enddo
 
    ! rearrange the order of family nodes
    do i=1, mfam
      n1=family(i)
      do j=i+1, mfam
        n2=family(j)
        if (n1.gt.n2) then
          n3=n1
          n1=n2
          n2=n3
        endif
      enddo
    enddo

    ! abort if the family array is not big enough.
 
    if(mfam > maxfam) then
      write(*,*) 'get_family: Too many nodes in family. maxfam=',maxfam, ' need ',mfam
      stop
    endif
 
    return
 
  end subroutine Get_family
}

void 
BondFamilyComputer::getFamilyDeformed(const Node* node,
                              NodeArray& family)
{
}

void 
BondFamilyComputer::sortNodesReference()
{
}

void 
BondFamilyComputer::sortNodesDeformed()
{
}


    void getDefFamily(const Node* node,
                      NodeArray& family);
  subroutine get_def_family(m)
 
    ! Finds the deformed family of node m.
    ! This is all the nodes currently within a distance rad_search of node m.
 
    ! input:
    !   m...                          node number (local)
    !   rad_search...                 search distance
    !   y1_lo,y1_hi,...               upper and lower coordinates of search region
 
    ! output:
    !   def_family...                 array containing the family of node m
    !   m_def_fam...                  number of nodes in the family of node m
 
    ! parameters
    int m
 
    ! local variables
    int cell1, cell2
    int bin_loc, bin_loc_start, bin_loc_end
    double dist_sq
    double dy1_cell, dy2_cell
    int search_lo_1, search_lo_2
    int search_hi_1, search_hi_2
    int search_cell_1, search_cell_2
    int search_node
    int i
    double yt1_search_node, yt2_search_node
    real(8), dimension(max_global_nodes) :: yt1,yt2
 
    ! initialize
    family = 0

    ! find the deformed position of node m.

    ! find the deformed node positions.
    do i=1, nnodes
      yt1(i) = nodes(i)%pos(1)+nodes(i)%disp(1)
      yt2(i) = nodes(i)%pos(2)+nodes(i)%disp(2)
    enddo

    ! set the cell limits and widths.
    dy1_cell = (y1_hi-y1_lo)/ncell_def1
    dy2_cell = (y2_hi-y2_lo)/ncell_def2
 
    ! find the cell number for node m.
    cell1 = 1+(yt1(m)-y1_lo)/dy1_cell
    cell2 = 1+(yt2(m)-y2_lo)/dy2_cell
 
    ! find the search distance for cells that might contain deformed family members.
    search_lo_1 = 1+(yt1(m)-rad_search-y1_lo)/dy1_cell
    search_hi_1 = 1+(yt1(m)+rad_search-y1_lo)/dy1_cell
    search_lo_2 = 1+(yt2(m)-rad_search-y2_lo)/dy2_cell
    search_hi_2 = 1+(yt2(m)+rad_search-y2_lo)/dy2_cell
 
    search_lo_1 = max(1, search_lo_1)
    search_lo_2 = max(1, search_lo_2)
    search_hi_1 = min(ncell_def1, search_hi_1)
    search_hi_2 = min(ncell_def2, search_hi_2)
 
    ! do the search.
    m_def_fam = 0
    do search_cell_1 = search_lo_1, search_hi_1
      do search_cell_2 = search_lo_2, search_hi_2
        if(number_of_nodes_in_cell_def(search_cell_1,search_cell_2)>0) then
          bin_loc_start = cell_location_in_bin_def(search_cell_1,search_cell_2)
          bin_loc_end = bin_loc_start + number_of_nodes_in_cell_def(search_cell_1,search_cell_2)-1
          do bin_loc = bin_loc_start, bin_loc_end
            search_node = nodes_in_bin_def(bin_loc)
            yt1_search_node = nodes(search_node)%pos(1)+nodes(search_node)%disp(1)
            yt2_search_node = nodes(search_node)%pos(2)+nodes(search_node)%disp(2)
            dist_sq = (yt1_search_node-yt1(m))**2 + (yt2_search_node-yt2(m))**2 
            if(m/=search_node .and. dist_sq < rad_search**2) then
              m_def_fam = m_def_fam+1
              if(m_def_fam <= maxfam) then
                def_family(m_def_fam) = search_node
              endif
            endif
          enddo
        endif
      enddo
    enddo
 
    ! abort if the family array is not big enough.
    if(m_def_fam > maxfam) then
      write(*,*) 'get_def_family: Too many nodes in family. maxfam=',maxfam, ' need ',m_def_fam
      stop
    endif
 
    return
 
  end subroutine get_def_family

  subroutine Sort_ref()
    ! Sorts nodes according to position in the reference configuration.

    ! A "cell" is one element in a 3D array representing a rectangular region of space.
    ! There are ncell1*ncell2*ncell3 cells. Each node is placed into one of these cells.
    ! Any number of nodes can occupy a cell. Cells can possibly contain no nodes.

    ! The node numbers in each cell are compressed into the "bin", which is a
    ! 1D array. No space is wasted in the bin, saving a lot of memory.

    ! input:
    !   xt1_global,...                node positions
    !   max_global_nodes...           dimension of the node position arrays
    !   total_global_nodes...         number of nodes in the node position arrays
    !   ncell1,ncell2,ncell3...       number of cells in each direction for sorting

    ! output:
    !   cell_location_in_bin...       pointer to the element in the bin corresponding to a cell
    !   number_of_nodes_in_cell...    number of nodes whose position in each cell
    !   nodes_in_bin...               node numbers sorted by cell compressed into the bin
    !   x1_lo,x1_hi,...               upper and lower coordinates of search region

    ! local variables
    int cell1, cell2
    int nodes_in_last_cell, last_bin, number_of_nodes_in_this_cell
    int m, i
    int bin_loc, bin_loc_start, bin_loc_rel
    double range
    double dx1_cell, dx2_cell
    ! double xt1(nnodes), xt2(nnodes)
 
    save msg1
    logical :: msg1
    data msg1/.false./
 
    ! initialize
 
    cell_location_in_bin = 0
    number_of_nodes_in_cell = 0
    nodes_in_bin = 0

    ! find the spatial limits of the nodes.
    do i = 1, nnodes
      xt1(i)=nodes(i)%pos(1)
      xt2(i)=nodes(i)%pos(2)
    enddo
    x1_lo = minval(xt1(1:nnodes))
    x1_hi = maxval(xt1(1:nnodes))
    x2_lo = minval(xt2(1:nnodes))
    x2_hi = maxval(xt2(1:nnodes))

    ! x1_lo = 1.d16
    ! x1_hi = -1.d16
    ! x2_lo = 1.d16
    ! x2_hi = -1.d16
    ! do i = 1, nnodes
    !   if (x1_lo .lt. nodes(i)%pos(1)) x1_lo = nodes(i)%pos(1)
    !   if (x1_hi .gt. nodes(i)%pos(1)) x1_hi = nodes(i)%pos(1)
    !   if (x2_lo .lt. nodes(i)%pos(2)) x2_lo = nodes(i)%pos(2)
    !   if (x2_hi .gt. nodes(i)%pos(2)) x2_hi = nodes(i)%pos(2)
    ! enddo
 
    ! adjust limits to avoid zero width in any direction and to avoid roundoff problems.
    range = dmax1(x1_hi-x1_lo, x2_hi-x2_lo)
    if(range==zero) then
      write(*,*) 'sort_ref: range=0.0 fatal'
      stop
    endif
    x1_lo = x1_lo-range*0.01d0
    x1_hi = x1_hi+range*0.01d0
    x2_lo = x2_lo-range*0.01d0
    x2_hi = x2_hi+range*0.01d0
 
    ! set the cell limits and widths.
    dx1_cell = (x1_hi-x1_lo)/ncell1
    dx2_cell = (x2_hi-x2_lo)/ncell2
 
    ! count the number of nodes in each cell.
    do m=1,nnodes
      cell1 = 1+(nodes(m)%pos(1)-x1_lo)/dx1_cell
      cell2 = 1+(nodes(m)%pos(2)-x2_lo)/dx2_cell
 
      if( (cell1<1.or.cell1>ncell1) .or. (cell2<1.or.cell2>ncell2)) then
        if (.not.msg1) then
          write(*,*) 'sort_ref: There are nodes outside the reference mesh.'
        end if
        msg1 = .true.
      else
        number_of_nodes_in_cell(cell1,cell2) = number_of_nodes_in_cell(cell1,cell2)+1
      endif
    end do
 
    ! find the starting location in the bin corresponding to each cell.
    nodes_in_last_cell = 0
    last_bin = 1
    do cell2=1,ncell2
      do cell1=1,ncell1
        cell_location_in_bin(cell1,cell2) = last_bin+nodes_in_last_cell
        last_bin = cell_location_in_bin(cell1,cell2)
        nodes_in_last_cell = number_of_nodes_in_cell(cell1,cell2)
      end do
    end do
 
    ! put the node numbers in the bin.
    do m=1,nnodes
      cell1 = 1+(nodes(m)%pos(1)-x1_lo)/dx1_cell
      cell2 = 1+(nodes(m)%pos(2)-x2_lo)/dx2_cell
 
      if( (cell1<1.or.cell1>ncell1) .or. (cell2<1.or.cell2>ncell2) ) then
      else
        bin_loc_start = cell_location_in_bin(cell1,cell2)
        number_of_nodes_in_this_cell = number_of_nodes_in_cell(cell1,cell2)
 
        do bin_loc_rel=0,number_of_nodes_in_this_cell-1
          bin_loc = bin_loc_start+bin_loc_rel
 
          if(nodes_in_bin(bin_loc) == 0) then
            nodes_in_bin(bin_loc) = m
            go to 11
          endif
        enddo
        write(*,*) 'sort_ref: no room in bin'
11      continue
      endif
    enddo

    return

  end subroutine Sort_ref

  subroutine sort_def()
 
    ! Sorts nodes according to position in the deformed configuration.
    ! This is done only for the nodes currently on each processor.
 
    ! A "cell" is one element in a 3D array representing a rectangular region of space.
    ! There are n_def_cell1*n_def_cell2*n_def_cell3 cells. Each node is placed into one of these cells.
    ! Any number of nodes can occupy a cell. Cells can possibly contain no nodes.
 
    ! The node numbers in each cell are compressed into the "bin", which is a
    ! 1D array. No space is wasted in the bin, saving a lot of memory.
 
    ! input:
    !   xt1,...                       node reference positions
    !   ut1,...                       node displacements
 
    ! output:
    !   cell_location_in_bin_def...   pointer to the element in the bin corresponding to a cell
    !   number_of_nodes_in_cell_def...number of nodes whose position in each cell
    !   nodes_in_bin_def...           node numbers sorted by cell compressed into the bin
    !   y1_lo,y1_hi,...               upper and lower coordinates of search region
 
    ! local variables
    int cell1, cell2
    int nodes_in_last_cell, last_bin, number_of_nodes_in_this_cell
    int m, i
    int bin_loc, bin_loc_start, bin_loc_rel
    double range
    double dy1_cell, dy2_cell

    ! double yt1(nnodes), yt2(nnodes)

    ! initialize
    cell_location_in_bin_def = 0
    number_of_nodes_in_cell_def = 0
    nodes_in_bin_def = 0

    ! find the deformed node positions.
    do i=1, nnodes
      yt1(i) = nodes(i)%pos(1)+nodes(i)%disp(1)
      yt2(i) = nodes(i)%pos(2)+nodes(i)%disp(2)
    enddo
 
    ! find the spatial limits of the deformed nodes.
    y1_lo = minval(yt1)
    y1_hi = maxval(yt1)
    y2_lo = minval(yt2)
    y2_hi = maxval(yt2)

    ! y1_lo = 1.d16
    ! y1_hi = -1.d16
    ! y2_lo = 1.d16
    ! y2_hi = -1.d16
    ! do i = 1, nnodes
    !   if (y1_lo .lt. nodes(i)%pos(1)+nodes(i)%disp(1)) x1_lo = nodes(i)%pos(1)+nodes(i)%disp(1)
    !   if (y1_hi .gt. nodes(i)%pos(1)+nodes(i)%disp(1)) x1_hi = nodes(i)%pos(1)+nodes(i)%disp(1)
    !   if (y2_lo .lt. nodes(i)%pos(2)+nodes(i)%disp(2)) x2_lo = nodes(i)%pos(2)+nodes(i)%disp(2)
    !   if (y2_hi .gt. nodes(i)%pos(2)+nodes(i)%disp(2)) x2_hi = nodes(i)%pos(2)+nodes(i)%disp(2)
    ! enddo
 
    ! adjust limits to avoid zero width in any direction and to avoid roundoff problems.
    range = dmax1(y1_hi-y1_lo, y2_hi-y2_lo)
    if(range==0.0) then
      range = 1.0
    endif
    y1_lo = y1_lo-range/100.0
    y1_hi = y1_hi+range/100.0
    y2_lo = y2_lo-range/100.0
    y2_hi = y2_hi+range/100.0
 
    ! set the cell limits and widths.
    dy1_cell = (y1_hi-y1_lo)/ncell_def1
    dy2_cell = (y2_hi-y2_lo)/ncell_def2
 
    ! count the number of nodes in each cell.
    do m=1,nnodes
      cell1 = 1+(yt1(m)-y1_lo)/dy1_cell
      cell2 = 1+(yt2(m)-y2_lo)/dy2_cell
 
      if( (cell1<1.or.cell1>ncell_def1) .or. (cell2<1.or.cell2>ncell_def2) ) then
        write(*,*) 'sort_def: There are nodes outside the deformed mesh.'
      else
        number_of_nodes_in_cell_def(cell1,cell2) = number_of_nodes_in_cell_def(cell1,cell2)+1
      endif
    end do
 
    ! find the starting location in the bin corresponding to each cell.
    nodes_in_last_cell = 0
    last_bin = 1
    do cell2=1,ncell_def2
      do cell1=1,ncell_def1
        cell_location_in_bin_def(cell1,cell2) = last_bin+nodes_in_last_cell
        last_bin = cell_location_in_bin_def(cell1,cell2)
        nodes_in_last_cell = number_of_nodes_in_cell_def(cell1,cell2)
      end do
    end do
 
    ! put the node numbers in the bin.
    do m=1,nnodes
      cell1 = 1+(yt1(m)-y1_lo)/dy1_cell
      cell2 = 1+(yt2(m)-y2_lo)/dy2_cell
 
      if( (cell1<1.or.cell1>ncell_def1) .or. (cell2<1.or.cell2>ncell_def2) ) then
      else
        bin_loc_start = cell_location_in_bin_def(cell1,cell2)
        number_of_nodes_in_this_cell = number_of_nodes_in_cell_def(cell1,cell2)
 
        do bin_loc_rel=0,number_of_nodes_in_this_cell-1
          bin_loc = bin_loc_start+bin_loc_rel
 
          if(nodes_in_bin_def(bin_loc) == 0) then
            nodes_in_bin_def(bin_loc) = m
            go to 11
          endif
        enddo
        write(*,*) 'sort_def: no room in bin'
11      continue
      endif
    enddo
 
    ! open(4,file='debug')
    ! write(4,*) 'cell_location_in_bin_def & number_of_nodes_in_cell_def:'
    !! write(*,*) 'cell_location_in_bin_def & number_of_nodes_in_cell_def:'
    ! do cell3=1,ncell_def3
    !   do cell2=1,ncell_def2
    !     do cell1=1,ncell_def1
    !       write(4,80000) cell1,cell2,cell3,cell_location_in_bin_def(cell1,cell2,cell3),number_of_nodes_in_cell_def(cell1,cell2,cell3)
    !!       write(*,80000) cell1,cell2,cell3,cell_location_in_bin_def(cell1,cell2,cell3),number_of_nodes_in_cell_def(cell1,cell2,cell3)
    !     enddo
    !   enddo
    ! enddo
    ! write(4,*) 'nodes sorted by position:'
    !! write(*,*) 'nodes sorted by position:'
    ! do cell3=1,ncell_def3
    !   do cell2=1,ncell_def2
    !     do cell1=1,ncell_def1
    !       write(4,80000) cell1,cell2,cell3
    !!      write(*,80000) cell1,cell2,cell3
    !       do bin_loc=cell_location_in_bin_def(cell1,cell2,cell3),  &
    !            cell_location_in_bin_def(cell1,cell2,cell3)+number_of_nodes_in_cell_def(cell1,cell2,cell3)-1
    !         m = nodes_in_bin_def(bin_loc)
    !         write(4,80001) m,yt1(m),yt2(m),yt3(m)
    !!        write(*,80001) m,yt1(m),yt2(m),yt3(m)
    !       enddo
    !     enddo
    !   enddo
    ! enddo
    ! close(4)
    !80000     format(3i5,2x,i8,2x,i8)
    !80001       format(i15,2x,3f12.6)
 
    return
 
  end subroutine sort_def

!****************************************************************************
end module bond_family_computations
!****************************************************************************

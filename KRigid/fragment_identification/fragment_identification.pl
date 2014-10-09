# Written by Kumar Mithraratne, University of Auckland




# 1. Grid geometry
##############################

$xcells = 10;                      # no of cells in x-direction
$ycells = 10;                      # no of cells in y-direction
$zcells = 10;                      # no of cells in z-direction
$xlength = 1.0;                    # length of side parallel to x                             
$ylength = 1.0;                    # length of side parallel to y   
$zlength = 1.0;                    # length of side parallel to z  

$delx = $xlength/$xcells;          # cell size in x
$dely = $ylength/$ycells;          # cell size in y
$delz = $zlength/$zcells;          # cell size in z
$total_cells = $xcells*$ycells*$zcells;
##############################

# $cell_type[$c] = 1 : corner cells - 8 corners
#                = 2 : edge cells   - 12 edges
#                = 3 : surface cells - 6 surfaces
#                = 4 : internal cells

for ($c = 1; $c < $total_cells + 1; $c++)
    {
     $cell_type[$c] = 4;
    }


# 2. Boundary cell information
##############################

# Corner cells (8 corners)
$corner[1] = 1;
$corner[2] = $xcells;
$corner[3] = $xcells*($ycells - 1) + 1;
$corner[4] = $xcells*$ycells;
$corner[5] = 1 + $xcells*$ycells*($zcells - 1);
$corner[6] = $xcells + $xcells*$ycells*($zcells - 1);
$corner[7] = $xcells*($ycells - 1) + 1 + $xcells*$ycells*($zcells - 1);
$corner[8] = $xcells*$ycells + $xcells*$ycells*($zcells - 1);
#for ($c = 1; $c < 9; $c++) {print "$corner[$c] ";} print "\n";

for ($cc = 1; $cc < 9; $cc++)
    {
     $cell_type[$corner[$cc]] = 1;
    }
   
# Edge cells (12 edges)

$cells_per_edge_x = $xcells - 2;
$cells_per_edge_y = $ycells - 2; 
$cells_per_edge_z = $zcells - 2; 

for ($c = 1; $c < $cells_per_edge_x + 1; $c++)
    {
     $edge[1][$c] = $c + 1;                                                                  
     $edge[2][$c] = $c + 1 + $xcells*($ycells - 1);                                          
     $edge[3][$c] = $c + 1 + $xcells*$ycells*($zcells - 1);                                  
     $edge[4][$c] = $c + 1 + $xcells*($ycells - 1) + $xcells*$ycells*($zcells - 1);          

     $cell_type[$edge[1][$c]] = 2;
     $cell_type[$edge[2][$c]] = 2;   
     $cell_type[$edge[3][$c]] = 2;
     $cell_type[$edge[4][$c]] = 2;             
    }
#for ($e = 1; $e < 5; $e++) {for ($c = 1; $c < $cells_per_edge_x + 1; $c++) {print "$edge[$e][$c] ";} print "\n";}

for ($c = 1; $c < $cells_per_edge_y + 1; $c++)
    {
     $edge[5][$c] = $c*$xcells + 1;                                                          
     $edge[6][$c] = $c*$xcells + $xcells;                                                    
     $edge[7][$c] = $c*$xcells + 1 + $xcells*$ycells*($zcells - 1);                          
     $edge[8][$c] = $c*$xcells + $xcells + $xcells*$ycells*($zcells - 1);   

     $cell_type[$edge[5][$c]] = 2;
     $cell_type[$edge[6][$c]] = 2;   
     $cell_type[$edge[7][$c]] = 2;
     $cell_type[$edge[8][$c]] = 2;             
    }
#for ($e = 5; $e < 9; $e++) {for ($c = 1; $c < $cells_per_edge_y + 1; $c++) {print "$edge[$e][$c] ";} print "\n";}

for ($c = 1; $c < $cells_per_edge_z + 1; $c++)
    {
     $edge[9][$c] = $c*$xcells*$ycells + 1;                                                  
     $edge[10][$c] = $c*$xcells*$ycells + $xcells;                                           
     $edge[11][$c] = $c*$xcells*$ycells + 1 + $xcells*($ycells - 1);                         
     $edge[12][$c] = $c*$xcells*$ycells + $xcells*$ycells; 

     $cell_type[$edge[9][$c]] = 2;
     $cell_type[$edge[10][$c]] = 2;   
     $cell_type[$edge[11][$c]] = 2;
     $cell_type[$edge[12][$c]] = 2;             
    }
#for ($e = 9; $e < 13; $e++) {for ($c = 1; $c < $cells_per_edge_z + 1; $c++) {print "$edge[$e][$c] ";} print "\n";}

# Surface cells (6 surfaces)

$c = 0;
for ($cz = 1; $cz < $zcells - 2 + 1; $cz++)
    {
     for ($cy = 1; $cy < $ycells - 2 + 1; $cy++)
         {
          $c++;
	   $surface[1][$c] = $cz*$xcells*$ycells + $cy*$xcells + 1;
	   $cell_type[$surface[1][$c]] = 3;
	  }
    }
#$c = 0; for ($cz = 1; $cz < $zcells - 2 + 1; $cz++) {for ($cy = 1; $cy < $ycells - 2 + 1; $cy++){$c++; print "$surface[1][$c] "} print "\n";}

$c = 0;
for ($cz = 1; $cz < $zcells - 2 + 1; $cz++)
    {
     for ($cy = 1; $cy < $ycells - 2 + 1; $cy++)
         {
          $c++;
	   $surface[2][$c] = $cz*$xcells*$ycells + $cy*$xcells + $xcells;
	   $cell_type[$surface[2][$c]] = 3;
	  }
    }
#$c = 0; for ($cz = 1; $cz < $zcells - 2 + 1; $cz++) {for ($cy = 1; $cy < $ycells - 2 + 1; $cy++){$c++; print "$surface[2][$c] "} print "\n";}

$c = 0;
for ($cx = 1; $cx < $xcells - 2 + 1; $cx++)
    {
     for ($cz = 1; $cz < $zcells - 2 + 1; $cz++)
         {
          $c++;
	   $surface[3][$c] = $cz*$xcells*$ycells + $cx + 1;
	   $cell_type[$surface[3][$c]] = 3;
	  }
    }
#$c = 0; for ($cx = 1; $cx < $xcells - 2 + 1; $cx++) {for ($cz = 1; $cz < $zcells - 2 + 1; $cz++){$c++; print "$surface[3][$c] "} print "\n";}

$c = 0;
for ($cx = 1; $cx < $xcells - 2 + 1; $cx++)
    {
     for ($cz = 1; $cz < $zcells - 2 + 1; $cz++)
         {
          $c++;
	   $surface[4][$c] = $cz*$xcells*$ycells + $cx + $xcells*($ycells - 1) + 1;
	   $cell_type[$surface[4][$c]] = 3;
	  }
    }
#$c = 0; for ($cx = 1; $cx < $xcells - 2 + 1; $cx++) {for ($cz = 1; $cz < $zcells - 2 + 1; $cz++){$c++; print "$surface[4][$c] "} print "\n";}

$c = 0;
for ($cy = 1; $cy < $ycells - 2 + 1; $cy++)
    {
     for ($cx = 1; $cx < $xcells - 2 + 1; $cx++)
         {
          $c++;
	   $surface[5][$c] = $cy*$xcells + $cx + 1;
	   $cell_type[$surface[5][$c]] = 3;
	  }
    }
#$c = 0; for ($cy = 1; $cy < $ycells - 2 + 1; $cy++) {for ($cx = 1; $cx < $xcells - 2 + 1; $cx++){$c++; print "$surface[5][$c] "} print "\n";}

$c = 0;
for ($cy = 1; $cy < $ycells - 2 + 1; $cy++)
    {
     for ($cx = 1; $cx < $xcells - 2 + 1; $cx++)
         {
          $c++;
	   $surface[6][$c] = $cy*$xcells + $cx + $xcells*$ycells*($zcells - 1) + 1;
	   $cell_type[$surface[6][$c]] = 3;
	  }
    }
#$c = 0; for ($cy = 1; $cy < $ycells - 2 + 1; $cy++) {for ($cx = 1; $cx < $xcells - 2 + 1; $cx++){$c++; print "$surface[6][$c] "} print "\n";}



# Interior cells 

$ic = 0;
for ($c = 1; $c < $total_cells + 1; $c++)
    {
     if ($cell_type[$c] == 4)  # other (interior) cells	
	 {
	  $ic++;
	  $interior[$ic] = $c;	  		  
        }	  
    }
$total_interior_cells = $ic;


# 2. Adjacent cell information
##############################

# Corner cells (each has 7 adjacent cells)

$boundary_cells_corner[1][1] =  $corner[1] + 1;
$boundary_cells_corner[1][2] =  $corner[1] + $xcells;
$boundary_cells_corner[1][3] =  $corner[1] + $xcells + 1;   
$boundary_cells_corner[1][4] =  $corner[1] + $xcells*$ycells;    
$boundary_cells_corner[1][5] =  $corner[1] + $xcells*$ycells + 1;        
$boundary_cells_corner[1][6] =  $corner[1] + $xcells*$ycells + $xcells;    
$boundary_cells_corner[1][7] =  $corner[1] + $xcells*$ycells + $xcells + 1;      
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[1][$c] ";} print "\n";
    
$boundary_cells_corner[2][1] =  $corner[2] - 1; 
$boundary_cells_corner[2][2] =  $corner[2] + $xcells;
$boundary_cells_corner[2][3] =  $corner[2] + $xcells - 1;  
$boundary_cells_corner[2][4] =  $corner[2] + $xcells*$ycells; 
$boundary_cells_corner[2][5] =  $corner[2] + $xcells*$ycells - 1; 
$boundary_cells_corner[2][6] =  $corner[2] + $xcells*$ycells + $xcells;  
$boundary_cells_corner[2][7] =  $corner[2] + $xcells*$ycells + $xcells - 1;       
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[2][$c] ";} print "\n";
    
$boundary_cells_corner[3][1] =  $corner[3] + 1; 
$boundary_cells_corner[3][2] =  $corner[3] - $xcells;  
$boundary_cells_corner[3][3] =  $corner[3] - $xcells + 1; 
$boundary_cells_corner[3][4] =  $corner[3] + $xcells*$ycells;   
$boundary_cells_corner[3][5] =  $corner[3] + $xcells*$ycells + 1;     
$boundary_cells_corner[3][6] =  $corner[3] + $xcells*$ycells - $xcells;   
$boundary_cells_corner[3][7] =  $corner[3] + $xcells*$ycells - $xcells + 1;        
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[3][$c] ";} print "\n"; 
    
$boundary_cells_corner[4][1] =  $corner[4] - 1; 
$boundary_cells_corner[4][2] =  $corner[4] - $xcells;  
$boundary_cells_corner[4][3] =  $corner[4] - $xcells - 1; 
$boundary_cells_corner[4][4] =  $corner[4] + $xcells*$ycells;   
$boundary_cells_corner[4][5] =  $corner[4] + $xcells*$ycells - 1;     
$boundary_cells_corner[4][6] =  $corner[4] + $xcells*$ycells - $xcells;   
$boundary_cells_corner[4][7] =  $corner[4] + $xcells*$ycells - $xcells - 1;       
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[4][$c] ";} print "\n";  

$boundary_cells_corner[5][1] =  $corner[5] + 1;
$boundary_cells_corner[5][2] =  $corner[5] + $xcells;
$boundary_cells_corner[5][3] =  $corner[5] + $xcells + 1;   
$boundary_cells_corner[5][4] =  $corner[5] - $xcells*$ycells;    
$boundary_cells_corner[5][5] =  $corner[5] - $xcells*$ycells + 1;        
$boundary_cells_corner[5][6] =  $corner[5] - $xcells*$ycells + $xcells;    
$boundary_cells_corner[5][7] =  $corner[5] - $xcells*$ycells + $xcells + 1;      
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[5][$c] ";} print "\n";
     
$boundary_cells_corner[6][1] =  $corner[6] - 1; 
$boundary_cells_corner[6][2] =  $corner[6] + $xcells;
$boundary_cells_corner[6][3] =  $corner[6] + $xcells - 1;  
$boundary_cells_corner[6][4] =  $corner[6] - $xcells*$ycells; 
$boundary_cells_corner[6][5] =  $corner[6] - $xcells*$ycells - 1; 
$boundary_cells_corner[6][6] =  $corner[6] - $xcells*$ycells + $xcells;  
$boundary_cells_corner[6][7] =  $corner[6] - $xcells*$ycells + $xcells - 1;       
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[6][$c] ";} print "\n";

$boundary_cells_corner[7][1] =  $corner[7] + 1; 
$boundary_cells_corner[7][2] =  $corner[7] - $xcells;  
$boundary_cells_corner[7][3] =  $corner[7] - $xcells + 1; 
$boundary_cells_corner[7][4] =  $corner[7] - $xcells*$ycells;   
$boundary_cells_corner[7][5] =  $corner[7] - $xcells*$ycells + 1;     
$boundary_cells_corner[7][6] =  $corner[7] - $xcells*$ycells - $xcells;   
$boundary_cells_corner[7][7] =  $corner[7] - $xcells*$ycells - $xcells + 1;       
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[7][$c] ";} print "\n";   

$boundary_cells_corner[8][1] =  $corner[8] - 1; 
$boundary_cells_corner[8][2] =  $corner[8] - $xcells;  
$boundary_cells_corner[8][3] =  $corner[8] - $xcells - 1; 
$boundary_cells_corner[8][4] =  $corner[8] - $xcells*$ycells;   
$boundary_cells_corner[8][5] =  $corner[8] - $xcells*$ycells - 1;     
$boundary_cells_corner[8][6] =  $corner[8] - $xcells*$ycells - $xcells;   
$boundary_cells_corner[8][7] =  $corner[8] - $xcells*$ycells - $xcells - 1;       
#for ($c = 1; $c < 8; $c++) {print "$boundary_cells_corner[8][$c] ";} print "\n";   

for ($cc = 1; $cc < 9; $cc++)
    {
     $cell_type[$corner[$c]] = 1;
     $no_of_neighbours[$corner[$cc]] = 7;
     for ($nc = 1; $nc < 8; $nc++)
         {
	   $neighbour_cells[$corner[$cc]][$nc] = $boundary_cells_corner[$cc][$nc];
         }
    }
 

# Edge cells (each has 11 adjacent cells)

for ($ec = 1; $ec < $cells_per_edge_x + 1; $ec++)
    {
     $boundary_cells_edge[1][1][$ec] =  $edge[1][$ec] - 1; 
     $boundary_cells_edge[1][2][$ec] =  $edge[1][$ec] + 1; 
     $boundary_cells_edge[1][3][$ec] =  $edge[1][$ec] + $xcells; 
     $boundary_cells_edge[1][4][$ec] =  $edge[1][$ec] + $xcells - 1; 
     $boundary_cells_edge[1][5][$ec] =  $edge[1][$ec] + $xcells + 1; 
     $boundary_cells_edge[1][6][$ec] =  $edge[1][$ec] + $xcells*$ycells;
     $boundary_cells_edge[1][7][$ec] =  $edge[1][$ec] + $xcells*$ycells - 1;
     $boundary_cells_edge[1][8][$ec] =  $edge[1][$ec] + $xcells*$ycells + 1;
     $boundary_cells_edge[1][9][$ec] =  $edge[1][$ec] + $xcells*$ycells + $xcells;
     $boundary_cells_edge[1][10][$ec] =  $edge[1][$ec] + $xcells*$ycells + $xcells - 1;
     $boundary_cells_edge[1][11][$ec] =  $edge[1][$ec] + $xcells*$ycells + $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[1][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_x + 1; $ec++)
    {
     $boundary_cells_edge[2][1][$ec] =  $edge[2][$ec] - 1; 
     $boundary_cells_edge[2][2][$ec] =  $edge[2][$ec] + 1; 
     $boundary_cells_edge[2][3][$ec] =  $edge[2][$ec] - $xcells; 
     $boundary_cells_edge[2][4][$ec] =  $edge[2][$ec] - $xcells - 1; 
     $boundary_cells_edge[2][5][$ec] =  $edge[2][$ec] - $xcells + 1; 
     $boundary_cells_edge[2][6][$ec] =  $edge[2][$ec] + $xcells*$ycells;
     $boundary_cells_edge[2][7][$ec] =  $edge[2][$ec] + $xcells*$ycells - 1;
     $boundary_cells_edge[2][8][$ec] =  $edge[2][$ec] + $xcells*$ycells + 1;
     $boundary_cells_edge[2][9][$ec] =  $edge[2][$ec] + $xcells*$ycells - $xcells;
     $boundary_cells_edge[2][10][$ec] =  $edge[2][$ec] + $xcells*$ycells - $xcells - 1;
     $boundary_cells_edge[2][11][$ec] =  $edge[2][$ec] + $xcells*$ycells - $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[2][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_x + 1; $ec++)
    {
     $boundary_cells_edge[3][1][$ec] =  $edge[3][$ec] - 1; 
     $boundary_cells_edge[3][2][$ec] =  $edge[3][$ec] + 1; 
     $boundary_cells_edge[3][3][$ec] =  $edge[3][$ec] + $xcells; 
     $boundary_cells_edge[3][4][$ec] =  $edge[3][$ec] + $xcells - 1; 
     $boundary_cells_edge[3][5][$ec] =  $edge[3][$ec] + $xcells + 1; 
     $boundary_cells_edge[3][6][$ec] =  $edge[3][$ec] - $xcells*$ycells;
     $boundary_cells_edge[3][7][$ec] =  $edge[3][$ec] - $xcells*$ycells - 1;
     $boundary_cells_edge[3][8][$ec] =  $edge[3][$ec] - $xcells*$ycells + 1;
     $boundary_cells_edge[3][9][$ec] =  $edge[3][$ec] - $xcells*$ycells + $xcells;
     $boundary_cells_edge[3][10][$ec] =  $edge[3][$ec] - $xcells*$ycells + $xcells - 1;
     $boundary_cells_edge[3][11][$ec] =  $edge[3][$ec] - $xcells*$ycells + $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[3][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_x + 1; $ec++)
    {
     $boundary_cells_edge[4][1][$ec] =  $edge[4][$ec] - 1; 
     $boundary_cells_edge[4][2][$ec] =  $edge[4][$ec] + 1; 
     $boundary_cells_edge[4][3][$ec] =  $edge[4][$ec] - $xcells; 
     $boundary_cells_edge[4][4][$ec] =  $edge[4][$ec] - $xcells - 1; 
     $boundary_cells_edge[4][5][$ec] =  $edge[4][$ec] - $xcells + 1; 
     $boundary_cells_edge[4][6][$ec] =  $edge[4][$ec] - $xcells*$ycells;
     $boundary_cells_edge[4][7][$ec] =  $edge[4][$ec] - $xcells*$ycells - 1;
     $boundary_cells_edge[4][8][$ec] =  $edge[4][$ec] - $xcells*$ycells + 1;
     $boundary_cells_edge[4][9][$ec] =  $edge[4][$ec] - $xcells*$ycells - $xcells;
     $boundary_cells_edge[4][10][$ec] =  $edge[4][$ec] - $xcells*$ycells - $xcells - 1;
     $boundary_cells_edge[4][11][$ec] =  $edge[4][$ec] - $xcells*$ycells - $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[4][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_y + 1; $ec++)
    {
     $boundary_cells_edge[5][1][$ec] =  $edge[5][$ec] + 1; 
     $boundary_cells_edge[5][2][$ec] =  $edge[5][$ec] - $xcells; 
     $boundary_cells_edge[5][3][$ec] =  $edge[5][$ec] + $xcells; 
     $boundary_cells_edge[5][4][$ec] =  $edge[5][$ec] - $xcells + 1; 
     $boundary_cells_edge[5][5][$ec] =  $edge[5][$ec] + $xcells + 1; 
     $boundary_cells_edge[5][6][$ec] =  $edge[5][$ec] + $xcells*$ycells;
     $boundary_cells_edge[5][7][$ec] =  $edge[5][$ec] + $xcells*$ycells + 1;
     $boundary_cells_edge[5][8][$ec] =  $edge[5][$ec] + $xcells*$ycells - $xcells;
     $boundary_cells_edge[5][9][$ec] =  $edge[5][$ec] + $xcells*$ycells + $xcells;
     $boundary_cells_edge[5][10][$ec] =  $edge[5][$ec] + $xcells*$ycells - $xcells + 1;
     $boundary_cells_edge[5][11][$ec] =  $edge[5][$ec] + $xcells*$ycells + $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[5][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_y + 1; $ec++)
    {
     $boundary_cells_edge[6][1][$ec] =  $edge[6][$ec] - 1; 
     $boundary_cells_edge[6][2][$ec] =  $edge[6][$ec] - $xcells; 
     $boundary_cells_edge[6][3][$ec] =  $edge[6][$ec] + $xcells; 
     $boundary_cells_edge[6][4][$ec] =  $edge[6][$ec] - $xcells - 1; 
     $boundary_cells_edge[6][5][$ec] =  $edge[6][$ec] + $xcells - 1; 
     $boundary_cells_edge[6][6][$ec] =  $edge[6][$ec] + $xcells*$ycells;
     $boundary_cells_edge[6][7][$ec] =  $edge[6][$ec] + $xcells*$ycells - 1;
     $boundary_cells_edge[6][8][$ec] =  $edge[6][$ec] + $xcells*$ycells - $xcells;
     $boundary_cells_edge[6][9][$ec] =  $edge[6][$ec] + $xcells*$ycells + $xcells;
     $boundary_cells_edge[6][10][$ec] =  $edge[6][$ec] + $xcells*$ycells - $xcells - 1;
     $boundary_cells_edge[6][11][$ec] =  $edge[6][$ec] + $xcells*$ycells + $xcells - 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[6][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_y + 1; $ec++)
    {
     $boundary_cells_edge[7][1][$ec] =  $edge[7][$ec] + 1; 
     $boundary_cells_edge[7][2][$ec] =  $edge[7][$ec] - $xcells; 
     $boundary_cells_edge[7][3][$ec] =  $edge[7][$ec] + $xcells; 
     $boundary_cells_edge[7][4][$ec] =  $edge[7][$ec] - $xcells + 1; 
     $boundary_cells_edge[7][5][$ec] =  $edge[7][$ec] + $xcells + 1; 
     $boundary_cells_edge[7][6][$ec] =  $edge[7][$ec] - $xcells*$ycells;
     $boundary_cells_edge[7][7][$ec] =  $edge[7][$ec] - $xcells*$ycells + 1;
     $boundary_cells_edge[7][8][$ec] =  $edge[7][$ec] - $xcells*$ycells - $xcells;
     $boundary_cells_edge[7][9][$ec] =  $edge[7][$ec] - $xcells*$ycells + $xcells;
     $boundary_cells_edge[7][10][$ec] =  $edge[7][$ec] - $xcells*$ycells - $xcells + 1;
     $boundary_cells_edge[7][11][$ec] =  $edge[7][$ec] - $xcells*$ycells + $xcells + 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[7][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_y + 1; $ec++)
    {
     $boundary_cells_edge[8][1][$ec] =  $edge[8][$ec] - 1; 
     $boundary_cells_edge[8][2][$ec] =  $edge[8][$ec] - $xcells; 
     $boundary_cells_edge[8][3][$ec] =  $edge[8][$ec] + $xcells; 
     $boundary_cells_edge[8][4][$ec] =  $edge[8][$ec] - $xcells - 1; 
     $boundary_cells_edge[8][5][$ec] =  $edge[8][$ec] + $xcells - 1; 
     $boundary_cells_edge[8][6][$ec] =  $edge[8][$ec] - $xcells*$ycells;
     $boundary_cells_edge[8][7][$ec] =  $edge[8][$ec] - $xcells*$ycells - 1;
     $boundary_cells_edge[8][8][$ec] =  $edge[8][$ec] - $xcells*$ycells - $xcells;
     $boundary_cells_edge[8][9][$ec] =  $edge[8][$ec] - $xcells*$ycells + $xcells;
     $boundary_cells_edge[8][10][$ec] =  $edge[8][$ec] - $xcells*$ycells - $xcells - 1;
     $boundary_cells_edge[8][11][$ec] =  $edge[8][$ec] - $xcells*$ycells + $xcells - 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[8][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_z + 1; $ec++)
    {
     $boundary_cells_edge[9][1][$ec] =  $edge[9][$ec] + 1; 
     $boundary_cells_edge[9][2][$ec] =  $edge[9][$ec] - $xcells*$ycells; 
     $boundary_cells_edge[9][3][$ec] =  $edge[9][$ec] + $xcells*$ycells; 
     $boundary_cells_edge[9][4][$ec] =  $edge[9][$ec] - $xcells*$ycells + 1; 
     $boundary_cells_edge[9][5][$ec] =  $edge[9][$ec] + $xcells*$ycells + 1;      
     $boundary_cells_edge[9][6][$ec] =  $edge[9][$ec] + $xcells;
     $boundary_cells_edge[9][7][$ec] =  $edge[9][$ec] + 1 + $xcells;
     $boundary_cells_edge[9][8][$ec] =  $edge[9][$ec] - $xcells*$ycells + $xcells;
     $boundary_cells_edge[9][9][$ec] =  $edge[9][$ec] + $xcells*$ycells + $xcells;
     $boundary_cells_edge[9][10][$ec] =  $edge[9][$ec] - $xcells*$ycells + 1 + $xcells;
     $boundary_cells_edge[9][11][$ec] =  $edge[9][$ec] + $xcells*$ycells + 1 + $xcells;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[9][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_z + 1; $ec++)
    {
     $boundary_cells_edge[10][1][$ec] =  $edge[10][$ec] - 1; 
     $boundary_cells_edge[10][2][$ec] =  $edge[10][$ec] - $xcells*$ycells; 
     $boundary_cells_edge[10][3][$ec] =  $edge[10][$ec] + $xcells*$ycells; 
     $boundary_cells_edge[10][4][$ec] =  $edge[10][$ec] - $xcells*$ycells - 1; 
     $boundary_cells_edge[10][5][$ec] =  $edge[10][$ec] + $xcells*$ycells - 1;      
     $boundary_cells_edge[10][6][$ec] =  $edge[10][$ec] + $xcells;
     $boundary_cells_edge[10][7][$ec] =  $edge[10][$ec] + $xcells - 1;
     $boundary_cells_edge[10][8][$ec] =  $edge[10][$ec] - $xcells*$ycells + $xcells;
     $boundary_cells_edge[10][9][$ec] =  $edge[10][$ec] + $xcells*$ycells + $xcells;
     $boundary_cells_edge[10][10][$ec] =  $edge[10][$ec] - $xcells*$ycells + $xcells - 1;
     $boundary_cells_edge[10][11][$ec] =  $edge[10][$ec] + $xcells*$ycells + $xcells - 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[10][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_z + 1; $ec++)
    {
     $boundary_cells_edge[11][1][$ec] =  $edge[11][$ec] + 1; 
     $boundary_cells_edge[11][2][$ec] =  $edge[11][$ec] - $xcells*$ycells; 
     $boundary_cells_edge[11][3][$ec] =  $edge[11][$ec] + $xcells*$ycells; 
     $boundary_cells_edge[11][4][$ec] =  $edge[11][$ec] - $xcells*$ycells + 1; 
     $boundary_cells_edge[11][5][$ec] =  $edge[11][$ec] + $xcells*$ycells + 1;      
     $boundary_cells_edge[11][6][$ec] =  $edge[11][$ec] - $xcells;
     $boundary_cells_edge[11][7][$ec] =  $edge[11][$ec] + 1 - $xcells;
     $boundary_cells_edge[11][8][$ec] =  $edge[11][$ec] - $xcells*$ycells - $xcells;
     $boundary_cells_edge[11][9][$ec] =  $edge[11][$ec] + $xcells*$ycells - $xcells;
     $boundary_cells_edge[11][10][$ec] =  $edge[11][$ec] - $xcells*$ycells + 1 - $xcells;
     $boundary_cells_edge[11][11][$ec] =  $edge[11][$ec] + $xcells*$ycells + 1 - $xcells;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[11][$c][$ec] ";} print "\n";  
    }

for ($ec = 1; $ec < $cells_per_edge_z + 1; $ec++)
    {
     $boundary_cells_edge[12][1][$ec] =  $edge[12][$ec] - 1; 
     $boundary_cells_edge[12][2][$ec] =  $edge[12][$ec] - $xcells*$ycells; 
     $boundary_cells_edge[12][3][$ec] =  $edge[12][$ec] + $xcells*$ycells; 
     $boundary_cells_edge[12][4][$ec] =  $edge[12][$ec] - $xcells*$ycells - 1; 
     $boundary_cells_edge[12][5][$ec] =  $edge[12][$ec] + $xcells*$ycells - 1;      
     $boundary_cells_edge[12][6][$ec] =  $edge[12][$ec] - $xcells;
     $boundary_cells_edge[12][7][$ec] =  $edge[12][$ec] - $xcells - 1;
     $boundary_cells_edge[12][8][$ec] =  $edge[12][$ec] - $xcells*$ycells - $xcells;
     $boundary_cells_edge[12][9][$ec] =  $edge[12][$ec] + $xcells*$ycells - $xcells;
     $boundary_cells_edge[12][10][$ec] =  $edge[12][$ec] - $xcells*$ycells - $xcells - 1;
     $boundary_cells_edge[12][11][$ec] =  $edge[12][$ec] + $xcells*$ycells - $xcells - 1;
     #for ($c = 1; $c < 12; $c++) {print "$boundary_cells_edge[12][$c][$ec] ";} print "\n";  
    }


for ($e = 1; $e < 5; $e++)
    {
     for ($ec = 1; $ec < $cells_per_edge_x + 1; $ec++)
         {
	   $cell_type[$edge[$e][$ec]] = 2;
	   $no_of_neighbours[$edge[$e][$ec]] = 11;
	   for ($nc = 1; $nc < 12; $nc++)
	       {
	        $neighbour_cells[$edge[$e][$ec]][$nc] = $boundary_cells_edge[$e][$nc][$ec]; 
	       } 
	  } 
    }

for ($e = 5; $e < 9; $e++)
    {
     for ($ec = 1; $ec < $cells_per_edge_y + 1; $ec++)
         {
	   $cell_type[$edge[$e][$ec]] = 2;
	   $no_of_neighbours[$edge[$e][$ec]] = 11;
	   for ($nc = 1; $nc < 12; $nc++)
	       {
	        $neighbour_cells[$edge[$e][$ec]][$nc] = $boundary_cells_edge[$e][$nc][$ec]; 
	       } 
	  } 
    }

for ($e = 9; $e < 13; $e++)
    {
     for ($ec = 1; $ec < $cells_per_edge_z + 1; $ec++)
         {
	   $cell_type[$edge[$e][$ec]] = 2;
	   $no_of_neighbours[$edge[$e][$ec]] = 11;
	   for ($nc = 1; $nc < 12; $nc++)
	       {
	        $neighbour_cells[$edge[$e][$ec]][$nc] = $boundary_cells_edge[$e][$nc][$ec]; 
	       } 
	  } 
    }


# Surface cells (each has 17 adjacent cells)

for ($sc = 1; $sc < ($ycells - 2)*($zcells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[1][1][$sc] = $surface[1][$sc] + $xcells;
     $boundary_cells_surface[1][2][$sc] = $surface[1][$sc] - $xcells;  
     $boundary_cells_surface[1][3][$sc] = $surface[1][$sc] + $xcells*$ycells;  
     $boundary_cells_surface[1][4][$sc] = $surface[1][$sc] - $xcells*$ycells;  
     $boundary_cells_surface[1][5][$sc] = $surface[1][$sc] + $xcells*$ycells + $xcells;  
     $boundary_cells_surface[1][6][$sc] = $surface[1][$sc] + $xcells*$ycells - $xcells;    
     $boundary_cells_surface[1][7][$sc] = $surface[1][$sc] - $xcells*$ycells + $xcells;  
     $boundary_cells_surface[1][8][$sc] = $surface[1][$sc] - $xcells*$ycells - $xcells;   
     $boundary_cells_surface[1][9][$sc] = $surface[1][$sc] + 1; 
     $boundary_cells_surface[1][10][$sc] = $surface[1][$sc] + $xcells + 1;
     $boundary_cells_surface[1][11][$sc] = $surface[1][$sc] - $xcells + 1;  
     $boundary_cells_surface[1][12][$sc] = $surface[1][$sc] + $xcells*$ycells + 1;  
     $boundary_cells_surface[1][13][$sc] = $surface[1][$sc] - $xcells*$ycells + 1;  
     $boundary_cells_surface[1][14][$sc] = $surface[1][$sc] + $xcells*$ycells + $xcells + 1; 
     $boundary_cells_surface[1][15][$sc] = $surface[1][$sc] + $xcells*$ycells - $xcells + 1; 
     $boundary_cells_surface[1][16][$sc] = $surface[1][$sc] - $xcells*$ycells + $xcells + 1; 
     $boundary_cells_surface[1][17][$sc] = $surface[1][$sc] - $xcells*$ycells - $xcells + 1; 
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[1][$c][$sc] ";} print "\n";   
    }

for ($sc = 1; $sc < ($ycells - 2)*($zcells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[2][1][$sc] = $surface[2][$sc] + $xcells;
     $boundary_cells_surface[2][2][$sc] = $surface[2][$sc] - $xcells;  
     $boundary_cells_surface[2][3][$sc] = $surface[2][$sc] + $xcells*$ycells;  
     $boundary_cells_surface[2][4][$sc] = $surface[2][$sc] - $xcells*$ycells;  
     $boundary_cells_surface[2][5][$sc] = $surface[2][$sc] + $xcells*$ycells + $xcells;  
     $boundary_cells_surface[2][6][$sc] = $surface[2][$sc] + $xcells*$ycells - $xcells;    
     $boundary_cells_surface[2][7][$sc] = $surface[2][$sc] - $xcells*$ycells + $xcells;  
     $boundary_cells_surface[2][8][$sc] = $surface[2][$sc] - $xcells*$ycells - $xcells;   
     $boundary_cells_surface[2][9][$sc] = $surface[2][$sc] - 1; 
     $boundary_cells_surface[2][10][$sc] = $surface[2][$sc] + $xcells - 1;
     $boundary_cells_surface[2][11][$sc] = $surface[2][$sc] - $xcells - 1;  
     $boundary_cells_surface[2][12][$sc] = $surface[2][$sc] + $xcells*$ycells - 1;  
     $boundary_cells_surface[2][13][$sc] = $surface[2][$sc] - $xcells*$ycells - 1;  
     $boundary_cells_surface[2][14][$sc] = $surface[2][$sc] + $xcells*$ycells + $xcells - 1; 
     $boundary_cells_surface[2][15][$sc] = $surface[2][$sc] + $xcells*$ycells - $xcells - 1; 
     $boundary_cells_surface[2][16][$sc] = $surface[2][$sc] - $xcells*$ycells + $xcells - 1; 
     $boundary_cells_surface[2][17][$sc] = $surface[2][$sc] - $xcells*$ycells - $xcells - 1; 
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[2][$c][$sc] ";} print "\n";   
    }

for ($sc = 1; $sc < ($zcells - 2)*($xcells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[3][1][$sc] = $surface[3][$sc] + 1;
     $boundary_cells_surface[3][2][$sc] = $surface[3][$sc] - 1;  
     $boundary_cells_surface[3][3][$sc] = $surface[3][$sc] + $xcells*$ycells;  
     $boundary_cells_surface[3][4][$sc] = $surface[3][$sc] - $xcells*$ycells;  
     $boundary_cells_surface[3][5][$sc] = $surface[3][$sc] + $xcells*$ycells + 1;  
     $boundary_cells_surface[3][6][$sc] = $surface[3][$sc] + $xcells*$ycells - 1;    
     $boundary_cells_surface[3][7][$sc] = $surface[3][$sc] - $xcells*$ycells + 1;  
     $boundary_cells_surface[3][8][$sc] = $surface[3][$sc] - $xcells*$ycells - 1;   
     $boundary_cells_surface[3][9][$sc] = $surface[3][$sc] + $xcells; 
     $boundary_cells_surface[3][10][$sc] = $surface[3][$sc] + 1 + $xcells;
     $boundary_cells_surface[3][11][$sc] = $surface[3][$sc] - 1 + $xcells;  
     $boundary_cells_surface[3][12][$sc] = $surface[3][$sc] + $xcells*$ycells + $xcells;  
     $boundary_cells_surface[3][13][$sc] = $surface[3][$sc] - $xcells*$ycells + $xcells;  
     $boundary_cells_surface[3][14][$sc] = $surface[3][$sc] + $xcells*$ycells + 1 + $xcells;  
     $boundary_cells_surface[3][15][$sc] = $surface[3][$sc] + $xcells*$ycells - 1 + $xcells;  
     $boundary_cells_surface[3][16][$sc] = $surface[3][$sc] - $xcells*$ycells + 1 + $xcells;  
     $boundary_cells_surface[3][17][$sc] = $surface[3][$sc] - $xcells*$ycells - 1 + $xcells;  
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[3][$c][$sc] ";} print "\n";   
    }

for ($sc = 1; $sc < ($zcells - 2)*($xcells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[4][1][$sc] = $surface[4][$sc] + 1;
     $boundary_cells_surface[4][2][$sc] = $surface[4][$sc] - 1;  
     $boundary_cells_surface[4][3][$sc] = $surface[4][$sc] + $xcells*$ycells;  
     $boundary_cells_surface[4][4][$sc] = $surface[4][$sc] - $xcells*$ycells;  
     $boundary_cells_surface[4][5][$sc] = $surface[4][$sc] + $xcells*$ycells + 1;  
     $boundary_cells_surface[4][6][$sc] = $surface[4][$sc] + $xcells*$ycells - 1;    
     $boundary_cells_surface[4][7][$sc] = $surface[4][$sc] - $xcells*$ycells + 1;  
     $boundary_cells_surface[4][8][$sc] = $surface[4][$sc] - $xcells*$ycells - 1;   
     $boundary_cells_surface[4][9][$sc] = $surface[4][$sc] - $xcells; 
     $boundary_cells_surface[4][10][$sc] = $surface[4][$sc] + 1 - $xcells;
     $boundary_cells_surface[4][11][$sc] = $surface[4][$sc] - 1 - $xcells;  
     $boundary_cells_surface[4][12][$sc] = $surface[4][$sc] + $xcells*$ycells - $xcells;  
     $boundary_cells_surface[4][13][$sc] = $surface[4][$sc] - $xcells*$ycells - $xcells;  
     $boundary_cells_surface[4][14][$sc] = $surface[4][$sc] + $xcells*$ycells + 1 - $xcells;  
     $boundary_cells_surface[4][15][$sc] = $surface[4][$sc] + $xcells*$ycells - 1 - $xcells;  
     $boundary_cells_surface[4][16][$sc] = $surface[4][$sc] - $xcells*$ycells + 1 - $xcells;  
     $boundary_cells_surface[4][17][$sc] = $surface[4][$sc] - $xcells*$ycells - 1 - $xcells;  
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[4][$c][$sc] ";} print "\n";   
    }

for ($sc = 1; $sc < ($xcells - 2)*($ycells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[5][1][$sc] = $surface[5][$sc] + 1;
     $boundary_cells_surface[5][2][$sc] = $surface[5][$sc] - 1;
     $boundary_cells_surface[5][3][$sc] = $surface[5][$sc] + $xcells; 
     $boundary_cells_surface[5][4][$sc] = $surface[5][$sc] - $xcells;  
     $boundary_cells_surface[5][5][$sc] = $surface[5][$sc] + 1 + $xcells;
     $boundary_cells_surface[5][6][$sc] = $surface[5][$sc] + 1 - $xcells;
     $boundary_cells_surface[5][7][$sc] = $surface[5][$sc] - 1 + $xcells;
     $boundary_cells_surface[5][8][$sc] = $surface[5][$sc] - 1 - $xcells;
     $boundary_cells_surface[5][9][$sc] = $surface[5][$sc] + $xcells*$ycells;
     $boundary_cells_surface[5][10][$sc] = $surface[5][$sc] + 1 + $xcells*$ycells;
     $boundary_cells_surface[5][11][$sc] = $surface[5][$sc] - 1 + $xcells*$ycells;
     $boundary_cells_surface[5][12][$sc] = $surface[5][$sc] + $xcells + $xcells*$ycells; 
     $boundary_cells_surface[5][13][$sc] = $surface[5][$sc] - $xcells + $xcells*$ycells;  
     $boundary_cells_surface[5][14][$sc] = $surface[5][$sc] + 1 + $xcells + $xcells*$ycells;
     $boundary_cells_surface[5][15][$sc] = $surface[5][$sc] + 1 - $xcells + $xcells*$ycells;
     $boundary_cells_surface[5][16][$sc] = $surface[5][$sc] - 1 + $xcells + $xcells*$ycells;
     $boundary_cells_surface[5][17][$sc] = $surface[5][$sc] - 1 - $xcells + $xcells*$ycells;
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[5][$c][$sc] ";} print "\n";   
    }

for ($sc = 1; $sc < ($xcells - 2)*($ycells - 2) + 1; $sc++)
    {
     $boundary_cells_surface[6][1][$sc] = $surface[6][$sc] + 1;
     $boundary_cells_surface[6][2][$sc] = $surface[6][$sc] - 1;
     $boundary_cells_surface[6][3][$sc] = $surface[6][$sc] + $xcells; 
     $boundary_cells_surface[6][4][$sc] = $surface[6][$sc] - $xcells;  
     $boundary_cells_surface[6][5][$sc] = $surface[6][$sc] + 1 + $xcells;
     $boundary_cells_surface[6][6][$sc] = $surface[6][$sc] + 1 - $xcells;
     $boundary_cells_surface[6][7][$sc] = $surface[6][$sc] - 1 + $xcells;
     $boundary_cells_surface[6][8][$sc] = $surface[6][$sc] - 1 - $xcells;
     $boundary_cells_surface[6][9][$sc] = $surface[6][$sc] - $xcells*$ycells;
     $boundary_cells_surface[6][10][$sc] = $surface[6][$sc] + 1 - $xcells*$ycells;
     $boundary_cells_surface[6][11][$sc] = $surface[6][$sc] - 1 - $xcells*$ycells;
     $boundary_cells_surface[6][12][$sc] = $surface[6][$sc] + $xcells - $xcells*$ycells; 
     $boundary_cells_surface[6][13][$sc] = $surface[6][$sc] - $xcells - $xcells*$ycells;  
     $boundary_cells_surface[6][14][$sc] = $surface[6][$sc] + 1 + $xcells - $xcells*$ycells;
     $boundary_cells_surface[6][15][$sc] = $surface[6][$sc] + 1 - $xcells - $xcells*$ycells;
     $boundary_cells_surface[6][16][$sc] = $surface[6][$sc] - 1 + $xcells - $xcells*$ycells;
     $boundary_cells_surface[6][17][$sc] = $surface[6][$sc] - 1 - $xcells - $xcells*$ycells;
     #for ($c = 1; $c < 18; $c++) {print "$boundary_cells_surface[6][$c][$sc] ";} print "\n";   
    }

for ($s = 1; $s < 3; $s++)
    {
     for ($sc = 1; $sc < ($ycells - 2)*($zcells - 2) + 1; $sc++)
         {
	   $cell_type[$surface[$s][$sc]] = 3;
	   $no_of_neighbours[$surface[$s][$sc]] = 17;
	   for ($nc = 1; $nc < 18; $nc++)
	       {
	        $neighbour_cells[$surface[$s][$sc]][$nc] = $boundary_cells_surface[$s][$nc][$sc]; 
	       } 
	  } 
    }
for ($s = 3; $s < 5; $s++)
    {
     for ($sc = 1; $sc < ($zcells - 2)*($xcells - 2) + 1; $sc++)
         {
	   $cell_type[$surface[$s][$sc]] = 3;
	   $no_of_neighbours[$surface[$s][$sc]] = 17;
	   for ($nc = 1; $nc < 18; $nc++)
	       {
	        $neighbour_cells[$surface[$s][$sc]][$nc] = $boundary_cells_surface[$s][$nc][$sc]; 
	       } 
	  } 
    }
for ($s = 5; $s < 7; $s++)
    {
     for ($sc = 1; $sc < ($xcells - 2)*($ycells - 2) + 1; $sc++)
         {
	   $cell_type[$surface[$s][$sc]] = 3;
	   $no_of_neighbours[$surface[$s][$sc]] = 17;
	   for ($nc = 1; $nc < 18; $nc++)
	       {
	        $neighbour_cells[$surface[$s][$sc]][$nc] = $boundary_cells_surface[$s][$nc][$sc]; 
	       } 
	  } 
    }


# Other (interior) cells (each has 26 adjacent cells)

for ($ic = 1; $ic < $total_interior_cells + 1; $ic++)
    {
     $interior_cells[1][$ic] = $interior[$ic] + 1;	  
     $interior_cells[2][$ic] = $interior[$ic] - 1;	  
     $interior_cells[3][$ic] = $interior[$ic] + $xcells;
     $interior_cells[4][$ic] = $interior[$ic] - $xcells;
     $interior_cells[5][$ic] = $interior[$ic] + $xcells + 1;
     $interior_cells[6][$ic] = $interior[$ic] + $xcells - 1;
     $interior_cells[7][$ic] = $interior[$ic] - $xcells + 1;
     $interior_cells[8][$ic] =  $interior[$ic] - $xcells - 1;
     $interior_cells[9][$ic] = $interior[$ic] + $xcells*$ycells;  
     $interior_cells[10][$ic] = $interior[$ic] + 1 + $xcells*$ycells;	   
     $interior_cells[11][$ic] = $interior[$ic] - 1 + $xcells*$ycells;	   
     $interior_cells[12][$ic] = $interior[$ic] + $xcells + $xcells*$ycells;
     $interior_cells[13][$ic] = $interior[$ic] - $xcells + $xcells*$ycells;
     $interior_cells[14][$ic] = $interior[$ic] + $xcells + 1 + $xcells*$ycells;
     $interior_cells[15][$ic] = $interior[$ic] + $xcells - 1 + $xcells*$ycells;
     $interior_cells[16][$ic] = $interior[$ic] - $xcells + 1 + $xcells*$ycells;
     $interior_cells[17][$ic] = $interior[$ic] - $xcells - 1 + $xcells*$ycells;
     $interior_cells[18][$ic] = $interior[$ic] - $xcells*$ycells;
     $interior_cells[19][$ic] = $interior[$ic] + 1 - $xcells*$ycells;    
     $interior_cells[20][$ic] = $interior[$ic] - 1 - $xcells*$ycells;    
     $interior_cells[21][$ic] = $interior[$ic] + $xcells - $xcells*$ycells;
     $interior_cells[22][$ic] = $interior[$ic] - $xcells - $xcells*$ycells;
     $interior_cells[23][$ic] = $interior[$ic] + $xcells + 1 - $xcells*$ycells;
     $interior_cells[24][$ic] = $interior[$ic] + $xcells - 1 - $xcells*$ycells;
     $interior_cells[25][$ic] = $interior[$ic] - $xcells + 1 - $xcells*$ycells;
     $interior_cells[26][$ic] = $interior[$ic] - $xcells - 1 - $xcells*$ycells;		
     #for ($c = 1; $c < 27; $c++) {print "$interior_cells[$c][$ic] ";} print "\n";   
    }	  

for ($ic = 1; $ic < $total_interior_cells + 1; $ic++)
    {
     $cell_type[$interior[$ic]] = 4;
     $no_of_neighbours[$interior[$ic]] = 26;
     for ($nc = 1; $nc < 27; $nc++)
  	  {
  	   $neighbour_cells[$interior[$ic]][$nc] = $interior_cells[$nc][$ic]; 
  	  } 
    } 

##############################

    
# Read spatial coordinates of particles and determine what cell each particle belongs to
# ======================================================================================

for ($c = 1; $c < $total_cells + 1; $c++)
    {
     $cell_has_particles[$c] = 0;
    }

open (PARTICLE, "particle.dat") || die (" NO particle.dat EXISTS");

@file_data=<PARTICLE>;                              #print "@file_data";
$total_lines = $#file_data;                     #print "$total_lines\n";

$p = 0;
for ($line = 0; $line < $total_lines + 1; $line = $line + 1)
#for ($line = 0; $line < 5; $line = $line + 1)    
    { 
     $p++;
     chomp $file_data[$line];
     @words = split(/\s+/,$file_data[$line]);    #print "$words[1]\n";
     
     @words1 = split(/,/,$words[1]);
     
     $x[$p] = $words1[0]; $y[$p] = $words1[1]; $z[$p] = $words1[2];   
     
     $cellsinz = int($z[$p]/$delz);
     $cellsiny = int($y[$p]/$dely);
     $cellsinx = int($x[$p]/$delx);
     
     # cell no of the data point
     $cell[$p] = $cellsinz*$xcells*$ycells + $cellsiny*$xcells + $cellsinx + 1;     
     
     $cell_has_particles[$cell[$p]] = 1;
    }

close (PARTICLE);

$total_particles = $p;
#print " total data points = $total_data\n";   

for ($c = 1; $c < $total_cells + 1; $c++)
    {
     $particles_in_cell[$c] = 0;
    } 

for ($p = 1; $p < $total_particles + 1; $p++)
    {
     $c = $cell[$p];
     $particles_in_cell[$c]++;
     $cell_particle_id[$c][$particles_in_cell[$c]] = $p;
    }


#for ($c = 1; $c < $total_cells + 1; $c++)
#    {
#     print " $c - $particles_in_cell[$c]\n";
#    } 
#$cell = 766;
#for ($pc = 1; $pc < $particles_in_cell[$cell] + 1; $pc++)
#    {
#     print "$pc - $cell_particle_id[$cell][$pc]\n";
#    }

###########################################################################

    
# Fragment identifcation
# ======================


$fragment_no = 0;
for ($c = 1; $c < $total_cells + 1; $c++)
    {
     $cell_included_in_a_fragment[$c] = 0;
     $cell_fragment_number[$c] = 0;
    }
    

for ($c = 1; $c < $total_cells + 1; $c++)
#for ($c = 1; $c < 501 + 1; $c++)    
    {
     
     #print " $c : $cell_type[$c] : $no_of_neighbours[$c] : "; 
     #for ($nc = 1; $nc < $no_of_neighbours[$c] + 1; $nc++) 
     #    {
     #	   print "$neighbour_cells[$c][$nc],";
     #	  }
     #print "\n"; 
	  
     
     if ($cell_has_particles[$c] == 1)
        {
	  if ($cell_included_in_a_fragment[$c] == 1)
	     {
	      #search neighbours who do not belong to a fragment
             for ($nc = 1; $nc < $no_of_neighbours[$c] + 1; $nc++) 
                 {
	           $neighbour = $neighbour_cells[$c][$nc];
		    if ($cell_has_particles[$neighbour] == 1)
		       {
			 if ($cell_included_in_a_fragment[$neighbour] == 0)
			    {
			     $cell_fragment_number[$neighbour] = $cell_fragment_number[$c];
			     $cell_included_in_a_fragment[$neighbour] = 1;			    
			    }
			 
			}
                 }                
	     }
	  elsif ($cell_included_in_a_fragment[$c] == 0)
	     {
	      #search neighbours to see whether any of them belong to a fragment
             for ($nc = 1; $nc < $no_of_neighbours[$c] + 1; $nc++) 
                 {
	           $neighbour = $neighbour_cells[$c][$nc];    #print "$nc - $neighbour\n"; 

		    if ($cell_has_particles[$neighbour] == 1)
		       {
			 if ($cell_included_in_a_fragment[$neighbour] == 1)
			    {
			     $cell_fragment_number[$c] = $cell_fragment_number[$neighbour];
			     $cell_included_in_a_fragment[$c] = 1;
			     #print "$neighbour - $cell_fragment_number[$neighbour]\n";
			    } 
			 elsif ($cell_included_in_a_fragment[$neighbour] == 0)
			    {   
			     
			    
			    }
			}
		    			
                 }
	      
	      if ($cell_included_in_a_fragment[$c] == 0)
	         {                  
                 $fragment_no++;   
                 $cell_fragment_number[$c] = $fragment_no;
	          $cell_included_in_a_fragment[$c] = 1;	   
		   
		   #search neighbours who do not belong to a fragment		   
                 for ($nc = 1; $nc < $no_of_neighbours[$c] + 1; $nc++) 
                     {
	               $neighbour = $neighbour_cells[$c][$nc];    #print "$nc - $neighbour\n"; 

		        if ($cell_has_particles[$neighbour] == 1)
		           {
		            if ($cell_included_in_a_fragment[$neighbour] == 0)
			        {
				  $cell_fragment_number[$neighbour] = $cell_fragment_number[$c];
				  $cell_included_in_a_fragment[$neighbour] = 1;
				 }
			    }
			 } 		   
	         }
 
	     }

	  # do once more to ensure all neiighbours who have particles included   
	  if ($cell_included_in_a_fragment[$c] == 1)
	     {
	      #search neighbours who do not belong to a fragment
             for ($nc = 1; $nc < $no_of_neighbours[$c] + 1; $nc++) 
                 {
	           $neighbour = $neighbour_cells[$c][$nc];
		    if ($cell_has_particles[$neighbour] == 1)
		       {
			 if ($cell_included_in_a_fragment[$neighbour] == 0)
			    {
			     $cell_fragment_number[$neighbour] = $cell_fragment_number[$c];
			     $cell_included_in_a_fragment[$neighbour] = 1;			    
			    }
			 
			}
                 }                
	     }
	     	  
	 }
     
     else # cell has no particles
        {
	  # do nothing
	 } 
    
    }

$total_fragments =  $fragment_no; #print " $total_fragments\n";   
###########################################################################


###########################################################################

for ($f = 1; $f < $total_fragments + 1; $f++)
    {
     $counter = 0;
     for ($c = 1; $c < $total_cells + 1; $c++)
         {
	   if ($cell_fragment_number[$c] == $f)
             {
	       $counter++;
		$cells_in_fragment[$f][$counter] = $c;
	      }	   
         }	   
    }

for ($f = 1; $f < $total_fragments + 1; $f++)
    {
     $pf = 0;
     print " fragment $f (cells)\n";
     print " ==================\n";     
     for ($c = 1; $c < $#{$cells_in_fragment[$f]} + 1; $c++)
         {
	   $cell = $cells_in_fragment[$f][$c];
	   print " $cell";
	   
	   for ($pc = 1; $pc < $particles_in_cell[$cell] + 1; $pc++)
	       {
		 $pf++;
		 $particles_in_fragment[$f][$pf] = $cell_particle_id[$cell][$pc]; 
		}
	   
	  }
     print "\n";  print "\n";
    }

for ($f = 1; $f < $total_fragments + 1; $f++)
    {
     print " fragment $f (particles)\n";
     print " ============\n";
     
     print " Too many particles to display !\n";
     
     for ($pf = 1; $pf < $#{$particles_in_fragment[$f]} + 1; $pf++)
         {
	   #print " $particles_in_fragment[$f][$pf]"; 
	  }
     print "\n";  print "\n"; 
    }

        SUBROUTINE output(rbname,translation,rotation,
     '             output_interval)
                
!***********************************************************************
!        Written by Kumar Mithraratne
!        © Auckland Bioengineering Institute
!        V1.1 July 2013
!***********************************************************************
        
        IMPLICIT NONE
        
        INCLUDE 'rigidbody1.cmn'
                
        REAL*8 translation(nbm,0:time_steps,9),
     '         rotation(nbm,0:time_steps,18),output_interval
        !INTEGER
        CHARACTER rbname(nbm)*100,filename*10
        
        REAL*8 vel_mag(0:time_steps),vel_min,vel_max
        INTEGER len_rbname,nc,nf,len_filename
        CHARACTER comma*1
        
        
        print *, 'rbname=', rbname(1) 
        comma='"'
        DO nf=0,out_time_steps_reqd,1
           vel_mag(nf)=0.0d0
        ENDDO
        
        CALL system("rm -f *.exdata")
        
        DO nc=1,100,1
           IF(rbname(1)(nc:nc+1).EQ.' ') THEN
              len_rbname=nc-1
              EXIT
            ENDIF
        ENDDO
        
        print *, 'writing file_names.txt'
        print *, 'out_time_steps_reqd', out_time_steps_reqd
        OPEN(unit=20,file='file_names.txt')
        DO nf=0,out_time_steps_reqd,1
           WRITE(20,'(I6.6)') nf           
        ENDDO   
        CLOSE(unit=20)
        
        OPEN(unit=20,file='file_names.txt')
        OPEN(unit=25,file='view.com')        
        
        vel_min=10000000.0d0
        vel_max=-10000000.0d0
        DO nf=0,out_time_steps_reqd,1
           READ(20,'(A)') filename           
           DO nc=1,10,1
              IF(filename(nc:nc+1).EQ.' ') THEN
                 len_filename=nc-1
                 EXIT
              ENDIF
           ENDDO

           WRITE(25,'("gfx read data ",A," time",f8.4";")') 
     '           filename(2:len_filename)//'.exdata',
     '           output_interval*nf
           
           print *, 'filename = ', filename
           print *, 'Creating .exdata file ', filename(2:len_filename)
           OPEN(unit=30,file=filename(2:len_filename)//'.exdata')

           WRITE(30,'(" Group name: ",A)') rbname(1)(1:len_rbname)
           WRITE(30,'(" #Fields=6")')
           WRITE(30,'(" 1) coordinates, coordinate, ",
     '                "rectangular cartesian, #Components=3")')
           WRITE(30,'("  x.  Value index=1, #Derivatives=0, ", 
     '                "#Versions=1")')
           WRITE(30,'("  y.  Value index=2, #Derivatives=0, ",
     '                "#Versions=1")')
           WRITE(30,'("  z.  Value index=3, #Derivatives=0, ",
     '                "#Versions=1")')
           WRITE(30,'(" 2) velocity, coordinate, ",
     '                "rectangular cartesian, #Components=3")')
           WRITE(30,'("  1.  Value index=4, #Derivatives=0, ", 
     '                "#Versions=1")')
           WRITE(30,'("  2.  Value index=5, #Derivatives=0, ",
     '                "#Versions=1")')
           WRITE(30,'("  3.  Value index=6, #Derivatives=0, ",
     '                "#Versions=1")')
           WRITE(30,'(" 3) velocity_magnitude, coordinate, ",
     '                "rectangular cartesian, #Components=1")')
           WRITE(30,'("  1.  Value index=7, #Derivatives=0, ", 
     '                "#Versions=1")')        
           WRITE(30,'(" 4) principal_axis1, coordinate, ",
     '                "rectangular cartesian, #Components=3")');
           WRITE(30,'("  1.  Value index=8, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  2.  Value index=9, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  3.  Value index=10, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'(" 5) principal_axis2, coordinate, ",
     '                "rectangular cartesian, #Components=3")');  
           WRITE(30,'("  1.  Value index=11, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  2.  Value index=12, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  3.  Value index=13, #Derivatives=0, ",
     '                "#Versions=1")');        
           WRITE(30,'(" 6) principal_axis3, coordinate, ",
     '                "rectangular cartesian, #Components=3")');     
           WRITE(30,'("  1.  Value index=14, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  2.  Value index=15, #Derivatives=0, ",
     '                "#Versions=1")');
           WRITE(30,'("  3.  Value index=16, #Derivatives=0, ",
     '                "#Versions=1")');        

           WRITE(30,'(" Node: 1")') 
           
           WRITE(30,'("",3(f10.3))') translation(1,nf,1)*mm2m,
     '           translation(1,nf,2)*mm2m,translation(1,nf,3)*mm2m        
           WRITE(30,'("",3(f10.3))') translation(1,nf,4),
     '           translation(1,nf,5),translation(1,nf,6)        
           
           vel_mag(nf)=(translation(1,nf,4)**2.0d0+
     '                  translation(1,nf,5)**2.0d0+
     '                  translation(1,nf,6)**2.0d0)**0.5d0
           IF(vel_mag(nf).LT.vel_min) vel_min=vel_mag(nf)
           IF(vel_mag(nf).GT.vel_max) vel_max=vel_mag(nf)
                              
           WRITE(30,'("",1(f10.3))') vel_mag(nf)        
           WRITE(30,'("",3(f10.3))') rotation(1,nf,1),
     '           rotation(1,nf,2),rotation(1,nf,3)        
           WRITE(30,'("",3(f10.3))') rotation(1,nf,4),
     '           rotation(1,nf,5),rotation(1,nf,6)        
           WRITE(30,'("",3(f10.3))') rotation(1,nf,7),
     '           rotation(1,nf,8),rotation(1,nf,9)        
     
           CLOSE(unit=30)
        ENDDO   
        
        WRITE(25,'("")');        
        WRITE(25,'("gfx create axes length 500;")');
        WRITE(25,'("gfx draw axes;")');
        
        WRITE(25,'("")');
        WRITE(25,'("gfx modify spectrum default linear reverse range",
     '             2(1x,f8.4)," extend_above extend_below rainbow ",
     '             "colour_range 0 1 component 1;")') vel_min,vel_max
                
        WRITE(25,'("")');
        WRITE(25,'("gfx modify g_element projectile general clear ",
     '             "circle_discretization 6 default_coordinate ",
     '             "coordinates element_discretization ",A,"4*4*4",A,
     '             " native_discretization none;")')comma,comma
        WRITE(25,'("gfx modify g_element projectile data_points ",
     '                    "glyph sphere general size ",A,"50*50*50",A,
     '             " centre 0,0,0 font default select_on material ",
     '             "bone selected_material default_selected;")')
     '        comma,comma
        WRITE(25,'("gfx modify g_element projectile data_points ",
     '             "glyph arrow_solid general size ",A,"20*20*20",A,
     '             " centre 0,0,0 font default orientation velocity ",
     '             "variable_scale velocity_magnitude scale_factors ",
     '             A,"1*1*1",A," select_on material bone data ",
     '             "velocity_magnitude spectrum default ",
     '             "selected_material default_selected;")')
     '        comma,comma,comma,comma 
        WRITE(25,'("gfx modify g_element projectile data_points ",
     '             "glyph arrow_line general size ",A,"100*100*100",A,
     '             " centre 0,0,0 font default orientation ",
     '             "principal_axis1 scale_factors ",A,"1*1*1",A,
     '             " select_on material bone selected_material ",
     '             "default_selected;")') comma,comma,comma,comma
        WRITE(25,'("gfx modify g_element projectile data_points glyph",
     '             " arrow_line general size ",A,"100*100*100",A,
     '             " centre 0,0,0 font default orientation ",
     '             "principal_axis2 scale_factors ",A,"1*1*1",A,
     '             " select_on material bone selected_material ",
     '             "default_selected;")')comma,comma,comma,comma 
        WRITE(25,'("gfx modify g_element projectile data_points ",
     '             "glyph arrow_line general size ",A,"100*100*100",A,
     '             " centre 0,0,0 font default orientation ",
     '             "principal_axis3 scale_factors ",A,"1*1*1",A,
     '             " select_on material bone selected_material", 
     '             " default_selected;")')comma,comma,comma,comma 

        WRITE(25,'("gfx create window 1 double_buffer;")')
        WRITE(25,'("gfx modify window 1 image scene default ",
     '             "light_model default;")')
        WRITE(25,'("gfx modify window 1 image add_light default;")')
        WRITE(25,'("gfx modify window 1 layout simple ortho_axes ",
     '             "z -y eye_spacing 0.25 width 1083 height 643;")')
        WRITE(25,'("gfx modify window 1 set current_pane 1;")')
        WRITE(25,'("gfx modify window 1 background colour 0 0 0 ",
     '             "texture none;")')
        WRITE(25,'("gfx modify window 1 view parallel eye_point ",
     '             "3020.05 -1370.6 1309.04 interest_point ",
     '             "1537.73 1228.84 702.155 up_vector ",
     '             "0.0176225 0.236839 0.971389 view_angle 62.3129 ",
     '             "near_clipping_plane 30.533 far_clipping_plane ",
     '             "10911.5 relative_viewport ndc_placement -1 1 2 2 ",
     '             "viewport_coordinates 0 0 1 1;")')
        WRITE(25,'("gfx modify window 1 overlay scene none;")')
        WRITE(25,'("gfx modify window 1 set transform_tool ",
     '             "current_pane 1 std_view_angle 40 normal_lines ",
     '             "no_antialias depth_of_field 0.0 ",
     '             "fast_transparency blend_normal;")')


        CLOSE(unit=25)
        CLOSE(unit=20)        
        
        
        




             
        
        RETURN
        END



! This include block defines parameters used in the GeoModel.
!
! If you find that direct use of this include impedes your installation
! of the GeoModel, please contact the GeoModel developers 
!         brannon@mech.utah.edu
!         oestrac@sandia.gov
! and we will work to resolve the problem.
!
! Altering the GeoModel source code (or its includes) will result in
! loss of technical support.

!
! Just numbers
      PARAMETER (&
     &   PZERO  = 0.0D0        ,&
     &   P1MEPS = 1.0D0 - 1.0D-10,&
     &   PFORTH = 0.25D0       ,&
     &   PTHIRD = 0.3333333333333333333333333333333333333333333333D0,&
     &   P2THIRD = 0.6666666666666666666666666666666666666666666666D0,&
     &   PHALF  = 0.5D0        ,&
     &   PONE   = 1.0D0        ,&
     &   P4THIRD= 1.3333333333333333333333333333333333333333333333D0,&
     &   P10SEVENTH=10.0D0/7   ,&
     &   P3HALF = 1.5D0        ,&
     &   PTWO   = 2.0D0        ,&
     &   PTHREE = 3.0D0        ,&
     &   PFOUR  = 4.0D0        ,&
     &   PFIVE  = 5.0D0        ,&
     &   PSIX   = 6.0d0        ,&
     &   P5TEEN = 15.0D0&
     &   )
      PARAMETER (&
     &   P20    = 2.0D1        ,&
     &   BIGNUM = 1.0D30       ,  &
     &   P14TH  = 0.25D0       ,&
     &   P34TH  = 0.75D0       ,       &
     &   P1P05  = 1.05D0       ,&
     &   P1P1   = 1.1D0        ,&
     &   P1P14  = 1.14D0       ,&
     &   P1P2   = 1.2D0        ,&
     &   P1P5   = 1.5D0        ,&
     &   P2P1   = 2.1D0        ,&
     &   P2P5   = 2.5D0        ,&
     &   P3P5   = 3.5D0        ,&
     &   P3P524 = 3.524D0&
     &   )
      PARAMETER (&
     &   P0005  = 0.0005D0     ,&
     &   P0032  = 0.0032D0     ,&
     &   P001   = 0.001D0      ,&
     &   P002   = 0.002D0      ,&
     &   P007   = 0.007D0      ,&
     &   POINT01= 0.01D0       ,&
     &   POINT02= 0.02D0       ,&
     &   POINT03= 0.03D0       ,&
     &   POINT05= 0.05D0       ,&
     &   POINT06= 0.06D0       ,&
     &   POINT1 = 0.1D0        ,&
     &   P125   = 0.125D0      ,&
     &   POINT15= 0.15D0&
     &   )
      PARAMETER (&
     &   POINT33= 0.33D0       ,&
     &   POINT4 = 0.4D0        ,&
     &   POINT49= 0.49D0       ,&
     &   POINT501= 0.501D0     ,&
     &   POINT55= 0.55D0       ,&
     &   POINT9 = 0.9D0        ,&
     &   POINT95= 0.95D0       ,&
     &   POINT99= 0.99D0       ,&
     &   POINT999= 0.999D0     ,&
     &   POINT9999=0.9999D0    ,&
     &   PSIXTH = 1.0D0/6      ,&
     &   TWOTHD = 2.0D0/3&
     &   )
      PARAMETER (&
     &   P5SIXTHS= 5.0D0/6     ,&
     &   S9THS  = 16.0D0/9     ,&
     &   SIXTIETH=1.0D0/60     ,&
     &   P1M30  = 1.0D-30      ,&
     &   P1M20  = 1.0D-20      ,&
     &   P1M15  = 1.0D-15      ,&
     &   P1M9   = 1.0D-9       ,&
     &   P1M7   = 1.0D-7       ,&
     &   P1E3   = 1.0D3        ,&
     &   P1E4   = 1.0D4        ,&
     &   P1E6   = 1.0D6        ,&
     &   P1E10  = 1.0D10       ,&
     &   P1E12  = 1.0D12       ,&
     &   P1E15  = 1.0D15       ,&
     &   P1E20  = 1.0D20&
     &   )
!
!---.----1----.----2----.----3----.----4----.----5----.----6----.----7--
! Particular numbers
      PARAMETER (&
     &   ROOT2 = 1.414213562373095048801688724209698078569671875377D0,&
     &   ROOTO2 = PONE/ROOT2,&
     &   SQRT3  = 1.73205080756887729352744634150587236694280525381D0,&
     &   TWOORT3= PTWO/SQRT3,&
     &   SQRT2O100 = ROOT2/1.0d-2,&
     &   ROOT23 = 0.81649658092772603273242802490196379732198249355D0,&
     &   ROOT32 = PONE/ROOT23,&
     &   COS120 = -PHALF,&
     &   SIN120 = SQRT3/PTWO&
     &   )
!
! All kinds of pi
      PARAMETER (&
     &   PI     = 3.1415926535897932384626433832795028841971693993D0,&
     &   PIO2   = PI/PTWO,&
     &   TWOPI  = PTWO*PI,&
     &   RADDEG = PI/180,&
     &   DEGRAD = 180/PI,&
     &   TWOTHDPI = TWOPI/PTHREE&
     &   )
!
! Tolerances
      PARAMETER (&
     &   TOL1M3 = 1.0D-3,&
     &   TOL1M4 = 1.0D-4,&
     &   TOL1M5 = 1.0D-5,&
     &   TOL3M6 = 3.0D-6,&
     &   TOL1M6 = 1.0D-6,&
     &   TOL1M7 = 1.0D-7,&
     &   TOL1M8 = 1.0D-8,&
     &   TOL1M9 = 1.0D-9,&
     &   TOL1M10= 1.0D-10,&
     &   TOL1M12= 1.0D-12,&
     &   TOL1M14= 1.0D-14,&
     &   TOL1M20= 1.0D-20&
     &   )
!
! Flags
      PARAMETER (&
     &   UNDEF  = 1.23456D-7,&
     &   NOTDEF = -654321&
     &   )

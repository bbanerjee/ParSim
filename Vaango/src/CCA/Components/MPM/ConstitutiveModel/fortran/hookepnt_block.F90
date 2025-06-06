! This include block defines pointers and dimensioning information
! for the Hookes law model.
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!! important !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CODE OWNERS: This file is your primary reference for how things should
!  be ordered in calling arguments!  
!
!  Search this file for the string "CODE OWNERS" for further instructions.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@   d i m e n s i o n i n g    p a r a m e t e r s @@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  NBASICINPUTS: length of PROP array 
!  NHOOKEPROP: length of the PROP array including extra "add-on" properties
!  -----------------------------------------------------------------------
#ifdef BRANNON_IMPLNONE
      INTEGER NBASICINPUTS,NEOSINPUT,NHOOKEPROP
#endif
!----------------------------------------------------------------------------
!
      PARAMETER (NBASICINPUTS=2)
#ifdef HOOKE_EOS
      PARAMETER (NUIEOSMG=22,NDCEOSMG=13,NVIEOSMG=5)
      PARAMETER (NEOSINPUTS=NUIEOSMG+NDCEOSMG+NVIEOSMG)
#endif
!
!
!     Total number of properties
#ifdef HOOKE_EOS
      PARAMETER (NHOOKEPROP=NBASICINPUTS+NEOSINPUTS)
#else
      PARAMETER (NHOOKEPROP=NBASICINPUTS)
#endif
!
!
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@   p o i n t e r s   t o  p r o p e r t i e s @@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  IP[propname]: pointers to property array
!                Examples: the property B0 is in PROP(IPB0)
!                          the property FSLOPEI is in PROP(IPFSLOPI)
!                          and so on...
!-------------------------------------------------------------------------
#ifdef BRANNON_IMPLNONE
      INTEGER IPB0, IPG0
#ifdef HOOKE_EOS&
     & ,IPRHO0,IPTMPR0,IPSNDSP0,IPS1MG,IPGRPAR,IPCV,IPESFT,&
     & IPRP,IPPS,IPPE,IPCE,IPRNSUB,IPS2MG,IPTYP,IPRO,IPTO,IPSMG,&
     & IPGRPARO,IPBMG,IPXB,IPRNBMG,IPRPWR,IPA1MG,IPA2MG,IPA3MG,&
     & IPA4MG,IPA5MG,IPA0MG,IPAEMG,IPFK0,IPAF,IPPF,IPXF,IPCF,IPRMX,&
     & IPVI1MG,IPVI2MG,IPVI3MG,IPVI4MG,IPVI5MG
#endif
#endif

!-------------------------------------------------------------------------
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!! important !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CODE OWNERS: The host code is expected to read user inputs and
!     save them into a single array using the ordering indicated below.
!     The list below defines pointers to the property array. Every entry
!     is of the form "IPname", where 'name' is the string that we suggest
!     should be used as the property name. For example, the user input
!     that is called "B0" in the Hookes law model  should be read
!     from the user with the keyword "B0" and, as indicated below,
!     this parameter should be saved in 14th spot of the property array.
!     Note that some pointers farther down this list are wrapped in
!     ifdefs. Therefore your input parser will need similar ifdefs.
!
!   NBASICINPUTS points to last regular user input. 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER(IPB0       =  1) !Initial intact elastic bulk modulus (stress)
      PARAMETER(IPG0       =  2) !Initial intact elastic shear modulus (stress)
#ifdef HOOKE_EOS
      PARAMETER(IEOSMGCT   = NBASICINPUTS)
      PARAMETER(IPRHO0     = IEOSMGCT +  1)
      PARAMETER(IPTMPR0    = IEOSMGCT +  2)
      PARAMETER(IPSNDSP0   = IEOSMGCT +  3)
      PARAMETER(IPS1MG     = IEOSMGCT +  4)
      PARAMETER(IPGRPAR    = IEOSMGCT +  5)
      PARAMETER(IPCV       = IEOSMGCT +  6)
      PARAMETER(IPESFT     = IEOSMGCT +  7)
      PARAMETER(IPRP       = IEOSMGCT +  8)
      PARAMETER(IPPS       = IEOSMGCT +  9)
      PARAMETER(IPPE       = IEOSMGCT + 10)
      PARAMETER(IPCE       = IEOSMGCT + 11)
      PARAMETER(IPNSUB     = IEOSMGCT + 12)
      PARAMETER(IPS2MG     = IEOSMGCT + 13)
      PARAMETER(IPTYP      = IEOSMGCT + 14)
      PARAMETER(IPRO       = IEOSMGCT + 15)
      PARAMETER(IPTO       = IEOSMGCT + 16)
      PARAMETER(IPS        = IEOSMGCT + 17)
      PARAMETER(IPGRPARO   = IEOSMGCT + 18)
      PARAMETER(IPB        = IEOSMGCT + 19)
      PARAMETER(IPXB       = IEOSMGCT + 20)
      PARAMETER(IPNB       = IEOSMGCT + 21)
      PARAMETER(IPPWR      = IEOSMGCT + NUIEOSMG)
!     ________________________________________________________________________
!     EOSMG derived constants
      PARAMETER(IDCEOSMGCT=IEOSMGCT+NUIEOSMG)
      PARAMETER(IPA1MG     = IDCEOSMGCT +  1)
      PARAMETER(IPA2MG     = IDCEOSMGCT +  2)
      PARAMETER(IPA3MG     = IDCEOSMGCT +  3)
      PARAMETER(IPA4MG     = IDCEOSMGCT +  4)
      PARAMETER(IPA5MG     = IDCEOSMGCT +  5)
      PARAMETER(IPA0MG     = IDCEOSMGCT +  6)
      PARAMETER(IPAEMG     = IDCEOSMGCT +  7)
      PARAMETER(IPFK0      = IDCEOSMGCT +  8)
      PARAMETER(IPAF       = IDCEOSMGCT +  9)
      PARAMETER(IPPF       = IDCEOSMGCT + 10)
      PARAMETER(IPXF       = IDCEOSMGCT + 11)
      PARAMETER(IPCF       = IDCEOSMGCT + 12)
      PARAMETER(IPRMX      = IDCEOSMGCT + NDCEOSMG)
!     ________________________________________________________________________
!     EOSMG VI 
      PARAMETER(IVIEOSMGCT=IDCEOSMGCT+NDCEOSMG)
      PARAMETER(IPVI1MG     = IVIEOSMGCT +  1)
      PARAMETER(IPVI2MG     = IVIEOSMGCT +  2)
      PARAMETER(IPVI3MG     = IVIEOSMGCT +  3)
      PARAMETER(IPVI4MG     = IVIEOSMGCT +  4)
      PARAMETER(IPVI5MG     = IVIEOSMGCT +  NVIEOSMG)
#endif
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@   p o i n t e r s   t o  s t a t e   v a r i a b l e s  @@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  K[isvname]: pointers to the state variable array
!---.----1----.----2----.----3----.----4----.----5----.----6----.----7--
#ifdef BRANNON_IMPLNONE
      INTEGER NHOOKEISV,KMSG
#endif

#ifdef BRANNON_IMPLNONE
      INTEGER KTMPR,KENRGY,KRHO,KSNDSP,KPRES,KALPHA
#endif
#ifdef HOOKE_EOS
      PARAMETER (NHOOKEISV=8) !hardwired for SQA
#else
      PARAMETER (NHOOKEISV=1) !hardwired for SQA
#endif
!
!   selected pointers to the state variable array
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!! important !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CODE OWNERS:
!   Cross-check against subroutine HOOKERXV.   Even if you don't use HOOKERXV,
!   that's the routine to examine if you want to understand these values.
!   The following list shows pointers to internal state variables (ISVs).
!   A pointer with the name "Kname" should have a plot keyword "name",
!   and the value of the pointer shows where the ISV should be stored
!   in the ISV calling argument, called SVARG in the main physics routine.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER (KMSG     =1)  !Fake ISV for tracking installation
#ifdef HOOKE_EOS
      PARAMETER (KTMPR   =2 )  !Temperature
      PARAMETER (KENRGY  =3 )  !Energy
      PARAMETER (KDISTENRGY  =4 )  !Energy
      PARAMETER (KRHO    =5 )  !Density
      PARAMETER (KSNDSP  =6 )  !Soundspeed
      PARAMETER (KPRES   =7 )  !Pressure
      PARAMETER (KALPHA  =8 )  !Porosity parameter
#endif

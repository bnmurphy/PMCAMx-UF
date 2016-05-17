c-----CAMx v4.02 030709
c  
c     CHMSTRY.COM contains all chemistry variables 
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c 
c-----------------------------------------------------------------------
c     Parameters for some of the switches:
c
c     CDCMC  -- string for requesting the CMC (standard) chemistry solver
c     CDIEH  -- string for requesting the IEH chemistry solver
c     IDCMC  -- code for using the CMC (standard) chemistry solver
c     IDIEH  -- code for using the IEH chemistry solver
c-----------------------------------------------------------------------
c
      character*10 CDCMC
      character*10 CDIEH
      integer      IDCMC
      integer      IDIEH
c
      parameter( CDCMC = "CMC       " )
      parameter( CDIEH = "IEH       " )
      parameter( IDCMC = 1 )
      parameter( IDIEH = 2 )
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    ngas   --  number of gas species being modeled
c    naero  --  number of aersol species being modeled
c    nspec  --  total number of modeled species
c    nrad   --  number of radical species being modeled
c    nreact --  number of chemical reactions
c    nspfst --  number of "fast" species -- handled by the fast solver
c    iessrad--  number of radicals in steady state for IEH solver
c    idmech --  the code which determines which chemical mechanism is used
c    idsolv --  the code which determines which chemstry solver to use
c    navspc --  number of species to write to output average file
c    nicspc --  number of species in the initial conditions file
c    nbcspc --  number of species in the boundary conditions file
c    nptspc --  number of species in the point source emissions file
c    narspc --  number of species in the surface emissions file
c-----------------------------------------------------------------------
c
       integer   ngas
       integer   naero
       integer   nspec
       integer   nrad
       integer   nreact
       integer   nspfst
       integer   iessrad
       integer   idmech
       integer   idsolv
       integer   navspc
       integer   nicspc
       integer   nbcspc
       integer   nptspc
       integer   narspc
c
      common /chm1/ ngas, naero, nspec, nrad, nreact, nspfst, iessrad,
     &              idmech, idsolv, navspc, nicspc, nbcspc, nptspc, 
     &              narspc
c
c-----------------------------------------------------------------------
c     Variables for keeping track of where chmistry is being performed:
c     NOTE:  Used for diagnostic and error messages.
c
c     igrdchm  --  grid number of current chemistry step
c     ichm     --  column for the current chemistry step
c     jchm     --  row for the current chemistry step
c     kchm     --  layer for the current chemistry step
c-----------------------------------------------------------------------
c
      integer   igrdchm
      integer   ichm
      integer   jchm
      integer   kchm
c
      common /ijkgrd/ igrdchm, ichm, jchm, kchm
c
c-----------------------------------------------------------------------
c     Variables for storing chemical reaction data:
c
c     rk     -- reaction rate constant (ppm/hr)
c     ltdep  -- flag to determine if rate constant is temperature dependent
c     lpdep  -- flag to determine if rate constant is pressure dependent
c     bdnl   -- lower vound value for each modeled species (ppm)
c     bdlrad -- lower bound value for each radical species (ppm)
c-----------------------------------------------------------------------
c
      real    rk(MXRXN)
      logical ltdep(MXRXN)
      logical lpdep(MXRXN)
      real    bdnl(MXSPEC)
      real    bdlrad
      real    nflag  ! turns hno3 production off at night (tmg,06/19/04)
c
      common /chmratep/ rk, nflag
      common /chmrate/ ltdep, lpdep, bdnl, bdlrad
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     nphot1   -- number of primary photolysis reactions
c     nphot2   -- number of secondary (scaled) photolysis reactions
c     idphot1  -- ID of primary photolysis reactions
c     idphot2  -- ID of secondary (scaled) photolysis reactions 
c     idphot3  -- ID of primary photolysis reaction to scale to obtain
c                 the secondary photolysis reaction
c     phtscl   -- photolysis reaction scaling factor
c-----------------------------------------------------------------------
c
      integer   nphot1
      integer   nphot2
      integer   idphot1(MXPHT1)
      integer   idphot2(MXPHT2)
      integer   idphot3(MXPHT2)
      real      phtscl(MXPHT2)
c
      common /photmap/ nphot1, nphot2, idphot1, idphot2, idphot3, phtscl
c 
c-----------------------------------------------------------------------
c     Variables for species names:
c
c     spname  --  name of each modeled species
c     spavg   --  name of each species to be written to the output file
c     nmrad   --  name of each radical species
c     depsp   --  name of each deposition species output to file
c-----------------------------------------------------------------------
c
      character*10 spname(MXSPEC+1)
      character*10 spavg(MXSPEC)
      character*10 nmrad(MXRADCL+1)
      character*10 depsp(4*MXSPEC)
c
      common /cname/ spname, spavg, nmrad, depsp
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     krad     -- mapping of radical species to specific mechanism order
c     kmap     -- mapping of species on chemistry parameters file to
c                 internal order
c     lbcmap   -- mapping of species in the boundary condition file
c     lavmap   -- mapping of species written to average file
c     lptmap   -- mapping of species in the point source emissions file
c     lptrdmap -- mapping of species in the NetCDF point source file
c     larmap   -- mapping of species in the surface emissions file
c     licmap   -- mapping of species in the initial conditions file
c     lbcmapn  -- mapping of species on north edge of NetCDF boundary file
c     lbcmaps  -- mapping of species on south edge of NetCDF boundary file
c     lbcmape  -- mapping of species on east edge of NetCDF boundary file
c     lbcmapw  -- mapping of species on west edge of NetCDF boundary file
c     lgenrmap -- mapping of species in NetCDF general area emission file
c     lbiomap  -- mapping of species in NetCDF biogenic area emission file
c     lmoblmap -- mapping of species in NetCDF mobile area emission file
c     lspmap   -- mapping of species in NetCDF instant concentration file
c     lavwrmap -- mapping of species in NetCDF average concentration file
c-----------------------------------------------------------------------
c
      integer   krad(NRADNM)
      integer   kmap(NSPNAM)
      integer   lbcmap(MXSPEC)
      integer   lavmap(MXSPEC)
      integer   lptmap(MXSPEC)
      integer   lptrdmap(MXSPEC)
      integer   larmap(MXSPEC,MXGRID) 
      integer   licmap(MXSPEC,MXGRID) 
c
      common /kname/ krad, kmap, lbcmap,
     &               lavmap, lptmap, lptrdmap,
     &               larmap, licmap,
     &               lbcmapn(MXSPEC),lbcmaps(MXSPEC),lbcmape(MXSPEC),
     &               lbcmapw(MXSPEC),lgenrmap(MXSPEC,MXGRID),
     &               lbiomap(MXSPEC,MXGRID),lmoblmap(MXSPEC,MXGRID),
     &               lspmap(MXSPEC,MXGRID),lavwrmap(MXSPEC,MXGRID)
c
c      Gases
      integer   kno   ,kno2  ,ko3  
      integer   kpan  ,kcres ,kpan2
      integer   kmpan ,kpbzn ,knphe
      integer   krno3 ,kdcb2 ,kdcb3
      integer   khno4 ,kacet ,kald2
      integer   kalk1 ,kalk2 ,kalk3
      integer   kalk4 ,kalk5 ,karo1
      integer   karo2 ,kbacl ,kbald
      integer   kbcl1 ,kbcl2 ,kbuta
      integer   kccho ,kccrs ,kcg1 
      integer   kcg2  ,kcg3  ,kcg4 
      integer   kcl2  ,kco   ,kco2h
      integer   kco3h ,kcooh ,kcprm
      integer   kdcb1 ,keth  ,kethe
      integer   ketoh ,kfcrs ,kfmcl
      integer   kform ,kfprm ,kgly 
      integer   kh2o2 ,khc2h ,khcho
      integer   khcl  ,khono ,khno3
      integer   kho2h ,khocl ,kicl1
      integer   kicl2 ,kisop ,kispd
      integer   kmek  ,kmeoh ,kmeth
      integer   kmgly ,kmvk  ,kna  
      integer   knh3  ,kntr  ,knxoy
      integer   kole  ,kole1 ,kole2
      integer   kopen ,kpar  ,kpcl 
      integer   kpec  ,kphen ,kpna 
      integer   kpnh4 ,kpno3 ,kpoa 
      integer   kprod ,kpso4 ,krc2h
      integer   krc3h ,krcho ,krooh
      integer   kso2  ,ksoa1 ,ksoa2
      integer   ksoa3 ,ksoa4 ,ksulf
      integer   kterp ,ktol  ,kxn  
      integer   kxyl  ,kcg5  ,kamine     
      
c      Aerosols
      integer   ksoa1_1  ,ksoa1_2
      integer   ksoa1_3  ,ksoa1_4  ,ksoa1_5
      integer   ksoa1_6  ,ksoa1_7  ,ksoa1_8
      integer   ksoa1_9  ,ksoa1_10 ,ksoa2_1
      integer   ksoa2_2  ,ksoa2_3  ,ksoa2_4
      integer   ksoa2_5  ,ksoa2_6  ,ksoa2_7
      integer   ksoa2_8  ,ksoa2_9  ,ksoa2_10
      integer   ksoa3_1  ,ksoa3_2  ,ksoa3_3
      integer   ksoa3_4  ,ksoa3_5  ,ksoa3_6
      integer   ksoa3_7  ,ksoa3_8  ,ksoa3_9
      integer   ksoa3_10 ,ksoa4_1  ,ksoa4_2
      integer   ksoa4_3  ,ksoa4_4  ,ksoa4_5
      integer   ksoa4_6  ,ksoa4_7  ,ksoa4_8
      integer   ksoa4_9  ,ksoa4_10 
      integer   ksoa5_1  ,ksoa5_2  ,ksoa5_3
      integer   ksoa5_4  ,ksoa5_5  ,ksoa5_6
      integer   ksoa5_7  ,ksoa5_8  ,ksoa5_9
      integer   ksoa5_10 ,kpoc_1
      integer   kpoc_2   ,kpoc_3   ,kpoc_4
      integer   kpoc_5   ,kpoc_6   ,kpoc_7
      integer   kpoc_8   ,kpoc_9   ,kpoc_10
      integer   kpec_1   ,kpec_2   ,kpec_3
      integer   kpec_4   ,kpec_5   ,kpec_6
      integer   kpec_7   ,kpec_8   ,kpec_9
      integer   kpec_10  ,kcrust_1 ,kcrust_2
      integer   kcrust_3 ,kcrust_4 ,kcrust_5
      integer   kcrust_6 ,kcrust_7 ,kcrust_8
      integer   kcrust_9 ,kcrust_10,kph2o_1
      integer   kph2o_2  ,kph2o_3  ,kph2o_4
      integer   kph2o_5  ,kph2o_6  ,kph2o_7
      integer   kph2o_8  ,kph2o_9  ,kph2o_10
      integer   kpcl_1   ,kpcl_2   ,kpcl_3
      integer   kpcl_4   ,kpcl_5   ,kpcl_6
      integer   kpcl_7   ,kpcl_8   ,kpcl_9
      integer   kpcl_10  ,kna_1    ,kna_2
      integer   kna_3    ,kna_4    ,kna_5
      integer   kna_6    ,kna_7    ,kna_8
      integer   kna_9    ,kna_10   ,kpnh4_1
      integer   kpnh4_2  ,kpnh4_3  ,kpnh4_4
      integer   kpnh4_5  ,kpnh4_6  ,kpnh4_7
      integer   kpnh4_8  ,kpnh4_9  ,kpnh4_10
      integer   kpno3_1  ,kpno3_2  ,kpno3_3
      integer   kpno3_4  ,kpno3_5  ,kpno3_6
      integer   kpno3_7  ,kpno3_8  ,kpno3_9
      integer   kpno3_10 ,kpso4_1  ,kpso4_2
      integer   kpso4_3  ,kpso4_4  ,kpso4_5
      integer   kpso4_6  ,kpso4_7  ,kpso4_8
      integer   kpso4_9  ,kpso4_10 ,kpamine_1
      integer   kpamine_2,kpamine_3,kpamine_4
      integer   kpamine_5,kpamine_6,kpamine_7
      integer   kpamine_8,kpamine_9,kpamine_10
      integer   kph2o
c
      equivalence (kmap(1), kno  ), (kmap(2), kno2 ), (kmap(3), ko3  ),
     &            (kmap(4), kpan ), (kmap(5), kcres), (kmap(6), kpan2),
     &            (kmap(7), kmpan), (kmap(8), kpbzn), (kmap(9), knphe),
     &            (kmap(10),krno3), (kmap(11),kdcb2), (kmap(12),kdcb3),
     &            (kmap(13),khno4), (kmap(14),kacet), (kmap(15),kald2),
     &            (kmap(16),kalk1), (kmap(17),kalk2), (kmap(18),kalk3),
     &            (kmap(19),kalk4), (kmap(20),kalk5), (kmap(21),karo1),
     &            (kmap(22),karo2), (kmap(23),kbacl), (kmap(24),kbald),
     &            (kmap(25),kbcl1), (kmap(26),kbcl2), (kmap(27),kbuta),
     &            (kmap(28),kccho), (kmap(29),kccrs), (kmap(30),kcg1 ),
     &            (kmap(31),kcg2 ), (kmap(32),kcg3 ), (kmap(33),kcg4 ),
     &            (kmap(34),kcl2 ), (kmap(35),kco  ), (kmap(36),kco2h),
     &            (kmap(37),kco3h), (kmap(38),kcooh), (kmap(39),kcprm),
     &            (kmap(40),kdcb1), (kmap(41),keth ), (kmap(42),kethe),
     &            (kmap(43),ketoh), (kmap(44),kfcrs), (kmap(45),kfmcl),
     &            (kmap(46),kform), (kmap(47),kfprm), (kmap(48),kgly ),
     &            (kmap(49),kh2o2), (kmap(50),khc2h), (kmap(51),khcho),
     &            (kmap(52),khcl ), (kmap(53),khono), (kmap(54),khno3),
     &            (kmap(55),kho2h), (kmap(56),khocl), (kmap(57),kicl1),
     &            (kmap(58),kicl2), (kmap(59),kisop), (kmap(60),kispd),
     &            (kmap(61),kmek ), (kmap(62),kmeoh), (kmap(63),kmeth),
     &            (kmap(64),kmgly), (kmap(65),kmvk ), (kmap(66),kna  ),
     &            (kmap(67),knh3 ), (kmap(68),kntr ), (kmap(69),knxoy),
     &            (kmap(70),kole ), (kmap(71),kole1), (kmap(72),kole2),
     &            (kmap(73),kopen), (kmap(74),kpar ), (kmap(75),kpcl ),
     &            (kmap(76),kpec ), (kmap(77),kphen), (kmap(78),kpna ),
     &            (kmap(79),kpnh4), (kmap(80),kpno3), (kmap(81),kpoa ),
     &            (kmap(82),kprod), (kmap(83),kpso4), (kmap(84),krc2h),
     &            (kmap(85),krc3h), (kmap(86),krcho), (kmap(87),krooh),
     &            (kmap(88),kso2 ), (kmap(89),ksoa1), (kmap(90),ksoa2),
     &            (kmap(91),ksoa3), (kmap(92),ksoa4), (kmap(93),ksoa5),
     &            (kmap(94),ksulf), (kmap(95),kterp), (kmap(96),ktol ),
     &            (kmap(97),kxn  ), (kmap(98),kxyl ), (kmap(99),kamine)
      equivalence (kmap(100), ksoa1_1  ),(kmap(101),ksoa1_2  ),  
     &(kmap(102),ksoa1_3  ),(kmap(103),ksoa1_4  ),(kmap(104),ksoa1_5  ),
     &(kmap(105),ksoa1_6  ),(kmap(106),ksoa1_7  ),(kmap(107),ksoa1_8  ),
     &(kmap(108),ksoa1_9  ),(kmap(109),ksoa1_10 ),(kmap(110),ksoa1_11 ),
     &(kmap(111),ksoa1_12 ),(kmap(112),ksoa1_13 ),(kmap(113),ksoa1_14 ),
     &(kmap(114),ksoa1_15 ),(kmap(115),ksoa1_16 ),(kmap(116),ksoa1_17 ),
     &(kmap(117),ksoa1_18 ),(kmap(118),ksoa1_19 ),(kmap(119),ksoa1_20 ),
     &(kmap(120),ksoa1_21 ),(kmap(121),ksoa1_22 ),(kmap(122),ksoa1_23 ),
     &(kmap(123),ksoa1_24 ),(kmap(124),ksoa1_25 ),(kmap(125),ksoa1_26 ),
     &(kmap(126),ksoa1_27 ),(kmap(127),ksoa1_28 ),(kmap(128),ksoa1_29 ),
     &(kmap(129),ksoa1_30 ),(kmap(130),ksoa1_31 ),(kmap(131),ksoa1_32 ),
     &(kmap(132),ksoa1_33 ),(kmap(133),ksoa1_34 ),(kmap(134),ksoa1_35 ),
     &(kmap(135),ksoa1_36 ),(kmap(136),ksoa1_37 ),(kmap(137),ksoa1_38 ),
     &(kmap(138),ksoa1_39 ),(kmap(139),ksoa1_40 ),(kmap(140),ksoa1_41 ),
     &(kmap(141),ksoa1_42 ),(kmap(142),ksoa1_43 ),
     &(kmap(143),ksoa2_1  ),
     &(kmap(144),ksoa2_2  ),(kmap(145),ksoa2_3  ),(kmap(146),ksoa2_4  ),
     &(kmap(147),ksoa2_5  ),(kmap(148),ksoa2_6  ),(kmap(149),ksoa2_7  ),
     &(kmap(150),ksoa2_8  ),(kmap(151),ksoa2_9  ),(kmap(152),ksoa2_10 ),
     &(kmap(153),ksoa2_11 ),(kmap(154),ksoa2_12 ),(kmap(155),ksoa2_13 ),
     &(kmap(156),ksoa2_14 ),(kmap(157),ksoa2_15 ),(kmap(158),ksoa2_16 ),
     &(kmap(159),ksoa2_17 ),(kmap(160),ksoa2_18 ),(kmap(161),ksoa2_19 ),
     &(kmap(162),ksoa2_20 ),(kmap(163),ksoa2_21 ),(kmap(164),ksoa2_22 ),
     &(kmap(165),ksoa2_23 ),(kmap(166),ksoa2_24 ),(kmap(167),ksoa2_25 ),
     &(kmap(168),ksoa2_26 ),(kmap(169),ksoa2_27 ),(kmap(170),ksoa2_28 ),
     &(kmap(171),ksoa2_29 ),(kmap(172),ksoa2_30 ),(kmap(173),ksoa2_31 ),
     &(kmap(174),ksoa2_32 ),(kmap(175),ksoa2_33 ),(kmap(176),ksoa2_34 ),
     &(kmap(177),ksoa2_35 ),(kmap(178),ksoa2_36 ),(kmap(179),ksoa2_37 ),
     &(kmap(180),ksoa2_38 ),(kmap(181),ksoa2_39 ),(kmap(182),ksoa2_40 ),
     &(kmap(183),ksoa2_41 ),(kmap(184),ksoa2_42 ),(kmap(185),ksoa2_43 ),
     &(kmap(186),ksoa3_1  ),(kmap(187),ksoa3_2  ),(kmap(188),ksoa3_3  ),
     &(kmap(189),ksoa3_4  ),(kmap(190),ksoa3_5  ),(kmap(191),ksoa3_6  ),
     &(kmap(192),ksoa3_7  ),(kmap(193),ksoa3_8  ),(kmap(194),ksoa3_9  ),
     &(kmap(195),ksoa3_10 ),(kmap(196),ksoa3_11 ),(kmap(197),ksoa3_12 ),
     &(kmap(198),ksoa3_13 ),(kmap(199),ksoa3_14 ),(kmap(200),ksoa3_15 ),
     &(kmap(201),ksoa3_16 ),(kmap(202),ksoa3_17 ),(kmap(203),ksoa3_18 ),
     &(kmap(204),ksoa3_19 ),(kmap(205),ksoa3_20 ),(kmap(206),ksoa3_21 ),
     &(kmap(207),ksoa3_22 ),(kmap(208),ksoa3_23 ),(kmap(209),ksoa3_24 ),
     &(kmap(210),ksoa3_25 ),(kmap(211),ksoa3_26 ),(kmap(212),ksoa3_27 ),
     &(kmap(213),ksoa3_28 ),(kmap(214),ksoa3_29 ),(kmap(215),ksoa3_30 ),
     &(kmap(216),ksoa3_31 ),(kmap(217),ksoa3_32 ),(kmap(218),ksoa3_33 ),
     &(kmap(219),ksoa3_34 ),(kmap(220),ksoa3_35 ),(kmap(221),ksoa3_36 ),
     &(kmap(222),ksoa3_37 ),(kmap(223),ksoa3_38 ),(kmap(224),ksoa3_39 ),
     &(kmap(225),ksoa3_40 ),(kmap(226),ksoa3_41 ),(kmap(227),ksoa3_42 ),
     &(kmap(228),ksoa3_43 ),
     &(kmap(229),ksoa4_1  ),(kmap(230),ksoa4_2  ),
     &(kmap(231),ksoa4_3  ),(kmap(232),ksoa4_4  ),(kmap(233),ksoa4_5  ),
     &(kmap(234),ksoa4_6  ),(kmap(235),ksoa4_7  ),(kmap(236),ksoa4_8  ),
     &(kmap(237),ksoa4_9  ),(kmap(238),ksoa4_10 ),(kmap(239),ksoa4_11 ),
     &(kmap(240),ksoa4_12 ),(kmap(241),ksoa4_13 ),(kmap(242),ksoa4_14 ),
     &(kmap(243),ksoa4_15 ),(kmap(244),ksoa4_16 ),(kmap(245),ksoa4_17 ),
     &(kmap(246),ksoa4_18 ),(kmap(247),ksoa4_19 ),(kmap(248),ksoa4_20 ),
     &(kmap(249),ksoa4_21 ),(kmap(250),ksoa4_22 ),(kmap(251),ksoa4_23 ),
     &(kmap(252),ksoa4_24 ),(kmap(253),ksoa4_25 ),(kmap(254),ksoa4_26 ),
     &(kmap(255),ksoa4_27 ),(kmap(256),ksoa4_28 ),(kmap(257),ksoa4_29 ),
     &(kmap(258),ksoa4_30 ),(kmap(259),ksoa4_31 ),(kmap(260),ksoa4_32 ),
     &(kmap(261),ksoa4_33 ),(kmap(262),ksoa4_34 ),(kmap(263),ksoa4_35 ),
     &(kmap(264),ksoa4_36 ),(kmap(265),ksoa4_37 ),(kmap(266),ksoa4_38 ),
     &(kmap(267),ksoa4_39 ),(kmap(268),ksoa4_40 ),(kmap(269),ksoa4_41 ),
     &(kmap(270),ksoa4_42 ),(kmap(271),ksoa4_43 ),
     &(kmap(272),ksoa5_1  ),(kmap(273),ksoa5_2  ),
     &(kmap(274),ksoa5_3  ),(kmap(275),ksoa5_4  ),(kmap(276),ksoa5_5  ),
     &(kmap(277),ksoa5_6  ),(kmap(278),ksoa5_7  ),(kmap(279),ksoa5_8  ),
     &(kmap(280),ksoa5_9  ),(kmap(281),ksoa5_10 ),(kmap(282),ksoa5_11 ),
     &(kmap(283),ksoa5_12 ),(kmap(284),ksoa5_13 ),(kmap(285),ksoa5_14 ),
     &(kmap(286),ksoa5_15 ),(kmap(287),ksoa5_16 ),(kmap(288),ksoa5_17 ),
     &(kmap(289),ksoa5_18 ),(kmap(290),ksoa5_19 ),(kmap(291),ksoa5_20 ),
     &(kmap(292),ksoa5_21 ),(kmap(293),ksoa5_22 ),(kmap(294),ksoa5_23 ),
     &(kmap(295),ksoa5_24 ),(kmap(296),ksoa5_25 ),(kmap(297),ksoa5_26 ),
     &(kmap(298),ksoa5_27 ),(kmap(299),ksoa5_28 ),(kmap(300),ksoa5_29 ),
     &(kmap(301),ksoa5_30 ),(kmap(302),ksoa5_31 ),(kmap(303),ksoa5_32 ),
     &(kmap(304),ksoa5_33 ),(kmap(305),ksoa5_34 ),(kmap(306),ksoa5_35 ),
     &(kmap(307),ksoa5_36 ),(kmap(308),ksoa5_37 ),(kmap(309),ksoa5_38 ),
     &(kmap(310),ksoa5_39 ),(kmap(311),ksoa5_40 ),(kmap(312),ksoa5_41 ),
     &(kmap(313),ksoa5_42 ),(kmap(314),ksoa5_43 ),
     &(kmap(315),kpoc_1   ),
     &(kmap(316),kpoc_2   ),(kmap(317),kpoc_3   ),(kmap(318),kpoc_4   ),
     &(kmap(319),kpoc_5   ),(kmap(320),kpoc_6   ),(kmap(321),kpoc_7   ),
     &(kmap(322),kpoc_8   ),(kmap(323),kpoc_9   ),(kmap(324),kpoc_10  ),
     &(kmap(325),kpoc_11  ),(kmap(326),kpoc_12  ),(kmap(327),kpoc_13  ),
     &(kmap(328),kpoc_14  ),(kmap(329),kpoc_15  ),(kmap(330),kpoc_16  ),
     &(kmap(331),kpoc_17  ),(kmap(332),kpoc_18  ),(kmap(333),kpoc_19  ),
     &(kmap(334),kpoc_20  ),(kmap(335),kpoc_21  ),(kmap(336),kpoc_22  ),
     &(kmap(337),kpoc_23  ),(kmap(338),kpoc_24  ),(kmap(339),kpoc_25  ),
     &(kmap(340),kpoc_26  ),(kmap(341),kpoc_27  ),(kmap(342),kpoc_28  ),
     &(kmap(343),kpoc_29  ),(kmap(344),kpoc_30  ),(kmap(345),kpoc_31  ),
     &(kmap(346),kpoc_32  ),(kmap(347),kpoc_33  ),(kmap(348),kpoc_34  ),
     &(kmap(349),kpoc_35  ),(kmap(350),kpoc_36  ),(kmap(351),kpoc_37  ),
     &(kmap(352),kpoc_38  ),(kmap(353),kpoc_39  ),(kmap(354),kpoc_40  ),
     &(kmap(355),kpoc_41  ),(kmap(356),kpoc_42  ),(kmap(357),kpoc_43  ),
     &(kmap(358),kpec_1   ),(kmap(359),kpec_2   ),(kmap(360),kpec_3   ),
     &(kmap(361),kpec_4   ),(kmap(362),kpec_5   ),(kmap(363),kpec_6   ),
     &(kmap(364),kpec_7   ),(kmap(365),kpec_8   ),(kmap(366),kpec_9   ),
     &(kmap(367),kpec_10  ),(kmap(368),kpec_11  ),(kmap(369),kpec_12  ),
     &(kmap(370),kpec_13  ),(kmap(371),kpec_14  ),(kmap(372),kpec_15  ),
     &(kmap(373),kpec_16  ),(kmap(374),kpec_17  ),(kmap(375),kpec_18  ),
     &(kmap(376),kpec_19  ),(kmap(377),kpec_20  ),(kmap(378),kpec_21  ),
     &(kmap(379),kpec_22  ),(kmap(380),kpec_23  ),(kmap(381),kpec_24  ),
     &(kmap(382),kpec_25  ),(kmap(383),kpec_26  ),(kmap(384),kpec_27  ),
     &(kmap(385),kpec_28  ),(kmap(386),kpec_29  ),(kmap(387),kpec_30  ),
     &(kmap(388),kpec_31  ),(kmap(389),kpec_32  ),(kmap(390),kpec_33  ),
     &(kmap(391),kpec_34  ),(kmap(392),kpec_35  ),(kmap(393),kpec_36  ),
     &(kmap(394),kpec_37  ),(kmap(395),kpec_38  ),(kmap(396),kpec_39  ),
     &(kmap(397),kpec_40  ),(kmap(398),kpec_41  ),(kmap(399),kpec_42  ),
     &(kmap(400),kpec_43  ),
     &(kmap(401),kcrust_1 ),(kmap(402),kcrust_2 ),
     &(kmap(403),kcrust_3 ),(kmap(404),kcrust_4 ),(kmap(405),kcrust_5 ),
     &(kmap(406),kcrust_6 ),(kmap(407),kcrust_7 ),(kmap(408),kcrust_8 ),
     &(kmap(409),kcrust_9 ),(kmap(410),kcrust_10),(kmap(411),kcrust_11),
     &(kmap(412),kcrust_12),(kmap(413),kcrust_13),(kmap(414),kcrust_14),
     &(kmap(415),kcrust_15),(kmap(416),kcrust_16),(kmap(417),kcrust_17),
     &(kmap(418),kcrust_18),(kmap(419),kcrust_19),(kmap(420),kcrust_20),
     &(kmap(421),kcrust_21),(kmap(422),kcrust_22),(kmap(423),kcrust_23),
     &(kmap(424),kcrust_24),(kmap(425),kcrust_25),(kmap(426),kcrust_26),
     &(kmap(427),kcrust_27),(kmap(428),kcrust_28),(kmap(429),kcrust_29),
     &(kmap(430),kcrust_30),(kmap(431),kcrust_31),(kmap(432),kcrust_32),
     &(kmap(433),kcrust_33),(kmap(434),kcrust_34),(kmap(435),kcrust_35),
     &(kmap(436),kcrust_36),(kmap(437),kcrust_37),(kmap(438),kcrust_38),
     &(kmap(439),kcrust_39),(kmap(440),kcrust_40),(kmap(441),kcrust_41),
     &(kmap(442),kcrust_42),(kmap(443),kcrust_43),
     &(kmap(444),kph2o_1  ),
     &(kmap(445),kph2o_2  ),(kmap(446),kph2o_3  ),(kmap(447),kph2o_4  ),
     &(kmap(448),kph2o_5  ),(kmap(449),kph2o_6  ),(kmap(450),kph2o_7  ),
     &(kmap(451),kph2o_8  ),(kmap(452),kph2o_9  ),(kmap(453),kph2o_10 ),
     &(kmap(454),kph2o_11 ),(kmap(455),kph2o_12 ),(kmap(456),kph2o_13 ),
     &(kmap(457),kph2o_14 ),(kmap(458),kph2o_15 ),(kmap(459),kph2o_16 ),
     &(kmap(460),kph2o_17 ),(kmap(461),kph2o_18 ),(kmap(462),kph2o_19 ),
     &(kmap(463),kph2o_20 ),(kmap(464),kph2o_21 ),(kmap(465),kph2o_22 ),
     &(kmap(466),kph2o_23 ),(kmap(467),kph2o_24 ),(kmap(468),kph2o_25 ),
     &(kmap(469),kph2o_26 ),(kmap(470),kph2o_27 ),(kmap(471),kph2o_28 ),
     &(kmap(472),kph2o_29 ),(kmap(473),kph2o_30 ),(kmap(474),kph2o_31 ),
     &(kmap(475),kph2o_32 ),(kmap(476),kph2o_33 ),(kmap(477),kph2o_34 ),
     &(kmap(478),kph2o_35 ),(kmap(479),kph2o_36 ),(kmap(480),kph2o_37 ),
     &(kmap(481),kph2o_38 ),(kmap(482),kph2o_39 ),(kmap(483),kph2o_40 ),
     &(kmap(484),kph2o_41 ),(kmap(485),kph2o_42 ),(kmap(486),kph2o_43 ),
     &(kmap(487),kpcl_1   ),(kmap(488),kpcl_2   ),(kmap(489),kpcl_3   ),
     &(kmap(490),kpcl_4   ),(kmap(491),kpcl_5   ),(kmap(492),kpcl_6   ),
     &(kmap(493),kpcl_7   ),(kmap(494),kpcl_8   ),(kmap(495),kpcl_9   ),
     &(kmap(496),kpcl_10  ),(kmap(497),kpcl_11  ),(kmap(498),kpc1_12  ),
     &(kmap(499),kpcl_13  ),(kmap(500),kpcl_14  ),(kmap(501),kpcl_15  ),
     &(kmap(502),kpcl_16  ),(kmap(503),kpcl_17  ),(kmap(504),kpcl_18  ),
     &(kmap(505),kpcl_19  ),(kmap(506),kpcl_20  ),(kmap(507),kpcl_21  ),
     &(kmap(508),kpcl_22  ),(kmap(509),kpcl_23  ),(kmap(510),kpc1_24  ),
     &(kmap(511),kpcl_25  ),(kmap(512),kpcl_26  ),(kmap(513),kpc1_27  ),
     &(kmap(514),kpcl_28  ),(kmap(515),kpcl_29  ),(kmap(516),kpcl_30  ),
     &(kmap(517),kpcl_31  ),(kmap(518),kpcl_32  ),(kmap(519),kpcl_33  ),
     &(kmap(520),kpcl_34  ),(kmap(521),kpcl_35  ),(kmap(522),kpcl_36  ),
     &(kmap(523),kpcl_37  ),(kmap(524),kpcl_38  ),(kmap(525),kpc1_39  ),
     &(kmap(526),kpcl_40  ),(kmap(527),kpcl_41  ),(kmap(528),kpc1_42  ),
     &(kmap(529),kpcl_43  ),
     &(kmap(530),kna_1    ),(kmap(531),kna_2    ),
     &(kmap(532),kna_3    ),(kmap(533),kna_4    ),(kmap(534),kna_5    ),
     &(kmap(535),kna_6    ),(kmap(536),kna_7    ),(kmap(537),kna_8    ),
     &(kmap(538),kna_9    ),(kmap(539),kna_10   ),(kmap(540),kna_11   ),
     &(kmap(541),kna_12   ),(kmap(542),kna_13   ),(kmap(543),kna_14   ),
     &(kmap(544),kna_15   ),(kmap(545),kna_16   ),(kmap(546),kna_17   ),
     &(kmap(547),kna_18   ),(kmap(548),kna_19   ),(kmap(549),kna_20   ),
     &(kmap(550),kna_21   ),(kmap(551),kna_22   ),(kmap(552),kna_23   ),
     &(kmap(553),kna_24   ),(kmap(554),kna_25   ),(kmap(555),kna_26   ),
     &(kmap(556),kna_27   ),(kmap(557),kna_28   ),(kmap(558),kna_29   ),
     &(kmap(559),kna_30   ),(kmap(560),kna_31   ),(kmap(561),kna_32   ),
     &(kmap(562),kna_33   ),(kmap(563),kna_34   ),(kmap(564),kna_35   ),
     &(kmap(565),kna_36   ),(kmap(566),kna_37   ),(kmap(567),kna_38   ),
     &(kmap(568),kna_39   ),(kmap(569),kna_40   ),(kmap(570),kna_41   ),
     &(kmap(571),kna_42   ),(kmap(572),kna_43   ),
     &(kmap(573),kpnh4_1  ),
     &(kmap(574),kpnh4_2  ),(kmap(575),kpnh4_3  ),(kmap(576),kpnh4_4  ),
     &(kmap(577),kpnh4_5  ),(kmap(578),kpnh4_6  ),(kmap(579),kpnh4_7  ),
     &(kmap(580),kpnh4_8  ),(kmap(581),kpnh4_9  ),(kmap(582),kpnh4_10 ),
     &(kmap(583),kpnh4_11 ),(kmap(584),kpnh4_12 ),(kmap(585),kpnh4_13 ),
     &(kmap(586),kpnh4_14 ),(kmap(587),kpnh4_15 ),(kmap(588),kpnh4_16 ),
     &(kmap(589),kpnh4_17 ),(kmap(590),kpnh4_18 ),(kmap(591),kpnh4_19 ),
     &(kmap(592),kpnh4_20 ),(kmap(593),kpnh4_21 ),(kmap(594),kpnh4_22 ),
     &(kmap(595),kpnh4_23 ),(kmap(596),kpnh4_24 ),(kmap(597),kpnh4_25 ),
     &(kmap(598),kpnh4_26 ),(kmap(599),kpnh4_27 ),(kmap(600),kpnh4_28 ),
     &(kmap(601),kpnh4_29 ),(kmap(602),kpnh4_30 ),(kmap(603),kpnh4_31 ),
     &(kmap(604),kpnh4_32 ),(kmap(605),kpnh4_33 ),(kmap(606),kpnh4_34 ),
     &(kmap(607),kpnh4_35 ),(kmap(608),kpnh4_36 ),(kmap(609),kpnh4_37 ),
     &(kmap(610),kpnh4_38 ),(kmap(611),kpnh4_39 ),(kmap(612),kpnh4_40 ),
     &(kmap(613),kpnh4_41 ),(kmap(614),kpnh4_42 ),(kmap(615),kpnh4_43 ),
     &(kmap(616),kpno3_1  ),(kmap(617),kpno3_2  ),(kmap(618),kpno3_3  ),
     &(kmap(619),kpno3_4  ),(kmap(620),kpno3_5  ),(kmap(621),kpno3_6  ),
     &(kmap(622),kpno3_7  ),(kmap(623),kpno3_8  ),(kmap(624),kpno3_9  ),
     &(kmap(625),kpno3_10 ),(kmap(626),kpno3_11 ),(kmap(627),kpno3_12 ),
     &(kmap(628),kpno3_13 ),(kmap(629),kpno3_14 ),(kmap(630),kpno3_15 ),
     &(kmap(631),kpno3_16 ),(kmap(632),kpno3_17 ),(kmap(633),kpno3_18 ),
     &(kmap(634),kpno3_19 ),(kmap(635),kpno3_20 ),(kmap(636),kpno3_21 ),
     &(kmap(637),kpno3_22 ),(kmap(638),kpno3_23 ),(kmap(639),kpno3_24 ),
     &(kmap(640),kpno3_25 ),(kmap(641),kpno3_26 ),(kmap(642),kpno3_27 ),
     &(kmap(643),kpno3_28 ),(kmap(644),kpno3_29 ),(kmap(645),kpno3_30 ),
     &(kmap(646),kpno3_31 ),(kmap(647),kpno3_32 ),(kmap(648),kpno3_33 ),
     &(kmap(649),kpno3_34 ),(kmap(650),kpno3_35 ),(kmap(651),kpno3_36 ),
     &(kmap(652),kpno3_37 ),(kmap(653),kpno3_38 ),(kmap(654),kpno3_39 ),
     &(kmap(655),kpno3_40 ),(kmap(656),kpno3_41 ),(kmap(657),kpno3_42 ),
     &(kmap(658),kpno3_43),
     &(kmap(659),kpso4_1  ),(kmap(660),kpso4_2  ),
     &(kmap(661),kpso4_3  ),(kmap(662),kpso4_4  ),(kmap(663),kpso4_5  ),
     &(kmap(664),kpso4_6  ),(kmap(665),kpso4_7  ),(kmap(666),kpso4_8  ),
     &(kmap(667),kpso4_9  ),(kmap(668),kpso4_10 ),(kmap(669),kpso4_11 ),
     &(kmap(670),kpso4_12 ),(kmap(671),kpso4_13 ),(kmap(672),kpso4_14 ),
     &(kmap(673),kpso4_15 ),(kmap(674),kpso4_16 ),(kmap(675),kpso4_17 ),
     &(kmap(676),kpso4_18 ),(kmap(677),kpso4_19 ),(kmap(678),kpso4_20 ),
     &(kmap(679),kpso4_21 ),(kmap(680),kpso4_22 ),(kmap(681),kpso4_23 ),
     &(kmap(682),kpso4_24 ),(kmap(683),kpso4_25 ),(kmap(684),kpso4_26 ),
     &(kmap(685),kpso4_27 ),(kmap(686),kpso4_28 ),(kmap(687),kpso4_29 ),
     &(kmap(688),kpso4_30 ),(kmap(689),kpso4_31 ),(kmap(690),kpso4_32 ),
     &(kmap(691),kpso4_33 ),(kmap(692),kpso4_34 ),(kmap(693),kpso4_35 ),
     &(kmap(694),kpso4_36 ),(kmap(695),kpso4_37 ),(kmap(696),kpso4_38 ),
     &(kmap(697),kpso4_39 ),(kmap(698),kpso4_40 ),(kmap(699),kpso4_41 ),
     &(kmap(700),kpso4_42 ),(kmap(701),kpso4_43 ),
     &(kmap(702),kpamine_1 ),(kmap(703),kpamine_2 ),(kmap(704),kpamine_3 ),
     &(kmap(705),kpamine_4 ),(kmap(706),kpamine_5 ),(kmap(707),kpamine_6 ),
     &(kmap(708),kpamine_7 ),(kmap(709),kpamine_8 ),(kmap(710),kpamine_9 ),
     &(kmap(711),kpamine_10),(kmap(712),kpamine_11),(kmap(713),kpamine_12),
     &(kmap(714),kpamine_13),(kmap(715),kpamine_14),(kmap(716),kpamine_15),
     &(kmap(717),kpamine_16),(kmap(718),kpamine_17),(kmap(719),kpamine_18),
     &(kmap(720),kpamine_19),(kmap(721),kpamine_20),(kmap(722),kpamine_21),
     &(kmap(723),kpamine_22),(kmap(724),kpamine_23),(kmap(725),kpamine_24),
     &(kmap(726),kpamine_25),(kmap(727),kpamine_26),(kmap(728),kpamine_27),
     &(kmap(729),kpamine_28),(kmap(730),kpamine_29),(kmap(731),kpamine_30),
     &(kmap(732),kpamine_31),(kmap(733),kpamine_32),(kmap(734),kpamine_33),
     &(kmap(735),kpamine_34),(kmap(736),kpamine_35),(kmap(737),kpamine_36),
     &(kmap(738),kpamine_37),(kmap(739),kpamine_38),(kmap(740),kpamine_39),
     &(kmap(741),kpamine_40),(kmap(742),kpamine_41),(kmap(743),kpamine_42),
     &(kmap(744),kpamine_43),
     &(kmap(745),knum_1   ),(kmap(746),knum_2   ),(kmap(747),knum_3   ),
     &(kmap(748),knum_4   ),(kmap(749),knum_5   ),(kmap(750),knum_6   ),
     &(kmap(751),knum_7   ),(kmap(752),knum_8   ),(kmap(753),knum_9   ),
     &(kmap(754),knum_10  ),(kmap(755),knum_11  ),(kmap(756),knum_12  ),
     &(kmap(757),knum_13  ),(kmap(758),knum_14  ),(kmap(759),knum_15  ),
     &(kmap(760),knum_16  ),(kmap(761),knum_17  ),(kmap(762),knum_18  ),
     &(kmap(763),knum_19  ),(kmap(764),knum_20  ),(kmap(765),knum_21  ),
     &(kmap(766),knum_22  ),(kmap(767),knum_23  ),(kmap(768),knum_24  ),
     &(kmap(769),knum_25  ),(kmap(770),knum_26  ),(kmap(771),knum_27  ),
     &(kmap(772),knum_28  ),(kmap(773),knum_29  ),(kmap(774),knum_30  ),
     &(kmap(775),knum_31  ),(kmap(776),knum_32  ),(kmap(777),knum_33  ),
     &(kmap(778),knum_34  ),(kmap(779),knum_35  ),(kmap(780),knum_36  ),
     &(kmap(781),knum_37  ),(kmap(782),knum_38  ),(kmap(783),knum_39  ),
     &(kmap(784),knum_40  ),(kmap(785),knum_41  ),(kmap(786),knum_42  ),
     &(kmap(787),knum_43  ),
     &(kmap(788),kph2o    ),

     &(kmap(789),kbpin),(kmap(790),klimo),(kmap(791),kmono),
     &(kmap(792),ksesq),(kmap(793),kcg5)

c
      integer   ko1d  ,ko    ,kclo 
      integer   kcl   ,kn2o5 ,kno3 
      integer   koh   ,kho2  ,kc2o3
      integer   kxo2  ,kxo2n ,kto2 
      integer   kror  ,kcro  ,kro2r
      integer   kr2o2 ,kro2n ,kcco3
      integer   krco3 ,kmco3 ,kbzco
      integer   kcxo2 ,khco3 ,ktbuo
      integer   kbzo  ,kbzno
c
      equivalence (krad(1), ko1d ), (krad(2), ko   ), (krad(3), kclo ),
     &            (krad(4), kcl  ), (krad(5), kn2o5), (krad(6), kno3 ),
     &            (krad(7), koh  ), (krad(8), kho2 ), (krad(9), kc2o3),
     &            (krad(10),kxo2 ), (krad(11),kxo2n), (krad(12),kto2 ),
     &            (krad(13),kror ), (krad(14),kcro ), (krad(15),kro2r),
     &            (krad(16),kr2o2), (krad(17),kro2n), (krad(18),kcco3),
     &            (krad(19),krco3), (krad(20),kmco3), (krad(21),kbzco),
     &            (krad(22),kcxo2), (krad(23),khco3), (krad(24),ktbuo),
     &            (krad(25),kbzo ), (krad(26),kbzno)
c
c-----------------------------------------------------------------------
c     Variables for chemistry lookup tables:
c
c     tempr  -- temperature table
c     presr  -- pressure table
c     rktbl  -- temperature/pressure-dependent rate constant table
c     htint  -- height AGL table
c     zenint -- zenith angle table
c     prkn   -- reaction rate table
c-----------------------------------------------------------------------
c      
      common /tables/ tempr(NTEMPR), presr(NPRESR),
     &                rktbl(MXRXN,NTEMPR,NPRESR),
     &                htint(NHGHT), zenint(NZEN),
     &                prkn(NZEN,MXPHT1,NHGHT,NHAZE,NALB,NOZN)
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henry0   -- Henry's Law constant at STP (molar/atm)
c     tfact    -- Temperature dependence of Henry's Law constant (1/K)
c     diffrat  -- Species diffusivity
c     f0       -- Species reactivity parameter
c     rscale   -- Species scaling factor for surface resistance
c     henso20  -- Henry's Law constant at STP for SO2 (molar/atm)
c     tfactso2 -- Temperature dependence of SO2 Henry's Law constant (1/K)
c     nbin     -- Number of aerosol size bins
c     roprt    -- Aerosol density (g/m3)
c     dcut     -- Aerosol size bin cut points (um)
c     cwmin    -- Minimum cloud water threshold (g/m3)
c     tamin    -- Cloud water freezing threshold (K)
c-----------------------------------------------------------------------
c
      real cwmin,tamin
      common /depchm/ henry0(MXSPEC),tfact(MXSPEC),diffrat(MXSPEC),
     &                f0(MXSPEC),rscale(MXSPEC),henso20,tfactso2,cwmin,
     &                tamin
      common /aerochm/ nbin,roprt(MXSPEC),dcut(MXSPEC,2)
c
c-----------------------------------------------------------------------
c     Pointers used to lookup pig chemistry rate constants
c
c     ipigrxn  -- pointers to the nine reactions
c                 (1)   NO2 + O3 -> NO3
c                 (2)         O3 -> O(1D)
c                 (3)      O(1D) -> O(3P)
c                 (4)      O(1D) -> 2 OH
c                 (5)  NO3 + NO2 -> NO + NO2
c                 (6)  NO3 + NO2 -> N2O5
c                 (7) N2O5 + H2O -> 2 HNO3
c                 (8)       N2O5 -> NO3 + NO2
c                 (9)    NO + NO -> 2 NO2
c
      common /pigrxn/ ipigrxn(9)
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines 
c    
c     grd_time   -- time for calls to aerosol routines for each grid
c     date_aer   -- Julian date of current grd_time for each grid
c     dtaero     -- time interval between calls to aerosol routines
c     aero_dt    -- time between calls to aerosol routines for each grid
c     dt_aero    -- user input time interval between calls to aerosol routines
c
      real grd_time(MXGRID)
      real aero_dt(MXGRID)
      real dtaero
      integer date_aer(MXGRID)
c
      common /aero_t/ grd_time,aero_dt,date_aer,dtaero,dt_aero


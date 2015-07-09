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
      integer   kxyl  ,kamine     
      
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
      integer   ksoa4_9  ,ksoa4_10 ,kpoc_1
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
     &            (kmap(91),ksoa3), (kmap(92),ksoa4), (kmap(93),ksulf),
     &            (kmap(94),kterp), (kmap(95),ktol ), (kmap(96),kxn  ),
     &            (kmap(97),kxyl ), (kmap(98),kamine)
      equivalence (kmap(99), ksoa1_1  ),(kmap(100),ksoa1_2  ),  
     &(kmap(101),ksoa1_3  ),(kmap(102),ksoa1_4  ),(kmap(103),ksoa1_5  ),
     &(kmap(104),ksoa1_6  ),(kmap(105),ksoa1_7  ),(kmap(106),ksoa1_8  ),
     &(kmap(107),ksoa1_9  ),(kmap(108),ksoa1_10 ),(kmap(109),ksoa1_11 ),
     &(kmap(110),ksoa1_12 ),(kmap(111),ksoa1_13 ),(kmap(112),ksoa1_14 ),
     &(kmap(113),ksoa1_15 ),(kmap(114),ksoa1_16 ),(kmap(115),ksoa1_17 ),
     &(kmap(116),ksoa1_18 ),(kmap(117),ksoa1_19 ),(kmap(118),ksoa1_20 ),
     &(kmap(119),ksoa1_21 ),(kmap(120),ksoa1_22 ),(kmap(121),ksoa1_23 ),
     &(kmap(122),ksoa1_24 ),(kmap(123),ksoa1_25 ),(kmap(124),ksoa1_26 ),
     &(kmap(125),ksoa1_27 ),(kmap(126),ksoa1_28 ),(kmap(127),ksoa1_29 ),
     &(kmap(128),ksoa1_30 ),(kmap(129),ksoa1_31 ),(kmap(130),ksoa1_32 ),
     &(kmap(131),ksoa1_33 ),(kmap(132),ksoa1_34 ),(kmap(133),ksoa1_35 ),
     &(kmap(134),ksoa1_36 ),(kmap(135),ksoa1_37 ),(kmap(136),ksoa1_38 ),
     &(kmap(137),ksoa1_39 ),(kmap(138),ksoa1_40 ),(kmap(139),ksoa1_41 ),
     &(kmap(140),ksoa1_42 ),(kmap(141),ksoa1_43 ),
     &(kmap(142),ksoa2_1  ),
     &(kmap(143),ksoa2_2  ),(kmap(144),ksoa2_3  ),(kmap(145),ksoa2_4  ),
     &(kmap(146),ksoa2_5  ),(kmap(147),ksoa2_6  ),(kmap(148),ksoa2_7  ),
     &(kmap(149),ksoa2_8  ),(kmap(150),ksoa2_9  ),(kmap(151),ksoa2_10 ),
     &(kmap(152),ksoa2_11 ),(kmap(153),ksoa2_12 ),(kmap(154),ksoa2_13 ),
     &(kmap(155),ksoa2_14 ),(kmap(156),ksoa2_15 ),(kmap(157),ksoa2_16 ),
     &(kmap(158),ksoa2_17 ),(kmap(159),ksoa2_18 ),(kmap(160),ksoa2_19 ),
     &(kmap(161),ksoa2_20 ),(kmap(162),ksoa2_21 ),(kmap(163),ksoa2_22 ),
     &(kmap(164),ksoa2_23 ),(kmap(165),ksoa2_24 ),(kmap(166),ksoa2_25 ),
     &(kmap(167),ksoa2_26 ),(kmap(168),ksoa2_27 ),(kmap(169),ksoa2_28 ),
     &(kmap(170),ksoa2_29 ),(kmap(171),ksoa2_30 ),(kmap(172),ksoa2_31 ),
     &(kmap(173),ksoa2_32 ),(kmap(174),ksoa2_33 ),(kmap(175),ksoa2_34 ),
     &(kmap(176),ksoa2_35 ),(kmap(177),ksoa2_36 ),(kmap(178),ksoa2_37 ),
     &(kmap(179),ksoa2_38 ),(kmap(180),ksoa2_39 ),(kmap(181),ksoa2_40 ),
     &(kmap(182),ksoa2_41 ),(kmap(183),ksoa2_42 ),(kmap(184),ksoa2_43 ),
     &(kmap(185),ksoa3_1  ),(kmap(186),ksoa3_2  ),(kmap(187),ksoa3_3  ),
     &(kmap(188),ksoa3_4  ),(kmap(189),ksoa3_5  ),(kmap(190),ksoa3_6  ),
     &(kmap(191),ksoa3_7  ),(kmap(192),ksoa3_8  ),(kmap(193),ksoa3_9  ),
     &(kmap(194),ksoa3_10 ),(kmap(195),ksoa3_11 ),(kmap(196),ksoa3_12 ),
     &(kmap(197),ksoa3_13 ),(kmap(198),ksoa3_14 ),(kmap(199),ksoa3_15 ),
     &(kmap(200),ksoa3_16 ),(kmap(201),ksoa3_17 ),(kmap(202),ksoa3_18 ),
     &(kmap(203),ksoa3_19 ),(kmap(204),ksoa3_20 ),(kmap(205),ksoa3_21 ),
     &(kmap(206),ksoa3_22 ),(kmap(207),ksoa3_23 ),(kmap(208),ksoa3_24 ),
     &(kmap(209),ksoa3_25 ),(kmap(210),ksoa3_26 ),(kmap(211),ksoa3_27 ),
     &(kmap(212),ksoa3_28 ),(kmap(213),ksoa3_29 ),(kmap(214),ksoa3_30 ),
     &(kmap(215),ksoa3_31 ),(kmap(216),ksoa3_32 ),(kmap(217),ksoa3_33 ),
     &(kmap(218),ksoa3_34 ),(kmap(219),ksoa3_35 ),(kmap(220),ksoa3_36 ),
     &(kmap(221),ksoa3_37 ),(kmap(222),ksoa3_38 ),(kmap(223),ksoa3_39 ),
     &(kmap(224),ksoa3_40 ),(kmap(225),ksoa3_41 ),(kmap(226),ksoa3_42 ),
     &(kmap(227),ksoa3_43 ),
     &(kmap(228),ksoa4_1  ),(kmap(229),ksoa4_2  ),
     &(kmap(230),ksoa4_3  ),(kmap(231),ksoa4_4  ),(kmap(232),ksoa4_5  ),
     &(kmap(233),ksoa4_6  ),(kmap(234),ksoa4_7  ),(kmap(235),ksoa4_8  ),
     &(kmap(236),ksoa4_9  ),(kmap(237),ksoa4_10 ),(kmap(238),ksoa4_11 ),
     &(kmap(239),ksoa4_12 ),(kmap(240),ksoa4_13 ),(kmap(241),ksoa4_14 ),
     &(kmap(242),ksoa4_15 ),(kmap(243),ksoa4_16 ),(kmap(244),ksoa4_17 ),
     &(kmap(245),ksoa4_18 ),(kmap(246),ksoa4_19 ),(kmap(247),ksoa4_20 ),
     &(kmap(248),ksoa4_21 ),(kmap(249),ksoa4_22 ),(kmap(250),ksoa4_23 ),
     &(kmap(251),ksoa4_24 ),(kmap(252),ksoa4_25 ),(kmap(253),ksoa4_26 ),
     &(kmap(254),ksoa4_27 ),(kmap(255),ksoa4_28 ),(kmap(256),ksoa4_29 ),
     &(kmap(257),ksoa4_30 ),(kmap(258),ksoa4_31 ),(kmap(259),ksoa4_32 ),
     &(kmap(260),ksoa4_33 ),(kmap(261),ksoa4_34 ),(kmap(262),ksoa4_35 ),
     &(kmap(263),ksoa4_36 ),(kmap(264),ksoa4_37 ),(kmap(265),ksoa4_38 ),
     &(kmap(266),ksoa4_39 ),(kmap(267),ksoa4_40 ),(kmap(268),ksoa4_41 ),
     &(kmap(269),ksoa4_42 ),(kmap(270),ksoa4_43 ),
     &(kmap(271),kpoc_1   ),
     &(kmap(272),kpoc_2   ),(kmap(273),kpoc_3   ),(kmap(274),kpoc_4   ),
     &(kmap(275),kpoc_5   ),(kmap(276),kpoc_6   ),(kmap(277),kpoc_7   ),
     &(kmap(278),kpoc_8   ),(kmap(279),kpoc_9   ),(kmap(280),kpoc_10  ),
     &(kmap(281),kpoc_11  ),(kmap(282),kpoc_12  ),(kmap(283),kpoc_13  ),
     &(kmap(284),kpoc_14  ),(kmap(285),kpoc_15  ),(kmap(286),kpoc_16  ),
     &(kmap(287),kpoc_17  ),(kmap(288),kpoc_18  ),(kmap(289),kpoc_19  ),
     &(kmap(290),kpoc_20  ),(kmap(291),kpoc_21  ),(kmap(292),kpoc_22  ),
     &(kmap(293),kpoc_23  ),(kmap(294),kpoc_24  ),(kmap(295),kpoc_25  ),
     &(kmap(296),kpoc_26  ),(kmap(297),kpoc_27  ),(kmap(298),kpoc_28  ),
     &(kmap(299),kpoc_29  ),(kmap(300),kpoc_30  ),(kmap(301),kpoc_31  ),
     &(kmap(302),kpoc_32  ),(kmap(303),kpoc_33  ),(kmap(304),kpoc_34  ),
     &(kmap(305),kpoc_35  ),(kmap(306),kpoc_36  ),(kmap(307),kpoc_37  ),
     &(kmap(308),kpoc_38  ),(kmap(309),kpoc_39  ),(kmap(310),kpoc_40  ),
     &(kmap(311),kpoc_41  ),(kmap(312),kpoc_42  ),(kmap(313),kpoc_43  ),
     &(kmap(314),kpec_1   ),(kmap(315),kpec_2   ),(kmap(316),kpec_3   ),
     &(kmap(317),kpec_4   ),(kmap(318),kpec_5   ),(kmap(319),kpec_6   ),
     &(kmap(320),kpec_7   ),(kmap(321),kpec_8   ),(kmap(322),kpec_9   ),
     &(kmap(323),kpec_10  ),(kmap(324),kpec_11  ),(kmap(325),kpec_12  ),
     &(kmap(326),kpec_13  ),(kmap(327),kpec_14  ),(kmap(328),kpec_15  ),
     &(kmap(329),kpec_16  ),(kmap(330),kpec_17  ),(kmap(331),kpec_18  ),
     &(kmap(332),kpec_19  ),(kmap(333),kpec_20  ),(kmap(334),kpec_21  ),
     &(kmap(335),kpec_22  ),(kmap(336),kpec_23  ),(kmap(337),kpec_24  ),
     &(kmap(338),kpec_25  ),(kmap(339),kpec_26  ),(kmap(340),kpec_27  ),
     &(kmap(341),kpec_28  ),(kmap(342),kpec_29  ),(kmap(343),kpec_30  ),
     &(kmap(344),kpec_31  ),(kmap(345),kpec_32  ),(kmap(346),kpec_33  ),
     &(kmap(347),kpec_34  ),(kmap(348),kpec_35  ),(kmap(349),kpec_36  ),
     &(kmap(350),kpec_37  ),(kmap(351),kpec_38  ),(kmap(352),kpec_39  ),
     &(kmap(353),kpec_40  ),(kmap(354),kpec_41  ),(kmap(355),kpec_42  ),
     &(kmap(356),kpec_43  ),
     &(kmap(357),kcrust_1 ),(kmap(358),kcrust_2 ),
     &(kmap(359),kcrust_3 ),(kmap(360),kcrust_4 ),(kmap(361),kcrust_5 ),
     &(kmap(362),kcrust_6 ),(kmap(363),kcrust_7 ),(kmap(364),kcrust_8 ),
     &(kmap(365),kcrust_9 ),(kmap(366),kcrust_10),(kmap(367),kcrust_11),
     &(kmap(368),kcrust_12),(kmap(369),kcrust_13),(kmap(370),kcrust_14),
     &(kmap(371),kcrust_15),(kmap(372),kcrust_16),(kmap(373),kcrust_17),
     &(kmap(374),kcrust_18),(kmap(375),kcrust_19),(kmap(376),kcrust_20),
     &(kmap(377),kcrust_21),(kmap(378),kcrust_22),(kmap(379),kcrust_23),
     &(kmap(380),kcrust_24),(kmap(381),kcrust_25),(kmap(382),kcrust_26),
     &(kmap(383),kcrust_27),(kmap(384),kcrust_28),(kmap(385),kcrust_29),
     &(kmap(386),kcrust_30),(kmap(387),kcrust_31),(kmap(388),kcrust_32),
     &(kmap(389),kcrust_33),(kmap(390),kcrust_34),(kmap(391),kcrust_35),
     &(kmap(392),kcrust_36),(kmap(393),kcrust_37),(kmap(394),kcrust_38),
     &(kmap(395),kcrust_39),(kmap(396),kcrust_40),(kmap(397),kcrust_41),
     &(kmap(398),kcrust_42),(kmap(399),kcrust_43),
     &(kmap(400),kph2o_1  ),
     &(kmap(401),kph2o_2  ),(kmap(402),kph2o_3  ),(kmap(403),kph2o_4  ),
     &(kmap(404),kph2o_5  ),(kmap(405),kph2o_6  ),(kmap(406),kph2o_7  ),
     &(kmap(407),kph2o_8  ),(kmap(408),kph2o_9  ),(kmap(409),kph2o_10 ),
     &(kmap(410),kph2o_11 ),(kmap(411),kph2o_12 ),(kmap(412),kph2o_13 ),
     &(kmap(413),kph2o_14 ),(kmap(414),kph2o_15 ),(kmap(415),kph2o_16 ),
     &(kmap(416),kph2o_17 ),(kmap(417),kph2o_18 ),(kmap(418),kph2o_19 ),
     &(kmap(419),kph2o_20 ),(kmap(420),kph2o_21 ),(kmap(421),kph2o_22 ),
     &(kmap(422),kph2o_23 ),(kmap(423),kph2o_24 ),(kmap(424),kph2o_25 ),
     &(kmap(425),kph2o_26 ),(kmap(426),kph2o_27 ),(kmap(427),kph2o_28 ),
     &(kmap(428),kph2o_29 ),(kmap(429),kph2o_30 ),(kmap(430),kph2o_31 ),
     &(kmap(431),kph2o_32 ),(kmap(432),kph2o_33 ),(kmap(433),kph2o_34 ),
     &(kmap(434),kph2o_35 ),(kmap(435),kph2o_36 ),(kmap(436),kph2o_37 ),
     &(kmap(437),kph2o_38 ),(kmap(438),kph2o_39 ),(kmap(439),kph2o_40 ),
     &(kmap(440),kph2o_41 ),(kmap(441),kph2o_42 ),(kmap(442),kph2o_43 ),
     &(kmap(443),kpcl_1   ),(kmap(444),kpcl_2   ),(kmap(445),kpcl_3   ),
     &(kmap(446),kpcl_4   ),(kmap(447),kpcl_5   ),(kmap(448),kpcl_6   ),
     &(kmap(449),kpcl_7   ),(kmap(450),kpcl_8   ),(kmap(451),kpcl_9   ),
     &(kmap(452),kpcl_10  ),(kmap(453),kpcl_11  ),(kmap(454),kpc1_12  ),
     &(kmap(455),kpcl_13  ),(kmap(456),kpcl_14  ),(kmap(457),kpcl_15  ),
     &(kmap(458),kpcl_16  ),(kmap(459),kpcl_17  ),(kmap(460),kpcl_18  ),
     &(kmap(461),kpcl_19  ),(kmap(462),kpcl_20  ),(kmap(463),kpcl_21  ),
     &(kmap(464),kpcl_22  ),(kmap(465),kpcl_23  ),(kmap(466),kpc1_24  ),
     &(kmap(467),kpcl_25  ),(kmap(468),kpcl_26  ),(kmap(469),kpc1_27  ),
     &(kmap(470),kpcl_28  ),(kmap(471),kpcl_29  ),(kmap(472),kpcl_30  ),
     &(kmap(473),kpcl_31  ),(kmap(474),kpcl_32  ),(kmap(475),kpcl_33  ),
     &(kmap(476),kpcl_34  ),(kmap(477),kpcl_35  ),(kmap(478),kpcl_36  ),
     &(kmap(479),kpcl_37  ),(kmap(480),kpcl_38  ),(kmap(481),kpc1_39  ),
     &(kmap(482),kpcl_40  ),(kmap(483),kpcl_41  ),(kmap(484),kpc1_42  ),
     &(kmap(485),kpcl_43  ),
     &(kmap(486),kna_1    ),(kmap(487),kna_2    ),
     &(kmap(488),kna_3    ),(kmap(489),kna_4    ),(kmap(490),kna_5    ),
     &(kmap(491),kna_6    ),(kmap(492),kna_7    ),(kmap(493),kna_8    ),
     &(kmap(494),kna_9    ),(kmap(495),kna_10   ),(kmap(496),kna_11   ),
     &(kmap(497),kna_12   ),(kmap(498),kna_13   ),(kmap(499),kna_14   ),
     &(kmap(500),kna_15   ),(kmap(501),kna_16   ),(kmap(502),kna_17   ),
     &(kmap(503),kna_18   ),(kmap(504),kna_19   ),(kmap(505),kna_20   ),
     &(kmap(506),kna_21   ),(kmap(507),kna_22   ),(kmap(508),kna_23   ),
     &(kmap(509),kna_24   ),(kmap(510),kna_25   ),(kmap(511),kna_26   ),
     &(kmap(512),kna_27   ),(kmap(513),kna_28   ),(kmap(514),kna_29   ),
     &(kmap(515),kna_30   ),(kmap(516),kna_31   ),(kmap(517),kna_32   ),
     &(kmap(518),kna_33   ),(kmap(519),kna_34   ),(kmap(520),kna_35   ),
     &(kmap(521),kna_36   ),(kmap(522),kna_37   ),(kmap(523),kna_38   ),
     &(kmap(524),kna_39   ),(kmap(525),kna_40   ),(kmap(526),kna_41   ),
     &(kmap(527),kna_42   ),(kmap(528),kna_43   ),
     &(kmap(529),kpnh4_1  ),
     &(kmap(530),kpnh4_2  ),(kmap(531),kpnh4_3  ),(kmap(532),kpnh4_4  ),
     &(kmap(533),kpnh4_5  ),(kmap(534),kpnh4_6  ),(kmap(535),kpnh4_7  ),
     &(kmap(536),kpnh4_8  ),(kmap(537),kpnh4_9  ),(kmap(538),kpnh4_10 ),
     &(kmap(539),kpnh4_11 ),(kmap(540),kpnh4_12 ),(kmap(541),kpnh4_13 ),
     &(kmap(542),kpnh4_14 ),(kmap(543),kpnh4_15 ),(kmap(544),kpnh4_16 ),
     &(kmap(545),kpnh4_17 ),(kmap(546),kpnh4_18 ),(kmap(547),kpnh4_19 ),
     &(kmap(548),kpnh4_20 ),(kmap(549),kpnh4_21 ),(kmap(550),kpnh4_22 ),
     &(kmap(551),kpnh4_23 ),(kmap(552),kpnh4_24 ),(kmap(553),kpnh4_25 ),
     &(kmap(554),kpnh4_26 ),(kmap(555),kpnh4_27 ),(kmap(556),kpnh4_28 ),
     &(kmap(557),kpnh4_29 ),(kmap(558),kpnh4_30 ),(kmap(559),kpnh4_31 ),
     &(kmap(560),kpnh4_32 ),(kmap(561),kpnh4_33 ),(kmap(562),kpnh4_34 ),
     &(kmap(563),kpnh4_35 ),(kmap(564),kpnh4_36 ),(kmap(565),kpnh4_37 ),
     &(kmap(566),kpnh4_38 ),(kmap(567),kpnh4_39 ),(kmap(568),kpnh4_40 ),
     &(kmap(569),kpnh4_41 ),(kmap(570),kpnh4_42 ),(kmap(571),kpnh4_43 ),
     &(kmap(572),kpno3_1  ),(kmap(573),kpno3_2  ),(kmap(574),kpno3_3  ),
     &(kmap(575),kpno3_4  ),(kmap(576),kpno3_5  ),(kmap(577),kpno3_6  ),
     &(kmap(578),kpno3_7  ),(kmap(579),kpno3_8  ),(kmap(580),kpno3_9  ),
     &(kmap(581),kpno3_10 ),(kmap(582),kpno3_11 ),(kmap(583),kpno3_12 ),
     &(kmap(584),kpno3_13 ),(kmap(585),kpno3_14 ),(kmap(586),kpno3_15 ),
     &(kmap(587),kpno3_16 ),(kmap(588),kpno3_17 ),(kmap(589),kpno3_18 ),
     &(kmap(590),kpno3_19 ),(kmap(591),kpno3_20 ),(kmap(592),kpno3_21 ),
     &(kmap(593),kpno3_22 ),(kmap(594),kpno3_23 ),(kmap(595),kpno3_24 ),
     &(kmap(596),kpno3_25 ),(kmap(597),kpno3_26 ),(kmap(598),kpno3_27 ),
     &(kmap(599),kpno3_28 ),(kmap(600),kpno3_29 ),(kmap(601),kpno3_30 ),
     &(kmap(602),kpno3_31 ),(kmap(603),kpno3_32 ),(kmap(604),kpno3_33 ),
     &(kmap(605),kpno3_34 ),(kmap(606),kpno3_35 ),(kmap(607),kpno3_36 ),
     &(kmap(608),kpno3_37 ),(kmap(609),kpno3_38 ),(kmap(610),kpno3_39 ),
     &(kmap(611),kpno3_40 ),(kmap(612),kpno3_41 ),(kmap(613),kpno3_42 ),
     &(kmap(614),kpno3_43),
     &(kmap(615),kpso4_1  ),(kmap(616),kpso4_2  ),
     &(kmap(617),kpso4_3  ),(kmap(618),kpso4_4  ),(kmap(619),kpso4_5  ),
     &(kmap(620),kpso4_6  ),(kmap(621),kpso4_7  ),(kmap(622),kpso4_8  ),
     &(kmap(623),kpso4_9  ),(kmap(624),kpso4_10 ),(kmap(625),kpso4_11 ),
     &(kmap(626),kpso4_12 ),(kmap(627),kpso4_13 ),(kmap(628),kpso4_14 ),
     &(kmap(629),kpso4_15 ),(kmap(630),kpso4_16 ),(kmap(631),kpso4_17 ),
     &(kmap(632),kpso4_18 ),(kmap(633),kpso4_19 ),(kmap(634),kpso4_20 ),
     &(kmap(635),kpso4_21 ),(kmap(636),kpso4_22 ),(kmap(637),kpso4_23 ),
     &(kmap(638),kpso4_24 ),(kmap(639),kpso4_25 ),(kmap(640),kpso4_26 ),
     &(kmap(641),kpso4_27 ),(kmap(642),kpso4_28 ),(kmap(643),kpso4_29 ),
     &(kmap(644),kpso4_30 ),(kmap(645),kpso4_31 ),(kmap(646),kpso4_32 ),
     &(kmap(647),kpso4_33 ),(kmap(648),kpso4_34 ),(kmap(649),kpso4_35 ),
     &(kmap(650),kpso4_36 ),(kmap(651),kpso4_37 ),(kmap(652),kpso4_38 ),
     &(kmap(653),kpso4_39 ),(kmap(654),kpso4_40 ),(kmap(655),kpso4_41 ),
     &(kmap(656),kpso4_42 ),(kmap(657),kpso4_43 ),
     &(kmap(658),kpamine_1 ),(kmap(659),kpamine_2 ),(kmap(660),kpamine_3 ),
     &(kmap(661),kpamine_4 ),(kmap(662),kpamine_5 ),(kmap(663),kpamine_6 ),
     &(kmap(664),kpamine_7 ),(kmap(665),kpamine_8 ),(kmap(666),kpamine_9 ),
     &(kmap(667),kpamine_10),(kmap(668),kpamine_11),(kmap(669),kpamine_12),
     &(kmap(670),kpamine_13),(kmap(671),kpamine_14),(kmap(672),kpamine_15),
     &(kmap(673),kpamine_16),(kmap(674),kpamine_17),(kmap(675),kpamine_18),
     &(kmap(676),kpamine_19),(kmap(677),kpamine_20),(kmap(678),kpamine_21),
     &(kmap(679),kpamine_22),(kmap(680),kpamine_23),(kmap(681),kpamine_24),
     &(kmap(682),kpamine_25),(kmap(683),kpamine_26),(kmap(684),kpamine_27),
     &(kmap(685),kpamine_28),(kmap(686),kpamine_29),(kmap(687),kpamine_30),
     &(kmap(688),kpamine_31),(kmap(689),kpamine_32),(kmap(690),kpamine_33),
     &(kmap(691),kpamine_34),(kmap(692),kpamine_35),(kmap(693),kpamine_36),
     &(kmap(694),kpamine_37),(kmap(695),kpamine_38),(kmap(696),kpamine_39),
     &(kmap(697),kpamine_40),(kmap(698),kpamine_41),(kmap(699),kpamine_42),
     &(kmap(700),kpamine_43),
     &(kmap(701),knum_1   ),(kmap(702),knum_2   ),(kmap(703),knum_3   ),
     &(kmap(704),knum_4   ),(kmap(705),knum_5   ),(kmap(706),knum_6   ),
     &(kmap(707),knum_7   ),(kmap(708),knum_8   ),(kmap(709),knum_9   ),
     &(kmap(710),knum_10  ),(kmap(711),knum_11  ),(kmap(712),knum_12  ),
     &(kmap(713),knum_13  ),(kmap(714),knum_14  ),(kmap(715),knum_15  ),
     &(kmap(716),knum_16  ),(kmap(717),knum_17  ),(kmap(718),knum_18  ),
     &(kmap(719),knum_19  ),(kmap(720),knum_20  ),(kmap(721),knum_21  ),
     &(kmap(722),knum_22  ),(kmap(723),knum_23  ),(kmap(724),knum_24  ),
     &(kmap(725),knum_25  ),(kmap(726),knum_26  ),(kmap(727),knum_27  ),
     &(kmap(728),knum_28  ),(kmap(729),knum_29  ),(kmap(730),knum_30  ),
     &(kmap(731),knum_31  ),(kmap(732),knum_32  ),(kmap(733),knum_33  ),
     &(kmap(734),knum_34  ),(kmap(735),knum_35  ),(kmap(736),knum_36  ),
     &(kmap(737),knum_37  ),(kmap(738),knum_38  ),(kmap(739),knum_39  ),
     &(kmap(740),knum_40  ),(kmap(741),knum_41  ),(kmap(742),knum_42  ),
     &(kmap(743),knum_43  ),(kmap(744),kph2o    )
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


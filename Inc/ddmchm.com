c-----CAMx v4.02 030709
c
c     DDMCHM.COM sets species pointers for the DDM chemistry
c     These equivalences must be consistent with the internal
c     species lists defined in data statements in READCHM
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c               
c 
      common /lname/ lmap(NSPNAM), lrad(NRADNM)
c
      integer   lno   ,lno2  ,lo3  
      integer   lpan  ,lcres ,lpan2
      integer   lmpan ,lpbzn ,lnphe
      integer   lrno3 ,ldcb2 ,ldcb3
      integer   lhno4 ,lacet ,lald2
      integer   lalk1 ,lalk2 ,lalk3
      integer   lalk4 ,lalk5 ,laro1
      integer   laro2 ,lbacl ,lbald
      integer   lbcl1 ,lbcl2 ,lbuta
      integer   lccho ,lccrs ,lcg1 
      integer   lcg2  ,lcg3  ,lcg4 
      integer   lcl2  ,lco   ,lco2h
      integer   lco3h ,lcooh ,lcprm
      integer   ldcb1 ,leth  ,lethe
      integer   letoh ,lfcrs ,lfmcl
      integer   lform ,lfprm ,lgly 
      integer   lh2o2 ,lhc2h ,lhcho
      integer   lhcl  ,lhono ,lhno3
      integer   lho2h ,lhocl ,licl1
      integer   licl2 ,lisop ,lispd
      integer   lmek  ,lmeoh ,lmeth
      integer   lmgly ,lmvk  ,lna  
      integer   lnh3  ,lntr  ,lnxoy
      integer   lole  ,lole1 ,lole2
      integer   lopen ,lpar  ,lpcl 
      integer   lpec  ,lphen ,lpna 
      integer   lpnh4 ,lpno3 ,lpoa 
      integer   lprod ,lpso4 ,lrc2h
      integer   lrc3h ,lrcho ,lrooh
      integer   lso2  ,lsoa1 ,lsoa2
      integer   lsoa3 ,lsoa4 ,lsulf
      integer   lterp ,ltol  ,lxn  
      integer   lxyl,  lamine
      integer   lsoa1_1  ,lsoa1_2
      integer   lsoa1_3  ,lsoa1_4  ,lsoa1_5
      integer   lsoa1_6  ,lsoa1_7  ,lsoa1_8
      integer   lsoa1_9  ,lsoa1_10 ,lsoa2_1
      integer   lsoa2_2  ,lsoa2_3  ,lsoa2_4
      integer   lsoa2_5  ,lsoa2_6  ,lsoa2_7
      integer   lsoa2_8  ,lsoa2_9  ,lsoa2_10
      integer   lsoa3_1  ,lsoa3_2  ,lsoa3_3
      integer   lsoa3_4  ,lsoa3_5  ,lsoa3_6
      integer   lsoa3_7  ,lsoa3_8  ,lsoa3_9
      integer   lsoa3_10 ,lsoa4_1  ,lsoa4_2
      integer   lsoa4_3  ,lsoa4_4  ,lsoa4_5
      integer   lsoa4_6  ,lsoa4_7  ,lsoa4_8
      integer   lsoa4_9  ,lsoa4_10 ,lpoc_1
      integer   lpoc_2   ,lpoc_3   ,lpoc_4
      integer   lpoc_5   ,lpoc_6   ,lpoc_7
      integer   lpoc_8   ,lpoc_9   ,lpoc_10
      integer   lpec_1   ,lpec_2   ,lpec_3
      integer   lpec_4   ,lpec_5   ,lpec_6
      integer   lpec_7   ,lpec_8   ,lpec_9
      integer   lpec_10  ,lcrust_1 ,lcrust_2
      integer   lcrust_3 ,lcrust_4 ,lcrust_5
      integer   lcrust_6 ,lcrust_7 ,lcrust_8
      integer   lcrust_9 ,lcrust_10,lph2o_1
      integer   lph2o_2  ,lph2o_3  ,lph2o_4
      integer   lph2o_5  ,lph2o_6  ,lph2o_7
      integer   lph2o_8  ,lph2o_9  ,lph2o_10
      integer   lpcl_1   ,lpcl_2   ,lpcl_3
      integer   lpcl_4   ,lpcl_5   ,lpcl_6
      integer   lpcl_7   ,lpcl_8   ,lpcl_9
      integer   lpcl_10  ,lna_1    ,lna_2
      integer   lna_3    ,lna_4    ,lna_5
      integer   lna_6    ,lna_7    ,lna_8
      integer   lna_9    ,lna_10   ,lpnh4_1
      integer   lpnh4_2  ,lpnh4_3  ,lpnh4_4
      integer   lpnh4_5  ,lpnh4_6  ,lpnh4_7
      integer   lpnh4_8  ,lpnh4_9  ,lpnh4_10
      integer   lpno3_1  ,lpno3_2  ,lpno3_3
      integer   lpno3_4  ,lpno3_5  ,lpno3_6
      integer   lpno3_7  ,lpno3_8  ,lpno3_9
      integer   lpno3_10 ,lpso4_1  ,lpso4_2
      integer   lpso4_3  ,lpso4_4  ,lpso4_5
      integer   lpso4_6  ,lpso4_7  ,lpso4_8
      integer   lpso4_9  ,lpso4_10 ,lph2o
c
      equivalence (lmap(1), lno  ), (lmap(2), lno2 ), (lmap(3), lo3  ),
     &            (lmap(4), lpan ), (lmap(5), lcres), (lmap(6), lpan2),
     &            (lmap(7), lmpan), (lmap(8), lpbzn), (lmap(9), lnphe),
     &            (lmap(10),lrno3), (lmap(11),ldcb2), (lmap(12),ldcb3),
     &            (lmap(13),lhno4), (lmap(14),lacet), (lmap(15),lald2),
     &            (lmap(16),lalk1), (lmap(17),lalk2), (lmap(18),lalk3),
     &            (lmap(19),lalk4), (lmap(20),lalk5), (lmap(21),laro1),
     &            (lmap(22),laro2), (lmap(23),lbacl), (lmap(24),lbald),
     &            (lmap(25),lbcl1), (lmap(26),lbcl2), (lmap(27),lbuta),
     &            (lmap(28),lccho), (lmap(29),lccrs), (lmap(30),lcg1 ),
     &            (lmap(31),lcg2 ), (lmap(32),lcg3 ), (lmap(33),lcg4 ),
     &            (lmap(34),lcl2 ), (lmap(35),lco  ), (lmap(36),lco2h),
     &            (lmap(37),lco3h), (lmap(38),lcooh), (lmap(39),lcprm),
     &            (lmap(40),ldcb1), (lmap(41),leth ), (lmap(42),lethe),
     &            (lmap(43),letoh), (lmap(44),lfcrs), (lmap(45),lfmcl),
     &            (lmap(46),lform), (lmap(47),lfprm), (lmap(48),lgly ),
     &            (lmap(49),lh2o2), (lmap(50),lhc2h), (lmap(51),lhcho),
     &            (lmap(52),lhcl ), (lmap(53),lhono), (lmap(54),lhno3),
     &            (lmap(55),lho2h), (lmap(56),lhocl), (lmap(57),licl1),
     &            (lmap(58),licl2), (lmap(59),lisop), (lmap(60),lispd),
     &            (lmap(61),lmek ), (lmap(62),lmeoh), (lmap(63),lmeth),
     &            (lmap(64),lmgly), (lmap(65),lmvk ), (lmap(66),lna  ),
     &            (lmap(67),lnh3 ), (lmap(68),lntr ), (lmap(69),lnxoy),
     &            (lmap(70),lole ), (lmap(71),lole1), (lmap(72),lole2),
     &            (lmap(73),lopen), (lmap(74),lpar ), (lmap(75),lpcl ),
     &            (lmap(76),lpec ), (lmap(77),lphen), (lmap(78),lpna ),
     &            (lmap(79),lpnh4), (lmap(80),lpno3), (lmap(81),lpoa ),
     &            (lmap(82),lprod), (lmap(83),lpso4), (lmap(84),lrc2h),
     &            (lmap(85),lrc3h), (lmap(86),lrcho), (lmap(87),lrooh),
     &            (lmap(88),lso2 ), (lmap(89),lsoa1), (lmap(90),lsoa2),
     &            (lmap(91),lsoa3), (lmap(92),lsoa4), (lmap(93),lsulf),
     &            (lmap(94),lterp), (lmap(95),ltol ), (lmap(96),lxn  ),
     &            (lmap(97),lxyl ), (lmap(98),lamine)
      equivalence (lmap( 99),lsoa1_1  ),(lmap(100),lsoa1_2  ),
     &(lmap(101),lsoa1_3  ),(lmap(102),lsoa1_4  ),(lmap(103),lsoa1_5  ),
     &(lmap(104),lsoa1_6  ),(lmap(105),lsoa1_7  ),(lmap(106),lsoa1_8  ),
     &(lmap(107),lsoa1_9  ),(lmap(108),lsoa1_10 ),(lmap(109),lsoa2_1  ),
     &(lmap(110),lsoa2_2  ),(lmap(111),lsoa2_3  ),(lmap(112),lsoa2_4  ),
     &(lmap(113),lsoa2_5  ),(lmap(114),lsoa2_6  ),(lmap(115),lsoa2_7  ),
     &(lmap(116),lsoa2_8  ),(lmap(117),lsoa2_9  ),(lmap(118),lsoa2_10 ),
     &(lmap(119),lsoa3_1  ),(lmap(120),lsoa3_2  ),(lmap(121),lsoa3_3  ),
     &(lmap(122),lsoa3_4  ),(lmap(123),lsoa3_5  ),(lmap(124),lsoa3_6  ),
     &(lmap(125),lsoa3_7  ),(lmap(126),lsoa3_8  ),(lmap(127),lsoa3_9  ),
     &(lmap(128),lsoa3_10 ),(lmap(129),lsoa4_1  ),(lmap(130),lsoa4_2  ),
     &(lmap(131),lsoa4_3  ),(lmap(132),lsoa4_4  ),(lmap(133),lsoa4_5  ),
     &(lmap(134),lsoa4_6  ),(lmap(135),lsoa4_7  ),(lmap(136),lsoa4_8  ),
     &(lmap(137),lsoa4_9  ),(lmap(138),lsoa4_10 ),(lmap(139),lpoc_1   ),
     &(lmap(140),lpoc_2   ),(lmap(141),lpoc_3   ),(lmap(142),lpoc_4   ),
     &(lmap(143),lpoc_5   ),(lmap(144),lpoc_6   ),(lmap(145),lpoc_7   ),
     &(lmap(146),lpoc_8   ),(lmap(147),lpoc_9   ),(lmap(148),lpoc_10  ),
     &(lmap(149),lpec_1   ),(lmap(150),lpec_2   ),(lmap(151),lpec_3   ),
     &(lmap(152),lpec_4   ),(lmap(153),lpec_5   ),(lmap(154),lpec_6   ),
     &(lmap(155),lpec_7   ),(lmap(156),lpec_8   ),(lmap(157),lpec_9   ),
     &(lmap(158),lpec_10  ),(lmap(159),lcrust_1 ),(lmap(160),lcrust_2 ),
     &(lmap(161),lcrust_3 ),(lmap(162),lcrust_4 ),(lmap(163),lcrust_5 ),
     &(lmap(164),lcrust_6 ),(lmap(165),lcrust_7 ),(lmap(166),lcrust_8 ),
     &(lmap(167),lcrust_9 ),(lmap(168),lcrust_10),(lmap(169),lph2o_1  ),
     &(lmap(170),lph2o_2  ),(lmap(171),lph2o_3  ),(lmap(172),lph2o_4  ),
     &(lmap(173),lph2o_5  ),(lmap(174),lph2o_6  ),(lmap(175),lph2o_7  ),
     &(lmap(176),lph2o_8  ),(lmap(177),lph2o_9  ),(lmap(178),lph2o_10 ),
     &(lmap(179),lpcl_1   ),(lmap(180),lpcl_2   ),(lmap(181),lpcl_3   ),
     &(lmap(182),lpcl_4   ),(lmap(183),lpcl_5   ),(lmap(184),lpcl_6   ),
     &(lmap(185),lpcl_7   ),(lmap(186),lpcl_8   ),(lmap(187),lpcl_9   ),
     &(lmap(188),lpcl_10  ),(lmap(189),lna_1    ),(lmap(190),lna_2    ),
     &(lmap(191),lna_3    ),(lmap(192),lna_4    ),(lmap(193),lna_5    ),
     &(lmap(194),lna_6    ),(lmap(195),lna_7    ),(lmap(196),lna_8    ),
     &(lmap(197),lna_9    ),(lmap(198),lna_10   ),(lmap(199),lpnh4_1  ),
     &(lmap(200),lpnh4_2  ),(lmap(201),lpnh4_3  ),(lmap(202),lpnh4_4  ),
     &(lmap(203),lpnh4_5  ),(lmap(204),lpnh4_6  ),(lmap(205),lpnh4_7  ),
     &(lmap(206),lpnh4_8  ),(lmap(207),lpnh4_9  ),(lmap(208),lpnh4_10 ),
     &(lmap(209),lpno3_1  ),(lmap(210),lpno3_2  ),(lmap(211),lpno3_3  ),
     &(lmap(212),lpno3_4  ),(lmap(213),lpno3_5  ),(lmap(214),lpno3_6  ),
     &(lmap(215),lpno3_7  ),(lmap(216),lpno3_8  ),(lmap(217),lpno3_9  ),
     &(lmap(218),lpno3_10 ),(lmap(219),lpso4_1  ),(lmap(220),lpso4_2  ),
     &(lmap(221),lpso4_3  ),(lmap(222),lpso4_4  ),(lmap(223),lpso4_5  ),
     &(lmap(224),lpso4_6  ),(lmap(225),lpso4_7  ),(lmap(226),lpso4_8  ),
     &(lmap(227),lpso4_9  ),(lmap(228),lpso4_10 ),(lmap(229),lph2o    )
c
      integer   lo1d  ,lo    ,lclo 
      integer   lcl   ,ln2o5 ,lno3 
      integer   loh   ,lho2  ,lc2o3
      integer   lxo2  ,lxo2n ,lto2 
      integer   lror  ,lcro  ,lro2r
      integer   lr2o2 ,lro2n ,lcco3
      integer   lrco3 ,lmco3 ,lbzco
      integer   lcxo2 ,lhco3 ,ltbuo
      integer   lbzo  ,lbzno
c
      equivalence (lrad(1), lo1d ), (lrad(2), lo   ), (lrad(3), lclo ),
     &            (lrad(4), lcl  ), (lrad(5), ln2o5), (lrad(6), lno3 ),
     &            (lrad(7), loh  ), (lrad(8), lho2 ), (lrad(9), lc2o3),
     &            (lrad(10),lxo2 ), (lrad(11),lxo2n), (lrad(12),lto2 ),
     &            (lrad(13),lror ), (lrad(14),lcro ), (lrad(15),lro2r),
     &            (lrad(16),lr2o2), (lrad(17),lro2n), (lrad(18),lcco3),
     &            (lrad(19),lrco3), (lrad(20),lmco3), (lrad(21),lbzco),
     &            (lrad(22),lcxo2), (lrad(23),lhco3), (lrad(24),ltbuo),
     &            (lrad(25),lbzo ), (lrad(26),lbzno)

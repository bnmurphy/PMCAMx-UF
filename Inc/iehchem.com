c-----CAMx v4.02 030709
c
c     IEHCHEM.COM sets species pointers for the IEH chemistry solver
c     These equivalences must be consistent with the internal   
c     species lists defined in data statements in READCHM
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
      common /iname/ imap(NSPNAM), irad(NRADNM)
c
      integer   ino   ,ino2  ,io3  
      integer   ipan  ,icres ,ipan2
      integer   impan ,ipbzn ,inphe
      integer   irno3 ,idcb2 ,idcb3
      integer   ihno4 ,iacet ,iald2
      integer   ialk1 ,ialk2 ,ialk3
      integer   ialk4 ,ialk5 ,iaro1
      integer   iaro2 ,ibacl ,ibald
      integer   ibcl1 ,ibcl2 ,ibuta
      integer   iccho ,iccrs ,icg1 
      integer   icg2  ,icg3  ,icg4 
      integer   icl2  ,ico   ,ico2h
      integer   ico3h ,icooh ,icprm
      integer   idcb1 ,ieth  ,iethe
      integer   ietoh ,ifcrs ,ifmcl
      integer   iform ,ifprm ,igly 
      integer   ih2o2 ,ihc2h ,ihcho
      integer   ihcl  ,ihono ,ihno3
      integer   iho2h ,ihocl ,iicl1
      integer   iicl2 ,iisop ,iispd
      integer   imek  ,imeoh ,imeth
      integer   imgly ,imvk  ,ina  
      integer   inh3  ,intr  ,inxoy
      integer   iole  ,iole1 ,iole2
      integer   iopen ,ipar  ,ipcl 
      integer   ipec  ,iphen ,ipna 
      integer   ipnh4 ,ipno3 ,ipoa 
      integer   iprod ,ipso4 ,irc2h
      integer   irc3h ,ircho ,irooh
      integer   iso2  ,isoa1 ,isoa2
      integer   isoa3 ,isoa4 ,isulf
      integer   iterp ,itol  ,ixn  
      integer   ixyl  ,iamine
      
      integer   isoa1_1  ,isoa1_2
      integer   isoa1_3  ,isoa1_4  ,isoa1_5
      integer   isoa1_6  ,isoa1_7  ,isoa1_8
      integer   isoa1_9  ,isoa1_10 ,isoa2_1
      integer   isoa2_2  ,isoa2_3  ,isoa2_4
      integer   isoa2_5  ,isoa2_6  ,isoa2_7
      integer   isoa2_8  ,isoa2_9  ,isoa2_10
      integer   isoa3_1  ,isoa3_2  ,isoa3_3
      integer   isoa3_4  ,isoa3_5  ,isoa3_6
      integer   isoa3_7  ,isoa3_8  ,isoa3_9
      integer   isoa3_10 ,isoa4_1  ,isoa4_2
      integer   isoa4_3  ,isoa4_4  ,isoa4_5
      integer   isoa4_6  ,isoa4_7  ,isoa4_8
      integer   isoa4_9  ,isoa4_10 ,ipoc_1
      integer   ipoc_2   ,ipoc_3   ,ipoc_4
      integer   ipoc_5   ,ipoc_6   ,ipoc_7
      integer   ipoc_8   ,ipoc_9   ,ipoc_10
      integer   ipec_1   ,ipec_2   ,ipec_3
      integer   ipec_4   ,ipec_5   ,ipec_6
      integer   ipec_7   ,ipec_8   ,ipec_9
      integer   ipec_10  ,icrust_1 ,icrust_2
      integer   icrust_3 ,icrust_4 ,icrust_5
      integer   icrust_6 ,icrust_7 ,icrust_8
      integer   icrust_9 ,icrust_10,iph2o_1
      integer   iph2o_2  ,iph2o_3  ,iph2o_4
      integer   iph2o_5  ,iph2o_6  ,iph2o_7
      integer   iph2o_8  ,iph2o_9  ,iph2o_10
      integer   ipcl_1   ,ipcl_2   ,ipcl_3
      integer   ipcl_4   ,ipcl_5   ,ipcl_6
      integer   ipcl_7   ,ipcl_8   ,ipcl_9
      integer   ipcl_10  ,ina_1    ,ina_2
      integer   ina_3    ,ina_4    ,ina_5
      integer   ina_6    ,ina_7    ,ina_8
      integer   ina_9    ,ina_10   ,ipnh4_1
      integer   ipnh4_2  ,ipnh4_3  ,ipnh4_4
      integer   ipnh4_5  ,ipnh4_6  ,ipnh4_7
      integer   ipnh4_8  ,ipnh4_9  ,ipnh4_10
      integer   ipno3_1  ,ipno3_2  ,ipno3_3
      integer   ipno3_4  ,ipno3_5  ,ipno3_6
      integer   ipno3_7  ,ipno3_8  ,ipno3_9
      integer   ipno3_10 ,ipso4_1  ,ipso4_2
      integer   ipso4_3  ,ipso4_4  ,ipso4_5
      integer   ipso4_6  ,ipso4_7  ,ipso4_8
      integer   ipso4_9  ,ipso4_10 ,iph2o
c
      equivalence (imap(1), ino  ), (imap(2), ino2 ), (imap(3), io3  ),
     &            (imap(4), ipan ), (imap(5), icres), (imap(6), ipan2),
     &            (imap(7), impan), (imap(8), ipbzn), (imap(9), inphe),
     &            (imap(10),irno3), (imap(11),idcb2), (imap(12),idcb3),
     &            (imap(13),ihno4), (imap(14),iacet), (imap(15),iald2),
     &            (imap(16),ialk1), (imap(17),ialk2), (imap(18),ialk3),
     &            (imap(19),ialk4), (imap(20),ialk5), (imap(21),iaro1),
     &            (imap(22),iaro2), (imap(23),ibacl), (imap(24),ibald),
     &            (imap(25),ibcl1), (imap(26),ibcl2), (imap(27),ibuta),
     &            (imap(28),iccho), (imap(29),iccrs), (imap(30),icg1 ),
     &            (imap(31),icg2 ), (imap(32),icg3 ), (imap(33),icg4 ),
     &            (imap(34),icl2 ), (imap(35),ico  ), (imap(36),ico2h),
     &            (imap(37),ico3h), (imap(38),icooh), (imap(39),icprm),
     &            (imap(40),idcb1), (imap(41),ieth ), (imap(42),iethe),
     &            (imap(43),ietoh), (imap(44),ifcrs), (imap(45),ifmcl),
     &            (imap(46),iform), (imap(47),ifprm), (imap(48),igly ),
     &            (imap(49),ih2o2), (imap(50),ihc2h), (imap(51),ihcho),
     &            (imap(52),ihcl ), (imap(53),ihono), (imap(54),ihno3),
     &            (imap(55),iho2h), (imap(56),ihocl), (imap(57),iicl1),
     &            (imap(58),iicl2), (imap(59),iisop), (imap(60),iispd),
     &            (imap(61),imek ), (imap(62),imeoh), (imap(63),imeth),
     &            (imap(64),imgly), (imap(65),imvk ), (imap(66),ina  ),
     &            (imap(67),inh3 ), (imap(68),intr ), (imap(69),inxoy),
     &            (imap(70),iole ), (imap(71),iole1), (imap(72),iole2),
     &            (imap(73),iopen), (imap(74),ipar ), (imap(75),ipcl ),
     &            (imap(76),ipec ), (imap(77),iphen), (imap(78),ipna ),
     &            (imap(79),ipnh4), (imap(80),ipno3), (imap(81),ipoa ),
     &            (imap(82),iprod), (imap(83),ipso4), (imap(84),irc2h),
     &            (imap(85),irc3h), (imap(86),ircho), (imap(87),irooh),
     &            (imap(88),iso2 ), (imap(89),isoa1), (imap(90),isoa2),
     &            (imap(91),isoa3), (imap(92),isoa4), (imap(93),isulf),
     &            (imap(94),iterp), (imap(95),itol ), (imap(96),ixn  ),
     &            (imap(97),ixyl ), (imap(98),iamine)
      equivalence (imap( 99),isoa1_1  ),(imap(100),isoa1_2  ),
     &(imap(101),isoa1_3  ),(imap(102),isoa1_4  ),(imap(103),isoa1_5  ),
     &(imap(104),isoa1_6  ),(imap(105),isoa1_7  ),(imap(106),isoa1_8  ),
     &(imap(107),isoa1_9  ),(imap(108),isoa1_10 ),(imap(109),isoa2_1  ),
     &(imap(110),isoa2_2  ),(imap(111),isoa2_3  ),(imap(112),isoa2_4  ),
     &(imap(113),isoa2_5  ),(imap(114),isoa2_6  ),(imap(115),isoa2_7  ),
     &(imap(116),isoa2_8  ),(imap(117),isoa2_9  ),(imap(118),isoa2_10 ),
     &(imap(119),isoa3_1  ),(imap(120),isoa3_2  ),(imap(121),isoa3_3  ),
     &(imap(122),isoa3_4  ),(imap(123),isoa3_5  ),(imap(124),isoa3_6  ),
     &(imap(125),isoa3_7  ),(imap(126),isoa3_8  ),(imap(127),isoa3_9  ),
     &(imap(128),isoa3_10 ),(imap(129),isoa4_1  ),(imap(130),isoa4_2  ),
     &(imap(131),isoa4_3  ),(imap(132),isoa4_4  ),(imap(133),isoa4_5  ),
     &(imap(134),isoa4_6  ),(imap(135),isoa4_7  ),(imap(136),isoa4_8  ),
     &(imap(137),isoa4_9  ),(imap(138),isoa4_10 ),(imap(139),ipoc_1   ),
     &(imap(140),ipoc_2   ),(imap(141),ipoc_3   ),(imap(142),ipoc_4   ),
     &(imap(143),ipoc_5   ),(imap(144),ipoc_6   ),(imap(145),ipoc_7   ),
     &(imap(146),ipoc_8   ),(imap(147),ipoc_9   ),(imap(148),ipoc_10  ),
     &(imap(149),ipec_1   ),(imap(150),ipec_2   ),(imap(151),ipec_3   ),
     &(imap(152),ipec_4   ),(imap(153),ipec_5   ),(imap(154),ipec_6   ),
     &(imap(155),ipec_7   ),(imap(156),ipec_8   ),(imap(157),ipec_9   ),
     &(imap(158),ipec_10  ),(imap(159),icrust_1 ),(imap(160),icrust_2 ),
     &(imap(161),icrust_3 ),(imap(162),icrust_4 ),(imap(163),icrust_5 ),
     &(imap(164),icrust_6 ),(imap(165),icrust_7 ),(imap(166),icrust_8 ),
     &(imap(167),icrust_9 ),(imap(168),icrust_10),(imap(169),iph2o_1  ),
     &(imap(170),iph2o_2  ),(imap(171),iph2o_3  ),(imap(172),iph2o_4  ),
     &(imap(173),iph2o_5  ),(imap(174),iph2o_6  ),(imap(175),iph2o_7  ),
     &(imap(176),iph2o_8  ),(imap(177),iph2o_9  ),(imap(178),iph2o_10 ),
     &(imap(179),ipcl_1   ),(imap(180),ipcl_2   ),(imap(181),ipcl_3   ),
     &(imap(182),ipcl_4   ),(imap(183),ipcl_5   ),(imap(184),ipcl_6   ),
     &(imap(185),ipcl_7   ),(imap(186),ipcl_8   ),(imap(187),ipcl_9   ),
     &(imap(188),ipcl_10  ),(imap(189),ina_1    ),(imap(190),ina_2    ),
     &(imap(191),ina_3    ),(imap(192),ina_4    ),(imap(193),ina_5    ),
     &(imap(194),ina_6    ),(imap(195),ina_7    ),(imap(196),ina_8    ),
     &(imap(197),ina_9    ),(imap(198),ina_10   ),(imap(199),ipnh4_1  ),
     &(imap(200),ipnh4_2  ),(imap(201),ipnh4_3  ),(imap(202),ipnh4_4  ),
     &(imap(203),ipnh4_5  ),(imap(204),ipnh4_6  ),(imap(205),ipnh4_7  ),
     &(imap(206),ipnh4_8  ),(imap(207),ipnh4_9  ),(imap(208),ipnh4_10 ),
     &(imap(209),ipno3_1  ),(imap(210),ipno3_2  ),(imap(211),ipno3_3  ),
     &(imap(212),ipno3_4  ),(imap(213),ipno3_5  ),(imap(214),ipno3_6  ),
     &(imap(215),ipno3_7  ),(imap(216),ipno3_8  ),(imap(217),ipno3_9  ),
     &(imap(218),ipno3_10 ),(imap(219),ipso4_1  ),(imap(220),ipso4_2  ),
     &(imap(221),ipso4_3  ),(imap(222),ipso4_4  ),(imap(223),ipso4_5  ),
     &(imap(224),ipso4_6  ),(imap(225),ipso4_7  ),(imap(226),ipso4_8  ),
     &(imap(227),ipso4_9  ),(imap(228),ipso4_10 ),(imap(229),iph2o    )
c
      integer   io1d  ,io    ,iclo 
      integer   icl   ,in2o5 ,ino3 
      integer   ioh   ,iho2  ,ic2o3
      integer   ixo2  ,ixo2n ,ito2 
      integer   iror  ,icro  ,iro2r
      integer   ir2o2 ,iro2n ,icco3
      integer   irco3 ,imco3 ,ibzco
      integer   icxo2 ,ihco3 ,itbuo
      integer   ibzo  ,ibzno
c
      equivalence (irad(1), io1d ), (irad(2), io   ), (irad(3), iclo ),
     &            (irad(4), icl  ), (irad(5), in2o5), (irad(6), ino3 ),
     &            (irad(7), ioh  ), (irad(8), iho2 ), (irad(9), ic2o3),
     &            (irad(10),ixo2 ), (irad(11),ixo2n), (irad(12),ito2 ),
     &            (irad(13),iror ), (irad(14),icro ), (irad(15),iro2r),
     &            (irad(16),ir2o2), (irad(17),iro2n), (irad(18),icco3),
     &            (irad(19),irco3), (irad(20),imco3), (irad(21),ibzco),
     &            (irad(22),icxo2), (irad(23),ihco3), (irad(24),itbuo),
     &            (irad(25),ibzo ), (irad(26),ibzno)

      subroutine readchm
c
c-----CAMx v4.02 030709
c
c     READCHM reads the CAMx chemistry parameter file, which defines
c     the chemical system to be simulated
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:
c        10/26/99  Added check for zero Henry's law constants to avoid a
c                  divide by 0 later in the code
c        10/20/00  Added CAMx version as first record on chemistry parameters
c                  file
c        1/9/02    Aerosol size cut points and density now defined on
c                  chemistry parameters file; removed conversion of aerosol
c                  BDNL values from ug/m3 to umol/m3
c        1/18/02   Added EC, PFIN, and PCRS to mechanism 4 species list
c        12/12/02  Expanded species list for Mechanism 4
c        3/26/03   Added surface resistance scaling factor to gas params
c
c     Input arguments: 
c        none 
c 
c     Output arguments: 
c            
c     Routines called: 
c        EXPTBL
c        KTHERM
c            
c     Called by:
c        CHMPREP
c
      include 'camx.prm'
      include 'chmstry.com'
      include 'filunit.com'
      include 'flags.com'
      include 'ddmchm.com'
      include 'iehchem.com'
      include 'diameters.inc'
c
      parameter(ncrsspc = 2)
cjgj      character*180 record
      character*380 record
      character*10 nametmp, splist(NSPNAM), radlist(NRADNM) 
      character*10 nmrad1(14), nmrad2(12), nmrad5(18)
      character*10 crsspc(ncrsspc),spcsz
      character*10 blank,camxv,camxvin
      character*10 tmpnam, tmpnam1
      integer mchgas(10), mchaero(10), mchrad(10), mchrxn(10)
      integer mchphot(10), mchfast(10), mchiessr(10)
      integer ipigcb4(9), ipigsap(9)
      integer rxntyp(MXRXN), rxnord(MXRXN), npar(7)
      real rxnpar(MXRXN,12)
      real kdum(MXRXN,3), tdum(3), pdum(3)
      integer omp_get_num_procs
c 
c-----Data that define the mechanism/solver options
c     The fast state species must come first in SPLIST
c
      data camxv  /'VERSION4.0'/
      data blank  /'BLANK     '/
c
      data splist /'NO        ','NO2       ','O3        ',
     &             'PAN       ','CRES      ','PAN2      ',
     &             'MPAN      ','PBZN      ','NPHE      ',
     &             'RNO3      ','DCB2      ','DCB3      ',
     &             'HNO4      ','ACET      ','ALD2      ',
     &             'ALK1      ','ALK2      ','ALK3      ',
     &             'ALK4      ','ALK5      ','ARO1      ',
     &             'ARO2      ','BACL      ','BALD      ',
     &             'BCL1      ','BCL2      ','BUTA      ',
     &             'CCHO      ','CCRS      ','CG1       ',
     &             'CG2       ','CG3       ','CG4       ',
     &             'CL2       ','CO        ','CO2H      ',
     &             'CO3H      ','COOH      ','CPRM      ',
     &             'DCB1      ','ETH       ','ETHE      ',
     &             'ETOH      ','FCRS      ','FMCL      ',
     &             'FORM      ','FPRM      ','GLY       ',
     &             'H2O2      ','HC2H      ','HCHO      ',
     &             'HCL       ','HONO      ','HNO3      ',
     &             'HO2H      ','HOCL      ','ICL1      ',
     &             'ICL2      ','ISOP      ','ISPD      ',
     &             'MEK       ','MEOH      ','METH      ',
     &             'MGLY      ','MVK       ','NA        ',
     &             'NH3       ','NTR       ','NXOY      ',
     &             'OLE       ','OLE1      ','OLE2      ',
     &             'OPEN      ','PAR       ','PCL       ',
     &             'PEC       ','PHEN      ','PNA       ',
     &             'PNH4      ','PNO3      ','POA       ',
     &             'PROD      ','PSO4      ','RC2H      ',
     &             'RC3H      ','RCHO      ','ROOH      ',
     &             'SO2       ','SOA1      ','SOA2      ',
     &             'SOA3      ','SOA4      ','SULF      ',
     &             'TERP      ','TOL       ','XN        ',
     &             'XYL       ','AMINE     ',
     &             'SOA1_1    ','SOA1_2    ',
     &             'SOA1_3    ','SOA1_4    ','SOA1_5    ',
     &             'SOA1_6    ','SOA1_7    ','SOA1_8    ',
     &             'SOA1_9    ','SOA1_10   ','SOA1_11   ',
     &             'SOA1_12   ','SOA1_13   ','SOA1_14   ',
     &             'SOA1_15   ','SOA1_16   ','SOA1_17   ',
     &             'SOA1_18   ','SOA1_19   ','SOA1_20   ',
     &             'SOA1_21   ','SOA1_22   ','SOA1_23   ',
     &             'SOA1_24   ','SOA1_25   ','SOA1_26   ',
     &             'SOA1_27   ','SOA1_28   ','SOA1_29   ',
     &             'SOA1_30   ','SOA1_31   ','SOA1_32   ',
     &             'SOA1_33   ','SOA1_34   ','SOA1_35   ',
     &             'SOA1_36   ','SOA1_37   ','SOA1_38   ',
     &             'SOA1_39   ','SOA1_40   ','SOA1_41   ',
     &             'SOA1_42   ','SOA1_43   ',
     &             'SOA2_1    ',
     &             'SOA2_2    ','SOA2_3    ','SOA2_4    ',
     &             'SOA2_5    ','SOA2_6    ','SOA2_7    ',
     &             'SOA2_8    ','SOA2_9    ','SOA2_10   ',
     &             'SOA2_11   ','SOA2_12   ','SOA2_13   ',
     &             'SOA2_14   ','SOA2_15   ','SOA2_16   ',
     &             'SOA2_17   ','SOA2_18   ','SOA2_19   ',
     &             'SOA2_20   ','SOA2_21   ','SOA2_22   ',
     &             'SOA2_23   ','SOA2_24   ','SOA2_25   ',
     &             'SOA2_26   ','SOA2_27   ','SOA2_28   ',
     &             'SOA2_29   ','SOA2_30   ','SOA2_31   ',
     &             'SOA2_32   ','SOA2_33   ','SOA2_34   ',
     &             'SOA2_35   ','SOA2_36   ','SOA2_37   ',
     &             'SOA2_38   ','SOA2_39   ','SOA2_40   ',
     &             'SOA2_41   ','SOA2_42   ','SOA2_43   ',
     &             'SOA3_1    ','SOA3_2    ','SOA3_3    ',
     &             'SOA3_4    ','SOA3_5    ','SOA3_6    ',
     &             'SOA3_7    ','SOA3_8    ','SOA3_9    ',
     &             'SOA3_10   ','SOA3_11   ','SOA3_12   ',
     &             'SOA3_13   ','SOA3_14   ','SOA3_15   ',
     &             'SOA3_16   ','SOA3_17   ','SOA3_18   ',
     &             'SOA3_19   ','SOA3_20   ','SOA3_21   ',
     &             'SOA3_22   ','SOA3_23   ','SOA3_24   ',
     &             'SOA3_25   ','SOA3_26   ','SOA3_27   ',
     &             'SOA3_28   ','SOA3_29   ','SOA3_30   ',
     &             'SOA3_31   ','SOA3_32   ','SOA3_33   ',
     &             'SOA3_34   ','SOA3_35   ','SOA3_36   ',
     &             'SOA3_37   ','SOA3_38   ','SOA3_39   ',
     &             'SOA3_40   ','SOA3_41   ','SOA3_42   ',
     &             'SOA3_43   ',
     &             'SOA4_1    ','SOA4_2    ',
     &             'SOA4_3    ','SOA4_4    ','SOA4_5    ',
     &             'SOA4_6    ','SOA4_7    ','SOA4_8    ',
     &             'SOA4_9    ','SOA4_10   ','SOA4_11   ',
     &             'SOA4_12   ','SOA4_13   ','SOA4_14   ',
     &             'SOA4_15   ','SOA4_16   ','SOA4_17   ',
     &             'SOA4_18   ','SOA4_19   ','SOA4_20   ',
     &             'SOA4_21   ','SOA4_22   ','SOA4_23   ',
     &             'SOA4_24   ','SOA4_25   ','SOA4_26   ',
     &             'SOA4_27   ','SOA4_28   ','SOA4_29   ',
     &             'SOA4_30   ','SOA4_31   ','SOA4_32   ',
     &             'SOA4_33   ','SOA4_34   ','SOA4_35   ',
     &             'SOA4_36   ','SOA4_37   ','SOA4_38   ',
     &             'SOA4_39   ','SOA4_40   ','SOA4_41   ',
     &             'SOA4_42   ','SOA4_43   ',
     &             'POC_1     ',
     &             'POC_2     ','POC_3     ','POC_4     ',
     &             'POC_5     ','POC_6     ','POC_7     ',
     &             'POC_8     ','POC_9     ','POC_10    ',
     &             'POC_11    ','POC_12    ','POC_13    ',
     &             'POC_14    ','POC_15    ','POC_16    ',
     &             'POC_17    ','POC_18    ','POC_19    ',
     &             'POC_20    ','POC_21    ','POC_22    ',
     &             'POC_23    ','POC_24    ','POC_25    ',
     &             'POC_26    ','POC_27    ','POC_28    ',
     &             'POC_29    ','POC_30    ','POC_31    ',
     &             'POC_32    ','POC_33    ','POC_34    ',
     &             'POC_35    ','POC_36    ','POC_37    ',
     &             'POC_38    ','POC_39    ','POC_40    ',
     &             'POC_41    ','POC_42    ','POC_43    ',
     &             'PEC_1     ','PEC_2     ','PEC_3     ',
     &             'PEC_4     ','PEC_5     ','PEC_6     ',
     &             'PEC_7     ','PEC_8     ','PEC_9     ',
     &             'PEC_10    ','PEC_11    ','PEC_12    ',
     &             'PEC_13    ','PEC_14    ','PEC_15    ',
     &             'PEC_16    ','PEC_17    ','PEC_18    ',
     &             'PEC_19    ','PEC_20    ','PEC_21    ',
     &             'PEC_22    ','PEC_23    ','PEC_24    ',
     &             'PEC_25    ','PEC_26    ','PEC_27    ',
     &             'PEC_28    ','PEC_29    ','PEC_30    ',
     &             'PEC_31    ','PEC_32    ','PEC_33    ',
     &             'PEC_34    ','PEC_35    ','PEC_36    ',
     &             'PEC_37    ','PEC_38    ','PEC_39    ',
     &             'PEC_40    ','PEC_41    ','PEC_42    ',
     &             'PEC_43    ',
     &             'CRST_1    ','CRST_2    ',
     &             'CRST_3    ','CRST_4    ','CRST_5    ',
     &             'CRST_6    ','CRST_7    ','CRST_8    ',
     &             'CRST_9    ','CRST_10   ','CRST_11   ',
     &             'CRST_12   ','CRST_13   ','CRST_14   ',
     &             'CRST_15   ','CRST_16   ','CRST_17   ',
     &             'CRST_18   ','CRST_19   ','CRST_20   ',
     &             'CRST_21   ','CRST_22   ','CRST_23   ',
     &             'CRST_24   ','CRST_25   ','CRST_26   ',
     &             'CRST_27   ','CRST_28   ','CRST_29   ',
     &             'CRST_30   ','CRST_31   ','CRST_32   ',
     &             'CRST_33   ','CRST_34   ','CRST_35   ',
     &             'CRST_36   ','CRST_37   ','CRST_38   ',
     &             'CRST_39   ','CRST_40   ','CRST_41   ',
     &             'CRST_42   ','CRST_43   ',
     &             'PH2O_1    ',
     &             'PH2O_2    ','PH2O_3    ','PH2O_4    ',
     &             'PH2O_5    ','PH2O_6    ','PH2O_7    ',
     &             'PH2O_8    ','PH2O_9    ','PH2O_10   ',
     &             'PH2O_11   ','PH2O_12   ','PH2O_13   ',
     &             'PH2O_14   ','PH2O_15   ','PH2O_16   ',
     &             'PH2O_17   ','PH2O_18   ','PH2O_19   ',
     &             'PH2O_20   ','PH2O_21   ','PH2O_22   ',
     &             'PH2O_23   ','PH2O_24   ','PH2O_25   ',
     &             'PH2O_26   ','PH2O_27   ','PH2O_28   ',
     &             'PH2O_29   ','PH2O_30   ','PH2O_31   ',
     &             'PH2O_32   ','PH2O_33   ','PH2O_34   ',
     &             'PH2O_35   ','PH2O_36   ','PH2O_37   ',
     &             'PH2O_38   ','PH2O_39   ','PH2O_40   ',
     &             'PH2O_41   ','PH2O_42   ','PH2O_43   ',
     &             'PCL_1     ','PCL_2     ','PCL_3     ',
     &             'PCL_4     ','PCL_5     ','PCL_6     ',
     &             'PCL_7     ','PCL_8     ','PCL_9     ',
     &             'PCL_10    ','PCL_11    ','PCL_12    ',
     &             'PCL_13    ','PCL_14    ','PCL_15    ',
     &             'PCL_16    ','PCL_17    ','PCL_18    ',
     &             'PCL_19    ','PCL_20    ','PCL_21    ',
     &             'PCL_22    ','PCL_23    ','PCL_24    ',
     &             'PCL_25    ','PCL_26    ','PCL_27    ',
     &             'PCL_28    ','PCL_29    ','PCL_30    ',
     &             'PCL_31    ','PCL_32    ','PCL_33    ',
     &             'PCL_34    ','PCL_35    ','PCL_36    ',
     &             'PCL_37    ','PCL_38    ','PCL_39    ',
     &             'PCL_40    ','PCL_41    ','PCL_42    ',
     &             'PCL_43    ',
     &             'NA_1      ','NA_2      ',
     &             'NA_3      ','NA_4      ','NA_5      ',
     &             'NA_6      ','NA_7      ','NA_8      ',
     &             'NA_9      ','NA_10     ','NA_11     ',
     &             'NA_12     ','NA_13     ','NA_14     ',
     &             'NA_15     ','NA_16     ','NA_17     ',
     &             'NA_18     ','NA_19     ','NA_20     ',
     &             'NA_21     ','NA_22     ','NA_23     ',
     &             'NA_24     ','NA_25     ','NA_26     ',
     &             'NA_27     ','NA_28     ','NA_29     ',
     &             'NA_30     ','NA_31     ','NA_32     ',
     &             'NA_33     ','NA_34     ','NA_35     ',
     &             'NA_36     ','NA_37     ','NA_38     ',
     &             'NA_39     ','NA_40     ','NA_41     ',
     &             'NA_42     ','NA_43     ',
     &             'PNH4_1    ',
     &             'PNH4_2    ','PNH4_3    ','PNH4_4    ',
     &             'PNH4_5    ','PNH4_6    ','PNH4_7    ',
     &             'PNH4_8    ','PNH4_9    ','PNH4_10   ',
     &             'PNH4_11   ','PNH4_12   ','PNH4_13   ',
     &             'PNH4_14   ','PNH4_15   ','PNH4_16   ',
     &             'PNH4_17   ','PNH4_18   ','PNH4_19   ',
     &             'PNH4_20   ','PNH4_21   ','PNH4_22   ',
     &             'PNH4_23   ','PNH4_24   ','PNH4_25   ',
     &             'PNH4_26   ','PNH4_27   ','PNH4_28   ',
     &             'PNH4_29   ','PNH4_30   ','PNH4_31   ',
     &             'PNH4_32   ','PNH4_33   ','PNH4_34   ',
     &             'PNH4_35   ','PNH4_36   ','PNH4_37   ',
     &             'PNH4_38   ','PNH4_39   ','PNH4_40   ',
     &             'PNH4_41   ','PNH4_42   ','PNH4_43   ',
     &             'PNO3_1    ','PNO3_2    ','PNO3_3    ',
     &             'PNO3_4    ','PNO3_5    ','PNO3_6    ',
     &             'PNO3_7    ','PNO3_8    ','PNO3_9    ',
     &             'PNO3_10   ','PNO3_11   ','PNO3_12   ',
     &             'PNO3_13   ','PNO3_14   ','PNO3_15   ',
     &             'PNO3_16   ','PNO3_17   ','PNO3_18   ',
     &             'PNO3_19   ','PNO3_20   ','PNO3_21   ',
     &             'PNO3_22   ','PNO3_23   ','PNO3_24   ',
     &             'PNO3_25   ','PNO3_26   ','PNO3_27   ',
     &             'PNO3_28   ','PNO3_29   ','PNO3_30   ',
     &             'PNO3_31   ','PNO3_32   ','PNO3_33   ',
     &             'PNO3_34   ','PNO3_35   ','PNO3_36   ',
     &             'PNO3_37   ','PNO3_38   ','PNO3_39   ',
     &             'PNO3_40   ','PNO3_41   ','PNO3_42   ',
     &             'PNO3_43   ',
     &             'PSO4_1    ','PSO4_2    ',
     &             'PSO4_3    ','PSO4_4    ','PSO4_5    ',
     &             'PSO4_6    ','PSO4_7    ','PSO4_8    ',
     &             'PSO4_9    ','PSO4_10   ','PSO4_11   ',
     &             'PSO4_12   ','PSO4_13   ','PSO4_14   ',
     &             'PSO4_15   ','PSO4_16   ','PSO4_17   ',
     &             'PSO4_18   ','PSO4_19   ','PSO4_20   ',
     &             'PSO4_21   ','PSO4_22   ','PSO4_23   ',
     &             'PSO4_24   ','PSO4_25   ','PSO4_26   ',
     &             'PSO4_27   ','PSO4_28   ','PSO4_29   ',
     &             'PSO4_30   ','PSO4_31   ','PSO4_32   ',
     &             'PSO4_33   ','PSO4_34   ','PSO4_35   ',
     &             'PSO4_36   ','PSO4_37   ','PSO4_38   ',
     &             'PSO4_39   ','PSO4_40   ','PSO4_41   ',
     &             'PSO4_42   ','PSO4_43   ',
     &             'NUM_1     ','NUM_2     ','NUM_3     ',
     &             'NUM_4     ','NUM_5     ','NUM_6     ',
     &             'NUM_7     ','NUM_8     ','NUM_9     ',
     &             'NUM_10    ','NUM_11    ','NUM_12    ',
     &             'NUM_13    ','NUM_14    ','NUM_15    ',
     &             'NUM_16    ','NUM_17    ','NUM_18    ',
     &             'NUM_19    ','NUM_20    ','NUM_21    ',
     &             'NUM_22    ','NUM_23    ','NUM_24    ',
     &             'NUM_25    ','NUM_26    ','NUM_27    ',
     &             'NUM_28    ','NUM_29    ','NUM_30    ',
     &             'NUM_31    ','NUM_32    ','NUM_33    ',
     &             'NUM_34    ','NUM_35    ','NUM_36    ',
     &             'NUM_37    ','NUM_38    ','NUM_39    ',
     &             'NUM_40    ','NUM_41    ','NUM_42    ',
     &             'NUM_43    ',
     &             'PH2O      '/
c
      data crsspc /'CCRS      ','CPRM      '/
c
      data radlist/'O1D       ','O         ','CLO       ',
     &             'CL        ','N2O5      ','NO3       ',
     &             'OH        ','HO2       ','C2O3      ',
     &             'XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       ','RO2R      ',
     &             'R2O2      ','RO2N      ','CCO3      ',
     &             'RCO3      ','MCO3      ','BZCO      ',
     &             'CXO2      ','HCO3      ','TBUO      ',
     &             'BZO       ','BZNO      '/
c
      data nmrad1 /'O1D       ','CL        ','CLO       ','O         ',
     &             'N2O5      ','NO3       ','OH        ','HO2       ',
     &             'C2O3      ','XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       '/
      data nmrad2 /'O1D       ','O         ',
     &             'N2O5      ','NO3       ','OH        ','HO2       ',
     &             'C2O3      ','XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       '/
      data nmrad5 /'O1D       ','O         ','N2O5      ',
     &             'NO3       ','OH        ','HO2       ',
     &             'RO2R      ','R2O2      ','RO2N      ',
     &             'CCO3      ','RCO3      ','MCO3      ','BZCO      ',
     &             'CXO2      ','HCO3      ','TBUO      ',
     &             'BZO       ','BZNO      '/
c
c     Chemical Mechanism Specs
      data mchgas   / 34, 24, 25, 34, 56, 34,  0, 0, 0, 0 /
cjgj      data mchaero  /  0,  0,  0, 16,  0, 13,  0, 0, 0, 0 /
      data mchaero  /  0,  0,  0, 16,  0, 14,  0, 0, 0, 0 /
      data mchrad   / 14, 12, 12, 12, 18, 12,  0, 0, 0, 0 /
      data mchiessr /  4,  2,  2,  2,  2,  2,  0, 0, 0, 0 /
      data mchrxn   /110, 91, 96,100,211,100,  0, 0, 0, 0 /
      data mchphot  / 14, 11, 12, 12, 30, 12,  0, 0, 0, 0 /
      data mchfast  /  4,  4,  4,  4, 13,  4,  0, 0, 0, 0 /
      data mchidmin / 1 /
      data mchidmax / 6 /
      data ipigcb4  /  7,  9, 10, 11, 16, 17, 18, 19, 20 /
      data ipigsap  /  8, 18, 20, 19, 14, 11, 13, 12, 10 /
      data npar     /  1,  2,  4, 10,  5, 12,  8 /
      data tdum     / 298.,  273.,  298. /
      data pdum     / 1013., 1013., 491. /
c
c     Many rules about names, number and order of species are 
c     enforced here unless LCHEM is false.
c
c     Arrays MCHGAS, MCHAERO, MCHRXN and MCHPHOT allow for up to 10
c     mechanisms to be called by RADDRIVR and CHEMDRIV. 
c     Currently:
c      1 = CB4 (mech 3) + Chlorine
c      2 = CB4 + Updated radical-radical reactions
c      3 = CB4 + Carter one product isoprene mechanism
c      4 = CB4 + "1-atmosphere" aerosol mechanism
c      5 = SAPRC99 (56 specs, 221 rxns, 30 phot rxns)
c
c     The state species solver (TRAP) expects the fast species to be
c     ordered first in SPLIST (NO, NO2, O3, PAN etc) so don't change these.
c     Only state species that are named in SPLIST will be allowed.
c     Update CHMSTRY.COM if new species are added to SPLIST
c
c     The radical solvers (RADINIT, RADSLVR) expect radical species
c     in specific orders as defined in NMRAD lists, RADLIST is the
c     union of these lists in order matching CHMSTRY.COM equivalence
c
c
c-----Entry point
c 
c
c-----Mechanism ID
c
      idmech = 0
      read(ichem,'(a)') record
      camxvin = record(21:30)
      call jstlft(camxvin)
      call toupper(camxvin)
      if (camxvin.ne.camxv) then
        write(iout,'(/,a)') ' CAMx version in CHEMPARAM file INVALID'
        write(iout,'(a,a)') ' Expecting: ',camxv
        write(iout,'(a,a)') '     Found: ',camxvin
        goto 910
      endif
      read(ichem,'(a)') record
      read(record(21:80),*) idmech
      write(idiag,'(a,i4)') 'Using CHEMPARAM mechanism id ',idmech
      if (lchem) then
        if (idmech.lt.mchidmin .or. idmech.gt.mchidmax) then
          write(iout,'(/,a)') ' CHEMPARAM mechanism id INVALID'
          write(iout,'(a,i3,a,i3)')' use an id between',mchidmin,
     &                             'and',mchidmax
          goto 910
        endif
c
c --- Mechanism 4 is not allowed with OMP,
c     For compilers with no OMP support you should link in the 
c     routine in the dummy.f file.
c
        ncpus = omp_get_num_procs()
        if( (idmech .EQ. 4 .OR. idmech .EQ. 6) .AND. ncpus .GT. 1 ) then
           write(iout,'(//,a)') 'ERROR in READCHM:'
           write(iout,'(1X,2A)') 'OMP (multiprocessing) is not ',
     &                           'supported for Mechanism 4 or 6.'
           write(iout,'(1X,2A)') 'You must either compile without ',
     &                            'OMP or choose another mechanism.'
           call camxerr()
        endif
      endif
      read(ichem,'(a)') record
      write(idiag,'(a)') record(:istrln(record))
c
c-----Set radical species order. Most mechanisms use same
c     order as 2.  Save final name order.
c
      if (lchem) then
        nrad = mchrad(idmech)
        iessrad = mchiessr(idmech)
        if (nrad.gt.MXRADCL) then
          write(iout,'(/,a,i5,a,i5,a)') ' The number of radicals ',
     &         nrad, ' exceeds the value of MXRADCL ', MXRADCL,
     &         ' set in camx.prm'
          goto 910
        endif
c
        do l = 1,NRADNM
          krad(l) = nrad + 1
        enddo
c
        if (idmech.eq.1) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad1(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        if (idmech.eq.2.or.idmech.eq.3.or.idmech.eq.4
     &                                .or.idmech.eq.6) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad2(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        if (idmech.eq.5) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad5(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        nmrad(nrad+1) = blank
        write(idiag,'(a)') 'The radicals are'
        write(idiag,'(6a10)') (nmrad(i),i=1,nrad)
      endif
c
c-----Species number and reaction number
c
      read(ichem,'(a)') record
      read(record(21:80),*) ngas
      read(ichem,'(a)') record
      read(record(21:80),*) naero
      nspec = ngas + naero 
      if (idmech.EQ.6 .AND. naero.GT.0) then
        read(record(21:80),*) naero,nsec_c,dt_aero
        read(ichem,'(a)') record
        read(record(21:),*) (dsec_i(i),i=1,nsec_c+1)
        nspec = ngas + naero * nsec_c
      endif
      read(ichem,'(a)') record
      read(record(21:80),*) nreact
      if (nspec.lt.1) then 
        write(iout,'(/,a,i5,a)') 
     &    ' number of GAS + AERO species on CHEMPARAM file =', nspec, 
     &    ' is less than 1' 
          goto 910 
      endif
c 
c-----Check number of species and number of reactions against chosen
c     mechanism
c 
      if (lchem) then
        if (ngas.gt.mchgas(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of GAS species on CHEMPARAM file =', ngas,
     &     ' greater than ngas for this mechanism =', mchgas(idmech)
          goto 910
        endif
c
        if (naero.gt.mchaero(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of AERO species on CHEMPARAM file =', naero,
     &     ' greater than naero for this mechanism =', mchaero(idmech)
          goto 910
        endif
c
        if (nreact.ne.mchrxn(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of reactions on CHEMPARAM file =',nreact,
     &     ' not equal to reactions for this mechanism =',mchrxn(idmech)
          goto 910
        endif
      endif
c
c-----Check dimensional limit of species and reaction numbers
c
      if (nspec.gt.MXSPEC) then
        write(iout,'(/,a,i5,a,i5,a)') ' The number of species ',
     &       nspec, ' exceeds the value of MXSPEC ', MXSPEC,
     &       ' set in camx.prm'
        goto 910
      endif
c
      if (nreact.gt.MXRXN) then
        write(iout,'(/,a,i5,a,i5,a)') ' The number of reactions ',
     &       nreact, ' exceeds the value of MXRXN ', MXRXN,
     &       ' set in camx.prm'
        goto 910
      endif
c
c-----Read primary photolysis ID record
c
      read(ichem,'(a)') record
      read(record(21:80),*) nphot1,(idphot1(n),n=1,nphot1) 
      if (nphot1.gt.0) then
        write(idiag,'(a)') 'The primary photolysis reactions are'
        write(idiag,'(i6)') (idphot1(n),n=1,nphot1) 
      endif
c
c-----Read secondary (scaled) photolysis ID records
c
      read(ichem,'(a)') record
      read(record(21:80),*) nphot2
      if (nphot2.gt.0) then
        if (nphot1.eq.0) then
          write(iout,'(/,a)') 
     &     'Need at least one primary photolysis reaction'
          goto 910
        endif
        do n = 1,nphot2
          read(ichem,'(a)') record
          read(record(21:80),*) idphot2(n),idphot3(n),phtscl(n)
        enddo
        write(idiag,'(a)') 'The secondary photolysis reactions are' 
        write(idiag,'(i6,a,i6,a,1pe10.3)')  
     &    (idphot2(n),' =',idphot3(n),' *',phtscl(n),n=1,nphot2)
      endif
      nphot = nphot1 + nphot2
c
c-----Check photolysis ID records
c
      if (lchem) then
c
        if (nphot1.lt.1) then 
          write(iout,'(/,a)') 
     &      ' Need at least one primary photolysis reaction' 
          goto 910 
        endif 
c
        if (nphot1.gt.MXPHT1) then
          write(iout,'(/,a,i5)') 
     &     ' Number of primary photolysis reactions exceeds max of ',
     &       MXPHT1
          goto 910
        endif
c
        if (nphot2.gt.MXPHT2) then
          write(iout,'(/,a,i5)') 
     &     ' Number of secondary photolysis reactions exceeds max of ',
     &       MXPHT2
          goto 910
        endif
c
        if (nphot.ne.mchphot(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' Chemistry mechanism requires', mchphot(idmech),
     &     ' photolysis reactions, but CHEMPARAM file has' , nphot
          goto 910
        endif
c
        nerr = 0
        do n = 1,nphot1
          if (idphot1(n).gt.nreact) nerr = 1
        enddo
c
        do n = 1,nphot2
          if (idphot2(n).gt.nreact) nerr = 1
          if (idphot3(n).gt.nreact) nerr = 1
          nhit = 0
          do nn = 1,nphot1
            if (idphot3(n).eq.idphot1(nn)) nhit = 1
          enddo
          if (nhit.eq.0) nerr = 1
        enddo
c
        if (nerr.eq.1) then
           write(iout,'(/,a)') 'ERROR in the CHEMPARAM file:'
           write(iout,'(a)') 
     &     ' Bad reaction number in one of the photolysis reaction IDs.'
          goto 910
        endif
      endif
c
c-----Species records, gases come first
c
      read(ichem,'(a)') record
      read(ichem,'(a)') record
      write(idiag,'(a)') 'The state species are'
      write(idiag,'(a)') record(:istrln(record))
      do l=1,ngas
        read(ichem,'(5x,a10,2e10.0,4f10.0)')
     &             spname(l),bdnl(l),henry0(l),tfact(l),
     &             diffrat(l),f0(l),rscale(l)
        rscale(l) = amin1(1.,rscale(l))
        rscale(l) = amax1(0.,rscale(l))
        write(idiag,'(i3,2x,a10,2e10.2,4f10.2)')
     &             l,spname(l),bdnl(l),henry0(l),tfact(l),
     &             diffrat(l),f0(l),rscale(l)
      enddo
c
c-----Check over the gas phase species
c
      if (lchem) then
        do l=1,ngas
          if (henry0(l).eq.0.) then
            write(iout,'(/,a,i5)')
     &      'The Henry0 value must be non-zero for species ',l
            goto 910
          endif
        enddo
c
c-----Reorder fast species (NO, NO2, O3, and PAN etc) to come first
c
        nspfst=mchfast(idmech)
        do i=1,nspfst
          do j=1,ngas
            if (splist(i).eq.spname(j)) then
              nametmp = spname(i)
              spname(i) = spname(j)
              spname(j) = nametmp
              tmp = bdnl(i)
              bdnl(i) = bdnl(j)
              bdnl(j) = tmp
              tmp = henry0(i)
              henry0(i) = henry0(j)
              henry0(j) = tmp
              tmp = tfact(i)
              tfact(i) = tfact(j)
              tfact(j) = tmp
              tmp = diffrat(i)
              diffrat(i) = diffrat(j)
              diffrat(j) = tmp
              tmp = f0(i)
              f0(i) = f0(j)
              f0(j) = tmp
              tmp = rscale(i)
              rscale(i) = rscale(j)
              rscale(j) = tmp
            endif
          enddo
        enddo
      endif
c
c-----Read the aero species, if any
c
      if (naero.ne.0) then
        read(ichem,'(a)') record
        write(idiag,'(a)') record(:istrln(record))
        if (idmech.EQ.6) then
          do iaero = 1, naero
            read(ichem,'(5x,a10,e10.0,f10.0)')
     &           tmpnam,bdnl_tmp,roprt_tmp
            nl=istrln(tmpnam)
            do isec = 1, nsec_c
              l = ngas + (iaero-1) * nsec_c + isec
              if ( isec .ge. 10 ) then
                write(tmpnam1,'(a1,i2)') '_',isec
                spname(l)=tmpnam(1:nl)//tmpnam1(1:3)
              else
                write(tmpnam1,'(a1,i1)') '_',isec
                spname(l)=tmpnam(1:nl)//tmpnam1(1:2)
              endif
              bdnl(l) = bdnl_tmp
              roprt(l) = roprt_tmp * 1.e6
              dcut(l,1) = dsec_i(isec)
              dcut(l,2) = dsec_i(isec+1)
              write(idiag,'(i3,2x,a10,e10.2,3f10.2)')
     &           l,spname(l),bdnl(l),roprt_tmp,(dcut(l,m),m=1,2)
            enddo
          enddo
        else
          do l = ngas+1,nspec
            read(ichem,'(5x,a10,e10.0,f10.0)')
     &         spname(l),bdnl(l),roprt(l)
            dcut(l,1) = 0.04
            dcut(l,2) = 2.50
            spcsz = 'FINE      '
            do j = 1, ncrsspc
              if (spname(l).eq.crsspc(j)) then
                dcut(l,1) = 4.30
                dcut(l,2) = 10.0
                spcsz = 'COARSE    '
              endif
            enddo
            write(idiag,'(i3,2x,a10,e10.2,f10.2,2x,a10)')
     &         l,spname(l),bdnl(l),roprt(l),spcsz
            roprt(l) = roprt(l)*1.e6
          enddo
        endif
      endif
c
c-----Map species names to the internal name list.  The default setting
c     to nspec+1 allows species that are in the chem solvers to be
c     omitted from the species list for this run.  Testing if a named
c     pointer is set to nspec+1 is used to identify species that are not
c     in the run.  Named pointers (e.g., kno) are set by equivalence
c     in CHMSTRY.COM
c
      if (lchem) then
        do j = 1,NSPNAM
          kmap(j) = nspec+1
        enddo
c
        do 10 j = 1,nspec
          do i = 1,NSPNAM
            if (splist(i).eq.spname(j)) then
              kmap(i) = j
              goto 10
            endif
          enddo
          write(iout,'(/,3a)') 'species ', spname(j),
     &                        ' is not in the internal list'
          goto 910
  10    continue
        spname(nspec+1) = blank
        write(idiag,'(a)')'Internal species order is:'
        write(idiag,'(6a10)') (spname(i),i=1,nspec) 
c
c-----Check for a consistent set of ammonium/nitrate/sulfate 
c     species
c
        if (idmech.EQ.6) then ! MECH 6
          if (kso2.eq.nspec+1    .or. kh2o2.eq.nspec+1
     &   .or. kform.eq.nspec+1   .or. khono.eq.nspec+1
     &   .or. ko3.eq.nspec+1     .or. koh.eq.nspec+1
     &   .or. kho2.eq.nspec+1    .or. kno3.eq.nspec+1
     &   .or. kno.eq.nspec+1     .or. kno2.eq.nspec+1
     &   .or. kpan.eq.nspec+1    .or. kcg1.eq.nspec+1
     &   .or. kcg2.eq.nspec+1    .or. kcg3.eq.nspec+1
     &   .or. kcg4.eq.nspec+1    .or. khno3.eq.nspec+1
     &   .or. knh3.eq.nspec+1    .or. ksulf.eq.nspec+1
     &   .or. khcl.eq.nspec+1    .or. ksoa1_1.eq.nspec+1
     &   .or. ksoa2_1.eq.nspec+1 .or. ksoa3_1.eq.nspec+1
     &   .or. ksoa4_1.eq.nspec+1 .or. kcrust_1.eq.nspec+1
     &   .or. kpoc_1.eq.nspec+1  .or. kpec_1.eq.nspec+1
     &   .or. kph2o_1.eq.nspec+1 .or. kpcl_1.eq.nspec+1
     &   .or. kna_1.eq.nspec+1   .or. kpnh4_1.eq.nspec+1
     &   .or. kpno3_1.eq.nspec+1 .or. kpso4_1.eq.nspec+1) then
            write(iout,'(/,a)') ' You must have all of the species'
            write(iout,'(a)')   ' SO2,H2O2,FORM,HONO,O3,OH,HO2,NO3,'
            write(iout,'(a)')   ' NO,NO2,PAN,CG1,CG2,CG3,CG4,HNO3,'
            write(iout,'(a)')   ' NH3,SULF,HCL,SOA1,SOA2,SOA3,SOA4,'
            write(iout,'(a)')   ' CRST,POC,PEC,PH2O,PCL,NA,PNH4,PNO3,'
            write(iout,'(a)')   ' PSO4 to use Chemistry Mechanism 6.'
            goto 910
          endif
        else                  ! OTHER THAN MECH 6
          if (knh3.lt.nspec+1 .and. kpnh4.lt.nspec+1 .and.
     &       kso2.lt.nspec+1 .and. ksulf.lt.nspec+1 .and.
     &        kpso4.lt.nspec+1 .and. khno3.lt.nspec+1 .and.
     &         kpno3.lt.nspec+1 ) then
            continue
          elseif (kpnh4.eq.nspec+1 .and. kpso4.eq.nspec+1 .and.
     &         kpno3.eq.nspec+1 ) then
            continue
          else
            write(iout,'(/,a)') ' You must have all of the '
            write(iout,'(a)')   ' ammonium/sulfate/nitrate species '
            write(iout,'(a)')   ' NH3,NH4,SO2,SULF,PSO4,HNO3,PNO3. '
            write(iout,'(a)')   ' Or: '
            write(iout,'(a)')   ' Any combination of the gas-phase'
            write(iout,'(a)')   ' species NH3,SO2,HNO3. '
            goto 910
          endif
c
c-----Check for a consistent set of sea salt species, or none
c
          if (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                               khcl.eq.nspec+1) then 
            continue
          elseif (kna.lt.nspec+1 .and. kpcl.lt.nspec+1 .and.
     &                               khcl.lt.nspec+1) then 
            continue
          elseif (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                               khcl.lt.nspec+1) then 
            continue
          else
            write(iout,'(/,a)') ' You must have all or none of the '
            write(iout,'(a)')   ' sea salt aerosol species '
            write(iout,'(a)')   ' NA, PCL, HCL.'
            write(iout,'(a)')   ' Or, just HCL.'
            goto 910
          endif
c
c-----Check for a consistent set of CG/SOA species, or none
c
          if (kcg1.eq.nspec+1 .and. kcg2.eq.nspec+1 .and.
     &       kcg3.eq.nspec+1 .and. kcg4.eq.nspec+1 .and.
     &        ksoa1.eq.nspec+1 .and. ksoa2.eq.nspec+1 .and.
     &         ksoa3.eq.nspec+1 .and. ksoa4.eq.nspec+1 ) then
            continue
          elseif (kcg1.lt.nspec+1 .and. kcg2.lt.nspec+1 .and.
     &       kcg3.lt.nspec+1 .and. kcg4.lt.nspec+1 .and.
     &        ksoa1.lt.nspec+1 .and. ksoa2.lt.nspec+1 .and.
     &         ksoa3.lt.nspec+1 .and. ksoa4.lt.nspec+1 ) then
            continue
          else
            write(iout,'(/,a)') ' You must have all or none of the '
            write(iout,'(a)')   ' secondary organic aerosol species '
            write(iout,'(a)')   ' CG1,CG2,CG3,CG4,SOA1,SOA2,SOA3,SOA4 '
            goto 910
          endif
        endif                 ! MECH 6 ?
      endif
c
c-----Set up section diameters and check parameters for AERO routines
c
      if ( lchem .AND. idmech.EQ.6 .AND. naero.GT.0 ) then
        ierr = 0 
        call aeroset(nsec_c,dsec_i,ierr)
        if ( ierr .ne. 0 ) goto 910
      endif
c
c-----Reaction records
c
      if (nreact.gt.0) then
        read(ichem,'(a)') record
        read(ichem,'(a)') record
        write(idiag,'(a)') 'The reaction rate parameters are'
        do i = 1,nreact
          read(ichem,'(a)') record
          read(record,*) num,rxntyp(i)
          if (rxntyp(i).lt.1 .or. rxntyp(i).gt.7) goto 902
          read(record,*,err=900) num,rxntyp(i),rxnord(i),
     &                           (rxnpar(i,j),j=1,npar(rxntyp(i)))
c
          if (rxntyp(i).eq.5) then
            iref = anint(rxnpar(i,1))
            if (i.le.iref) goto 901
          endif
c
          write(idiag,'(3i3,1p12e12.4)')
     &       num,rxntyp(i),rxnord(i),(rxnpar(i,j),j=1,npar(rxntyp(i)))
        enddo
c
c-----Populate rate constant lookup table
c
        call exptbl(rxntyp,rxnord,rxnpar)
c
c-----Provide diagnostic info for checking rate expressions
c
        write(idiag,'(/,a,/,/,a)') 
     &        'Diagnostic info for checking rate expressions',
     &        'Rates at three temps and pressures in ppm-n min-1'
        do i=1,3
          call ktherm(tdum(i), pdum(i))
          do j=1,mxrxn
            kdum(j,i)=rk(j)/60.
          enddo
        enddo
        write(idiag,'(a,3F12.1)') 'Temp= ', tdum
        write(idiag,'(a,3F12.1)') 'Pres= ', pdum
        write(idiag,'(a)') 'Rxn Type'
        write(idiag,'(2i3,1p3e12.4)')
     &       (j,rxntyp(j),kdum(j,1),
     &        kdum(j,2),kdum(j,3),j=1,nreact)
c
c-----Set reaction pointers for pig chemistry
c
        if(idmech.le.4 .OR. idmech.eq.6)then
          do i=1,9
            ipigrxn(i)=ipigcb4(i)
          enddo
        elseif(idmech.eq.5)then
          do i=1,9
            ipigrxn(i)=ipigsap(i)
          enddo
        else
          write(iout,'(/,a)') 
     &      'Dont know how to set pig rxns for mech #', idmech
          goto 910
        endif
      endif
c
c---  Set IEH solver species pointers via equivalences in iehchem.com
c
      if (lchem) then
        do i = 1,nradnm
          if (krad(i).gt.nrad) then
            irad(i)=ngas+nrad+1
          elseif (krad(i).gt.iessrad) then
            irad(i)=krad(i)-iessrad
          else
            irad(i)=krad(i)+nspfst+nrad-iessrad
          endif
        enddo
        do i=1, nspnam
          if (kmap(i).gt.nspfst) then
            imap(i)=kmap(i)+nrad
          else
            imap(i)=kmap(i)+nrad-iessrad
          endif
        enddo
      endif
c
c======================== DDM Begin =======================
c
c---  Set species pointers via equivalences in ddmchm.com
c
      if (lchem) then
        do i = 1,nradnm
          if (krad(i).le.nrad) then
            lrad(i)=krad(i)
          else
            lrad(i)=ngas+nrad+1
          endif
        enddo
        do i=1, nspnam
          lmap(i)=kmap(i)+nrad
        enddo
      endif
c
c======================== DDM End   =======================
c
      write(idiag,*)
      call flush(idiag)
c
      return
c
 900  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(4(a,i4))')
     &  'reaction number', num, ' of type', rxntyp(i),
     &  ' should have', npar(rxntyp(i)), ' parameters'
      goto 910
c
 901  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)')  'For reaction type 5'
      write(iout,'(a,i4)') 'reference reaction # ', iref
      write(iout,'(a,i4)') 'must come before this reaction', i
      goto 910
c
 902  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)') 'Reaction type out of bounds (1-5)'
      write(iout,'(a,i4)') 'for reaction ', num
      write(iout,'(a,i4)') 'value was', rxntyp(i)
      write(iout,'(a)') 'Check this is CAMx3 chemparam file'
c
 910  write(*,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(idiag,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(iout,'(//,a)') 'ERROR in READCHM:'
      call camxerr()
c
      end

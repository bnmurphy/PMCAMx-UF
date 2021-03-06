      subroutine rateslo4(neq1,r,rate)
c
c-----CAMx v4.02 030709
c
c     RATESLOW computes reaction rates for slow state species
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        none
c
c     Called by:
c        TRAP
c
      include "camx.prm"
      include "chmstry.com"
c
      real loss(MXSPEC+1),gain(MXSPEC+1),rate(MXSPEC+1),r(MXRXN)
c
      do l=neq1+1,ngas
        Loss(l) = 0.
        Gain(l) = 0.
      enddo
c
c   PAN  HNO3  HONO   PNA  H2O2    CO  FORM  ALD2   PAR   NTR
c   OLE   ETH   TOL  CRES  OPEN  MGLY   XYL  ISOP   SO2  SULF
c  MEOH  ETOH  ISPD  OLE2   CG1   CG2   CG3   CG4
c
        Loss(kPAN  )= +( 1.000)*r( 48)
        Gain(kPAN  )= +( 1.000)*r( 47)
        Loss(kHNO3 )= +( 1.000)*r( 27)
        Gain(kHNO3 )= +( 2.000)*r( 18)+( 1.000)*r( 26)+( 1.000)*r( 41)
     &                +( 1.000)*r( 44)+( 1.000)*r( 67)+( 0.150)*r( 94)
        Loss(kHONO )= +( 1.000)*r( 23)+( 1.000)*r( 24)+( 2.000)*r( 25)
        Gain(kHONO )= +( 2.000)*r( 21)+( 1.000)*r( 22)
        Loss(kPNA  )= +( 1.000)*r( 30)+( 1.000)*r( 31)
        Gain(kPNA  )= +( 1.000)*r( 29)
        Loss(kH2O2 )= +( 1.000)*r( 34)+( 1.000)*r( 35)
        Gain(kH2O2 )= +( 1.000)*r( 32)+( 1.000)*r( 33)
        Loss(kCO   )= +( 1.000)*r( 36)
        Gain(kCO   )= +( 1.000)*r( 37)+( 1.000)*r( 38)+( 1.000)*r( 39)
     &                +( 1.000)*r( 40)+( 1.000)*r( 41)+( 1.000)*r( 45)
     &                +( 0.300)*r( 56)+( 0.330)*r( 58)+( 1.000)*r( 60)
     &                +( 0.420)*r( 62)+( 1.000)*r( 69)+( 2.000)*r( 70)
     &                +( 0.690)*r( 71)+( 1.000)*r( 74)+( 0.066)*r( 77)
     &                +( 0.334)*r( 92)+( 0.225)*r( 93)+( 0.643)*r( 94)
     &                +( 0.333)*r( 95)+( 0.300)*r( 97)+( 0.330)*r( 99)
        Loss(kFORM )= +( 1.000)*r( 37)+( 1.000)*r( 38)+( 1.000)*r( 39)
     &                +( 1.000)*r( 40)+( 1.000)*r( 41)
        Gain(kFORM )= +( 1.000)*r( 45)+( 1.000)*r( 46)+( 2.000)*r( 49)
     &                +( 0.790)*r( 50)+( 1.000)*r( 51)+( 0.200)*r( 56)
     &                +( 1.000)*r( 57)+( 0.740)*r( 58)+( 1.000)*r( 59)
     &                +( 1.000)*r( 60)+( 1.560)*r( 61)+( 1.000)*r( 62)
     &                +( 1.000)*r( 70)+( 0.700)*r( 71)+( 0.500)*r( 75)
     &                +( 0.629)*r( 76)+( 0.600)*r( 77)+( 1.000)*r( 84)
     &                +( 0.167)*r( 92)+( 0.150)*r( 93)+( 0.282)*r( 94)
     &                +( 0.900)*r( 95)+( 0.200)*r( 97)+( 1.000)*r( 98)
     &                +( 0.740)*r( 99)+( 1.000)*r(100)
        Loss(kALD2 )= +( 1.000)*r( 42)+( 1.000)*r( 43)+( 1.000)*r( 44)
     &                +( 1.000)*r( 45)
        Gain(kALD2 )= +( 0.110)*r( 52)+( 1.100)*r( 53)+( 0.630)*r( 56)
     &                +( 1.000)*r( 57)+( 0.500)*r( 58)+( 1.000)*r( 59)
     &                +( 0.220)*r( 61)+( 0.030)*r( 71)+( 0.150)*r( 77)
     &                +( 0.800)*r( 78)+( 1.000)*r( 85)+( 0.273)*r( 92)
     &                +( 0.020)*r( 93)+( 0.357)*r( 94)+( 0.067)*r( 95)
     &                +( 0.800)*r( 96)+( 0.630)*r( 97)+( 1.000)*r( 98)
     &                +( 0.500)*r( 99)+( 1.000)*r(100)
        Loss(kPAR  )= +( 1.000)*r( 52)
        Gain(kPAR  )= +(-0.110)*r( 52)+(-2.100)*r( 53)+( 0.220)*r( 56)
     &                +(-1.000)*r( 57)+(-1.000)*r( 58)+(-1.000)*r( 59)
     &                +( 1.100)*r( 72)+( 0.250)*r( 75)+( 0.350)*r( 77)
     &                +( 2.400)*r( 78)+( 1.565)*r( 92)+( 0.360)*r( 93)
     &                +( 1.282)*r( 94)+( 0.832)*r( 95)+( 2.400)*r( 96)
     &                +( 0.220)*r( 97)+(-1.000)*r( 98)+(-1.000)*r( 99)
     &                +(-1.000)*r(100)

        Gain(kNTR  )= +( 1.000)*r( 55)+( 0.100)*r( 64)+( 1.000)*r( 68)
     &                +( 0.800)*r( 78)+( 1.000)*r( 81)+( 0.850)*r( 94)
     &                +( 0.800)*r( 96)
        Loss(kOLE  )= +( 1.000)*r( 56)+( 1.000)*r( 57)+( 1.000)*r( 58)
     &                +( 1.000)*r( 59)

        Loss(kETH  )= +( 1.000)*r( 60)+( 1.000)*r( 61)+( 1.000)*r( 62)

        Loss(kTOL  )= +( 1.000)*r( 63)

        Loss(kCRES )= +( 1.000)*r( 66)+( 1.000)*r( 67)
        Gain(kCRES )= +( 0.360)*r( 63)+( 1.000)*r( 65)+( 0.200)*r( 72)
        Loss(kOPEN )= +( 1.000)*r( 69)+( 1.000)*r( 70)+( 1.000)*r( 71)
        Gain(kOPEN )= +( 0.900)*r( 64)+( 0.300)*r( 66)
        Loss(kMGLY )= +( 1.000)*r( 73)+( 1.000)*r( 74)
        Gain(kMGLY )= +( 0.200)*r( 71)+( 0.800)*r( 72)+( 0.168)*r( 92)
     &                +( 0.850)*r( 93)
        Loss(kXYL  )= +( 1.000)*r( 72)

        Loss(kISOP )= +( 1.000)*r( 75)+( 1.000)*r( 76)+( 1.000)*r( 77)
     &                +( 1.000)*r( 78)+( 1.000)*r( 96)

        Loss(kSO2  )= +( 1.000)*r( 82)+( 1.000)*r( 83)


        Gain(kSULF )= +( 1.000)*r( 82)+( 1.000)*r( 83)
        Loss(kMEOH )= +( 1.000)*r( 84)

        Loss(kETOH )= +( 1.000)*r( 85)

        Loss(kISPD )= +( 1.000)*r( 92)+( 1.000)*r( 93)+( 1.000)*r( 94)
     &                +( 1.000)*r( 95)
        Gain(kISPD )= +( 0.750)*r( 75)+( 0.912)*r( 76)+( 0.650)*r( 77)
     &                +( 0.200)*r( 78)+( 0.200)*r( 96)
        Loss(kOLE2 )= +( 1.000)*r( 97)+( 1.000)*r( 98)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)


        Gain(kCG1  )= +( 0.070)*r( 63)+( 0.044)*r( 72)

        Gain(kCG2  )= +( 0.137)*r( 63)+( 0.192)*r( 72)

        Gain(kCG3  )= +( .0024)*r( 52)+( .0024)*r( 56)+( .0024)*r( 57)
     &                +( .0024)*r( 58)+( .0024)*r( 59)+( 0.036)*r( 66)
     &                +( 0.036)*r( 67)

        Gain(kCG4  )= +( 0.136)*r( 97)+( 0.136)*r( 98)+( 0.136)*r( 99)
     &                +( 0.136)*r(100)
      do l=neq1+1,ngas
        rate(l) = gain(l) -loss(l)
      enddo
c
c
      return
      end

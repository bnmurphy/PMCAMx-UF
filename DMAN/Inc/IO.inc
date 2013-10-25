c-----------------------------------------------------------------------
c                 Inclusion file for IO parameters
c-----------------------------------------------------------------------
      integer idat, iout, itrk
c                 input, output and tracking file unit number
      character*8 dir(3)
c
      parameter (idat = 10, iout = 40, itrk = 80)
c
      data dir/'./Input/', './Outpt/', './ErrTk/'/

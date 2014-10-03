      Block Data soapdat
c
c-----------------------------------------------------------------------
c     Initialize common block variables for SOAP
c
c     mwsoap   - molecular weights of CG/SOA species (g/mol)
c     csat     - saturation concentrations of SOA species (ug/m3)
c     cstemp   - temperatures corresponding to saturation concentrations
c                of CG/SOA species (K)
c     deltah   - enthalpy of vaporization of CG/SOA species (kJ/mol)
c     flagsoap - 1 if SOA species forms solutions; 0 if not
c     lae3     - flag to set emulation of the CMAQ AE3 algorithm
c                which allows no evaporation of SOA (recommend false)
c-----------------------------------------------------------------------
c
      include 'soap.com'
c
      data mwsoap   /150., 150., 150., 180./
      data csat     /0.023, 0.674, 0.007, 0.008/
      data cstemp   /4*281.5/
      data deltah   /156250., 156250., 0., 0./
      data flagsoap /4*1/
      data lae3     /.false./
c
      end

#
MAIN  =  ./
INC  =  ./Inc
DAT = ./Input
OUT = ./Outpt
TRK = ./ErrTk
PSSA = ./cond_nuc_PSSA

FC   =  pgf77
#FLGS =  -I $(INC) -pc 64 -O1 # flags for the compiling of *.o files
FLGS =  -I $(INC) -pc 64 -g # flags for the compiling of *.o files

nucbox:

TARGT = $(MAIN)nucbox

#derivnucl.o \#
#experiment.o \#
#gas.o \#
#initconv.o \#
#ion_nucl.o \#
#jacobn.o \#
#lognorm.o \#
#lubksb.o \#
#ludcmp.o \#
#Mtot_sub.o \#
#neutral.o \#
#Nk2pn.o \#
#nucl_sub.o \#
#odeint.o \#
#pn2Nk.o \#
#printf.o \#
#readinput.o \#
#rkqs.o \#
#startup.o \#
#stiff.o \#

#derivnucl.o               : derivnucl.f                                        \
                          $(INC)/aervaria.inc $(INC)/nucleation.inc
#experiment.o 		  : experiment.f                                       \
                          $(INC)/aervaria.inc
#gas.o  			  : gas.f                                              \
                          $(INC)/aervaria.inc 
#initconv.o		  : initconv.f                                         \
                          $(INC)/aervaria.inc $(INC)/IO.inc
#ion_nucl.o 		  : ion_nucl.f
#jacobn.o                  : jacobn.f                                           \
                          $(INC)/aervaria.inc $(INC)/nucleation.inc
#lognorm.o                 : lognorm.f                                          \
                          $(INC)/aervaria.inc
#lubksb.o                  : lubksb.f                          
#ludcmp.o                  : ludcmp.f
#Mtot_sub.o                : Mtot_sub.f                                         \
                          $(INC)/aervaria.inc

#Nk2pn.o 		  : Nk2pn.f                                            \
	                  $(INC)/aervaria.inc
#pn2Nk.o 		  : pn2Nk.f                                            \
	                  $(INC)/aervaria.inc
#printf.o   		  : printf.f                                           \
	                  $(INC)/aervaria.inc
#readinput.o               : readinput.f                                        \
			  $(INC)/aervaria.inc $(INC)/IO.inc
#rkqs.o 			  : rkqs.f                                             \
                          $(INC)/aervaria.inc                                 
#neutral.o                 : neutral.f                                          \
                          $(INC)/aervaria.inc
#nucl_sub.o 		  : nucl_sub.f                                         \
	                  $(INC)/aervaria.inc
#odeint.o                  : odeint.f                                           \
                          $(INC)/aervaria.inc
#startup.o  		  : startup.f                                          \
	                  $(INC)/aervaria.inc
#stiff.o                   : stiff.f                                            \
                          $(INC)/aervaria.inc
OBJCTS = \
alpha_nh3.o \
dman.o \
dmdt_int.o \
gasdiff.o \
initbounds.o \
multicoag.o \
so4cond.o \
tmconds.o \
test.o \
$(PSSA)/aerodens_PSSA.o \
$(PSSA)/cond_nuc.o \
$(PSSA)/getCondSink.o \
$(PSSA)/getH2SO4conc.o \
$(PSSA)/nucleation.o \
$(PSSA)/ezcond.o \
$(PSSA)/eznh3eqm.o \
$(PSSA)/ezwatereqm.o \
$(PSSA)/getNucRate.o \
$(PSSA)/napa_nucl.o \
$(PSSA)/vehk_nucl.o \
$(PSSA)/mnfix_PSSA.o \
$(PSSA)/waterso4.o \
$(PSSA)/waternacl.o

$(TARGT):  $(OBJCTS)
	   $(FC) -o $(TARGT) $(FLGS) $(OBJCTS)
.f.o    :
	$(FC) -c -o $@ $(FLGS) $<
alpha_nh3.o               : alpha_nh3.f
dmdt_int.o 		  : dmdt_int.f
dman.o 		  	  : dman.f                                             \
                          $(INC)/aervaria.inc $(INC)/IO.inc
gasdiff.o  		  : gasdiff.f
initbounds.o	          : initbounds.f                                       \
			  $(INC)/aervaria.inc
multicoag.o  		  : multicoag.f                                        \
	                  $(INC)/aervaria.inc
so4cond.o  		  : so4cond.f                                          \
                          $(INC)/aervaria.inc
tmconds.o                 : tmconds.f                                          \
                          $(INC)/aervaria.inc
test.o                    : test.f                                             \
                          $(INC)/aervaria.inc
$(PSSA)/aerodens_PSSA.o   : $(PSSA)/aerodens_PSSA.f                            \
			  $(PSSA)/aerodens_PSSA.f
$(PSSA)/cond_nuc.o 	  : $(PSSA)/cond_nuc.f                                 \
			  $(INC)/sizecode.COM
$(PSSA)/ezcond.o  	  : $(PSSA)/ezcond.f                                   \
	                  $(INC)/sizecode.COM
$(PSSA)/eznh3eqm.o	  : $(PSSA)/eznh3eqm.f                                 \
			  $(INC)/sizecode.COM
$(PSSA)/ezwatereqm.o      : $(PSSA)/ezwatereqm.f                               \
	                  $(INC)/sizecode.COM
$(PSSA)/getCondSink.o     : $(PSSA)/getCondSink.f                              \
			  $(INC)/sizecode.COM
$(PSSA)/getH2SO4conc.o    : $(PSSA)/getH2SO4conc.f                             \
			  $(INC)/sizecode.COM
$(PSSA)/getNucRate.o	  : $(PSSA)/getNucRate.f                               \
			  $(INC)/sizecode.COM
$(PSSA)/mnfix_PSSA.o	  : $(PSSA)/mnfix_PSSA.f                               \
			  $(INC)/sizecode.COM
$(PSSA)/napa_nucl.o	  : $(PSSA)/napa_nucl.f                          
$(PSSA)/nucleation.o	  : $(PSSA)/nucleation.f                               \
			  $(INC)/sizecode.COM
$(PSSA)/vehk_nucl.o	  : $(PSSA)/vehk_nucl.f                      
$(PSSA)/waterso4.o	  : $(PSSA)/waterso4.f	
$(PSSA)/waternacl.o	  : $(PSSA)/waternacl.f 
clean:
	rm $(OBJCTS) $(OUT)/*.out $(OUT)/*.dat $(TRK)/*.trk fort.* core.* nucbox

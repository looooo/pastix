###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Bilel Hadri
# @author Hatem Ltaief
# @date 2010-11-15
#
###

PLASMA_DIR = ../..
include ../../Makefile.internal

ALINTST = \
   aladhd.o alaerh.o alaesm.o alahd.o alareq.o \
   alasum.o alasvm.o chkxer.o lsamen.o xlaenv.o xerbla.o \

SCLNTST= slaord.o

DZLNTST= dlaord.o

SLINTST = schkaa.o \
   schkge.o serrge.o sdrvge.o \
   serrvx.o \
   sget02.o sget04.o sget06.o \
   schkpo.o serrpo.o sdrvpo.o \
   sposvx.o spotri.o sporfs.o \
   spot01.o spot02.o spot03.o spot05.o \
   sdrvls.o serrls.o \
   schkqr.o serrqr.o \
   sqrt01.o sqrt02.o sqrt03.o sqrt13.o sqrt14.o sqrt15.o sqrt16.o sqrt17.o \
   schklq.o serrlq.o \
   slqt01.o slqt02.o slqt03.o \
   slarhs.o slatb4.o sgennd.o \
   sgeequ.o slaqge.o \
   spoequ.o slaqsy.o spocon.o slacn2.o slatrs.o slabad.o srscl.o \
   slascl.o slaset.o slatms.o slartg.o slatm1.o slarnd.o slarot.o \
   slaror.o slagsy.o slagge.o slaran.o \
   slauum.o strtri.o \
   slauu2.o strti2.o \

DLINTST = dchkaa.o \
   dchkge.o derrge.o ddrvge.o \
   derrvx.o \
   dget02.o dget04.o dget06.o \
   dchkpo.o derrpo.o ddrvpo.o \
   dposvx.o dpotri.o dporfs.o \
   dpot01.o dpot02.o dpot03.o dpot05.o \
   ddrvls.o derrls.o \
   dchkqr.o derrqr.o \
   dqrt01.o dqrt02.o dqrt03.o dqrt13.o dqrt14.o dqrt15.o dqrt16.o dqrt17.o \
   dchklq.o derrlq.o \
   dlqt01.o dlqt02.o dlqt03.o \
   dlarhs.o dlatb4.o dgennd.o \
   dgeequ.o dlaqge.o \
   dpoequ.o dlaqsy.o dpocon.o dlacn2.o dlatrs.o dlabad.o drscl.o \
   dlascl.o dlaset.o dlatms.o dlartg.o dlatm1.o dlarnd.o dlarot.o \
   dlaror.o dlagsy.o dlagge.o dlaran.o \
   dlauum.o dtrtri.o \
   dlauu2.o dtrti2.o \

CLINTST = cchkaa.o \
   cchkge.o cerrge.o cdrvge.o \
   cerrvx.o \
   cget02.o cget04.o \
   cchkpo.o cerrpo.o cdrvpo.o \
   cposvx.o cpotri.o cporfs.o \
   cpot01.o cpot02.o cpot03.o cpot05.o \
   cdrvls.o cerrls.o \
   cchkqr.o cerrqr.o \
   cqrt01.o cqrt02.o cqrt03.o cqrt13.o cqrt14.o cqrt15.o cqrt16.o cqrt17.o \
   cchklq.o cerrlq.o \
   clqt01.o clqt02.o clqt03.o \
   clarhs.o clatb4.o cgennd.o \
   cgeequ.o claqge.o \
   cpoequ.o claqsy.o cpocon.o clacn2.o clatrs.o csrscl.o \
   clascl.o claset.o clatms.o clartg.o clatm1.o clarnd.o clarot.o clanhe.o claqhe.o \
   claror.o clagsy.o clagge.o  \
   claipd.o csbmv.o \
   clauum.o ctrtri.o \
   clauu2.o ctrti2.o \
   cspmv.o csymv.o \
   sget06.o slabad.o slaran.o slarnd.o slascl.o slatm1.o scsum1.o icmax1.o \

ZLINTST = zchkaa.o \
   zchkge.o zerrge.o zdrvge.o \
   zerrvx.o \
   zget02.o zget04.o \
   zchkpo.o zerrpo.o zdrvpo.o \
   zposvx.o zpotri.o zporfs.o \
   zpot01.o zpot02.o zpot03.o zpot05.o \
   zdrvls.o zerrls.o \
   zchkqr.o zerrqr.o \
   zqrt01.o zqrt02.o zqrt03.o zqrt13.o zqrt14.o zqrt15.o zqrt16.o zqrt17.o \
   zchklq.o zerrlq.o \
   zlqt01.o zlqt02.o zlqt03.o \
   zlarhs.o zlatb4.o zgennd.o \
   zgeequ.o zlaqge.o \
   zpoequ.o zlaqsy.o zpocon.o zlacn2.o zlatrs.o zdrscl.o \
   zlascl.o zlaset.o zlatms.o zlartg.o zlatm1.o zlarnd.o zlarot.o zlanhe.o zlaqhe.o \
   zlaror.o zlagsy.o zlagge.o  \
   zlaipd.o zsbmv.o \
   zlauum.o ztrtri.o \
   zlauu2.o ztrti2.o \
   zspmv.o zsymv.o \
   dget06.o dlabad.o dlaran.o dlarnd.o dlascl.o dlatm1.o dzsum1.o izmax1.o \

all:  	$(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 	\
	single double complex complex16

single: xlintsts
double: xlintstd
complex: xlintstc
complex16: xlintstz

xlintsts xlintstd xlintstc xlintstz : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

xlintsts : $(ALINTST) $(SLINTST) $(SCLNTST)
	       $(LOADER) $(LDFLAGS) $(ALINTST) $(SLINTST) $(SCLNTST) $(LIB) -o $@

xlintstc : $(ALINTST) $(CLINTST) $(SCLNTST)
	       $(LOADER) $(LDFLAGS) $(ALINTST) $(CLINTST) $(SCLNTST) $(LIB) -o $@

xlintstd : $(ALINTST) $(DLINTST) $(DZLNTST)
	       $(LOADER) $(LDFLAGS) $(ALINTST) $(DLINTST) $(DZLNTST) $(LIB) -o $@

xlintstz : $(ALINTST) $(ZLINTST) $(DZLNTST)
	       $(LOADER) $(LDFLAGS) $(ALINTST) $(ZLINTST) $(DZLNTST) $(LIB) -o $@

$(ALINTST): $(FRC)
$(SCLNTST): $(FRC)
$(DZLNTST): $(FRC)
$(SLINTST): $(FRC)
$(DLINTST): $(FRC)
$(CLINTST): $(FRC)
$(ZLINTST): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o

cleanall: clean
	rm -f xlintsts xlintstd xlintstc xlintstz

schkaa.o: schkaa.f
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
dchkaa.o: dchkaa.f
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
cchkaa.o: cchkaa.f
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
zchkaa.o: zchkaa.f
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
.o:  
	$(LOADER) $(LDFLAGS) $(INC) $< -o $@ $(LIB) 

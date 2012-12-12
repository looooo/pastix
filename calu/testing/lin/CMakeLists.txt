###
#
# @file CMakeLists.txt
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

set (ALINTST aladhd.f alaerh.f alaesm.f alahd.f alareq.f alasum.f
   alasvm.f chkxer.f lsamen.f xlaenv.f xerbla.f )

set (SCLNTST slaord.f )

set (DZLNTST dlaord.f )

set (SLINTST schkaa.f schkge.f serrge.f sdrvge.f serrvx.f sget02.f
   sget04.f sget06.f schkpo.f serrpo.f sdrvpo.f sposvx.f spotri.f
   sporfs.f spot01.f spot02.f spot03.f spot05.f sdrvls.f serrls.f
   schkqr.f serrqr.f sqrt01.f sqrt02.f sqrt03.f sqrt13.f sqrt14.f
   sqrt15.f sqrt16.f sqrt17.f schklq.f serrlq.f slqt01.f slqt02.f
   slqt03.f slarhs.f slatb4.f sgennd.f sgeequ.f slaqge.f spoequ.f
   slaqsy.f spocon.f slacn2.f slatrs.f slabad.f srscl.f slascl.f
   slaset.f slatms.f slartg.f slatm1.f slarnd.f slarot.f slaror.f
   slagsy.f slagge.f slaran.f slauum.f strtri.f slauu2.f strti2.f )

set( DLINTST dchkaa.f dchkge.f derrge.f ddrvge.f derrvx.f dget02.f
   dget04.f dget06.f dchkpo.f derrpo.f ddrvpo.f dposvx.f dpotri.f
   dporfs.f dpot01.f dpot02.f dpot03.f dpot05.f ddrvls.f derrls.f
   dchkqr.f derrqr.f dqrt01.f dqrt02.f dqrt03.f dqrt13.f dqrt14.f
   dqrt15.f dqrt16.f dqrt17.f dchklq.f derrlq.f dlqt01.f dlqt02.f
   dlqt03.f dlarhs.f dlatb4.f dgennd.f dgeequ.f dlaqge.f dpoequ.f
   dlaqsy.f dpocon.f dlacn2.f dlatrs.f dlabad.f drscl.f dlascl.f
   dlaset.f dlatms.f dlartg.f dlatm1.f dlarnd.f dlarot.f dlaror.f
   dlagsy.f dlagge.f dlaran.f dlauum.f dtrtri.f dlauu2.f dtrti2.f )

set (CLINTST cchkaa.f cchkge.f cerrge.f cdrvge.f cerrvx.f cget02.f
   cget04.f cchkpo.f cerrpo.f cdrvpo.f cposvx.f cpotri.f cporfs.f
   cpot01.f cpot02.f cpot03.f cpot05.f cdrvls.f cerrls.f cchkqr.f
   cerrqr.f cqrt01.f cqrt02.f cqrt03.f cqrt13.f cqrt14.f cqrt15.f
   cqrt16.f cqrt17.f cchklq.f cerrlq.f clqt01.f clqt02.f clqt03.f
   clarhs.f clatb4.f cgennd.f cgeequ.f claqge.f cpoequ.f claqsy.f
   cpocon.f clacn2.f clatrs.f csrscl.f clascl.f claset.f clatms.f
   clartg.f clatm1.f clarnd.f clarot.f clanhe.f claqhe.f claror.f
   clagsy.f clagge.f claipd.f csbmv.f clauum.f ctrtri.f clauu2.f
   ctrti2.f cspmv.f csymv.f sget06.f slabad.f slaran.f slarnd.f
   slascl.f slatm1.f scsum1.f icmax1.f )

set (ZLINTST zchkaa.f zchkge.f zerrge.f zdrvge.f zerrvx.f zget02.f
   zget04.f zchkpo.f zerrpo.f zdrvpo.f zposvx.f zpotri.f zporfs.f
   zpot01.f zpot02.f zpot03.f zpot05.f zdrvls.f zerrls.f zchkqr.f
   zerrqr.f zqrt01.f zqrt02.f zqrt03.f zqrt13.f zqrt14.f zqrt15.f
   zqrt16.f zqrt17.f zchklq.f zerrlq.f zlqt01.f zlqt02.f zlqt03.f
   zlarhs.f zlatb4.f zgennd.f zgeequ.f zlaqge.f zpoequ.f zlaqsy.f
   zpocon.f zlacn2.f zlatrs.f zdrscl.f zlascl.f zlaset.f zlatms.f
   zlartg.f zlatm1.f zlarnd.f zlarot.f zlanhe.f zlaqhe.f zlaror.f
   zlagsy.f zlagge.f zlaipd.f zsbmv.f zlauum.f ztrtri.f zlauu2.f
   ztrti2.f zspmv.f zsymv.f dget06.f dlabad.f dlaran.f dlarnd.f
   dlascl.f dlatm1.f dzsum1.f izmax1.f )


add_executable( xlintsts ${ALINTST} ${SLINTST} ${SCLNTST} )
PLASMA_SETUP_FORTRAN_MAIN( xlintsts )

add_executable( xlintstc ${ALINTST} ${CLINTST} ${SCLNTST} )
PLASMA_SETUP_FORTRAN_MAIN( xlintstc )

add_executable( xlintstd ${ALINTST} ${DLINTST} ${DZLNTST} )
PLASMA_SETUP_FORTRAN_MAIN( xlintstd )

add_executable( xlintstz ${ALINTST} ${ZLINTST} ${DZLNTST} )
PLASMA_SETUP_FORTRAN_MAIN( xlintstz )

INSTALL( TARGETS xlintsts xlintstc xlintstd xlintstz DESTINATION testing/lin )

INSTALL( FILES lapack_testing.py stest.in ctest.in dtest.in ztest.in DESTINATION testing/lin )
INSTALL( FILES stestdyn.in ctestdyn.in dtestdyn.in ztestdyn.in DESTINATION testing/lin )

INSTALL( FILES ${ALINTST} ${SLINTST} ${SCLNTST} ${CLINTST} 
 ${DLINTST} ${DZLNTST} ${ZLINTST} CMakeLists.txt 
 DESTINATION testing/lin )

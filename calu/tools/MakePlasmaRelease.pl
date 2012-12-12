#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

my $major; 
my $minor; 
my $micro; 

#Options par défaut
my $DIRNAME;
my $BASENAME;
my $svn    = "https://icl.cs.utk.edu/svn/plasma/trunk";
my $svninst= "https://icl.cs.utk.edu/svn/plasma/plasma-installer";
my $user   = "";

my @file2delete = (
    "tools/BuildbotMakegen.py",
    "tools/code_cleanup",
    "tools/convert2eztrace.pl",
    "tools/genf77interface.pl",
    "tools/genf90interface.pl",
    "tools/MakePlasmaRelease.pl",
    "tools/pl_tstsrc.py",
    "tools/UpdateVersionNumber.sh",
    "tools/subs.pyc",
    "Makefile.gen",
    "quark/examples",
    "quark/docs",
    "include/Makefile",
    "include/Makefile.src",
    "GenerateZDCS.cmake.optional",
    "CMakeBuildNotes",
    "CMakeLists.txt",
    "control/CMakeLists.txt",
    "include/CMakeLists.txt",
    "timing/CMakeLists.txt",
    "docs/asciidoc/CMakeBuildNotes.txt",
    "core_blas/CMakeLists.txt",
    "compute/CMakeLists.txt",
    "testing/lin/CMakeLists.txt",
    "testing/CMakeLists.txt",
    "examples/CMakeLists.txt",
    "makes/cmake32.bat",
    "makes/cmake64.bat",
    "timing/lin",
    "compute/pzhbrdt_v2.c",
    "Makefile.src"
    );

my $RELEASE_PATH;
my %opts;
my $NUMREL = "";

my @vars = ( 'SSRC', 'DSRC', 'CSRC', 'SHDR', 'DHDR', 'CHDR' );



sub myCmd {
    my ($cmd) = @_ ;
    my $err = 0;

    print $cmd;
    $err = system($cmd);
    if ($err != 0) {
    	print "Error during execution of the following command:\n$cmd\n";
    	exit;
    }    
}
    
sub CleanMakefile {
    my ($file) = @_ ;
    my $str;
    my $dir;

    # Replace the [SDC]SRC and [SDC]HDR variables by their values
    print $file;
    $dir = `dirname $file`; chop $dir;

    foreach my $v (@vars)
    {
        # The first command is to be sure that Makefile.gen is
        # uptodate because we are touching at each step Makefile, so
        # it will be generated during next cvall otherwise
        myCmd("cd $dir && make .Makefile.gen");
        $str = `cd $dir && make print$v`; chop $str;
        myCmd("sed -i 's/${v}[ ]*=.*/$v = $str/' $file");
    }
}

sub CleanF90 {
    my ($file) = @_ ;
    my $str;
    my $dir;
    my @vars = ( 'SF90', 'DF90', 'CF90' );

    # Replace the [SDC]SRC and [SDC]HDR variables by their values
    print $file;
    $dir = `dirname $file`; chop $dir;

    foreach my $v (@vars)
    {
        # The first command is to be sure that Makefile.gen is
        # uptodate because we are touching at each step Makefile, so
        # it will be generated during next cvall otherwise
        myCmd("cd $dir && make .Makefile.gen");
        $str = `cd $dir && make print$v`; chop $str;
        myCmd("sed -i 's/${v}[ ]*=.*/$v = $str/' $file");
    }
}

sub MakeRelease {

    my $numversion = $major.'.'.$minor.'.'.$micro;
    my $cmd;

    $RELEASE_PATH = $ENV{ PWD}."/plasma_".$numversion;

    # Sauvegarde du rep courant
    my $dir = `pwd`; chop $dir;

    $cmd = 'svn export --force '.$NUMREL.' '.$user.' '.$svn.' '.$RELEASE_PATH;
    myCmd($cmd);
    
    chdir $RELEASE_PATH;

    # Change version in plasma.h
    #myCmd("sed -i 's/PLASMA_VERSION_MAJOR[ ]*[0-9]/PLASMA_VERSION_MAJOR $major/' include/plasma.h CMakeLists.txt"); 
    #myCmd("sed -i 's/PLASMA_VERSION_MINOR[ ]*[0-9]/PLASMA_VERSION_MINOR $minor/' include/plasma.h CMakeLists.txt");
    #myCmd("sed -i 's/PLASMA_VERSION_MICRO[ ]*[0-9]/PLASMA_VERSION_MICRO $micro/' include/plasma.h CMakeLists.txt");

    # Change the version in comments 
    #myCmd("find -name '*.[ch]' -type f -exec sed -i 's/\@version[ ]*[.0-9]*/\@version $numversion/' {} \\;");
    #
    #find -name "*.[ch]" | xargs sed -i 's/@version 2.3.1/@version 2.4.0/'
    #find -name "CMakeLists.txt" | xargs sed -i 's/@version 2.3.1/@version 2.4.0/'
    #find -name "Makefile" | xargs sed -i 's/@version 2.3.1/@version 2.4.0/'
    #sed -i 's/@version 2.3.1/@version 2.4.0/' Makefile.gen Makefile.internal Makefile.tau 
    #sed -i 's/@version 2.3.1/@version 2.4.0/' makes/make.inc.*
    #sed -i 's/@version 2.3.1/@version 2.4.0/' make.inc.example                            
    #sed -i 's/@version 2.3.1/@version 2.4.0/' docs/latex/contributors_guide/comments.tex
    #sed -i 's/@version 2.3.1/@version 2.4.0/' makes/cmake32.bat makes/cmake64.bat       
    #sed -i 's/@version 2.3.1/@version 2.4.0/' tools/convert2eztrace.pl           
    #sed -i 's/version 2.3.1/version 2.4.0/' examples/*.f
    #

    #Precision Generation
    print "Generate the different precision\n"; 
    myCmd("echo \"MAKE = make -j 4\" > make.inc");
    myCmd("echo \"PLASMA_F90 = 1\" >> make.inc");
    myCmd("make generation");

    # Clean the Makefile form the precision generation process
    foreach my $file (`find $RELEASE_PATH -name Makefile | xargs grep -l "[SCD]SRC"`){
        CleanMakefile( $file );
    }
    CleanF90( "$RELEASE_PATH/control/Makefile" );

    my $file = "include/Makefile";
    print "Remove $file\n";
    myCmd("rm -rf $RELEASE_PATH/$file");

    #Compile the documentation
    print "Compile the documentation\n"; 
    system("make -C ./docs");
    myCmd("rm -f make.inc");
    
    #Remove non required files (Makefile.gen)
    foreach my $file (@file2delete){
	print "Remove $file\n";
 	myCmd("rm -rf $RELEASE_PATH/$file");
    }
 
    # Remove 'include Makefile.gen from Makefile'
    myCmd("find $RELEASE_PATH -name Makefile -exec sed -i '/Makefile.gen/ d' {} \\;");

    # Remove the lines relative to include directory in root Makefile
    myCmd("sed -i '/cd include/ d' $RELEASE_PATH/Makefile");

    # Remove '.Makefile.gen files'
    myCmd("find $RELEASE_PATH -name .Makefile.gen -exec rm -f {} \\;");

    chdir $dir;

    # Save the InstallationGuide if we want to do a plasma-installer release
    myCmd("cp $RELEASE_PATH/InstallationGuide InstallationGuide-v${numversion}");
    myCmd("cp $RELEASE_PATH/ReleaseNotes      ReleaseNotes-v${numversion}");

    #Create tarball
    print "Create the tarball\n";
    $DIRNAME=`dirname $RELEASE_PATH`;
    $BASENAME=`basename $RELEASE_PATH`;
    chop $DIRNAME;
    chop $BASENAME;
    myCmd("(cd $DIRNAME && tar czf ${BASENAME}.tar.gz $BASENAME)");

}

sub MakeInstallerRelease {

    my $numversion = $major.'.'.$minor.'.'.$micro;
    my $cmd;

    $RELEASE_PATH = $ENV{ PWD}."/plasma-installer_".$numversion;

    # Sauvegarde du rep courant
    my $dir = `pwd`; chop $dir;

    $cmd = 'svn export --force '.$NUMREL.' '.$user.' '.$svninst.' '.$RELEASE_PATH;
    myCmd($cmd);
    
    # Save the InstallationGuide if we want to do a plasma-installer release
    myCmd("cp InstallationGuide-v${numversion} $RELEASE_PATH/README");

    #Create tarball
    print "Create the tarball\n";
    $DIRNAME=`dirname $RELEASE_PATH`;
    $BASENAME=`basename $RELEASE_PATH`;
    chop $DIRNAME;
    chop $BASENAME;
    myCmd("(cd $DIRNAME && tar cvzf ${BASENAME}.tar.gz $BASENAME)");
}

sub Usage {
    
    print "MakeRelease.pl [ -h ][ -d Directory ] [ -u username ] [ -r numrelease ] Major Minor Micro\n";
    print "   -h   Print this help\n";
    print "   -d   Choose directory for release\n";
    print "   -r   Choose svn release number\n";
    print "   -s   Choose plasma directory for export\n";
    print "   -u   username\n";

}

getopts("hd:u:r:s:",\%opts);

if ( defined $opts{h} ){
    Usage();
    exit;
}

if (defined $opts{d}){
    $RELEASE_PATH = $opts{d};
}
if (defined $opts{u}){
    $user = "--username $opts{u}";
}

if (defined $opts{r}){
    $NUMREL = "-r $opts{r}";
}
if (defined $opts{s}){
    $svn = $opts{s};
}

if ( ($#ARGV + 1) < 3 ) {
    Usage();
    exit;
}

$major = $ARGV[0];
$minor = $ARGV[1];
$micro = $ARGV[2];

MakeRelease();
MakeInstallerRelease();

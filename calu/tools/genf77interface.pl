#!/usr/bin/perl

use strict;
use POSIX;


my %lapack_functions = ();
my %tile_functions = ();
my %tile_async_functions = ();
my %workspace_functions = ();

sub fillStructures {

    my ($filename) = @_;

    open FILE, $filename or die $!;

    my $complex = 0;

    while (<FILE>) {
        my $line = $_;
        chop $line;

        # Enter a complex section
        if ( $line =~ /#ifdef COMPLEX/ ) {
            $complex = 1;
            next;
        }

        # End a complex section
        if ( $line =~ /#endif/ && $complex ) {
            $complex = 0;
            next;
        }

        if ( $line =~ /^(\w*)\s*PLASMA_(.*)\((.*)\);/ ) {
            my $rettype  = $1;
            my $fullname = $2;
            my $funcname = $2;
            my $args     = $3;

            #print $rettype."\n";
            #print $funcname."\n";
            #print $3."\n";

            if ( $funcname =~ /Alloc_Workspace_(.*)/ ) {
                $funcname = $1;
                $workspace_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_ALLOC_WORKSPACE_%s", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args  };
            }
            elsif ( $funcname =~ /Lapack_to_Tile_Async/ ||
                    $funcname =~ /Tile_to_Lapack_Async/ ) {
                $lapack_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_%s", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args   };
            }
            elsif ( $funcname =~ /Lapack_to_Tile/ ||
                    $funcname =~ /Tile_to_Lapack/ ) {
                $lapack_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_%s", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args   };
            }
            elsif ( $funcname =~ /(.*)_Tile_Async/ ) {
                $funcname = $1;
                $tile_async_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_%s_TILE_ASYNC", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args   };
            }
            elsif ( $funcname =~ /(.*)_Tile/ ) {
                $funcname = $1;
                $tile_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_%s_TILE", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args   };
            } else {
                $lapack_functions{ $funcname } = {
                    'fullname'=> $fullname,
                    'full'    => sprintf("PLASMA_%s", uc( $funcname )),
                    'lower'   => lc( $funcname ),
                    'upper'   => uc( $funcname ),
                    'rettype' => $rettype,
                    'complex' => $complex,
                    'args'    => $args };
            }
        }
    }
    close FILE;
}

sub genHeader {

    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (simple interface)\n";
    print " **/\n";
    foreach my $key (sort keys %lapack_functions )
    {
        my $str = sprintf( "%s#define %-22s PLASMA_FNAME(%-13s, %-13s)\n%s",
                           $lapack_functions{ $key }{ 'complex' } ? "#ifdef COMPLEX\n" : "",
                           $lapack_functions{ $key }{ 'full' },
                           $lapack_functions{ $key }{ 'lower' },
                           $lapack_functions{ $key }{ 'upper' },
                           $lapack_functions{ $key }{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }

    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (native interface)\n";
    print " **/\n";
    foreach my $key (sort keys %tile_functions )
    {
        my $str = sprintf( "%s#define %-27s PLASMA_TILE_FNAME(%-13s, %-13s)\n%s",
                           $tile_functions{ $key }{ 'complex' } ? "#ifdef COMPLEX\n" : "",
                           $tile_functions{ $key }{ 'full' },
                           $tile_functions{ $key }{ 'lower' },
                           $tile_functions{ $key }{ 'upper' },
                           $tile_functions{ $key }{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }

    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (asynchronous interface)\n";
    print " **/\n";
    foreach my $key (sort keys %tile_async_functions )
    {
        my $str = sprintf( "%s#define %-33s PLASMA_ASYNC_FNAME(%-13s, %-13s)\n%s",
                           $tile_async_functions{ $key }{ 'complex' } ? "#ifdef COMPLEX\n" : "",
                           $tile_async_functions{ $key }{ 'full' },
                           $tile_async_functions{ $key }{ 'lower' },
                           $tile_async_functions{ $key }{ 'upper' },
                           $tile_async_functions{ $key }{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }

    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - workspace allocation\n";
    print " **/\n";
    foreach my $key (sort keys %workspace_functions )
    {
        my $str = sprintf( "%s#define %-42s PLASMA_WS_FNAME(%-18s, %-18s)\n%s",
                           $workspace_functions{ $key }{ 'complex' } ? "#ifdef COMPLEX\n" : "",
                           $workspace_functions{ $key }{ 'full' },
                           $workspace_functions{ $key }{ 'lower' },
                           $workspace_functions{ $key }{ 'upper' },
                           $workspace_functions{ $key }{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }
    print "\n";
    print "\n";
}

sub genOneGroupOfInterface {

    my ($work, %functions) = @_;

    foreach my $key (sort keys %functions )
    {
        my $retval;
        my $params;
        my $callargs;

        my %function = %{$functions{ $key }};
        my @arguments = split(/,/, $function{'args'} );

        if ( $function{ 'rettype' } =~ /int/ ) {
            $retval = "info";
        } else {
            $retval = "result";
        }

        foreach my $arg (@arguments) {
            my $type = $arg;
            my $varname;

            $type =~ s/ *(.*[ \*])([a-zA-Z0-9_]+)/\1/;
            $varname = $2;

            if ( !($type =~ /^.*\*$/) ||
                 ( !$work && ($varname =~ /desc[LT]/ ) ) ) {
                $params   .= $type.'*'.$varname.', ';
                $callargs .= '*'.$varname.', ';
            } else {
                $params   .= $type.$varname.', ';
                $callargs .= $varname.', ';
            }
        }

        $params .= $function{'rettype'}.' *'.$retval;
        chop $callargs; chop $callargs;

        my $str = sprintf( "%svoid %s(%s)\n{ *%s = PLASMA_%s(%s); }\n%s\n",
                           $function{ 'complex' } ? "#ifdef COMPLEX\n" : "",
                           $function{ 'full' },
                           $params, $retval,
                           $function{ 'fullname' }, $callargs,
                           $function{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }
}

sub genInterfaces {

    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (simple interface)\n";
    print " **/\n";
    genOneGroupOfInterface( 0, %lapack_functions );
    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (native interface)\n";
    print " **/\n";
    genOneGroupOfInterface( 0, %tile_functions );
    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - math functions (asynchronous interface)\n";
    print " **/\n";
    genOneGroupOfInterface( 0, %tile_async_functions );
    print "\n";
    print "/***************************************************************************//**\n";
    print " *  FORTRAN API - workspace allocation\n";
    print " **/\n";
    genOneGroupOfInterface( 1, %workspace_functions );
    print "\n";
    print "\n";
}


print '/**'."\n";
print ' *'."\n";
print ' * @file plasma_zf77.c'."\n";
print ' *'."\n";
print ' *  PLASMA computational routines'."\n";
print ' *  PLASMA is a software package provided by Univ. of Tennessee,'."\n";
print ' *  Univ. of California Berkeley and Univ. of Colorado Denver'."\n";
print ' *  This file is automatically generated by tools/genf77interface.pl'."\n";
print ' *'."\n";
print ' * @version 2.5.0'."\n";
print ' * @author Bilel Hadri'."\n";
print ' * @author Mathieu Faverge'."\n";
print ' * @date 2010-11-15'."\n";
print ' * @precisions normal z -> c d s'."\n";
print ' *'."\n";
print ' **/'."\n";
print '#include <stdlib.h>'."\n";
print '#include "common.h"'."\n";
print '#undef REAL'."\n";
print '#define COMPLEX'."\n";
print "\n";

fillStructures( "include/plasma_z.h" );
genHeader();
genInterfaces();

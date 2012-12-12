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

        if ( $line =~ /^(\w*)\s*PLASMA_([a-zA-Z0-9_]*) *\((.*)\);/ ) {
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

sub genOneGroupOfInterface {

    my ($work, %functions) = @_;

    foreach my $key (sort keys %functions )
    {
        my $retval;
        my $params;
        my $code;

        my %function = %{$functions{ $key }};
        my @arguments = split(/,/, $function{'args'} );

        if ( $function{ 'rettype' } =~ /int/ ) {
            $retval = "integer(kind=c_int)";
        } elsif ( $function{ 'rettype' } =~ /double/ ) {
            $retval = "real(kind=c_double)";
        } else {
            print 'Return type '.$function{ 'rettype' }." unknown\n";
            die;
        }

        foreach my $arg (@arguments) {
            my $type = $arg;
            my $varname;
            my $ftype;

            $type =~ s/ *(.*[ \*])([a-zA-Z0-9_]+)/\1/;
            $varname = $2;

            $type =~ s/(.*) +/$1/;

            if ( !($type =~ /^.*\*$/) ) {
                if ( $type =~ /unsigned long long int/ ) {
                    $ftype = "integer(kind=c_long_long)";
                } elsif ( $type =~ /int/ || $type =~ /PLASMA_enum/ ) {
                    $ftype = "integer(kind=c_int)";
                } elsif ( $type =~ /double/ ) {
                    $ftype = "real(kind=c_double)";
                } elsif ( $type =~ /PLASMA_Complex64_t/ ) {
                    $ftype = "complex(kind=c_double_complex)";
                } else {
                    print $type." unknown\n";;
                    die;
                }
                $code .= "            ".$ftype.", value :: ".$varname."\n";
            }
            # Pointers
            else {
                if ( !($type =~ /^.*\*\*$/) ) {
                    $code .= "            type(c_ptr), value :: ".$varname."\n";
                } else {
                    $code .= "            type(c_ptr) :: ".$varname." ! ".$varname." is **, so pass by reference\n";
                }
            }
            $params .= $varname.',';
        }

        chop $params;
        chop $code;
        
        my $str = sprintf( 
            "
%s      interface
         function PLASMA_%s_c(%s) &
          & bind(c, name='PLASMA_%s')
            use iso_c_binding
            implicit none
            %s :: PLASMA_%s_c
%s
          end function PLASMA_%s_c
      end interface
%s",
            $function{ 'complex' } ? "#if defined(PRECISION_z) || defined(PRECISION_c)\n" : "",
            $function{ 'fullname' }, $params,
            $function{ 'fullname' },
            $retval, $function{ 'fullname' },
            $code,
            $function{ 'fullname' },
            $function{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }
}

sub genInterfaces {

    print '
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (simple interface)
    !';
    genOneGroupOfInterface( 0, %lapack_functions );
    print '
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (native interface)
    !';
    genOneGroupOfInterface( 0, %tile_functions );
    print '
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - math functions (asynchronous interface)
    !';
    genOneGroupOfInterface( 0, %tile_async_functions );
    print '
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  FORTRAN API - workspace allocation
    !';
    genOneGroupOfInterface( 1, %workspace_functions );
    print "\n";
}

sub genOneGroupOfSubRoutine {

    my ($work, %functions) = @_;

    foreach my $key (sort keys %functions )
    {
        my $retval;
        my $params;
        my $code;

        my %function = %{$functions{ $key }};
        my @arguments = split(/,/, $function{'args'} );

        if ( $function{ 'rettype' } =~ /int/ ) {
            $retval = "integer(kind=c_int)";
        } elsif ( $function{ 'rettype' } =~ /double/ ) {
            $retval = "real(kind=c_double)";
        } else {
            print 'Return type '.$function{ 'rettype' }." unknown\n";
            die;
        }

        foreach my $arg (@arguments) {
            my $type = $arg;
            my $varname;
            my $ftype;

            $type =~ s/ *(.*[ \*])([a-zA-Z0-9_]+)/\1/;
            $varname = $2;

            $type =~ s/(.*) +/$1/;

            if ( !($type =~ /^.*\*$/) ) {
                if ( $type =~ /int/ || $type =~ /PLASMA_enum/ ) {
                    $ftype = "integer(kind=c_int)";
                } elsif ( $type =~ /double/ ) {
                    $ftype = "real(kind=c_double)";
                } elsif ( $type =~ /PLASMA_Complex64_t/ ) {
                    $ftype = "complex(kind=c_double_complex)";
                } else {
                    print $type." unknown\n";;
                    die;
                }
                $code .= "            ".$ftype.", value :: ".$varname."\n";
            }
            # Pointers
            else {
                $code .= "            type(c_ptr), value :: ".$varname."\n";
            }
            $params .= $varname.',';
        }

        chop $params;
        chop $code;
        
        my $str = sprintf( 
            "
%s
      subroutine PLASMA_%s_c(%s)
         use iso_c_binding
         implicit none
            %s :: PLASMA_%s_c
%s
          end function PLASMA_%s_c
      end interface
%s",
            $function{ 'complex' } ? "#if defined(PRECISION_z) || defined(PRECISION_c)\n" : "",
            $function{ 'fullname' }, $params,
            $function{ 'fullname' },
            $retval, $function{ 'fullname' },
            $code,
            $function{ 'fullname' },
            $function{ 'complex' } ? "#endif\n" : "" );
        print $str;
    }
}

print '!
!     Copyright © 2011 The Numerical Algorithms Group Ltd. All rights reserved.
!   
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are
!     met:
!     - Redistributions of source code must retain the above copyright notice,
!       this list of conditions, and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer listed in
!       this license in the documentation and/or other materials provided with
!       the distribution.
!     - Neither the name of the copyright holders nor the names of its
!       contributors may be used to endorse or promote products derived from
!       this software without specific prior written permission.
!     
!     This software is provided by the copyright holders and contributors "as
!     is" and any express or implied warranties, including, but not limited
!     to, the implied warranties of merchantability and fitness for a
!     particular purpose are disclaimed. in no event shall the copyright owner
!     or contributors be liable for any direct, indirect, incidental, special,
!     exemplary, or consequential damages (including, but not limited to,
!     procurement of substitute goods or services; loss of use, data, or
!     profits; or business interruption) however caused and on any theory of
!     liability, whether in contract, strict liability, or tort (including
!     negligence or otherwise) arising in any way out of the use of this
!     software, even if advised of the possibility of such damage.
!
!
! @file plasma_zf90.F90
!
!  PLASMA Fortran 90 interfaces using Fortran 2003 ISO C bindings
!  PLASMA is a software package provided by Univ. of Tennessee,
!  Univ. of California Berkeley and Univ. of Colorado Denver
!
! @version 2.4.6
! @author Numerical Algorithms Group
! @author Mathieu Faverge
! @date 2011-12-15
! @precisions normal z -> c d s
!
#define PRECISION_z

module plasma_z';

fillStructures( "include/plasma_z.h" );
genInterfaces();
print "  contains\n";
#genSubRoutine();

###
#
#  @file scotch64.rb
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for Scotch 6.*
#
#  @version 6.1.1
#  @author Pierre Ramet
#  @author Mathieu Faverge
#  @date 2023-04-06
#
###
class Scotch64 < Formula
  desc "Package for graph and mesh partitioning"
  homepage "https://gitlab.inria.fr/scotch/scotch"
  url "https://gitlab.inria.fr/scotch/scotch/-/archive/v6.1.1/scotch-v6.1.1.tar.gz"
  sha256 "14daf151399fc67f83fd3ff2933854f5e8d2207c7d35dd66a05660bf0bbd583c"
  revision 1

  option "without-test", "skip build-time tests (not recommended)"
  deprecated_option "without-check" => "without-test"

  depends_on "open-mpi"
  depends_on "xz" => :optional # Provides lzma compression.

  conflicts_with "scotch", because: "not compatible"

  def install
    cd "src" do
      ln_s "Make.inc/Makefile.inc.i686_mac_darwin10", "Makefile.inc"
      # MPI implementation is not threadsafe, do not use DSCOTCH_PTHREAD

      cflags   = %w[-O3 -fPIC -Drestrict=__restrict -DCOMMON_PTHREAD_BARRIER
                    -DIDXSIZE64 -DINTSIZE64
                    -DCOMMON_PTHREAD
                    -DSCOTCH_CHECK_AUTO -DCOMMON_RANDOM_FIXED_SEED
                    -DCOMMON_TIMING_OLD -DSCOTCH_RENAME
                    -DCOMMON_FILE_COMPRESS_BZ2 -DCOMMON_FILE_COMPRESS_GZ]
      ldflags  = %w[-lm -lz -lpthread -lbz2]

      cflags  += %w[-DCOMMON_FILE_COMPRESS_LZMA]   if build.with? "xz"
      ldflags += %W[-L#{Formula["xz"].lib} -llzma] if build.with? "xz"

      make_args = ["CCS=#{ENV["CC"]}",
                   "CCP=mpicc",
                   "CCD=mpicc",
                   "FC=gfortran",
                   "RANLIB=echo",
                   "CFLAGS=#{cflags.join(" ")}",
                   "LDFLAGS=#{ldflags.join(" ")}"]

      if OS.mac?
        make_args << "LIB=.dylib"
        make_args << "ARFLAGS=-dynamiclib -Wl,-single_module -install_name #{lib}/$(notdir $@) -undefined dynamic_lookup -o "
        if MacOS.version >= :mojave
          make_args << "AR=#{ENV["CC"]}"
        else
          make_args << "AR=libtool"
        end
      else
       make_args << "LIB=.so"
       make_args << "ARCH=ar"
       make_args << "ARCHFLAGS=-ruv"
      end

      system "make", "scotch", "VERBOSE=ON", *make_args
      system "make", "ptscotch", "VERBOSE=ON", *make_args
      system "make", "install", "prefix=#{prefix}", *make_args
      system "make", "check", "ptcheck", "EXECP=mpirun -np 2", *make_args unless build.without? "test"
    end

    # Install documentation + sample graphs and grids.
    doc.install Dir["doc/*"]
    (share / "scotch").install "grf", "tgt"
  end

  test do
    mktemp do
      system "echo cmplt 7 | #{bin}/gmap #{share}/scotch/grf/bump.grf.gz - bump.map"
      system "#{bin}/gmk_m2 32 32 | #{bin}/gmap - #{share}/scotch/tgt/h8.tgt brol.map"
      system "#{bin}/gout -Mn -Oi #{share}/scotch/grf/4elt.grf.gz #{share}/scotch/grf/4elt.xyz.gz - graph.iv"
    end
  end
end

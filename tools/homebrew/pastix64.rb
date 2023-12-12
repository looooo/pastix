###
#
#  @file pastix64.rb
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for PaStiX 6.*
#
#  @version 6.3.2
#  @author Pierre Ramet
#  @author Mathieu Faverge
#  @date 2023-07-21
#
###
class Pastix64 < Formula
  desc "Parallel solver for sparse linear systems based on direct methods"
  homepage "https://gitlab.inria.fr/solverstack/pastix"
  url "https://files.inria.fr/pastix/releases/v6/pastix-6.3.2.tar.gz"
  sha256 "c4da8802d1933eecf8c09d7e63c014c81ccf353fe623142e9f5c5fc65ed82ee0"
  head "https://gitlab.inria.fr/solverstack/pastix.git"
  license "LGPL"

  depends_on "openblas"
  depends_on "openmpi"
  depends_on "gfortran"
  depends_on "hwloc"               # Could be optional but strongly recommanded
  depends_on "scotch64"            # Could be optional but strongly recommanded
  depends_on "metis"  => :optional # Enable METIS ordering
  depends_on "starpu" => :optional # Enable STARPU runtime
  depends_on "gcc"    => :build    # GNU Fortran is now provided as part of GCC
  depends_on "cmake"  => :build

  conflicts_with "pastix", because: "not compatible"

  def install
    ENV.append "LDFLAGS", "-Wl,-ld_classic"
    args = ["-DCMAKE_INSTALL_PREFIX=#{prefix}",
            "-DBUILD_SHARED_LIBS=ON",
            "-DBUILD_DOCUMENTATION=OFF",
            "-DPASTIX_INT64=ON",
            "-DPASTIX_ORDERING_SCOTCH=ON",
            "-DPASTIX_WITH_FORTRAN=ON",
            "-DPASTIX_WITH_MPI=ON",
            "-DPASTIX_WITH_CUDA=OFF",
            "-DPASTIX_WITH_PARSEC=OFF"]
    args += ["-DPASTIX_ORDERING_METIS=ON"] if build.with? "metis"
    args += ["-DPASTIX_WITH_STARPU=ON"] if build.with? "starpu"
    mkdir "build" do
      system "cmake", "..", *args
      system "make"
      system "make", "install"
    end
    pkgshare.install "example" # Contains all test programs.
  end

  def caveats; <<~EOS
    Set the PYTHONPATH environment variable:
      export PYTHONPATH=#{prefix}/lib/python:$PYTHONPATH
    Try python example with:
      python #{prefix}/examples/python/simple.py
    Or simple example with:
      #{prefix}/examples/simple -9 10:10:10
    EOS
  end

  test do
    system "#{prefix}/examples/simple", "-9", "10:10:10"
    ohai "All test output is in #{HOMEBREW_LOGS}/pastix. Please check."
  end
end

###
#
#  @file pastix6.rb
#  @copyright 2013-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for PaStiX 6
#
#  @version 6.0.3
#  @author Pierre Ramet
#  @date 2019-11-12
#
###
class Pastix6 < Formula
  desc "Parallel solver for sparse linear systems based on direct methods"
  homepage "https://gitlab.inria.fr/solverstack/pastix"
  url "https://gitlab.inria.fr/solverstack/pastix//uploads/98ce9b47514d0ee2b73089fe7bcb3391/pastix-6.1.0.tar.gz"
  sha256 "3d8df8ac6663488e107b56a1a55a7b6cd25ae29a31e893fee27f4cddef72a7d0"
  head "https://gitlab.inria.fr/solverstack/pastix.git"

  bottle :disable, "needs to be rebuilt"

  depends_on "openblas"
  depends_on "hwloc"              # Could be optinal but strongly recommanded
  depends_on "scotch5"            # Could be optinal but strongly recommanded
  depends_on "metis" => :optional # Use METIS ordering.
  depends_on "gcc"   => :build    # GNU Fortran is now provided as part of GCC
  depends_on "cmake" => :build

  def install
    args = ["-DCMAKE_INSTALL_PREFIX=#{prefix}",
            "-DBUILD_SHARED_LIBS=ON",
            "-DBUILD_DOCUMENTATION=OFF",
            "-DBUILD_64bits=OFF",
            "-DPASTIX_INT64=OFF",
            "-DPASTIX_ORDERING_SCOTCH=ON",
            "-DPASTIX_WITH_FORTRAN=ON",
            "-DPASTIX_WITH_MPI=OFF",
            "-DPASTIX_WITH_CUDA=OFF",
            "-DPASTIX_WITH_STARPU=OFF",
            "-DPASTIX_WITH_PARSEC=OFF"]
    args += ["-DPASTIX_ORDERING_METIS=ON"] if build.with? "metis"
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
      python #{prefix}/examples/simple.py
    Or simple example with:
      #{prefix}/examples/simple -9 10:10:10
    EOS
  end

  test do
    system "#{prefix}/examples/simple", "-9", "10:10:10"
    ohai "All test output is in #{HOMEBREW_LOGS}/pastix. Please check."
  end
end

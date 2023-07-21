###
#
#  @file pastix64.rb
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for PaStiX 6.*
#
#  @version 6.3.0
#  @author Pierre Ramet
#  @author Mathieu Faverge
#  @date 2023-04-06
#
###
class Pastix64 < Formula
  desc "Parallel solver for sparse linear systems based on direct methods"
  homepage "https://gitlab.inria.fr/solverstack/pastix"
  url "https://gitlab.inria.fr/solverstack/pastix//uploads/3fc798eb3ca6282e21506349df7f9da2/pastix-6.2.1.tar.gz"
  sha256 "b680cbfc265df8cba18d3a7093fcc02e260198c4a2d6a86d1e684bb291e309dd"
  head "https://gitlab.inria.fr/solverstack/pastix.git"
  license "LGPL"

  depends_on "openblas"
  depends_on "openmpi"
  depends_on "hwloc"               # Could be optional but strongly recommanded
  depends_on "scotch64"            # Could be optional but strongly recommanded
  depends_on "metis"  => :optional # Enable METIS ordering
  depends_on "starpu" => :optional # Enable STARPU runtime
  depends_on "gcc"    => :build    # GNU Fortran is now provided as part of GCC
  depends_on "cmake"  => :build

  conflicts_with "pastix", because: "not compatible"

  def install
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

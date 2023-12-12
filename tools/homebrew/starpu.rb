###
#
#  @file starpu.rb
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for StarPU 1.3.*
#
#  @version 6.3.2
#  @author Pierre Ramet
#  @author Mathieu Faverge
#  @date 2023-07-21
#
###
class Starpu < Formula
  desc "StarPU is a task programming library for hybrid architectures"
  homepage "https://starpu.gitlabpages.inria.fr/"
  url "https://files.inria.fr/starpu/starpu-1.3.10/starpu-1.3.10.tar.gz"
  sha256 "757cd9a54f53751d37364965ac36102461a85df3a50b776447ac0acc0e1e2612"
  license "GNU GPL v2.1"

  depends_on "autoconf" => :build
  depends_on "automake" => :build
  depends_on "libtool" => :build
  depends_on "pkg-config" => [:build, :test]
  depends_on "hwloc"
  depends_on "openmpi"

  def install
    system "./autogen.sh" if build.head?
    system "./configure", *std_configure_args
    system "make", "install"
  end

  test do
    (testpath/"test.c").write <<~EOS
      #include <stdio.h>
      #include <stdlib.h>
      #include <starpu.h>

      struct starpu_codelet cl =
      {
        .where = STARPU_NOWHERE,
      };

      int main(int argc, char* argv[])
      {
        int ret = starpu_init(NULL);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");
        ret = starpu_task_insert(&cl, 0);
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
        ret = starpu_task_wait_for_all();
	STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_wait_for_all");
        starpu_shutdown();
        return 0;
      }
    EOS

    pkg_config_flags = `pkg-config --cflags --libs starpu-1.3`.chomp.split
    system ENV.cc, "test.c", *pkg_config_flags
    system "./a.out"
  end
end

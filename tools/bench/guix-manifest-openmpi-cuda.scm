;; What follows is a "manifest" equivalent to the command line you gave.
;; You can store it in a file that you may then pass to any 'guix' command
;; that accepts a '--manifest' (or '-m') option.

(concatenate-manifests
  (list (specifications->manifest
          (list "coreutils"
                "gawk"
                "grep"
                "jube"
                "mkl"
                "nss-certs"
                "openssh"
                "perl"
                "python-click"
                "python-certifi"
                "python-elasticsearch"
                "python-gitpython"
                "python-matplotlib"
                "python-pandas"
                "python-seaborn"
                "r-ggplot2"
                "r-plyr"
                "r-reshape2"
                "sed"
                "slurm"
                "zlib"))
        (package->development-manifest
          (specification->package "pastix-cuda"))))

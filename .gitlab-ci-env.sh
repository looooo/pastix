#export version_cppcheck=1.79
#export PATH=/builds/sonar/cppcheck-${version_cppcheck}/install/bin:$PATH
export version_drmemory=1.11.0-2
export PATH=/builds/sonar/DrMemory-Linux-${version_drmemory}/bin:$PATH
#export PATH=/builds/sonar/sonar-scanner-2.8/bin:$PATH
export PATH=/builds/sonar/sonar-scanner-2.9.0.670/bin:$PATH
export PATH=/builds/sonar/lcov-to-cobertura-xml-1.6/lcov_cobertura:$PATH

#StarPU
export LD_LIBRARY_PATH=/builds/install/starpu/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/builds/install/starpu/lib/pkgconfig:$PKG_CONFIG_PATH
#PaRSEC
export PATH=/builds/install/parsec/bin:$PATH
export LD_LIBRARY_PATH=/builds/install/parsec/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/builds/install/parsec/lib/pkgconfig:$PKG_CONFIG_PATH
#Scotc
export LD_LIBRARY_PATH=/builds/install/scotch/lib:$LD_LIBRARY_PATH
export SCOTCH_DIR=/builds/install/scotch/

#export version_cppcheck=1.79
#export PATH=/builds/sonar/cppcheck-${version_cppcheck}/install/bin:$PATH
export version_drmemory=1.11.0-2
export PATH=/builds/sonar/DrMemory-Linux-${version_drmemory}/bin:$PATH
#export PATH=/builds/sonar/sonar-scanner-2.8/bin:$PATH
export PATH=/builds/sonar/sonar-scanner-2.9.0.670/bin:$PATH
export COBERTURA_PATH=/builds/sonar/lcov-to-cobertura-xml-1.6/lcov_cobertura/

#StarPU
export LD_LIBRARY_PATH=/builds/tmp/starpu/install/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/builds/tmp/starpu/install/lib/pkgconfig:$PKG_CONFIG_PATH
#PaRSEC
export PATH=/builds/tmp/parsec/install/bin:$PATH
export LD_LIBRARY_PATH=/builds/tmp/parsec/install/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/builds/tmp/parsec/install/lib/pkgconfig:$PKG_CONFIG_PATH
#Scotc
export LD_LIBRARY_PATH=/builds/tmp/scotch/lib:$LD_LIBRARY_PATH
export SCOTCH_DIR=/builds/tmp/scotch/

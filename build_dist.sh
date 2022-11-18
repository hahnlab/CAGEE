mkdir distro

tar --transform "s,^,CAGEE/," --owner=0 --group=0 -czf distro/CAGEE0.1.0.0.tar.gz src test.cpp main.cpp diffmat_precalc.cpp CHANGELOG.md LICENSE CMakeLists.txt README.md config.h.in examples LBFGSpp/ docs/manual

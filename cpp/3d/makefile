all:
	mex -v GCC='/usr/bin/gcc-4.9' FLAGS='$CFLAGS -O3' MC3Dmex.cpp
	g++ MC3D.cpp -o MC3D.a -O3 
	cp MC3D.a ../..
	cp MC3D*.mex* ../..

debug:
	g++ MC3D.cpp -o MC3D.a -ggdb -fsanitize=address -fno-omit-frame-pointer

parallel:
	g++ MC3D.cpp -o MC3D.a -O3 -DUSE_OMP -fopenmp -O3
	cp MC3D.a ../..
	cp MC3D*.mex* ../..

parallelmex:
	mex -DUSE_OMP -v GCC='/usr/bin/gcc-4.9' CXXFLAGS='${CXXFLAGS} -fopenmp -fpic -O3' MC3Dmex.cpp LDFLAGS='${LDFLAGS} -fopenmp'
	cp MC3D*.mex* ../..


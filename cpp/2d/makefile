all:
	mex -v GCC='/usr/bin/gcc-4.9' FLAGS='$CFLAGS -O3' MC2Dmex.cpp
	g++ MC2D.cpp -o MC2D.a -O3
	cp MC2D.a ../..
	cp MC2D*.mex* ../..

debugsymbols:
	g++ MC2D.cpp -o MC2D.a -g -D_GLIBCXX_USE_CXX11_ABI=0 
	cp MC2D.a ../.. 


debugaddresses:
	g++ MC2D.cpp -o MC2D.a -fsanitize=address -fno-omit-frame-pointer
	cp MC2D.a ../.. 

parallel:
	g++ MC2D.cpp -o MC2D.a -O3 -DUSE_OMP -fopenmp -O3
	cp MC2D.a ../..

parallelmex:
	mex -DUSE_OMP -v GCC='/usr/bin/gcc-4.9' CXXFLAGS='${CXXFLAGS} -fopenmp -fpic -O3' MC2Dmex.cpp LDFLAGS='${LDFLAGS} -fopenmp'

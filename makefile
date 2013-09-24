bin/aban2: src/mesh.h src/domain.h src/main.cpp src/vector.h src/advection.h src/solver.h src/diffusion.h src/pressure.h src/tests.h
	g++ -std=c++11 src/main.cpp -ljsoncpp -o bin/aban2
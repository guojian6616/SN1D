#include "Solver.h"


#include <iostream>

int main() {

	Quadrature* quad = new Quadrature(6);
	// quad->initializeArray();
	Solver* solver = new Solver(10000,6);
	solver->setQuadrature(quad);
	solver->run();

	std::cout << "Residual = " << solver->getResidual() << std::endl;

	return 0;
}

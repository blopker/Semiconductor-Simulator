/*
 * Laplacian.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#include "Laplacian.h"
using namespace arma;
using namespace std;

void laplacian(int N, double dx, int lap_kind, mat* lap) {
	dx = pow(dx,2);
	if (lap_kind==1) {
		lap->diag() = ones<vec>(N)*-2/dx;
		lap->diag(-1) = ones<vec>(N-1)/dx;
		lap->diag(1) = ones<vec>(N-1)/dx;
	}
}

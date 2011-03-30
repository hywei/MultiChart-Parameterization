#include "solver.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Solver::Solver(int nb_variables) {
	nb_variables_ = nb_variables ;
	variable_ = new SolverVariable[nb_variables] ;

}

Solver::~Solver() {
	delete[] variable_ ;
	variable_ = 0 ;
}
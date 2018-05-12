#include <iostream>
#include "solver.h"

int main()
{
    structan::solver_t solver("models/test.json");
    solver.write("aa0.vtk");
    solver.step();
    // solver.step();
    // solver.step();
    // solver.step();
    // solver.step();
    // solver.step();
    // solver.step();
    // solver.step();
    solver.write("aa1.vtk");
    return 0;
}

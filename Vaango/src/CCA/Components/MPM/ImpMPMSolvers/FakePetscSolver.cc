/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/ImpMPMSolvers/PetscSolver.h>
#include <Core/Exceptions/InternalError.h>


using namespace Uintah;


//#define LOG
#undef LOG
#undef DEBUG_PETSC

MPMPetscSolver::MPMPetscSolver()
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::MPMPetscSolver(), specify simple instead of petsc in input file!", __FILE__, __LINE__ );

}

MPMPetscSolver::~MPMPetscSolver()
{
}


void MPMPetscSolver::initialize()
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::initialize()!", __FILE__, __LINE__ );

}

void 
MPMPetscSolver::createLocalToGlobalMapping([[maybe_unused]] const ProcessorGroup* d_myworld,
                                           [[maybe_unused]] const PatchSet* perproc_patches,
                                           [[maybe_unused]] const PatchSubset* patches,
                                           [[maybe_unused]] const int DOFsPerNode,
                                           [[maybe_unused]] const int n8or27)
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::createLocalToGlobalMapping()!", __FILE__, __LINE__ );
}

void MPMPetscSolver::solve([[maybe_unused]] vector<double>& guess)
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::solve()!", __FILE__, __LINE__ );
}

void MPMPetscSolver::createMatrix([[maybe_unused]] const ProcessorGroup* d_myworld,
                                  [[maybe_unused]] const  std::map< int,int>& dof_diag)
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::createMatrix()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::destroyMatrix([[maybe_unused]] bool recursion)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::destroyMatrix()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::flushMatrix()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::flushMatrix()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::fillVector([[maybe_unused]] int i,[[maybe_unused]] double v,[[maybe_unused]] bool add)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::fillVector()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::fillTemporaryVector([[maybe_unused]] int i,[[maybe_unused]] double v)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::fillTemporaryVector()!", __FILE__, __LINE__ );
}

void MPMPetscSolver::fillFluxVector([[maybe_unused]] int i,[[maybe_unused]] double v)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::fillFluxVector()!", __FILE__, __LINE__ );
}



void MPMPetscSolver::assembleVector()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::assembleVector()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::assembleTemporaryVector()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::assembleTemporaryVector()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::applyBCSToRHS()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::applyBCSToRHS()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::copyL2G([[maybe_unused]] Array3<int>& mapping,[[maybe_unused]] const Patch* patch)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::copyL2G()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::removeFixedDOF()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::removeFixedDOF()!", __FILE__, __LINE__ );
}

void MPMPetscSolver::removeFixedDOFHeat()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::removeFixedDOFHeat()!", __FILE__, __LINE__ );
}


void MPMPetscSolver::finalizeMatrix()
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::finalizeMatrix()!", __FILE__, __LINE__ );
}


int MPMPetscSolver::getSolution([[maybe_unused]] vector<double>& xPetsc)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::getSolution()!", __FILE__, __LINE__ );
}


int MPMPetscSolver::getRHS([[maybe_unused]] vector<double>& QPetsc)
{
 throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::getRHS()!", __FILE__, __LINE__ );
}

void MPMPetscSolver::fillMatrix([[maybe_unused]] int numi,[[maybe_unused]] int i[],[[maybe_unused]] int numj,
                                       [[maybe_unused]] int j[],[[maybe_unused]] double value[])
{
  throw InternalError( "Don't have PETSc so shouldn't be calling MPMPetscSolver::fillMatrix()!", __FILE__, __LINE__ );

}

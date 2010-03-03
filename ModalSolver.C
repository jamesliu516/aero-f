#include <IoData.h>
#include <Communicator.h>
#include <Domain.h>
#include <Modal.h>

//------------------------------------------------------------------------------

void startModalSolver(Communicator *com, IoData &ioData, Domain &domain)
{

  ModalSolver<5> mSolver(com, ioData, domain);
  domain.createVecPat(5, &ioData);
  mSolver.solve();
  
}

//------------------------------------------------------------------------------

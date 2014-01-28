#include <CCA/Components/Peridynamics/PeridynamicsCommon.h> 
#include <CCA/Components/Peridynamics/PeridynamicsSimulationState.h> 
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Vaango;

PeridynamicsCommon::PeridynamicsCommon(const Uintah::ProcessorGroup* myworld)
  : d_myworld(myworld)
{
}

PeridynamicsCommon::~PeridynamicsCommon()
{
}

void PeridynamicsCommon::materialProblemSetup(const Uintah::ProblemSpecP& prob_spec, 
                                     PeridynamicsSimulationStateP& sharedState,
                                     PeridynamicsFlags* flags)
{
  //Search for the MaterialProperties block and then get the Peridynamics section
  Uintah::ProblemSpecP mat_ps = prob_spec->findBlockWithOutAttribute("MaterialProperties");
  Uintah::ProblemSpecP peridynamics_mat_ps = mat_ps->findBlock("Peridynamics");

  for (Uintah::ProblemSpecP ps = peridynamics_mat_ps->findBlock("material"); ps != 0;
       ps = ps->findNextBlock("material") ) {

    std::string index("");
    ps->getAttribute("index",index);
    std::stringstream id(index);
    const int DEFAULT_VALUE = -1;
    int index_val = DEFAULT_VALUE;

    id >> index_val;

    if( !id ) {
      // stringstream parsing failed... on many (most) systems, the
      // original value assigned to index_val would be left
      // intact... but on some systems (redstorm) it inserts garbage,
      // so we have to manually restore the value.
      index_val = DEFAULT_VALUE;
    }
    // cout << "Material attribute = " << index_val << ", " << index << ", " << id << "\n";

    //Create and register as an Peridynamics material
    PeridynamicsMaterial *mat = scinew PeridynamicsMaterial(ps, sharedState, flags);
    // When doing restart, we need to make sure that we load the materials
    // in the same order that they were initially created.  Restarts will
    // ALWAYS have an index number as in <material index = "0">.
    // Index_val = -1 means that we don't register the material by its 
    // index number.
    if (index_val > -1){
      sharedState->registerPeridynamicsMaterial(mat, index_val);
    }
    else{
      sharedState->registerPeridynamicsMaterial(mat);
    }
  }
}

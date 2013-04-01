#include <Damage/BondDamageModel.h> 

using namespace Matiti;

BondDamageModel::BondDamageModel()
{
}

BondDamageModel::~BondDamageModel()
{
}

void 
BondDamageModel::initializeBrokenBonds(ProblemSpecP ps,
                                       DomainP domain,
                                       DataWarehouse* new_dw)
{
  // Read preexisting crack information from problem spec
  // the crack segments may be joined or disjointed
  ProblemSpecP crack_ps = ps->findBlock("Crack");
  int num_cracks = 0;
  crack_ps->get("num_cracks", num_cracks);
  vector<int> crack_id_set;
  vector<Point> crack_start_point_set;
  vector<Point> crack_end_point_set;
  for (int ii=0; ii < num_cracks; ii++) {
    int crack_id = 0;
    Point crack_start_point(0.0,0.0,0.0);
    Point crack_end_point(0.0,0.0,0.0);
    crack_ps->get("crack_id", crack_id);
    crack_ps->get("crack_start_point", crack_start_point);
    crack_ps->get("crack_end_point", crack_end_point);
    crack_id_set.push_back(crack_id);
    crack_start_point_set.push_back(crack_start_point);
    crack_end_point_set.push_back(crack_end_point);
  }

  // Get the variables needed
  MeshNodeSubset mset = domain->getMeshNodeSubset();
  constMeshNodeVariable<Point> mPosition;
  constMeshNodeVariable<double> mHorizon;
  constMeshNodeVariable<MeshNodeFamily> mNodeFamily;
  constMeshNodeVariable<MeshBondFamily> mBondFamily;

  old_dw->get(mPosition,     mPositionLabel, mset);
  old_dw->get(mHorizon,      mHorizonLabel,  mset);
  old_dw->get(mNodeFamily,   mNodeFamilyLabel,   mset);
  old_dw->get(mBondFamily,   mBondFamilyLabel,   mset);

  MeshNodeVariable<double> mDamageIndex;
  new_dw->allocateAndPut(mDamageIndex, mDamageIndexLabel, mset);
  
  // Loop over nodes
  for (MeshNodeIterator iter = mset->begin(); iter != mset->end(); iter++) {

    int numBrokenBonds = 0;
    int numBonds = 0;
    meshNodeIndex idx = *iter;

    // Loop through the bond family of the current mesh node
    MeshBondSubset mBondSet = mBondFamily[idx].getMeshBondSubset();
    for (MeshBondIterator bondIter = mBondSet->begin(); bondIter != mBondSet->end();
           bondIter++) {
      meshBondIndex bondIndex = *bondIter;
      numBonds++;
      for (jj = 0; jj < num_cracks; jj++) {
        if (mBondSet[bondIndex].intersect(crack_start_point_set[jj], crack_end_point_set[jj])) {
          numBrokenBonds++;
          mBondSet[bondIndex] = 0;
        }
      }
    }

    mDamageIndex[idx] = (double)numBrokenBonds/(double)numBonds;
  }

}

void BondDamageModel::computeCriticalStrain()
{
}

void BondDamageModel::updateBrokenBonds()
{
}


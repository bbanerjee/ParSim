#include <Solver/VelocityVerletIntegrator.h>

using namespace Matiti;

VelocityVerletIntegrator::VelocityVerletIntegrator()
{
}

VelocityVerletIntegrator::~VelocityVerletIntegrator()
{
}

VelocityVerletIntegrator::integrate()
{

  // Get the node list
  MeshNodeSet* nodes = new_dw->getMeshNodeSet(mesh);

  // Get the mesh node variables
  constMeshNodeVariable<bool> mOmitNode;
  constMeshNodeVariable<bool> mHangingNode;
  new_dw->get(mOmitNode, mOmitNodeLabel, nodes);
  new_dw->get(mHangingNode, mHangingNodeLabel, nodes);

  constMeshNodeVariable<Point> mPos_old;
  constMeshNodeVariable<Vector> mDisp_old;
  constMeshNodeVariable<Vector> mVel_old;
  constMeshNodeVariable<int> mNumFailedBonds;
  constMeshNodeVariable<double> mHorizonSize;
  new_dw->get(mPos_old, mPositionLabel, nodes);
  new_dw->get(mDisp_old, mDisplacementLabel, nodes);
  new_dw->get(mVel_old, mVelocityLabel, nodes);
  new_dw->get(mNumFailedBonds, mNumFailedBondsLabel, nodes);
  new_dw->get(mHorizonSize, mHorizonSizeLabel, nodes);

  MeshNodeVariable<Vector> mForce_new;
  MeshNodeVariable<double> mWeight_new;
  allocateTemporary(mForce_new, nodes);
  allocateTemporary(mWeight_new, nodes);

  // Iterate through the nodes
  MeshNodeSet::iterator = nodes->begin();
  for (; iter != nodes->end(); iter++) {
    meshNodeIndex idx = *iter;

    // If the node is to be omitted go to the next node
    if (mOmitNode[idx] || mHangingNode[idx]) continue;

    // Initialize
    mForce_new[idx] = Vector(0.0,0.0,0.0);
    mWeight_new[idx] = 0.0;

    // Bond broken?
    bool bondBroken = false;
    
  }
}

// Assume background uniform grid for search
VelocityVerletIntegrator::getMeshNodesInHorizon(const Point& pos,
                                                const double& horizon,
                                                const Domain* domain)
{
  // Find cell index in which the current node sits
  Point cellPos = domain->positionToIndex(pos);

  // Set up the search radius
  Point xPosMin = pos - Point(horizon, 0.0, 0.0); 
  Point zPosMin = pos - Point(0.0, 0.0, horizon); 
  Point yPosMin = pos - Point(0.0, horizon, 0.0); 
  Point xPosMax = pos + Point(horizon, 0.0, 0.0); 
  Point yPosMax = pos + Point(0.0, horizon, 0.0); 
  Point zPosMax = pos + Point(0.0, 0.0, horizon); 

  // Find cell index range for search in three directions
  Point xMinCellPos = domain->positionToIndex(xPosMin);
  Point yMinCellPos = domain->positionToIndex(yPosMin);
  Point zMinCellPos = domain->positionToIndex(zPosMin);
  Point xMaxCellPos = domain->positionToIndex(xPosMax);
  Point yMaxCellPos = domain->positionToIndex(yPosMax);
  Point zMaxCellPos = domain->positionToIndex(zPosMax);

  // Set up the range
  int ixmin = Floor(xMinCellPos.x());
  int iymin = Floor(yMinCellPos.y());
  int izmin = Floor(zMinCellPos.z());
  int ixmax = Floor(xMaxCellPos.x());
  int iymax = Floor(yMaxCellPos.y());
  int izmax = Floor(zMaxCellPos.z());

  ixmin = max(domain->getCellLowIndex.x(), ixmin);
  iymin = max(domain->getCellLowIndex.y(), iymin);
  izmin = max(domain->getCellLowIndex.z(), izmin);
  ixmax = min(ixmax, domain->getCellHighIndex.x());
  iymax = min(iymax, domain->getCellHighIndex.y());
  izmax = min(izmax, domain->getCellHighIndex.z());

  // Find the nodes within the range
  // (Assume straightforward cell-based storage of particles lists)
  for (int ii = ixmin; ii < ixmax, ii++) {
    for (int jj = iymin; jj < iymax, jj++) {
      for (int kk = izmin; kk < izmax, kk++) {
        
      }
    }
  }

}

#include <Peridynamics.h>
#include <FamilyComputer.h>

Peridynamics::Peridynamics(const GlobalFlags* flags) 
{
  d_num_broken_bonds = 0;
  d_flags = flags;
  d_damping = {{0.0, 0.0, 0.0}};

  // Hard-coded for now (dynamic relaxation?)
  if !(d_flags->isDynamic) {
    d_damping[0] = 0.2;
    d_damping[1] = 0.2;
    d_damping[2] = 0.2;
  }

  d_modified_mesh = true;
  d_use_canonical_micromodulus = true;
}

Peridynamics::~Peridynamics() 
{
}

void 
Peridynamics::updateDisplacementVelocityVerlet(const int& iteration,
		                               NodeArray& nodes)
{
  if (iteration == 1 || d_modified_mesh) {
    computeNodeFamily();
    computeBondFamily();
    computeInternalForce(nodes);
    d_modified_mesh = false;
  }

  integrateNodalAcceleration();
  breakBonds();
}

void 
Peridynamics::computeNodeFamily()
{
  FamilyComputer bfc;
  bfc.computeNodeFamily(d_node_family);
}

void
Peridynamics::computeBondFamily(const NodeArray& nodes)
{
  for (NodeIterator iter = nodes.begin(); iter != nodes.end(); iter++) {
    Node* cur_node = *iter;
    NodeArray family_nodes;
    getFamilyNodes(cur_node, family_nodes);
    for (NodeIterator family_iter = family_nodes.begin(); 
		    family_iter != family_nodes.end(); family_iter++) {
      Node* family_node = *fam_iter;
      Bond* bond = new Bond(cur_node, family_node);
      d_bond_family.insert(NodeBondPair(cur_node, bond));
    }
  }
}

void
Peridynamics::getFamilyNodes(const Node* node,
		             NodeArray& familyNodes) const
{
  std::pair<NodeFamilyIterator, NodeFamilyIterator> range = d_node_family.equal_range(node);
  for (NodeFamilyIterator iter=range.first; iter != range.second; iter++) {
    familyNodes->push_back((*iter).second); 
  }
}

void
Peridynamics::getFamilyBonds(const Node* node,
		             BondArray& familyBonds) const
{
  std::pair<BondFamilyIterator, BondFamilyIterator> range = d_bond_family.equal_range(node);
  for (BondFamilyIterator iter=range.first; iter != range.second; iter++) {
    familyBonds->push_back((*iter).second); 
  }
}

// force update at previous configuration (n)
void
Peridynamics::computeInternalForce(NodeArray& nodes)
{
  for (NodeIterator iter = nodes.begin(); iter != nodes.end(); iter++) {
    Node* cur_node = *iter;
    if (cur_node->omit()) continue;  // skip this node

    // Get the family of node mi (all the nodes within its horizon, delta).
    double delta = cur_node->getHorizonSize();
    BondArray family_bonds;
    getFamilyBonds(cur_node, family_bonds);

    // Initialize the nodal internal force
    Array3 int_force = {{0.0, 0.0, 0.0}};

    // Initialize strain energy and spsum
    double wt = 0.0;
    double spsum = 0.0;

    // Loop over bonds in the family of current node mi.
    for (BondIterator family_iter = family_bonds.begin(); 
		    family_iter != family_bonds.end(); family_iter++) {
      Bond* bond = *fam_iter;
      Node* family_node = bond->second();

      if (family_node->omit() || bond->isBroken()) continue;  // skip this node

      // Find the peridynamic interparticle force.
      Array3 bond_force = {{0.0, 0.0, 0.0}};
      double bond_length_ref = 0.0;
      double bond_length_new = 0.0;
      double bond_strain = 0.0;
      double bond_strain_energy = 0.0;
      double micro_modulus = 0.0;
      computeBondForce(bond, bond_force, bond_length_ref, bond_length_new,
		       bond_strain, bond_strain_energy, micro_modulus);
      bond->setStrain(bond_strain);

      // Find the volume contribution from the family node
      // reduce volume if it is not fully within the horizon of current node.
      double fam_volume = family_node->getVolume();
      double xi = bond_length_ref;
      Array3 fam_interval = {{0.0, 0.0, 0.0}};
      family_node->getInterval(fam_interval);
      double fam_radij = 0.5*std::max(fam_interval[0], fam_interval[1]);

      double volume_fac = 0.0;
      if (fam_radij > 0.0) {
        if (xi <= delta - fam_radij) {
          volume_fac = 1.0;
	} else if (xi <= delta + fam_radij) {
          volume_fac = (delta + fam_radij - xi)/(2.0*fam_radij);
	} else {
          volume_fac = 0.0;
	}
      } 

      // Sum up the force on node mi.
      // force at the current configuration (n+1)
      bond_force[0] *= (fam_volume*volume_fac);
      bond_force[1] *= (fam_volume*volume_fac);

      int_force[0] += bond_force[0];
      int_force[1] += bond_force[1];

      wt += 0.5*bond_strain_energy*fam_volume*volume_fac;
      spsum += fam_volume*micro_modulus/cur_node->getDensity();
    }
    cur_node->setInternalForce(int_force);
    cur_node->setStrainEnergy(wt);
    cur_node->setSpSum(spsum);
  }
}

// returns force density per unit volume due to peridynamic interaction between nodes
void
Peridynamics::computeBondForce(Bond* bond, 
		               Array3& bondForce,
			       double& bondLengthInit,
			       double& bondLengthNew,
			       double& bondStrain,
			       double& bondStrainEnergy,
			       double& micromodulus)
{
  // Get bond data
  Node* cur_node = bond->first();
  Node* family_node = bond->second();
  bool broken = bond->isBroken();

  // Get the reference position, displacement of node
  Array3 cur_pos, cur_disp;
  cur_node->getPosition(cur_pos);
  cur_node->getDisplacement(cur_disp);

  double x1_cur = cur_pos[0];
  double x2_cur = cur_pos[1];
  double u1_cur = cur_disp[0];
  double u2_cur = cur_disp[1];
		                        
  // Get the reference position and displacement of family node
  Array3 fam_pos, fam_disp;
  family_node->getPosition(fam_pos);
  family_node->getDisplacement(fam_disp);

  double x1_fam = fam_pos[0];
  double x2_fam = fam_pos[1];
  double u1_fam = fam_disp[0];
  double u2_fam = fam_disp[1];

  // Set up variables used in EMUNE
  double xi1 = x1_fam - x1_cur;
  double xi2 = x2_fam - x2_cur;
  double eta1 = u1_fam - u1_cur;
  double eta2 = u2_fam - u2_cur;

  // rr is the initial distance between the nodes
  double rr = std::sqrt(xi1*xi1 + xi2*xi2);

  // Compute updated nodal positions
  double p1 = xi1 + eta1;
  double p2 = xi2 + eta2;

  // pp is the new distance between the nodes
  double pp = std::sqrt(p1*p1 +p2*p2);
  if (pp > 1.0e16) {
    // diverges - break the bond
    bond->isBroken(true);
    return;
  }
  bondLengthInit = rr;
  bondLengthNew = pp;

  // Compute bond micromodulus
  double delta = cur_node->getHorizonSize();
  double young = cur_node->getModulus();
  double micromodulus = computeMicroModulus(rr, delta, young);

  // elasticity model.
  // find bond strain.
  // find the force
  if (pp > 0.0 && rr < delta) {

    // u is the displacement between the two nodes
    double uu = pp - rr; 
    if (rr > 0.0) {
      bondStrain  = uu/rr;
    } else {
      bondStrain = 0.0;
    }
 
    bondForce[0] = micromodulus*bondStrain*(p1/pp);
    bondForce[1] = micromodulus*bondStrain*(p2/pp);
    bondStrainEnergy = 0.5*micromodulus*(uu*uu);
  } else {
    bondForce[0] = 0.0;
    bondForce[1] = 0.0;
    bondStrainEnergy = 0.0;
    micromodulus = 0.0;
    bondStrain = 0.0;
  }

  return;
}

// Purpose : define the micromodulus function C
// Options :
// 2D constant micromodulus   
//    dmicroF =  6.0d0*young/(pi*thickness*(horizon**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> thickness effect will be vanished in volume integration
//    dmicroF =  13.5d0*young/(pi*(horizon**3)) ==> thickness = 1
// 2D canonical micromodulus
//    dmicroF = 24.0d0*young*(1.0d0-bondlength/horizon)/(pi*thickness*(horizon**3)*(1.d0/3.d0)*(4.d0/3.d0))  ==> thickness effect will be vanished in volume integration
//    dmicroF = 54.0d0*young*(1.0d0-bondlength/horizon)/(pi*(horizon**3))
double
Peridynamics::computeMicromodulus(const double& bondLengthInitial, 
		                  const double& horizonRadius,
				  const double& youngsModulus)
{
  double micromodulus = 0.0;
  double rad_cubed = horizonRadius*horizonRadius*horizonRadius;
  if (d_use_canonical_micromodulus) {
    micromodulus = 54.0*youngsModulus*(1.0-bondLengthInitial/horizonRadius)/(M_PI*rad_cubed);
  } else {
    micromodulus = 13.5*youngsModulus/(M_PI*rad_cubed);
  }
  return micromodulus;
}

// Integrates the node accelerations due to peridynamic ("structured") interaction.
// Update the node velocities
// [one-step Velocity-Verlet formulation]
// 1. v(n+1/2) = v(n) + dt/2m * f(q(n))
// 2. q(n+1) = q(n) + dt * v(n+1/2)
// 3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
void
Peridynamics::integrateNodalAcceleration(NodeArray& nodes)
{
  // 1st step & 2nd step (intermediate velocity update and position update)
  for (NodeIteration iter=nodes.begin(); iter != nodes.end(); ++iter) {
    Node* cur_node = *iter;

    Array3 vel_old = {{0.0, 0.0, 0.0}};
    Array3 disp_old = {{0.0, 0.0, 0.0}};
    cur_node->getDisplacement(disp_old);
    cur_node->getVelocity(vel_old);

    Array3 vel_new = {{0.0, 0.0, 0.0}};
    Array3 disp_new = disp_old;

    if (cur_node->omit() || cur_node->isHanging()) {
      cur_node->setVelocity(vel_new);
      cur_node->setDisplacement(disp_new);
    } else {
      // 1. v(n+1/2) = v(n) + (dt/2m) * f(u(n))
      double node_mass = cur_node->getDensity()*cur_node->getVolume();
      double delt = state->delT;
      double fac = delt/(2.0*node_mass);
      Array3 ext_force = {{0.0, 0.0, 0.0}};
      Array3 int_force = {{0.0, 0.0, 0.0}};
      cur_node->getExternalForce(ext_force);
      cur_node->getInternalForce(int_force);
      vel_new[0] = vel_old[0] + (ext_force[0]+int_force[0])*fac;
      vel_new[1] = vel_old[1] + (ext_force[1]+int_force[1])*fac;
      cur_node->setVelocity(vel_new);
      // 2. u(n+1) = u(n) + dt * v(n+1/2)
      disp_new[0] = disp_old[0] + vel_old[0]*delt;
      disp_new[1] = disp_old[1] + vel_old[1]*delt;
      cur_node->setDisplacement(disp_new);
    }
  }

  //   3rd step (force update and velocity update)
  //   3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
  // Update the internal force
  computeInternalForce(nodes);

  // Update the velocity using the updated internal force
  for (NodeIteration iter=nodes.begin(); iter != nodes.end(); ++iter) {
    Node* cur_node = *iter;
    Array3 vel_old = {{0.0, 0.0, 0.0}};
    cur_node->getVelocity(vel_old);

    Array3 vel_new = {{0.0, 0.0, 0.0}};

    if (cur_node->omit() || cur_node->isHanging()) {
      cur_node->setVelocity(vel_new);
    } else {
      //   3. v(n+1) = v(n+1/2) + dt/2m * f(q(n+1))
      double node_mass = cur_node->getDensity()*cur_node->getVolume();
      double delt = state->delT;
      double fac = delt/(2.0*node_mass);
      Array3 ext_force = {{0.0, 0.0, 0.0}};
      Array3 int_force = {{0.0, 0.0, 0.0}};
      cur_node->getExternalForce(ext_force);
      cur_node->getInternalForce(int_force);
      vel_new[0] = vel_old[0] + (ext_force[0]+int_force[0])*fac;
      vel_new[1] = vel_old[1] + (ext_force[1]+int_force[1])*fac;
      cur_node->setVelocity(vel_new);
    }
  }
}

void 
Peridynamics::breakBonds()
{
  for (NodeIterator iter = nodes.begin(); iter != nodes.end(); iter++) {
    Node* cur_node = *iter;
    if (cur_node->omit()) continue;  // skip this node

    // Get the family of node mi (all the nodes within its horizon, delta).
    double delta = cur_node->getHorizonSize();
    BondArray family_bonds;
    getFamilyBonds(cur_node, family_bonds);

    // Loop over bonds in the family of current node mi.
    int num_broken_bonds = 0;
    int num_family_bonds = 0;
    for (BondIterator family_iter = family_bonds.begin(); 
		    family_iter != family_bonds.end(); family_iter++) {
      
      num_family_bonds++;

      Bond* bond = *fam_iter;
      Node* fam_node = bond->second();

      if (fam_node->omit() || bond->isBroken()) {
	num_broken_bonds++;
	continue;  // skip this node
      }

      // Critical stretch as a function of damage.
      double damage_index_cur_node = cur_node->getDamageIndex();
      double damage_index_fam_node = fam_node->getDamageIndex();
      double dmgij = std::max(damage_index_cur_node, damage_index_family_node);

      double coef1 = 0.0, coef2 = 0.0, coef3 = 0.0;
      damageModel->getDamageCoeff(coef1, coef2, coef3);
      double damage_fac = 1.0;
      if (coef2 > 0.0 && coef3 > 1.0) {
        if (dmgij <= coef1) {
          damage_fac = 1.0;
	} else if (dmgij <0.9999) {
          damage_fac = 1.0 + coef2*(dmgij-coef1)/(1.0-dmgij);
          damage_fac = std::min(damage_fac, coef3);
        } else {
          damage_fac = coef3;
        }
      } 

      // Break bond if critical stretch exceeded.
      double critical_strain_cur = cur_node->getCriticalStrain();
      double ecr2 = critical_strain_cur*damage_fac;
      double str = bond->getStrain();
      if (str > ecr2) {
        bond->isBroken(true);
	num_broken_bonds++;
	cur_node->numBrokenBonds(num_broken_bonds);
      }
    }

    if (num_family_bonds > 0) {
      double damage_index = (double)num_broken_bonds/(double)num_family_bonds;
      cur_node->setDamageIndex(damage_index);
    } 
  }
}


  subroutine shortrange_motion
 
    ! Integrates the equation of motion for short-range forces ("structureless interaction") between nodes.
 
    ! Structureless interaction between target nodes:
    ! For each node on this processor, find the deformed family.
    ! Sum up the structureless force from each node in the deformed family.


    ! Sort the nodes on this processor according to deformed position.
    call sort_def()

    do mi=1,nnodes
      if(nodes(mi)%iflag.eqv. .TRUE.) continue

      x1i = nodes(mi)%pos(1)
      x2i = nodes(mi)%pos(2)
      y1i = nodes(mi)%pos(1)+nodes(mi)%disp(1)
      y2i = nodes(mi)%pos(2)+nodes(mi)%disp(2)
      v1i = nodes(mi)%veloc(1)
      v2i = nodes(mi)%veloc(2)          
      dvoli = nodes(mi)%volume
      nodti = node_type(mi)
      xmassi = denst*dvoli
      !<<
      ! radnodi = nodes(mi)%horizon_size/horizon_factor
      radnodi = radnod(mi)
      !>>02272009_YounDoh

      ! Find the deformed family for node m. This is all the nodes that
      ! could possibly have a structureless interaction with node m.

      rad_search = 3.d0*radnod_max
      ! rad_search = 3.0*radnodi
      ! call get_family(mi)
      call get_def_family(mi)


      ! Loop over nodes in the deformed family.
      !<<
      ! do mj=1,m_def_fam
      do js=1,m_def_fam
        mj = def_family(js)
        ! radnodj = nodes(mj)%horizon_size/horizon_factor
        radnodj = radnod(mj)
        !>>02272009_YounDoh
        x1j = nodes(mj)%pos(1)
        x2j = nodes(mj)%pos(2)
        y1j = nodes(mj)%pos(1)+nodes(mj)%disp(1)
        y2j = nodes(mj)%pos(2)+nodes(mj)%disp(2)
        v1j = nodes(mj)%veloc(1)
        v2j = nodes(mj)%veloc(2)
        dvolj = nodes(mj)%volume
        nodtj = node_type(mj)
        xmassj = denst*dvolj

        ! >>031223

        ! Set the interaction distance for this pair of nodes (mi and mj).
        ! This is the threshold distance beyond which there is no unstructured interaction.

        ! Dist_init is the distance between the 2 nodes in the reference configuration.
        ! Dist_fac_init is the coefficient that multiplies the above to get the shortrange interaction
        ! distance for nodes initially far apart. This is usually > 1.

        ! Dist_nom is the nominal contact distance based on the max of the node radii.
        ! (The max is used so small nodes cannot sneak between large nodes.)
        ! Dist_fac_nom is the coefficient that multiplies the above. This is usually < 1.
        ! This determines the shortrange interaction distance for nodes initially close together.

        dist_init = dsqrt( (x1j-x1i)**2+(x2j-x2i)**2 )
        ! dist_fac_init = max( shortrange_dist_fac_init(nodti), shortrange_dist_fac_init(nodtj) )
        dist_fac_init = 0.9d0
        dist_nom = 2.d0*dmax1(radnodi, radnodj)
        ! dist_fac_nom = max( shortrange_dist_fac_nom(nodti), shortrange_dist_fac_nom(nodtj) )
        dist_fac_nom = 1.35d0
        ! dist_fac_nom = 1.75d0

        deul = dmin1(dist_fac_nom*dist_nom, dist_fac_init*dist_init)

        !deul = 0.9*2.0*max(radnodi,radnodj)
        !deul = 1.5*deul
        !deul = min(deul, 0.9*dist_ref)

        ! Get the shortrange force factor (multiplies the peridynamic spring constant).
        ! amp = 0.5*( shortrange_force_fac(nodti) + shortrange_force_fac(nodtj) )
        ! amp = 15.0
        amp = 30.d0
        ! <<


        ! Find the current distance between this pair of nodes (mi and mj).
        distsq = (y1j-y1i)**2+(y2j-y2i)**2
        if(distsq<=deul**2.and.distsq>1.0d-6*deul**2) then

          ! Find the effective spring constant.
          ! Find the force on node mi due to unstructured interaction with node mj.
          fcof = One
          xmassij = 0.5d0*(xmassi+xmassj)
          call shortrange_force(y1i,y2i,y1j,y2j,deul,fcof,amp,v1i,v2i,v1j,v2j,xmassij,dvoli,dvolj,df1i,df2i,str,spcoef,dw)
          wt(mi) = wt(mi)+0.5d0*dw
          spsum(mi) = spsum(mi)+spcoef/denst
          nodes(mi)%new_veloc(1) = nodes(mi)%new_veloc(1) + df1i*dt/xmassi
          nodes(mi)%new_veloc(2) = nodes(mi)%new_veloc(2) + df2i*dt/xmassi


        end if
      end do
    end do

    return

  end subroutine shortrange_motion
 
 
  !-- shortrange_force
  subroutine shortrange_force(y1a,y2a,y1b,y2b,deul,fcof,amp,v1a,v2a,v1b,v2b,xmass,dvola,dvolb,df1a,df2a,str,spcoef,dw)
 
    ! returns force on point a due to structureless interaction
    ! with point b.
 
    ! p is the scalar distance from a to b.
    p1 = y1b-y1a
    p2 = y2b-y2a
    p = dsqrt(p1**2+p2**2)
 
    ! find dp/dt.
    if(p>zero) then
      pdot = ((v1b-v1a)*(y1b-y1a)+(v2b-v2a)*(y2b-y2a))/p
    else
      pdot = zero
    end if
 
    ! amp is the amplification factor relative to the
    ! structured interaction forces
    !amp = 15.

    ! find the critical damping factor.
    ccrit = 2.d0*dsqrt(amp*fcof*dvola*dvolb*xmass)
    cdamp = 0.05d0*ccrit
 
    ! tension
    if(p>deul) then
      df1a = zero
      df2a = zero
      str = zero
      spcoef = zero
      dw = zero
 
    ! compression
    else if(p>=1.0d-10) then
      felast = amp*fcof*(p-deul)*dvola*dvolb
      fdamp = cdamp*pdot
      df1a = (felast+fdamp)*p1/p
      df2a = (felast+fdamp)*p2/p
      str = p-deul
      spcoef = amp*fcof*dvola
      dw = 0.5d0*amp*fcof*(p-deul)**2*dvola*dvolb
    else
      df1a = zero
      df2a = zero
      str = zero
      spcoef = zero
      dw = zero
    end if
 
    return
 
  end subroutine shortrange_force

!****************************************************************************
end module peridynamic_computations
!****************************************************************************

import os
import math
import xml.etree.ElementTree as ET
import json
import matplotlib.cm as cm
import numpy as np
from scipy.spatial import ConvexHull
from scipy.special import comb
from shapely.geometry import LineString, Point


#-----------------------------------------------------------------------------
# Read the input xml file and save relevant data to a dictionary
#-----------------------------------------------------------------------------
def get_yield_surface_data(uda_path, PRINTOUT=False):

  # Parse xml
  try:
    ups_file = os.path.abspath(uda_path) + '/input.xml.orig'
    tree = ET.parse(ups_file)
  except:
    ups_file = os.path.abspath(uda_path) + '/input.xml'
    tree = ET.parse(ups_file)

  # Read the plasticity parameters
  material_dict = {}
  root = tree.getroot()

  for geom in root.iter("geom_object"):
    temperature = geom.find('temperature').text
    material_dict['temperature'] = float(temperature)

  for material in root.iter("material"):
    density = material.find('density').text
    specific_heat = material.find('specific_heat').text
    melt_temp = material.find('melt_temp').text
    material_dict['density'] = float(density)
    material_dict['specific_heat'] = float(specific_heat)
    material_dict['melt_temp'] = float(melt_temp)

  for model in root.iter('constitutive_model'):

    if (model.attrib['type'] == 'elastic_plastic'):

      for child in model:

        for eos_model in child.iter('equation_of_state'):
          if (eos_model.attrib['type'] == 'mie_gruneisen'):
            C0 = eos_model.find('C_0').text
            Gamma0 = eos_model.find('Gamma_0').text
            Salpha = eos_model.find('S_alpha').text
            rho0 = eos_model.find('rho_0').text
            material_dict['C0'] = float(C0)
            material_dict['Gamma0'] = float(Gamma0)
            material_dict['Salpha'] = float(Salpha)
            material_dict['rho0'] = float(rho0)

        for shear_model in child.iter('shear_modulus_model'):
          if (shear_model.attrib['type'] == 'np_shear'):
            mu0 = shear_model.find('mu_0').text
            dmup_dmu0 = shear_model.find('slope_mu_p_over_mu0').text
            zeta = shear_model.find('zeta').text
            C = shear_model.find('C').text
            m = shear_model.find('m').text
            material_dict['mu0'] = float(mu0)
            material_dict['dmup_dmu0'] = float(dmup_dmu0)
            material_dict['zeta'] = float(zeta)
            material_dict['C'] = float(C)
            material_dict['m'] = float(m)

        for flow_model in child.iter('flow_model'):
          if (flow_model.attrib['type'] == 'preston_tonks_wallace'):
            theta = flow_model.find('theta').text
            p = flow_model.find('p').text
            s0 = flow_model.find('s0').text
            sinf = flow_model.find('sinf').text
            kappa = flow_model.find('kappa').text
            gamma = flow_model.find('gamma').text
            y0 = flow_model.find('y0').text
            yinf = flow_model.find('yinf').text
            y1 = flow_model.find('y1').text
            y2 = flow_model.find('y2').text
            beta = flow_model.find('beta').text
            M = flow_model.find('M').text
            G0 = flow_model.find('G0').text
            alpha = flow_model.find('alpha').text
            alphap = flow_model.find('alphap').text
            material_dict['theta'] = float(theta)
            material_dict['p'] = float(p)
            material_dict['s0'] = float(s0)
            material_dict['sinf'] = float(sinf)
            material_dict['kappa'] = float(kappa)
            material_dict['gamma'] = float(gamma)
            material_dict['y0'] = float(y0)
            material_dict['yinf'] = float(yinf)
            material_dict['y1'] = float(y1)
            material_dict['y2'] = float(y2)
            material_dict['beta'] = float(beta)
            material_dict['M'] = float(M)
            material_dict['G0'] = float(G0)
            material_dict['alpha'] = float(alpha)
            material_dict['alphap'] = float(alphap)

        for melt_model in child.iter('melting_temp_model'):
          if (melt_model.attrib['type'] == 'scg_Tm'):
            Tm0 = melt_model.find('T_m0').text
            Gamma_0 = melt_model.find('Gamma_0').text
            a = melt_model.find('a').text
            material_dict['Tm0'] = float(Tm0)
            material_dict['Gamma_0'] = float(Gamma_0)
            material_dict['a'] = float(a)

  if PRINTOUT:
    print('--Material Specification--')
    for key in material_dict:
      print(key, ':', material_dict[key])

  material_dict = convertToBoldFont(material_dict, PRINTOUT)

  return material_dict


#-----------------------------------------------------------------------------
# Convert material dictionary data to bold font and add string to material_dict
#-----------------------------------------------------------------------------
def convertToBoldFont(material_dict, PRINTOUT=False):

  tmp_string = r'$\mathbf{\underline{Material\phantom{1}Properties:}}$' + '\n'
  key_list = list(material_dict.keys())
  key_list.sort()
  for key in key_list:
    if '_' in key:
      tmp = key.split('_')
      tmp = str_to_mathbf(tmp[0] + '_' + '{' + tmp[1] + '}')
      if (isinstance(material_dict[key], str)):
        tmp_string += tmp + str_to_mathbf(' = ') + str_to_mathbf(
            material_dict[key]) + '\n'
      else:
        tmp_string += tmp + str_to_mathbf(' = ') + str_to_mathbf(
            format(material_dict[key], '1.3e')) + '\n'
    else:
      tmp = key
      if key == 'subcycling char num':
        tmp_string += str_to_mathbf(tmp + ' = ' +
                                    format(material_dict[key], '4.1f')) + '\n'
      else:
        if (isinstance(material_dict[key], str)):
          tmp_string += str_to_mathbf(tmp + ' = ' + material_dict[key]) + '\n'
        else:
          tmp_string += str_to_mathbf(tmp + ' = ' +
                                      format(material_dict[key], '1.3e')) + '\n'

  material_dict['material string'] = tmp_string[0:-1]
  if PRINTOUT:
    print(tmp_string)

  return material_dict

#-----------------------------------------------------------------------------
# Make the text bold
# Only works with single spaces no leading space
#-----------------------------------------------------------------------------
def str_to_mathbf(string):
  string = string.split()
  return_string = ''
  for elem in string:
    #elem = r'$\mathbf{'+elem+'}$'
    elem = r'' + elem + ''
    return_string += elem + '  '
  return return_string[0:-1]


#-----------------------------------------------------------------------------
# Compute pressure
#-----------------------------------------------------------------------------
def computePressure(material_dict, density, temperature):
  C0 = material_dict['C0']
  Gamma0 = material_dict['Gamma0']
  Salpha = material_dict['Salpha']
  rho0 = material_dict['rho0']

  # Get the state data
  rho = density
  T = temperature
  T_0 = material_dict['temperature']

  # Get original density
  rho_0 = material_dict['density']

  # Get specific heat
  C_p = material_dict['specific_heat']

  # Calc. zeta
  zeta = (rho / rho_0 - 1.0)

  # Calculate internal energy E
  E = C_p * (T - T_0) * rho_0

  # Calculate the pressure
  p = Gamma0 * E
  if rho != rho_0:
    numer = rho_0 * (C0 * C0) * (1.0 / zeta + (1.0 - 0.5 * Gamma0))
    denom = 1.0 / zeta - (Salpha - 1.0)
    if denom == 0.0:
      print("rh0_0 = ", rho_0, " zeta = ", zeta, " numer = ", numer)
      denom = 1.0e-5

    p += numer / (denom * denom)

  return -p


#-----------------------------------------------------------------------------
# Compute melt temperature
#-----------------------------------------------------------------------------
def computeMeltTemp(material_dict, density):
  Tm0 = material_dict['Tm0']
  Gamma0 = material_dict['Gamma_0']
  a = material_dict['a']
  rho_0 = material_dict['density']

  eta = density / rho_0
  power = 2.0 * (Gamma0 - a - 1.0 / 3.0)
  Tm = Tm0 * np.exp(2.0 * a * (1.0 - 1.0 / eta)) * np.power(eta, power)

  #print(" SCG Melting Temp : ", Tm, " eta = ", eta)
  return Tm


#-----------------------------------------------------------------------------
# Compute shear modulus
#-----------------------------------------------------------------------------
def computeShearModulus(material_dict, temperature, density, pressure):
  mu0 = material_dict['mu0']
  dmup_dmu0 = material_dict['dmup_dmu0']
  zeta = material_dict['zeta']
  C = material_dict['C']
  m = material_dict['m']
  T_m = material_dict['melt_temp']
  rho_0 = material_dict['density']

  That = temperature / T_m
  if That <= 0:
    return mu0

  mu = 1.0e-8
  # Small value to keep the code from crashing
  if That > 1.0 + zeta:
    return mu

  j_denom = zeta * (1.0 - That / (1.0 + zeta))
  J = 1.0 + np.exp((That - 1.0) / j_denom)
  if not np.isfinite(J):
    return mu

  eta = density / rho_0
  eta = np.power(eta, 1.0 / 3.0)

  # Pressure is +ve in this calculation
  P = -pressure
  t1 = mu0 * (1.0 + dmup_dmu0 * P / eta)
  t2 = 1.0 - That
  k_amu = 1.3806503e4 / 1.6605402
  t3 = density * k_amu * temperature / (C * m)
  mu = 1.0 / J * (t1 * t2 + t3)

  if (mu < 1.0e-8):
    print("mu = ", mu, " T = ", temperature, " Tm = ", T_m, " T/Tm = ", That,
          " J = ", J, " rho/rho_0 = ", eta, " p = ", P, " t1 = ", t1, " t2 = ",
          t2, " t3 = ", t3)
    mu = 1.0e-8

  return mu


#-----------------------------------------------------------------------------
# Compute yield function
#-----------------------------------------------------------------------------
def computeYieldFunction(stress, backStress, yieldStress):

  xi = sigma_dev(stress) - sigma_dev(backStress)
  sigy = yieldStress
  xiNorm = sigma_mag(xi)
  f = np.sqrt(3.0 / 2.0) * xiNorm - sigy
  return f


#-----------------------------------------------------------------------------
# Compute flow stress
#-----------------------------------------------------------------------------
def computeFlowStress(material_dict, eqPlasticStrainRate, eqPlasticStrain,
                      temperature, shearModulus, density, meltTemp):

  theta = material_dict['theta']
  p = material_dict['p']
  s0 = material_dict['s0']
  sinf = material_dict['sinf']
  kappa = material_dict['kappa']
  gamma = material_dict['gamma']
  y0 = material_dict['y0']
  yinf = material_dict['yinf']
  y1 = material_dict['y1']
  y2 = material_dict['y2']
  beta = material_dict['beta']
  M = material_dict['M']
  G0 = material_dict['G0']

  # Retrieve plastic strain and strain rate
  epdot = eqPlasticStrainRate
  epdot = 1.0e-8 if epdot <= 0.0 else epdot

  ep = eqPlasticStrain

  # Check if temperature is correct
  T = temperature
  Tm = meltTemp

  # Check if shear modulus is correct
  mu = shearModulus

  # Get the current mass density
  rho = density

  # Convert the atomic mass to kg
  Mkg = M * 1.66043998903379e-27

  # Compute invxidot - the time required for a transverse wave to cross
  # an atom
  if (mu <= 0.0) or rho < 0.0 or (T > Tm) or (T <= 0.0):
    print("**ERROR** PTWFlow::computeFlowStress: mu = ", mu, " rho = ",
          rho << " T = ", T, " Tm = ", Tm)

  xidot = 0.5 * np.power(4.0 * np.pi * rho / (3.0 * Mkg),
                       (1.0 / 3.0)) * np.sqrt(mu / rho)

  # Compute the dimensionless plastic strain rate
  edot = epdot / xidot
  if not (xidot > 0.0) or not (edot > 0.0):
    print("**ERROR** PTWFlow::computeFlowStress: xidot = ", xidot, " edot = ",
          edot)

  # Compute the dimensionless temperature
  That = T / Tm

  # Calculate the dimensionless Arrhenius factor
  arrhen = kappa * That * np.log(gamma / edot)

  # Calculate the saturation hardening flow stress in the thermally
  # activated glide regime
  tauhat_s = s0 - (s0 - sinf) * math.erf(arrhen)

  # Calculate the yield stress in the thermally activated glide regime
  tauhat_y = y0 - (y0 - yinf) * math.erf(arrhen)

  # The overdriven shock regime
  if (epdot > 1.0e3):

    # Calculate the saturation hardening flow stress in the overdriven
    # shock regime
    shock_tauhat_s = s0 * np.power(edot / gamma, beta)

    # Calculate the yield stress in the overdriven shock regime
    shock_tauhat_y_jump = y1 * np.power(edot / gamma, y2)
    shock_tauhat_y = np.min(shock_tauhat_y_jump, shock_tauhat_s)

    # Calculate the saturation stress and yield stress
    tauhat_s = np.max(tauhat_s, shock_tauhat_s)
    tauhat_y = np.max(tauhat_y, shock_tauhat_y)

  # Compute the dimensionless flow stress
  tauhat = tauhat_s
  if (tauhat_s != tauhat_y):
    A = (s0 - tauhat_y) / p
    B = tauhat_s - tauhat_y
    D = np.exp(B / A)
    C = D - 1.0
    F = C / D
    E = theta / (A * C)
    exp_EEp = np.exp(-E * ep)
    tauhat = tauhat_s + A * np.log(1.0 - F * exp_EEp)

  sigma = 2.0 * tauhat * mu
  return sigma


#-----------------------------------------------------------------------------
# Compute points on the yield surface
#-----------------------------------------------------------------------------
def computeYieldSurfacePoints(material_dict, pbar_sim_list, epdot_sim, ep_sim,
                              T_sim, rho_sim):

  num_pts = 100
  pmin = min(pbar_sim_list)
  pmax = max(pbar_sim_list)
  pbar_yield = np.linspace(pmin, pmax, num=num_pts)

  pbar = computePressure(material_dict, rho_sim, T_sim)
  Tm = computeMeltTemp(material_dict, rho_sim)
  mu = computeShearModulus(material_dict, T_sim, rho_sim, pbar)
  sigy = computeFlowStress(material_dict, epdot_sim, ep_sim, T_sim, mu, rho_sim,
                           Tm)
  q_yield = [sigy] * num_pts

  return pbar_yield, q_yield


#-----------------------------------------------------------------------------
# Plot the yield surface (compression positive)
#-----------------------------------------------------------------------------
def plotPQYieldSurfaceSim(plt,
                          material_dict,
                          pbar_sim_list,
                          epdot,
                          ep,
                          T,
                          rho,
                          compression='negative',
                          plt_color='b'):

  # Compute the yield surface
  pbar_yield, q_yield = computeYieldSurfacePoints(material_dict, pbar_sim_list,
                                                  epdot, ep, T, rho)

  # Plot
  #print('Compression = ', compression)
  if (compression == 'negative'):
    ps = list(map(lambda pbar: -pbar, pbar_yield))
    qs = list(map(lambda q: q, q_yield))
    #print("yield surface: pmin = ", pmin, "pmax = ", pmax, "p = ", ps)
    #print("yield surface: qmax = ", qmax, "q = ", qs)

    line1 = plt.plot(ps, qs, '-b', linewidth=1)
    line2 = plt.plot(ps, list(map(lambda q: -q, qs)), '-b', linewidth=1)
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size': 8})
    #plt.axis('equal')

    axes = plt.gca()
    #axes.set_xlim([1.3*pmin, 1.3*pmax])
    #axes.set_ylim([-1.3*qmax, 1.3*qmax])
  else:
    ps = list(map(lambda pbar: pbar, pbar_yield))
    qs = list(map(lambda q: q, q_yield))
    #print("yield surface: pmin = ", pmin, "pmax = ", pmax, "p = ", ps)
    #print("yield surface: qmax = ", qmax, "q = ", qs)
    line1 = plt.plot(ps, qs, '-b', linewidth=1)
    line2 = plt.plot(ps, list(map(lambda q: -q, qs)), '-b', linewidth=1)
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size': 8})
    plt.axis('equal')

    axes = plt.gca()
    #axes.set_xlim([1.3*pmin, 1.3*pmax])
    #axes.set_ylim([-1.3*qmax, 1.3*qmax])

  return min(pbar_yield), max(q_yield)

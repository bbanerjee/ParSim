/*
 * The MIT License
 *
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_SupportVector.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>

#define DEBUG_SVR

using namespace Vaango;

/** Read inputs and JSON */
ElasticModuli_SupportVector::ElasticModuli_SupportVector(Uintah::ProblemSpecP& ps)
{
  ps->require("filename", d_bulk.filename);
  ps->require("G0", d_shear.G0);
  ps->require("nu", d_shear.nu);

  readJSONFile(d_bulk.filename);
  checkInputParameters();
}

/** Check that the input parameters are reasonable */
void
ElasticModuli_SupportVector::checkInputParameters()
{
  std::ostringstream warn;
  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

/** Construct a copy of the model */
ElasticModuli_SupportVector::ElasticModuli_SupportVector(const ElasticModuli_SupportVector* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

/** Write out info for restart files */
void
ElasticModuli_SupportVector::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "support_vector");

  elasticModuli_ps->appendElement("filename", d_bulk.filename);
  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("nu", d_shear.nu);
}

/** Compute initial elastic moduli */
ElasticModuli
ElasticModuli_SupportVector::getInitialElasticModuli() const
{
  double K = computeBulkModulus(1.0e-6, 0);
  double G = computeShearModulus(K);
  return ElasticModuli(K, G);
}

/** Compute current elastic moduli */
ElasticModuli
ElasticModuli_SupportVector::getCurrentElasticModuli(const ModelStateBase* state_input) const
{
  // This model requires the extra variables in ModelState_Tabular
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Make sure the quantities are positive in compression
  double ev_e_bar = -(state->elasticStrainTensor).Trace();
  double ev_p_bar = -(state->plasticStrainTensor).Trace();
  ev_p_bar = (ev_p_bar < 0) ? 0 : ev_p_bar;
  ev_e_bar = (ev_p_bar < 0) ? ev_e_bar + ev_p_bar : ev_e_bar;

  // Compute the elastic moduli
  double K = computeBulkModulus(ev_e_bar, ev_p_bar);
  double G = computeShearModulus(K);

  #ifdef DEBUG_SVR
    if (K < 0 || !std::isfinite(K)) {
      std::cout << "ev_e = " << ev_e_bar << " ev_p = " << ev_p_bar << "\n";
      std::cout << " K = " << K << " G = " << G << std::endl;
    }
  #endif

  return ElasticModuli(K, G);
}

/** Compute the bulk modulus */
double 
ElasticModuli_SupportVector::computeBulkModulus(const double& elasticVolStrain,
                                                const double& plasticVolStrain) const
{
  Eigen::Vector2d strain;
  strain << elasticVolStrain + plasticVolStrain, plasticVolStrain;
  strain = strain.array() / d_bulk.strain_conversion_factor.array();

  auto num_support_vec = d_bulk.svr_support_vectors.rows();

  Eigen::Vector2d strain_scaled = (strain - d_bulk.strain_min).array() * d_bulk.strain_scale_factor.array();
  auto strain_scaled_mat = (strain_scaled.transpose()).replicate(num_support_vec, 1);
  MatrixN2d strain_svr_diff = d_bulk.svr_support_vectors - strain_scaled_mat;
  MatrixN1d strain_svr_diff_normSq = strain_svr_diff.rowwise().squaredNorm() * (-d_bulk.rbf_kernel_gamma);
  MatrixN1d svr_kernel = strain_svr_diff_normSq.array().exp();
  MatrixN1d strain_svr_diff_gamma = strain_svr_diff.col(0) * (2.0 * d_bulk.rbf_kernel_gamma);
  MatrixN1d svr_kernel_deriv_totalstrain = svr_kernel.array() * strain_svr_diff_gamma.array();

  double K = d_bulk.svr_dual_coeffs.transpose() * svr_kernel_deriv_totalstrain;
  //std::cout << "strain = " << strain.transpose() 
  //          << " strain_scaled = " << strain_scaled.transpose() 
  //          << " K_scaled = " << K << "\n";
  K /= d_bulk.pressure_scale_factor;
  K *= (d_bulk.strain_scale_factor[0] * d_bulk.pressure_conversion_factor / d_bulk.strain_conversion_factor[0]);

  #ifdef DEBUG_SVR
    if (K < 0 || !std::isfinite(K)) {
      std::cout <<std::setprecision(16) 
                << "ee_v = " << elasticVolStrain
                << " ep_v = " << plasticVolStrain
                << " K = " << K << "\n";
    }
  #endif

  return K;
}

double 
ElasticModuli_SupportVector::computeShearModulus(const double& K) const
{
  double nu = d_shear.nu;
  double G = (nu > -1.0 && nu < 0.5) 
             ? 1.5*K*(1.0 - 2.0*nu)/(1.0 + nu) 
             : d_shear.G0;
  return G;
}

/* Get the elastic moduli and their derivatives with respect to a single
   plastic internal variable */
std::pair<ElasticModuli, ElasticModuli>
ElasticModuli_SupportVector::getElasticModuliAndDerivatives(const ModelStateBase* state_input) const
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  //if (!state) {
  //  std::ostringstream out;
  //  out << "**ERROR** The correct ModelState object has not been passed."
  //      << " Need ModelState_Tabular.";
  //  throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  //}

  // Make sure the quantities are positive in compression
  double ev_e_bar = -(state->elasticStrainTensor).Trace();
  double ev_p_bar = -(state->plasticStrainTensor).Trace();
  ev_p_bar = (ev_p_bar < 0) ? 0 : ev_p_bar;
  ev_e_bar = (ev_p_bar < 0) ? ev_e_bar + ev_p_bar : ev_e_bar;

  Eigen::Vector2d strain;
  strain << ev_e_bar+ev_p_bar, ev_p_bar;
  strain = strain.array() / d_bulk.strain_conversion_factor.array();

  auto num_support_vec = d_bulk.svr_support_vectors.rows();
  Eigen::Vector2d strain_scaled = (strain - d_bulk.strain_min).array() * d_bulk.strain_scale_factor.array();
  auto strain_scaled_mat = (strain_scaled.transpose()).replicate(num_support_vec, 1);
  MatrixN2d strain_svr_diff = d_bulk.svr_support_vectors - strain_scaled_mat;
  MatrixN1d strain_svr_diff_normSq = strain_svr_diff.rowwise().squaredNorm() * (-d_bulk.rbf_kernel_gamma);
  MatrixN1d svr_kernel = strain_svr_diff_normSq.array().exp();
  MatrixN1d strain_svr_diff_gamma = strain_svr_diff.col(0) * (2.0 * d_bulk.rbf_kernel_gamma);
  MatrixN1d svr_kernel_deriv_totalstrain = svr_kernel.array() * strain_svr_diff_gamma.array();

  MatrixN1d strain_svr_diffsq = strain_svr_diff.col(0).array() * strain_svr_diff.col(1).array();
  MatrixN1d strain_svr_diffsq_gamma = strain_svr_diffsq * (4.0 * d_bulk.rbf_kernel_gamma * d_bulk.rbf_kernel_gamma);
  MatrixN1d svr_kernel_dderiv_plasticstrain = svr_kernel.array() * strain_svr_diffsq_gamma.array();

  double K = d_bulk.svr_dual_coeffs.transpose() * svr_kernel_deriv_totalstrain;
  double dK_deps_p = d_bulk.svr_dual_coeffs.transpose() * svr_kernel_dderiv_plasticstrain;
  K /= d_bulk.pressure_scale_factor;
  dK_deps_p /= d_bulk.pressure_scale_factor;
  K *= (d_bulk.strain_scale_factor[0] * d_bulk.pressure_conversion_factor / 
        d_bulk.strain_conversion_factor[0]);
  dK_deps_p *= (d_bulk.strain_scale_factor[0] * d_bulk.strain_scale_factor[1] * 
    d_bulk.pressure_conversion_factor / 
    (d_bulk.strain_conversion_factor[0] * d_bulk.strain_conversion_factor[1]));

  /*
  std::cout << "num_support_vec = " << num_support_vec << "\n";
  std::cout << "strain_scaled = " << strain_scaled.transpose() << "\n";
  std::cout << "strain_scaled_mat = " << strain_scaled_mat << "\n";
  std::cout << "strain_svr_diff = " << strain_svr_diff << "\n";
  std::cout << "strain_svr_diff_normsq = " << strain_svr_diff_normSq.transpose() << "\n";
  std::cout << "svr_kernel = " << svr_kernel.transpose() << "\n";
  std::cout << "strain_svr_diff_gamma = " << strain_svr_diff_gamma.transpose() << "\n";
  std::cout << "svr_kernel_deriv_totalstrain = " << svr_kernel_deriv_totalstrain.transpose() << "\n";
  std::cout << "svr_support_vectors = " << d_bulk.svr_support_vectors << "\n";
  std::cout << "svr_dual_coeffs = " << d_bulk.svr_dual_coeffs.transpose() << "\n";

  std::cout << "strain_svr_diffsq = " << strain_svr_diffsq << "\n";
  std::cout << "strain_svr_diffsq_gamma = " << strain_svr_diffsq_gamma << "\n";
  std::cout << "svr_kernel_dderiv_plasticstrain = " << svr_kernel_dderiv_plasticstrain.transpose() << "\n";
  std::cout << "K = " << K << " dK_eps_p = " << dK_deps_p << "\n";
  */

  double nu = d_shear.nu;
  double G = (nu > -1.0 && nu < 0.5) 
             ? 1.5*K*(1.0 - 2.0*nu)/(1.0 + nu) 
             : d_shear.G0;
  double dG_deps_p = (nu > -1.0 && nu < 0.5) 
             ? 1.5*dK_deps_p*(1.0 - 2.0*nu)/(1.0 + nu) 
             : 0.0;

  return std::make_pair(ElasticModuli(K, G),
                        ElasticModuli(dK_deps_p, dG_deps_p));
}

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_SupportVector::computeDModuliDIntVar(const ModelStateBase* state) const
{
  auto data = computeModuliAndDModuliDIntVar(state);
  return data.second;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_SupportVector::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  auto data = getElasticModuliAndDerivatives(state);
  std::vector<ElasticModuli> derivs;
  derivs.push_back(data.second);
  return std::make_pair(data.first, derivs);
}

void
ElasticModuli_SupportVector::readJSONFile(const std::string& filename)
{
  std::ifstream inputFile(filename);
  if (!inputFile) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot read support vector regression JSON file " << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::stringstream inputStream;
  inputStream << inputFile.rdbuf();
  inputFile.close();

  nlohmann::json doc = loadJSON(inputStream, filename);
  
  d_bulk.strain_conversion_factor = getVector2dJSON(doc, "X_conversion_factor", filename); 
  d_bulk.strain_scale_factor = getVector2dJSON(doc, "X_scale", filename); 
  d_bulk.strain_min = getVector2dJSON(doc, "X_min", filename); 
  d_bulk.strain_max = getVector2dJSON(doc, "X_max", filename); 
  d_bulk.pressure_conversion_factor = getDoubleJSON(doc, "y_conversion_factor", filename); 
  d_bulk.pressure_scale_factor = getVector1dJSON(doc, "y_scale", filename); 
  d_bulk.pressure_min = getVector1dJSON(doc, "y_min", filename); 
  d_bulk.pressure_max = getVector1dJSON(doc, "y_max", filename); 
  d_bulk.rbf_kernel_gamma = getDoubleJSON(doc, "gamma", filename);
  d_bulk.svr_support_vectors = getMatrixN2dJSON(doc, "support_vectors", filename);
  d_bulk.svr_dual_coeffs = getMatrixN1dJSON(doc, "dual_coeffs", filename);
  d_bulk.svr_intercept = getVector1dJSON(doc, "intercept", filename); 
}

nlohmann::json 
ElasticModuli_SupportVector::loadJSON(std::stringstream& inputStream,
                                      const std::string& filename) const
{
  nlohmann::json doc;
  try {
    inputStream >> doc;
  } catch (std::invalid_argument err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot parse support vector regression JSON file " << filename << "\n"
        << " Please check that the file is valid JSON using a linter";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return doc;
}

double 
ElasticModuli_SupportVector::getDoubleJSON(const nlohmann::json& object, 
                                           const std::string key,
                                           const std::string& filename) const
{
  // First check if the key exists
  nlohmann::json data;
  try {
    data = object.at(key);
    //std::cout << "Data = " << data << std::endl;
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" not found in support vector regression JSON file "
        << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  double val;
  try {
    val = data.get<double>();
    //std::cout << key << ":";
    //std::cout << val << " " << "\n";
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains data in the wrong format in JSON file "
        << filename << "\n"
        << " Data expected to be of the form `key : val'.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return val;
}

double 
ElasticModuli_SupportVector::getVector1dJSON(const nlohmann::json& object, 
                                             const std::string key,
                                             const std::string& filename) const
{
  // First check if the key exists
  nlohmann::json data;
  try {
    data = object.at(key);
    //std::cout << "Data = " << data << std::endl;
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" not found in support vector regression JSON file "
        << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::vector<double> vec;
  try {
    vec = data.get<std::vector<double>>();
    //std::cout << key << ":";
    //for (const auto& val : vec) {
    //  std::cout << val << " ";
    //}
    //std::cout << std::endl;
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains data in the wrong format in JSON file "
        << filename << "\n"
        << " Data expected to be of the form `key : [val]'.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return vec[0];
}

Eigen::Vector2d 
ElasticModuli_SupportVector::getVector2dJSON(const nlohmann::json& object, 
                                             const std::string key,
                                             const std::string& filename) const
{
  // First check if the key exists
  nlohmann::json data;
  try {
    data = object.at(key);
    //std::cout << "Data = " << data << std::endl;
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" not found in support vector regression JSON file "
        << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::vector<double> vec;
  try {
    vec = data.get<std::vector<double>>();
    //std::cout << key << ":";
    //for (const auto& val : vec) {
    //  std::cout << val << " ";
    //}
    //std::cout << std::endl;
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains data in the wrong format in JSON file "
        << filename << "\n"
        << " Data expected to be of the form 'key : [x1, x2]'.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return Eigen::Vector2d(vec.data());
}

Eigen::Matrix<double, Eigen::Dynamic, 1>
ElasticModuli_SupportVector::getMatrixN1dJSON(const nlohmann::json& object,
                                              const std::string key,
                                              const std::string& filename) const
{
  // First check if the key exists
  const auto matrix = object.find(key);
  if (matrix == object.end()) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Matrix \"" << key << "\" not found in support vector regression JSON file "
        << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::vector<double> vec;
  try {
    for (const auto& val : *matrix) {
      vec = val.get<std::vector<double>>();
    }
    /*
    std::cout << key << ":";
    for (const auto& val : vec) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
    */
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains data in the wrong format in JSON file " 
        << filename << "\n"
        << " Data expected to be of the form 'key : [x1, x2, ..., x_n]^T'.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> mapped_vec(vec.data(), vec.size(), 1);
  return mapped_vec;
}

Eigen::Matrix<double, Eigen::Dynamic, 2>
ElasticModuli_SupportVector::getMatrixN2dJSON(const nlohmann::json& object,
                                              const std::string key,
                                              const std::string& filename) const
{
  // First check if the key exists
  const auto matrix = object.find(key);
  if (matrix == object.end()) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Matrix \"" << key << "\" not found in support vector regression JSON file "
        << filename;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::vector<double> col1, col2;
  try {
    //std::cout << key << ":";
    for (const auto& val : *matrix) {
      auto vec = val.get<std::vector<double>>();
      col1.push_back(vec[0]);
      col2.push_back(vec[1]);
    }
    /*
    for (const auto& val : col1) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
    for (const auto& val : col2) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
    */
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains data in the wrong format in JSON file "
        << filename << "\n"
        << " Data expected to be of the form 'key : [[x1, x2], [x3, x4], ..., [x_p, x_q]]^T'.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::move(col2.begin(), col2.end(), std::back_inserter(col1));

  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 2>> mapped_vec(col1.data(), col2.size(), 2);
  return mapped_vec;
}

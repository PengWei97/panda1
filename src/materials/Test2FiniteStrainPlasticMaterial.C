//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Test2FiniteStrainPlasticMaterial.h"
#include "libmesh/utility.h"
#include "UserObject.h"

registerMooseObject("TensorMechanicsApp", Test2FiniteStrainPlasticMaterial);

InputParameters
Test2FiniteStrainPlasticMaterial::validParams()
{
  InputParameters params = ComputeStressBase::validParams();

  params.addRequiredParam<std::vector<Real>>(
      "yield_stress",
      "Input data as pairs of equivalent plastic strain and yield stress: Should "
      "start with equivalent plastic strain 0");
  params.addParam<Real>("rtol", 1e-8, "Plastic strain NR tolerance");
  params.addParam<Real>("ftol", 1e-4, "Consistency condition NR tolerance");
  params.addParam<Real>("eptol", 1e-7, "Equivalent plastic strain NR tolerance");
  params.addClassDescription("Associative J2 plasticity with isotropic hardening.");
  params.addRequiredParam<UserObjectName>(
    "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

Test2FiniteStrainPlasticMaterial::Test2FiniteStrainPlasticMaterial(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _grain_boundary(declareProperty<Real>(_base_name + "grain_boundary")),
    _yield_stress_vector(getParam<std::vector<Real>>("yield_stress")), // Read from input file
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "elastic_strain")),
    _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "plastic_strain")),
    _eqv_plastic_strain(declareProperty<Real>(_base_name + "eqv_plastic_strain")),
    _eqv_plastic_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _eqven_plasticity_strain(declareProperty<Real>(_base_name + "eqven_plasticity_strain")),
    // _eqv_plasticity_strain(declareProperty<std::vector<Real>>(_base_name + "eqv_plasticity_strain")),
    _eqv_plastic_strain_op(declareProperty<RealVectorValue>(_base_name + "eqv_plastic_strain_op")),
    // _eqv_plasticity_strain_old(getMaterialPropertyOld<Real>(_base_name + "eqv_plastic_strain")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _strain_increment(getMaterialProperty<RankTwoTensor>(_base_name + "strain_increment")),
    _total_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "total_strain")),
    _rotation_increment(getMaterialProperty<RankTwoTensor>(_base_name + "rotation_increment")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _hard_factor(declareProperty<Real>(_base_name + "hard_factor")),
    _grain_tracker(getUserObject<GrainDataTracker<RankFourTensor>>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _rtol(getParam<Real>("rtol")),
    _ftol(getParam<Real>("ftol")),
    _eptol(getParam<Real>("eptol")),
    _deltaOuter(RankTwoTensor::Identity().outerProduct(RankTwoTensor::Identity())),
    _deltaMixed(RankTwoTensor::Identity().mixedProductIkJl(RankTwoTensor::Identity()))
{
}

void
Test2FiniteStrainPlasticMaterial::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _elastic_strain[_qp].zero();
  _plastic_strain[_qp].zero();
  _eqv_plastic_strain[_qp] = 0.0;
  
  // _von_mises_stress[_qp] = 0; 
  _hard_factor[_qp] = 0.0;
  _grain_boundary[_qp] = 0.0;
  _eqven_plasticity_strain[_qp] = 0.0;
}

// void
// Test2FiniteStrainPlasticMaterial::initJ2StatefulProperties()
// {
//   ComputeStressBase::initQpStatefulProperties();
//   _elastic_strain[_qp].zero();
//   _plastic_strain[_qp].zero();
//   _eqv_plastic_strain[_qp] = 0.0;
//   // _eqv_plasticity_strain[_qp] = 0.0;
//   // _von_mises_stress[_qp] = 0; 
// }

void
Test2FiniteStrainPlasticMaterial::computeQpStress()
{
  // perform the return-mapping algorithm
  returnMap(_stress_old[_qp],
            _grain_boundary[_qp],
            _eqv_plastic_strain_old[_qp],
            // _eqv_plasticity_strain_old[_qp],
            _plastic_strain_old[_qp],
            _strain_increment[_qp],
            _total_strain_old[_qp],
            _elasticity_tensor[_qp],
            _stress[_qp],
            _eqv_plastic_strain[_qp],
            _eqven_plasticity_strain[_qp],
            _hard_factor[_qp],
            _plastic_strain[_qp]);

  // Rotate the stress tensor to the current configuration
   

  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

  // _eqv_plasticity_strain[_qp] = computeIntegral();

  // Rotate plastic strain tensor to the current configuration
  _plastic_strain[_qp] =
      _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();

  // Calculate the elastic strain_increment
  _elastic_strain[_qp] = _mechanical_strain[_qp] - _plastic_strain[_qp];

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  // _von_mises_stress[_qp] = std::sqrt(3*_stress[_qp].secondInvariant());  
}

/**
 * Implements the return-map algorithm via a Newton-Raphson process.
 * This idea is fully explained in Simo and Hughes "Computational
 * Inelasticity" Springer 1997, for instance, as well as many other
 * books on plasticity.
 * The basic idea is as follows.
 * Given: sig_old - the stress at the start of the "time step"
 *        plastic_strain_old - the plastic strain at the start of the "time step"
 *        eqvpstrain_old - equivalent plastic strain at the start of the "time step"
 *                         (In general we would be given some number of internal
 *                         parameters that the yield function depends upon.)
 *        delta_d - the prescribed strain increment for this "time step"
 * we want to determine the following parameters at the end of this "time step":
 *        sig - the stress
 *        plastic_strain - the plastic strain
 *        eqvpstrain - the equivalent plastic strain (again, in general, we would
 *                     have an arbitrary number of internal parameters).
 *
 * To determine these parameters, introduce
 *    the "yield function", f
 *    the "consistency parameter", flow_incr
 *    the "flow potential", flow_dirn_ij
 *    the "internal potential", internalPotential (in general there are as many internalPotential
 *        functions as there are internal parameters).
 * All three of f, flow_dirn_ij, and internalPotential, are functions of
 * sig and eqvpstrain.
 * To find sig, plastic_strain and eqvpstrain, we need to solve the following
 *   resid_ij = 0
 *   f = 0
 *   rep = 0
 * This is done by using Newton-Raphson.
 * There are 8 equations here: six from resid_ij=0 (more generally there are nine
 * but in this case resid is symmetric); one from f=0; one from rep=0 (more generally, for N
 * internal parameters there are N of these equations).
 *
 * resid_ij = flow_incr*flow_dirn_ij - (plastic_strain - plastic_strain_old)_ij
 *          = flow_incr*flow_dirn_ij - (E^{-1}(trial_stress - sig))_ij
 * Here trial_stress = E*(strain - plastic_strain_old)
 *      sig = E*(strain - plastic_strain)
 * Note: flow_dirn_ij is evaluated at sig and eqvpstrain (not the old values).
 *
 * f is the yield function, evaluated at sig and eqvpstrain
 *
 * rep = -flow_incr*internalPotential - (eqvpstrain - eqvpstrain_old)
 * Here internalPotential are evaluated at sig and eqvpstrain (not the old values).
 *
 * The Newton-Raphson procedure has sig, flow_incr, and eqvpstrain as its
 * variables.  Therefore we need the derivatives of resid_ij, f, and rep
 * with respect to these parameters
 *
 * In this associative J2 with isotropic hardening, things are a little more specialised.
 * (1) f = sqrt(3*s_ij*s_ij/2) - K(eqvpstrain)    (this is called "isotropic hardening")
 * (2) associativity means that flow_dirn_ij = df/d(sig_ij) = s_ij*sqrt(3/2/(s_ij*s_ij)), and
 *     this means that flow_dirn_ij*flow_dirn_ij = 3/2, so when resid_ij=0, we get
 *     (plastic_strain_dot)_ij*(plastic_strain_dot)_ij = (3/2)*flow_incr^2, where
 *     plastic_strain_dot = plastic_strain - plastic_strain_old
 * (3) The definition of equivalent plastic strain is through
 *     eqvpstrain_dot = sqrt(2*plastic_strain_dot_ij*plastic_strain_dot_ij/3), so
 *     using (2), we obtain eqvpstrain_dot = flow_incr, and this yields
 *     internalPotential = -1 in the "rep" equation.
 */
void
Test2FiniteStrainPlasticMaterial::returnMap(const RankTwoTensor & sig_old,
                                       Real & grain_boundary, 
                                       const Real eqvpstrain_old,
                                      //  const Real eqvpTeststrain_old,
                                       const RankTwoTensor & plastic_strain_old,                          const RankTwoTensor & delta_d,
                                       const RankTwoTensor & epsilon_total_old,
                                       const RankFourTensor & E_ijkl,
                                       RankTwoTensor & sig,
                                       Real & eqvpstrain,
                                       Real & eqven_pstrain,
                                       Real & hard_factor,
                                       RankTwoTensor & plastic_strain)
{
  // the yield function, must be non-positive
  // Newton-Raphson sets this to zero if trial stress enters inadmissible region
  Real f;

  // the consistency parameter, must be non-negative
  // change in plastic strain in this timestep = flow_incr*flow_potential
  Real flow_incr = 0.0;

  // direction of flow defined by the potential
  RankTwoTensor flow_dirn;

  // Newton-Raphson sets this zero
  // resid_ij = flow_incr*flow_dirn_ij - (plastic_strain - plastic_strain_old)
  RankTwoTensor resid;

  // Newton-Raphson sets this zero
  // rep = -flow_incr*internalPotential - (eqvpstrain - eqvpstrain_old)
  Real rep;

  // change in the stress (sig) in a Newton-Raphson iteration
  RankTwoTensor ddsig;

  // change in the consistency parameter in a Newton-Raphson iteration
  Real dflow_incr = 0.0;

  // change in equivalent plastic strain in one Newton-Raphson iteration
  Real deqvpstrain = 0.0;

  // convenience variable that holds the change in plastic strain incurred during the return
  // delta_dp = plastic_strain - plastic_strain_old
  // delta_dp = E^{-1}*(trial_stress - sig), where trial_stress = E*(strain - plastic_strain_old)
  RankTwoTensor delta_dp;

  // d(yieldFunction)/d(stress)
  RankTwoTensor df_dsig;

  // d(resid_ij)/d(sigma_kl)
  RankFourTensor dr_dsig;

  // dr_dsig_inv_ijkl*dr_dsig_klmn = 0.5*(de_ij de_jn + de_ij + de_jm), where de_ij = 1 if i=j, but
  // zero otherwise
  RankFourTensor dr_dsig_inv;

  // d(yieldFunction)/d(eqvpstrain)
  Real fq;

  // yield stress at the start of this "time step" (ie, evaluated with
  // eqvpstrain_old).  It is held fixed during the Newton-Raphson return,
  // even if eqvpstrain != eqvpstrain_old.
  Real yield_stress;

  // measures of whether the Newton-Raphson process has converged
  Real err1, err2, err3;

  // number of Newton-Raphson iterations performed
  unsigned int iter = 0;

  // maximum number of Newton-Raphson iterations allowed
  unsigned int maxiter = 100;

  // plastic loading occurs if yieldFunction > toly
  Real toly = 1.0e-8;

  // Assume this strain increment does not induce any plasticity
  // This is the elastic-predictor
  // sig = sig_old + E_ijkl * delta_d; // the trial stress \sigma = \sigmna_{n} + E_{ijkl}*\delta_varsip

  sig = E_ijkl * (epsilon_total_old + delta_d - plastic_strain_old);
  eqvpstrain = eqvpstrain_old;
  // eqvpTeststrain = eqvpstrain; //eqvpTeststrain_old;

  

  plastic_strain = plastic_strain_old;

  computeHardFactor(eqvpstrain,plastic_strain);
  // eqven_pstrain = eqvpstrain;
  // computeIntegral(hard_factor);

  // eqvpstrain = eqven_pstrain;

  // if(grain_boundary >=0.45 && grain_boundary <=0.6 && )
  // {
  //   // eqvpTeststrain = 0.0;
  //   eqvpstrain = 0.0;
  //   // plastic_strain.zero();
  // }
  // if(grain_boundary <= 0.6 && grain_boundary <=)
  // {
  //   // Real 
  //   eqvpstrain = 0.0;
  //   plastic_strain.zero();
  // }
  // eqvpstrain = if(_grain_boundary[_qp]<=0.6,0,1)*eqvpstrain;

  yield_stress = getYieldStress(eqvpstrain, hard_factor); // yield stress at this equivalent plastic strain
  if (yieldFunction(sig, yield_stress) > toly)
  {
    // the sig just calculated is inadmissable.  We must return to the yield surface.
    // This is done iteratively, using a Newton-Raphson process.

    delta_dp.zero();

    // sig = E_ijkl * (epsilon_total+delta_d-plastic_strain_old);
    sig = E_ijkl * (epsilon_total_old + delta_d - plastic_strain);
    // sig = sig_old + E_ijkl * delta_d; // this is the elastic predictor

    flow_dirn = flowPotential(sig);

    resid = flow_dirn * flow_incr - delta_dp; // Residual 1 - refer Hughes Simo
    f = yieldFunction(sig, yield_stress);
    rep = -eqvpstrain + eqvpstrain_old - flow_incr * internalPotential(); // Residual 3 rep=0

    err1 = resid.L2norm();
    err2 = std::abs(f);
    err3 = std::abs(rep);

    while ((err1 > _rtol || err2 > _ftol || err3 > _eptol) &&
           iter < maxiter) // Stress update iteration (hardness fixed)
    {
      iter++;

      df_dsig = dyieldFunction_dstress(sig);
      getJac(sig, E_ijkl, flow_incr, dr_dsig);   // gets dr_dsig = d(resid_ij)/d(sig_kl)
      fq = dyieldFunction_dinternal(eqvpstrain,hard_factor); // d(f)/d(eqvpstrain)

      /**
       * The linear system is
       *   ( dr_dsig  flow_dirn  0  )( ddsig       )   ( - resid )
       *   ( df_dsig     0       fq )( dflow_incr  ) = ( - f     )
       *   (   0         1       -1 )( deqvpstrain )   ( - rep   )
       * The zeroes are: d(resid_ij)/d(eqvpstrain) = flow_dirn*d(df/d(sig_ij))/d(eqvpstrain) = 0
       * and df/d(flow_dirn) = 0  (this is always true, even for general hardening and
       * non-associative)
       * and d(rep)/d(sig_ij) = -flow_incr*d(internalPotential)/d(sig_ij) = 0
       */

      dr_dsig_inv = dr_dsig.invSymm();

      /**
       * Because of the zeroes and ones, the linear system is not impossible to
       * solve by hand.
       * NOTE: andy believes there was originally a sign-error in the next line.  The
       *       next line is unchanged, however andy's definition of fq is negative of
       *       the original definition of fq.  andy can't see any difference in any tests!
       */
      dflow_incr = (f - df_dsig.doubleContraction(dr_dsig_inv * resid) + fq * rep) /
                   (df_dsig.doubleContraction(dr_dsig_inv * flow_dirn) - fq);
      ddsig =
          dr_dsig_inv *
          (-resid -
           flow_dirn * dflow_incr); // from solving the top row of linear system, given dflow_incr
      deqvpstrain =
          rep + dflow_incr; // from solving the bottom row of linear system, given dflow_incr

      // update the variables
      flow_incr += dflow_incr;
      delta_dp -= E_ijkl.invSymm() * ddsig;
      sig += ddsig;
      eqvpstrain += deqvpstrain;
      // eqvpTeststrain += deqvpstrain;

      // if(grain_boundary >=0.45 && grain_boundary <=0.60)
      // {
      //   // eqvpTeststrain = 0.0;
      //   eqvpstrain = 0.0;
      //   // plastic_strain.zero();
      // }
      
      // evaluate the RHS equations ready for next Newton-Raphson iteration
      flow_dirn = flowPotential(sig);
      resid = flow_dirn * flow_incr - delta_dp;
      f = yieldFunction(sig, yield_stress);
      rep = -eqvpstrain + eqvpstrain_old - flow_incr * internalPotential();

      err1 = resid.L2norm();
      err2 = std::abs(f);
      err3 = std::abs(rep);
    }

    if (iter >= maxiter)
      mooseError("Constitutive failure");

    plastic_strain += delta_dp;
  }
}

Real
Test2FiniteStrainPlasticMaterial::yieldFunction(const RankTwoTensor & stress, const Real yield_stress)
{
  return getSigEqv(stress) - yield_stress;
}

RankTwoTensor
Test2FiniteStrainPlasticMaterial::dyieldFunction_dstress(const RankTwoTensor & sig)
{
  RankTwoTensor deriv = sig.dsecondInvariant();
  deriv *= std::sqrt(3.0 / sig.secondInvariant()) / 2.0;
  return deriv;
}

Real
Test2FiniteStrainPlasticMaterial::dyieldFunction_dinternal(const Real equivalent_plastic_strain, const Real hf)
{
  return -getdYieldStressdPlasticStrain(equivalent_plastic_strain,hf);
}

RankTwoTensor
Test2FiniteStrainPlasticMaterial::flowPotential(const RankTwoTensor & sig)
{
  return dyieldFunction_dstress(sig); // this plasticity model assumes associative flow
}

Real
Test2FiniteStrainPlasticMaterial::internalPotential()
{
  return -1;
}

Real
Test2FiniteStrainPlasticMaterial::getSigEqv(const RankTwoTensor & stress)
{
  return std::sqrt(3 * stress.secondInvariant());
}

// Jacobian for stress update algorithm
void
Test2FiniteStrainPlasticMaterial::getJac(const RankTwoTensor & sig,
                                    const RankFourTensor & E_ijkl,
                                    Real flow_incr,
                                    RankFourTensor & dresid_dsig)
{
  RankTwoTensor sig_dev, df_dsig, flow_dirn;
  RankTwoTensor dfi_dft, dfi_dsig;
  RankFourTensor dft_dsig, dfd_dft, dfd_dsig;
  Real sig_eqv;
  Real f1, f2, f3;
  RankFourTensor temp;

  sig_dev = sig.deviatoric();
  sig_eqv = getSigEqv(sig);
  df_dsig = dyieldFunction_dstress(sig);
  flow_dirn = flowPotential(sig);

  f1 = 3.0 / (2.0 * sig_eqv);
  f2 = f1 / 3.0;
  f3 = 9.0 / (4.0 * Utility::pow<3>(sig_eqv));

  dft_dsig = f1 * _deltaMixed - f2 * _deltaOuter - f3 * sig_dev.outerProduct(sig_dev);

  dfd_dsig = dft_dsig;
  dresid_dsig = E_ijkl.invSymm() + dfd_dsig * flow_incr;
}

// Obtain yield stress for a given equivalent plastic strain (input)
Real
Test2FiniteStrainPlasticMaterial::getYieldStress(const Real eqpe, const Real hf)
{
  unsigned nsize;

  nsize = _yield_stress_vector.size();

  if (_yield_stress_vector[0] > 0.0 || nsize % 2 > 0) // Error check for input inconsitency
    mooseError("Error in yield stress input: Should be a vector with eqv plastic strain and yield "
               "stress pair values.\n");

  unsigned int ind = 0;
  Real tol = 1e-8;

  // while (ind < nsize)
  // {
  //   if (std::abs(eqpe - _yield_stress_vector[ind]) < tol)
  //     return _yield_stress_vector[ind + 1];

  //   if (ind + 2 < nsize)
  //   {
  //     if (eqpe > _yield_stress_vector[ind] && eqpe < _yield_stress_vector[ind + 2])
  //       return _yield_stress_vector[ind + 1] +
  //              (eqpe - _yield_stress_vector[ind]) /
  //                  (_yield_stress_vector[ind + 2] - _yield_stress_vector[ind]) *
  //                  (_yield_stress_vector[ind + 3] - _yield_stress_vector[ind + 1]);
  //   }
  //   else
  //     return _yield_stress_vector[nsize - 1];

  //   ind += 2;
  // }

  Real value_0 = 4.0e3; //4.0e4; //3.0e3;

  return value_0 + hf*eqpe;
}

Real
Test2FiniteStrainPlasticMaterial::getdYieldStressdPlasticStrain(const Real eqpe, const Real hf)
{
  unsigned nsize;

  nsize = _yield_stress_vector.size();

  if (_yield_stress_vector[0] > 0.0 || nsize % 2 > 0) // Error check for input inconsitency
    mooseError("Error in yield stress input: Should be a vector with eqv plastic strain and yield "
               "stress pair values.\n");

  unsigned int ind = 0;

  Real tol = 1e-8;

  if (std::abs(eqpe-tol) < 0.0) // if (eqpe > 0)
    return hf;
    // return hard_factor;
  else
    return 0.0;

  return 0.0;
}

void
Test2FiniteStrainPlasticMaterial::computeHardFactor(Real & eqvpstrain,RankTwoTensor & plastic_strain)
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  _hard_factor[_qp] = 0.0;
  _grain_boundary[_qp] = 0.0;
  Real sum_h = 0.0;
  Real gb_bnd = 0.0;

  std::vector<Real> hardFactor = {0.9e5,1.1e5};

  // std::vector<Real> hardFactor = {2.0e4,2.0e4};

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;
    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all rotated elasticity tensors
    _hard_factor[_qp] += hardFactor[op_index] * h;
    sum_h += h;
    _grain_boundary[_qp] += (*_vals[op_index])[_qp]*(*_vals[op_index])[_qp];
  }
  
  if (_grain_boundary[_qp] >=0.47 && (*_vals[0])[_qp] >= 0.50) // 
  {
    eqvpstrain = 0;
    plastic_strain.zero();
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _hard_factor[_qp] /= sum_h;
}


// // std::vector<Real> 
// void
// Test2FiniteStrainPlasticMaterial::computeIntegral(Real & eqven_pstrain) // Real & eqv_pstrain, 
// {
//   const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
//   std::vector<Real> sum(2, 0); // 定义2个整数型元素的向量,且给出每个元素的初值为0
//   std::vector<Real> area(2, 0); // 
//   Real sum_hh = 0.0;
//   // Real eqv_strain = 0.0;
//   // std::vector<Real> average(2, 0);
//   // RealVectorValue average;
//   // _real_vec_prop[_qp](0) = 6;
//   _eqv_plastic_strain_op[_qp] = {0.0,0.0};
//   // if (_grain_boundary[_qp] <= 0.60)
//   //   {
//   //     eqven_pstrain = 0.0;
//   //     eqven_pstrain += _eqv_plasticity_strain[_qp](op_index)*hh;
//   //     sum_hh += hh;
//   //     _eqv_plasticity_strain[_qp] = 0;
//   //   }    
//   for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index) // 遍历 gr0 gr1
//   {   
//     Real hh = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0; // 插值函数

//     for (int qp = 0; qp < _qrule->n_points(); qp++)
//     {
//       Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[qp] - 0.5))) / 2.0; // 插值函数
//       sum[op_index] += _JxW[qp] * _coord[qp] * _eqv_plastic_strain[qp] * h; // 计算对于每个序参数gr_i每个网格上的等效塑性应变比重，$\sum_{i = 1}^{N}{h(\eta_i)\bar{\varepsilon}^p}$ 
//       area[op_index] += _JxW[qp] * _coord[qp] * h; // 每个序参数所占据的总面积； A(gr0) gr(gr1)
//     }

//     // _eqv_plastic_strain_op[_qp](op_index) = sum[op_index]/area[op_index]; // 每个晶粒内的平均等效塑性应变

//     eqven_pstrain += _eqv_plastic_strain_op[_qp](op_index)*hh; // 计算加权后的等效应变
//     sum_hh += hh;
//     // _eqv_plasticity_strain[op_index][qp] = sum[op_index]/area[op_index];
//     // if (_grain_boundary[_qp] <= 0.60)
//     // {
//     //   // eqven_pstrain = 0.0;
//     //   eqven_pstrain += _eqv_plasticity_strain[_qp](op_index)*hh;
//     //   sum_hh += hh;
//     //   // _eqv_plasticity_strain[_qp] = 0;
//     // }    
//   }
//   _eqv_plastic_strain_op[_qp](0) = sum[0]/area[0];
//   _eqv_plastic_strain_op[_qp](1) = sum[1]/area[1];
//   const Real tol = 1.0e-10;
//   sum_hh = std::max(sum_hh,tol);
//   eqven_pstrain /= sum_hh; // 计算加权后的等效应变
//   eqven_pstrain = 0.0;
//   // _eqv_plastic_strain_op[_qp] = {0.0,10.0};
//   // _eqv_plasticity_strain[]
//   // return average;
//   // eqv_pstrain = _eqv_plasticity_strain[_qp](1)*(*_vals[1])[_qp]+_eqv_plasticity_strain[_qp](1)*(*_vals[2])[_qp];
// }


// std::vector<Real> 
<<<<<<< HEAD
void
Test2FiniteStrainPlasticMaterial::computeIntegral(Real & eqven_pstrain) // Real & eqv_pstrain, 
{
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
  std::vector<Real> sum(2, 0); //
  std::vector<Real> area(2, 0); // 
  Real sum_hh = 0.0;
  // Real eqv_strain = 0.0;
  // std::vector<Real> average(2, 0);
  // RealVectorValue average;
  // _real_vec_prop[_qp](0) = 6;
  // if (_grain_boundary[_qp] <= 0.60)
  //   {
  //     eqven_pstrain = 0.0;
  //     eqven_pstrain += _eqv_plasticity_strain[_qp](op_index)*hh;
  //     sum_hh += hh;
  //     _eqv_plasticity_strain[_qp] = 0;
  //   }    
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index) // 遍历 gr0 gr1
  {   
    Real hh = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0; // 插值函数

    for (int qp = 0; qp < _qrule->n_points(); qp++)
    {
      Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[qp] - 0.5))) / 2.0; // 插值函数
      sum[op_index] += _JxW[qp] * _coord[qp] * _eqv_plastic_strain[qp] * h; // 计算对于每个序参数gr_i每个网格上的等效塑性应变比重，$\sum_{i = 1}^{N}{h(\eta_i)\bar{\varepsilon}^p}$ 
      area[op_index] += _JxW[qp] * _coord[qp] * h; // 每个序参数所占据的总面积；
    }
=======
// void
// Test2FiniteStrainPlasticMaterial::computeIntegral(Real & material_property) // Real & eqv_pstrain, 
// {
//   const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());
//   std::vector<Real> sum(2, 0); // 定义2个整数型元素的向量,且给出每个元素的初值为0
//   std::vector<Real> area(2, 0); // 定义2个整数型元素的向量,且给出每个元素的初值为0
//   Real sum_hh = 0.0; 
//   Real average = 0.0;
//   Real value = 1.0;

//   _eqv_plastic_strain_op[_qp] = {0.0,0.0};

//   // for (unsigned int op_index = 0; op_index < 2 ; ++op_index) // 0 1
//   // {   
//     unsigned int op_index = 0;
//     unsigned int size = _hard_factor.size();
//     for (unsigned int qp = 0; qp < _qrule->n_points(); qp++) // _qrule->n_points()--四节点单元网格
//     {
//       // Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[qp] - 0.5))) / 2.0; // 插值函数
//       // sum[op_index] += _JxW[qp] * _coord[qp] * _eqv_plastic_strain[qp] * h; // 计算对于每个序参数gr_i每个网格上的等效塑性应变比重，$\sum_{i = 1}^{N}{h(\eta_i)\bar{\varepsilon}^p}$ 
//       sum[op_index] += _JxW[qp] * _coord[qp] * value; //computeQpIntegral(); //  // _hard_factor[qp]; // 计算对于每个序参数gr_i每个网格上的等效塑性应变比重，$\sum_{i = 1}^{N}{h(\eta_i)\bar{\varepsilon}^p}$ 
//       // std::cout << "material_property = " << material_property << std::endl;
//       // std::cout << "sum[0] = " << sum[op_index] << std::endl;
//       area[op_index] += _JxW[qp] * _coord[qp]; //_coord[qp]; //* h _JxW[qp]* ; // 4, 每个序参数所占据的总面积； A(gr0) gr(gr1)
//       // std::cout << "area[0] = " << area[op_index] << std::endl;

//       // _JxW[qp] = 0.03;
//       // _coord[qp] = 1;
//       // _qrule->n_points() = 4;
//     }

    // _volume += this->_current_elem_volume; // 计算所的体积

    // _communicator.sum(area[op_index]);
    // std::cout << "area = " << area[0] << std::endl;
    // _communicator.sum(sum[op_index]);

    // std::cout << "sum[" << op_index << "] = " << sum[op_index];
    // std::cout << ", area[" << op_index << "] = " << area[op_index] << std::endl;
    // std::cout << ", size = " << size << std::endl;
    // std::cout << ", volume = " << _volume << std::endl;
    // _eqv_plastic_strain_op[_qp](op_index) = sum[op_index]/area[op_index]; // 每个晶粒内的平均等效塑性应变 
  // }
    // average = sum[0]/area[0];
    // std::cout << "average = " << average << std::endl;
  // _eqv_plastic_strain_op[_qp](1) = sum[1]/area[1];
>>>>>>> 44cc8e023973ccf6c4a443a690178e26204a8d48

    // eqven_pstrain = 0.0;
// }

// Real
// Test2FiniteStrainPlasticMaterial::computeQpIntegral()
// {
//   return MetaPhysicL::raw_value(_hard_factor[_qp]);
// }

// Real
// Test2FiniteStrainPlasticMaterial::getValue(Real & _integral_value)
// {
//   gatherSum(_integral_value);
//   return _integral_value;
// }
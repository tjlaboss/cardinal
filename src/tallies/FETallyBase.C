/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#ifdef ENABLE_OPENMC_COUPLING

#include "FETallyBase.h"

InputParameters
FETallyBase::validParams()
{
    auto params = TallyBase::validParams();
    params.addParam<std::string>("function_suffix","_function",
                            "The suffix to append to the score(s) as the name of the"
                            "function object(s) holding the functional expansion(s).");
    return params;
}

FETallyBase::FETallyBase(const InputParameters & parameters)
: TallyBase(parameters),
  _function_suffix(getParam<std::string>("function_suffix"))
{

    // FETs create FunctionSeries, not ElementalAuxVariables, so clear any names here
    // (in the base, so all FET types behave the same) and warn if the user set them.
    if (isParamValid("name"))
    {
      mooseWarning(this->_name + " does not have any ElementalAuxVariables "
                   "associated with it! "+this->_name+" creates functions, "
                   "the names of which are controllable by \"function_suffix\". "
                   "Clearing \"name\" parameter...");
    }
    _tally_name.clear();

    // OpenMC spatial FETs only support the collision estimator
    if (isParamValid("estimator"))
    {
        if (_estimator != openmc::TallyEstimator::COLLISION)
        paramError("estimator",
                    "Collision estimators are currently the only compatible "
                    "estimator type for Spatial Legendre expansion tallies!");
    }
    else
      _estimator = openmc::TallyEstimator::COLLISION;
}

void
FETallyBase::resetTally()
{
  // Erase the tally.
  openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index);

  for (int i = 0; i < getNumSptlFilters(); i++)
  {
    openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index + i);
  };
}

void
FETallyBase::computeSumAndMean()
{
  // Score-independent, so compute once and cache for storeResultsInner.
  _expansion_volume = computeVolume();

  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    // Read (do not mutate) the SUM slice: dividing the View in place would corrupt the
    // results_ tensor that relaxAndNormalizeTally re-reads. The zeroth coefficient is
    // the domain-integrated quantity, used directly as the undivided sum (as the base
    // computeSumAndMean sums raw SUM bins) so relaxation yields a shape with a unit
    // zeroth coefficient.
    auto coeffs_view = _local_tally->results_.slice(openmc::tensor::all,
                                                     score,
                                                     static_cast<int>(openmc::TallyResult::SUM));
    const Real zeroth_moment = coeffs_view[0];

    _local_sum_tally[score] = zeroth_moment;
    _local_mean_tally[score] = zeroth_moment / _local_tally->n_realizations_;
  }
}

Real
FETallyBase::storeResultsInner(const std::vector<unsigned int> & /* var_numbers */,
                               unsigned int local_score,
                               const std::vector<OMCTensor> & tally_vals,
                               bool norm_by_src_rate)
{
  // Relaxed, normalized coefficients for this score (unit zeroth coefficient).
  const auto & coeffs_tensor = tally_vals[local_score];

  // Scale the normalized shape to a physical density. tallyMultiplier is the same
  // power / source-rate factor a cell tally applies; (V_std / V_phys) is the volume
  // Jacobian MOOSE's own FX generation applies in FXIntegralBaseUserObject::finalize,
  // which we must reproduce here because we set the coefficients directly. The result
  // reconstructs to a density whose integral over the domain equals tallyMultiplier.
  Real scalar = 1.0;
  if (norm_by_src_rate)
  {
    const Real jacobian = _function->getStandardizedFunctionVolume() / _expansion_volume;
    scalar = _openmc_problem.tallyMultiplier(_tally_score[local_score],
                                             _local_mean_tally[local_score]) *
             jacobian;
  }

  std::vector<Real> coeffs(coeffs_tensor.begin(), coeffs_tensor.end());
  for (auto & c : coeffs)
    c *= scalar;
  _function->setCoefficients(coeffs);

  // Unit zeroth coefficient, mirroring a cell tally's sum-of-fractions for
  // checkNormalization.
  return coeffs_tensor.empty() ? 0.0 : coeffs_tensor[0];
}

Real
FETallyBase::computeVolume()
{
  Real volume = 0.0;
  for (const auto * elem :
       _openmc_problem.getMooseMesh().getMesh().active_local_element_ptr_range())
    volume += elem->volume();

  // Reduce the local partial volumes across all ranks.
  _openmc_problem.comm().sum(volume);

  return volume;
}
#endif

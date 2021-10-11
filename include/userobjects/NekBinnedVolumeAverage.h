#pragma once

#include "NekBinnedVolumeIntegral.h"

/**
 * Compute a volume average of the NekRS solution in spatial bins.
 */
class NekBinnedVolumeAverage : public NekBinnedVolumeIntegral
{
public:
  static InputParameters validParams();

  NekBinnedVolumeAverage(const InputParameters & parameters);

  virtual ~NekBinnedVolumeAverage();

  virtual void execute() override;

protected:
  /// Volumes of each bin
  double * _bin_volumes;
};
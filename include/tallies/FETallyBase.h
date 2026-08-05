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

#pragma once

#include "TallyBase.h"
#include "OpenMCCellAverageProblem.h"
#include "FunctionSeries.h"
#include <vector>
#include <string>
#include <utility>

class FETallyBase : public TallyBase
{
    public:

        static InputParameters validParams();

        FETallyBase(const InputParameters & parameters);

        virtual void resetTally() override;

        virtual void computeSumAndMean() override;

        // passing down to derived classes
        virtual std::pair<unsigned int, std::vector<openmc::Filter *>> spatialFilter() override
        {return std::pair<unsigned int, std::vector<openmc::Filter *>>(0, {nullptr});};

    protected:
        virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                       unsigned int local_score,
                                       const std::vector<OMCTensor> & tally_vals,
                                       bool norm_by_src_rate = true) override;

        /// Physical volume of the expansion domain (active local mesh, reduced across ranks).
        Real computeVolume();

        virtual int getNumBins() = 0;

        virtual int getNumSptlFilters() = 0;

        std::string _function_suffix;
        FunctionSeries* _function;

        /// Expansion volume cached by computeSumAndMean, reused by storeResultsInner.
        Real _expansion_volume = 0.0;
};

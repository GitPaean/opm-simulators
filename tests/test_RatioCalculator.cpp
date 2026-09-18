/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#define BOOST_TEST_MODULE RatioCalculatorTests
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/RatioCalculator.hpp>

#include <vector>

BOOST_AUTO_TEST_CASE(ThreePhaseProducerCrossflowPreservesOilGasSplit)
{
    Opm::RatioCalculator<double> calc(2, 1, 0, "producer");
    Opm::DeferredLogger logger;
    Opm::PerforationRates<double> split;
    const std::vector<double> rates{5.0, 10.3, 23.0};
    calc.perfRateInj(rates, split, 0.1, 2.0, 0.0, 0.0, 1e7, true, logger);
    BOOST_CHECK_CLOSE(split.free_gas, 3.0, 1e-10);
    BOOST_CHECK_CLOSE(split.free_oil, 10.0, 1e-10);
    BOOST_CHECK_CLOSE(split.free_gas + split.dis_gas, rates[2], 1e-10);
    BOOST_CHECK_CLOSE(split.free_oil + split.vap_oil, rates[1], 1e-10);
}

BOOST_AUTO_TEST_CASE(UnmixedConnectionsRetainBothFlowDirections)
{
    Opm::DeferredLogger logger;
    for (const auto& indices : {std::vector<int>{-1, 1, 0}, // oil/water
                               std::vector<int>{-1, 0, -1}, // oil only
                               std::vector<int>{0, -1, -1}}) { // gas only
        const int gas = indices[0], oil = indices[1], water = indices[2];
        const auto component = oil >= 0 ? oil : gas;
        Opm::RatioCalculator<double> calc(gas, oil, water, "unmixed");
        std::vector<double> rates(component + 1, 0.0);
        rates[component] = -10.0;
        Opm::PerforationRates<double> production, injection;
        calc.perfRateProd(rates, production, 0.0, 0.0, 0.0, 0.0);
        rates[component] = 3.0;
        calc.perfRateInj(rates, injection, 0.0, 0.0, 0.0, 0.0, 1e7, true, logger);
        const auto netFree = oil >= 0 ? production.free_oil + injection.free_oil
                                     : production.free_gas + injection.free_gas;
        BOOST_CHECK_EQUAL(netFree, -7.0);
    }
}

BOOST_AUTO_TEST_CASE(InjectorBackflowHasSolutionGas)
{
    Opm::RatioCalculator<double> calc(2, 1, 0, "injector");
    std::vector<double> rates{-5.0, -10.0, -3.0};
    Opm::PerforationRates<double> split;
    calc.perfRateProd(rates, split, 0.1, 2.0, 0.0, 0.0);
    BOOST_CHECK_EQUAL(split.free_gas, -3.0);
    BOOST_CHECK_EQUAL(split.dis_gas, -20.0);
    BOOST_CHECK_EQUAL(rates[2], -23.0);
    BOOST_CHECK_CLOSE(split.free_oil + split.vap_oil, rates[1], 1e-10);

    // Returning to injection must clear the previous solution contribution.
    Opm::DeferredLogger logger;
    const std::vector<double> injected{0.0, 0.0, 23.0};
    calc.perfRateInj(injected, split, 0.1, 2.0, 0.0, 0.0, 1e7, false, logger);
    BOOST_CHECK_EQUAL(split.free_gas, 23.0);
    BOOST_CHECK_EQUAL(split.dis_gas, 0.0);
    BOOST_CHECK_EQUAL(split.vap_oil, 0.0);
}

BOOST_AUTO_TEST_CASE(GasWaterSplitIncludesGasDissolvedInWater)
{
    Opm::RatioCalculator<double> calc(1, -1, 0, "gas-water");
    std::vector<double> rates{-10.0, -3.0};
    Opm::PerforationRates<double> produced, injected;
    calc.perfRateProd(rates, produced, 0.0, 0.0, 0.1, 2.0);
    BOOST_CHECK_EQUAL(produced.free_gas, -3.0);
    BOOST_CHECK_EQUAL(produced.dis_gas_in_water, -20.0);
    BOOST_CHECK_CLOSE(produced.free_gas + produced.dis_gas_in_water, rates[1], 1e-10);
    Opm::DeferredLogger logger;
    for (auto& rate : rates) { rate = -rate; }
    calc.perfRateInj(rates, injected, 0.0, 0.0, 0.1, 2.0, 1e7, true, logger);
    BOOST_CHECK_CLOSE(injected.free_gas, 3.0, 1e-10);
    BOOST_CHECK_CLOSE(injected.free_gas + injected.dis_gas_in_water, rates[1], 1e-10);
}

BOOST_AUTO_TEST_CASE(InvalidMixtureUsesDocumentedUnmixedFallback)
{
    Opm::RatioCalculator<double> calc(2, 1, 0, "invalid-mixture");
    Opm::DeferredLogger logger;
    Opm::PerforationRates<double> split;
    const std::vector<double> rates{5.0, 10.3, 23.0};
    calc.perfRateInj(rates, split, 0.5, 2.0, 0.1, 0.0, 1e7, true, logger);
    BOOST_CHECK_EQUAL(split.free_gas, rates[2]);
    BOOST_CHECK_EQUAL(split.free_oil, rates[1]);
    BOOST_CHECK_EQUAL(split.dis_gas, 0.0);
    BOOST_CHECK_EQUAL(split.vap_wat, 0.0);
}

BOOST_AUTO_TEST_CASE(ComponentFluxDerivativesArePreserved)
{
    using Eval = Opm::DenseAd::Evaluation<double, -1, 4u>;
    Opm::RatioCalculator<Eval> calc(1, 0, -1, "derivatives");
    std::vector<Eval> rates{Eval(4, -10.0, 0), Eval(4, -3.0, 1)};
    Opm::PerforationRates<double> split;
    calc.perfRateProd(rates, split, Eval(4, 0.1), Eval(4, 2.0), Eval(4, 0.0), Eval(4, 0.0));
    BOOST_CHECK_EQUAL(rates[1].derivative(0), 2.0);
    BOOST_CHECK_EQUAL(rates[1].derivative(1), 1.0);
    BOOST_CHECK_EQUAL(rates[0].derivative(0), 1.0);
    BOOST_CHECK_EQUAL(rates[0].derivative(1), 0.1);
    BOOST_CHECK_EQUAL(split.free_gas, -3.0);
}

BOOST_AUTO_TEST_CASE(MswComponentFluxDerivativesArePreserved)
{
    using Eval = Opm::DenseAd::Evaluation<double, 7>;
    Opm::RatioCalculator<Eval> calc(1, 0, -1, "msw-derivatives");
    std::vector<Eval> rates{Eval(-10.0), Eval(-3.0)};
    rates[0].setDerivative(0, 1.0);
    rates[1].setDerivative(1, 1.0);
    Opm::PerforationRates<double> split;
    calc.perfRateProd(rates, split, Eval(0.1), Eval(2.0), Eval(0.0), Eval(0.0));
    BOOST_CHECK_EQUAL(rates[1].derivative(0), 2.0);
    BOOST_CHECK_EQUAL(rates[1].derivative(1), 1.0);
    BOOST_CHECK_EQUAL(rates[0].derivative(0), 1.0);
    BOOST_CHECK_EQUAL(rates[0].derivative(1), 0.1);
    BOOST_CHECK_EQUAL(split.free_gas, -3.0);
}

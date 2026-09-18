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

#define BOOST_TEST_MODULE WellRateAllocationTests
#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/WellRateAllocation.hpp>

#include <limits>

BOOST_AUTO_TEST_CASE(PreservesFreeFractionAndAcceptedTotal)
{
    const auto split = Opm::allocateProductionRate(-22.0, -3.0, -20.0, 1e-12, false);
    BOOST_REQUIRE(split);
    BOOST_CHECK_CLOSE(split->free + split->solution, -22.0, 1e-10);
    BOOST_CHECK_CLOSE(split->free / split->solution, 3.0 / 20.0, 1e-10);
    BOOST_CHECK_LT(split->free, 0.0);
    BOOST_CHECK_LT(split->solution, 0.0);
}

BOOST_AUTO_TEST_CASE(ZeroFreeFluxDoesNotBecomeNegativeProduction)
{
    const auto split = Opm::allocateProductionRate(-19.99, 0.0, -20.0, 1e-12, false);
    BOOST_REQUIRE(split);
    BOOST_CHECK_EQUAL(split->free, 0.0);
    BOOST_CHECK_EQUAL(split->solution, -19.99);
}

BOOST_AUTO_TEST_CASE(StoppedWellHasZeroWellheadSplitDespiteConnectionCrossflow)
{
    const auto split = Opm::allocateProductionRate(0.0, -3.0, -20.0, 1e-12, true);
    BOOST_REQUIRE(split);
    BOOST_CHECK_EQUAL(split->free, 0.0);
    BOOST_CHECK_EQUAL(split->solution, 0.0);
}

BOOST_AUTO_TEST_CASE(RejectsUnsupportedOrUndefinedAllocation)
{
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, -3.0, -20.0, 1e-12, true));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, 3.0, -20.0, 1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, -3.0, 20.0, 1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(22.0, -3.0, -20.0, 1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, 0.0, 0.0, 1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, -1e-15, -1e-15, 1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, -3.0, -20.0, -1e-12, false));
    BOOST_CHECK(!Opm::allocateProductionRate(-22.0, std::numeric_limits<double>::quiet_NaN(), -20.0, 1e-12, false));
}

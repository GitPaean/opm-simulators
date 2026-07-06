// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \copydoc Opm::EQUIL::CompositionalInitializer
 */
#ifndef OPM_EQUIL_INITIALIZER_COMP_HPP
#define OPM_EQUIL_INITIALIZER_COMP_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>

#include <opm/simulators/flow/equil/InitStateEquil.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <fmt/format.h>

namespace Opm {
namespace EQUIL {

/*!
 * \brief Hydrostatic equilibration initializer for the compositional simulator.
 *
 * Given the EQUIL and ZMFVD keywords, this computes the initial per-cell state
 * needed by the compositional model: total composition \c z, phase pressures and
 * temperature.  The saturations and phase split are subsequently recomputed by the
 * flash in the intensive quantities, so they need not be provided here.
 *
 * The vertical pressure distribution is obtained by integrating the hydrostatic ODE
 *   dp/dz = rho(depth, p) * g
 * with the fourth-order Runge-Kutta integrator (\c Details::RK4IVP) that is shared
 * with the black-oil equilibration facility.  The density is evaluated from the
 * equation of state of the compositional fluid system.
 *
 * Two initialization types (EQUIL item 10) are supported:
 *  - type 1 (default): the ZMFVD composition is the total composition; the whole
 *    column is integrated as a single (liquid) phase.
 *  - type 3: the ZMFVD composition is the liquid composition.  The datum is taken at
 *    the gas-oil contact, where the pressure equals the bubble-point pressure of the
 *    liquid composition.  Below the contact the liquid composition is used with the
 *    liquid density; above the contact the equilibrium vapour composition at the
 *    contact is used with the gas density.
 */
template <class FluidSystem, class Scalar>
class CompositionalInitializer
{
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;

    using EOSType = CompositionalConfig::EOSType;
    using CompArray = std::array<Scalar, numComponents>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
    // Provides the total composition as a function of depth.
    using CompositionFunc = std::function<CompArray(Scalar)>;

public:
    using ScalarFluidState = FluidState;

    /// Compute the initial states.
    ///
    /// \param[in] eclState       Deck-derived state (EQUIL, ZMFVD, RTEMP, EQLNUM).
    /// \param[in] eosType        Equation of state used by the fluid system.
    /// \param[in] cellDepth      Cell-centre depth for every grid cell.
    /// \param[in] eqlnum         Zero-based equilibration region for every cell.
    /// \param[in] gravity        Gravitational acceleration (m/s^2).
    /// \param[in] numPressurePts Number of RK4 sub-steps per pressure table.
    CompositionalInitializer(const EclipseState& eclState,
                             const EOSType eosType,
                             const std::vector<Scalar>& cellDepth,
                             const std::vector<int>& eqlnum,
                             const Scalar gravity,
                             const int numPressurePts)
        : eosType_(eosType)
    {
        const auto& equilContainer = eclState.getInitConfig().getEquil();
        const std::vector<EquilRecord> records(equilContainer.begin(), equilContainer.end());
        const int numRegions = static_cast<int>(records.size());

        const auto& tables = eclState.getTableManager();
        const auto& zmfvdTables = tables.getZmfvdTables();
        if (zmfvdTables.empty()) {
            throw std::runtime_error("Compositional equilibration requires the ZMFVD keyword "
                                     "to specify the composition versus depth.");
        }

        // Build one region-equilibration description per equilibration region.
        std::vector<RegionEquil> regions;
        regions.reserve(numRegions);
        for (int r = 0; r < numRegions; ++r) {
            regions.push_back(makeRegion(records[r], zmfvdTables.template getTable<ZmfvdTable>(r),
                                         tables, gravity, numPressurePts,
                                         cellDepth, eqlnum, r));
        }

        // Assign the per-cell state.
        const std::size_t numCells = cellDepth.size();
        initialFluidStates_.resize(numCells);
        for (std::size_t c = 0; c < numCells; ++c) {
            const int reg = eqlnum[c];
            if (reg < 0 || reg >= numRegions) {
                throw std::runtime_error("EQLNUM value out of range in compositional equilibration.");
            }
            assignCell(initialFluidStates_[c], regions[reg], cellDepth[c]);
        }
    }

    const ScalarFluidState& initialFluidState(unsigned cellIdx) const
    { return initialFluidStates_[cellIdx]; }

    std::vector<ScalarFluidState>& initialFluidStates()
    { return initialFluidStates_; }

private:
    // Everything needed to evaluate the equilibrated state within one region.
    struct RegionEquil {
        Scalar gocDepth{0.0};
        bool twoPhase{false};       // gas zone above the contact (EQUIL item 10 == 3)
        CompArray vaporComp{};      // equilibrium vapour composition at the contact
        CompositionFunc oilComp;    // composition vs depth in the liquid (ZMFVD) zone
        TabulatedFunction tempVsDepth;
        // Phase-pressure distributions (Pa) as a function of depth.
        std::vector<Scalar> depthTab;
        std::vector<Scalar> oilPressTab;
        std::vector<Scalar> gasPressTab;
    };

    // Right-hand side of the hydrostatic ODE for a single (miscible) phase.
    class DensityOde {
    public:
        DensityOde(const CompositionFunc& comp,
                   const TabulatedFunction& tempVsDepth,
                   const int phaseIdx,
                   const EOSType eosType,
                   const Scalar gravity)
            : comp_(comp), tempVsDepth_(tempVsDepth)
            , phaseIdx_(phaseIdx), eosType_(eosType), g_(gravity)
        {}

        Scalar operator()(const Scalar depth, const Scalar press) const
        {
            const CompArray z = comp_(depth);
            const Scalar temp = tempVsDepth_.eval(depth, /*extrapolate=*/true);

            FluidState fs;
            fs.setTemperature(temp);
            fs.setPressure(oilPhaseIdx, press);
            fs.setPressure(gasPhaseIdx, press);
            for (int c = 0; c < numComponents; ++c) {
                fs.setMoleFraction(phaseIdx_, c, z[c]);
            }

            typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType_);
            paramCache.updatePhase(fs, phaseIdx_);
            const Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx_);
            return rho * g_;
        }

    private:
        const CompositionFunc& comp_;
        const TabulatedFunction& tempVsDepth_;
        int phaseIdx_;
        EOSType eosType_;
        Scalar g_;
    };

    RegionEquil makeRegion(const EquilRecord& rec,
                           const ZmfvdTable& zmfvd,
                           const TableManager& tables,
                           const Scalar gravity,
                           const int numPressurePts,
                           const std::vector<Scalar>& cellDepth,
                           const std::vector<int>& eqlnum,
                           const int regionIdx)
    {
        RegionEquil reg;

        // Composition versus depth from ZMFVD (interpreted as the liquid
        // composition for type 3, or as the total composition otherwise).
        std::vector<TabulatedFunction> zmfTab(numComponents);
        const auto& depthCol = zmfvd.getDepthColumn();
        std::vector<Scalar> depths(depthCol.begin(), depthCol.end());
        for (int c = 0; c < numComponents; ++c) {
            const auto& col = zmfvd.getMoleFractionColumn(c);
            std::vector<Scalar> vals(col.begin(), col.end());
            zmfTab[c].setXYContainers(depths, vals);
        }
        reg.oilComp = [zmfTab](Scalar depth) {
            CompArray z;
            Scalar sum = 0.0;
            for (int c = 0; c < numComponents; ++c) {
                z[c] = std::max(Scalar{0}, zmfTab[c].eval(depth, /*extrapolate=*/true));
                sum += z[c];
            }
            if (sum > 0.0) {
                for (int c = 0; c < numComponents; ++c) {
                    z[c] /= sum;
                }
            }
            return z;
        };

        // Temperature versus depth (RTEMPVD or the constant RTEMP).
        reg.tempVsDepth = makeTemperatureTable(tables, regionIdx);

        reg.gocDepth = static_cast<Scalar>(rec.gasOilContactDepth());
        const Scalar gocPc = static_cast<Scalar>(rec.gasOilContactCapillaryPressure());

        // EQUIL item 10: 1 (default) treats ZMFVD as the total composition, 3 treats
        // it as the liquid composition with a gas cap above the gas-oil contact.
        // Type 2 (vapour composition specified) would need the equilibrium liquid at
        // the contact instead and is not implemented.
        const int initType = rec.compositionalInitType();
        if (initType != 1 && initType != 3) {
            throw std::runtime_error(fmt::format("Compositional equilibration type {} (EQUIL item 10) "
                                                 "is not supported; only types 1 and 3 are.", initType));
        }
        reg.twoPhase = (initType == 3);

        // Depth span covered by the region's cells. For a two-phase region the span
        // is extended to the gas-oil contact, which anchors both pressure curves.
        // (The reference depth may lie outside the span; the integration handles
        // that, so the span is not extended to the datum for single-phase regions.)
        Scalar zmin = reg.gocDepth;
        Scalar zmax = reg.gocDepth;
        bool found = false;
        for (std::size_t c = 0; c < cellDepth.size(); ++c) {
            if (eqlnum[c] != regionIdx) {
                continue;
            }
            zmin = found ? std::min(zmin, cellDepth[c]) : cellDepth[c];
            zmax = found ? std::max(zmax, cellDepth[c]) : cellDepth[c];
            found = true;
        }
        if (reg.twoPhase) {
            zmin = std::min(zmin, reg.gocDepth);
            zmax = std::max(zmax, reg.gocDepth);
        }

        // Only the hydrocarbon column is equilibrated; a water zone within the
        // region cannot be represented yet.
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            const Scalar wocDepth = static_cast<Scalar>(rec.waterOilContactDepth());
            if (wocDepth < zmax) {
                throw std::runtime_error(fmt::format("Compositional equilibration does not support "
                                                     "a water zone within the region yet: the water-oil "
                                                     "contact ({} m) lies above the deepest cell ({} m) "
                                                     "of region {}.", wocDepth, zmax, regionIdx + 1));
            }
            OpmLog::warning(fmt::format("Compositional equilibration region {}: the water phase is "
                                        "initialized with zero saturation (connate water is ignored).",
                                        regionIdx + 1));
        }

        // Reference point of the oil-pressure integration.
        Scalar oilRefDepth;
        Scalar oilRefPress;
        if (reg.twoPhase) {
            // Datum is the gas-oil contact; the pressure there is the bubble-point
            // pressure of the liquid composition at the contact.
            if (std::abs(static_cast<Scalar>(rec.datumDepth()) - reg.gocDepth) > 0) {
                OpmLog::warning(fmt::format("Compositional equilibration region {}: the datum depth "
                                            "({} m) must be at the gas-oil contact when EQUIL item 10 "
                                            "is 3; it is reset to the contact depth ({} m).",
                                            regionIdx + 1, rec.datumDepth(), reg.gocDepth));
            }
            const CompArray zGoc = reg.oilComp(reg.gocDepth);
            const Scalar tGoc = reg.tempVsDepth.eval(reg.gocDepth, /*extrapolate=*/true);
            const Scalar psat = bubblePointPressure(zGoc, tGoc, reg.vaporComp);
            oilRefDepth = reg.gocDepth;
            oilRefPress = rec.setToSaturationPressure() ? psat
                        : static_cast<Scalar>(rec.datumDepthPressure());

            OpmLog::info(fmt::format("Compositional equilibration region {}: two-phase "
                                     "(liquid composition specified). Gas-oil contact at {:.4g} m, "
                                     "bubble-point pressure {:.6g} bar.",
                                     regionIdx + 1, reg.gocDepth, psat / 1e5));
        } else {
            oilRefDepth = static_cast<Scalar>(rec.datumDepth());
            oilRefPress = static_cast<Scalar>(rec.datumDepthPressure());

            OpmLog::info(fmt::format("Compositional equilibration region {}: single-phase "
                                     "(total composition specified). Datum {:.4g} m at {:.6g} bar.",
                                     regionIdx + 1, oilRefDepth, oilRefPress / 1e5));
        }

        // Integrate the phase pressures over a fine depth table.
        buildPressureTables(reg, gravity, numPressurePts, zmin, zmax,
                            oilRefDepth, oilRefPress, gocPc);

        return reg;
    }

    // Integrate the oil (and, for two-phase, gas) pressure distribution and store
    // it on an equidistant depth table for later interpolation.
    void buildPressureTables(RegionEquil& reg,
                             const Scalar gravity,
                             const int numPressurePts,
                             const Scalar zmin,
                             const Scalar zmax,
                             const Scalar oilRefDepth,
                             const Scalar oilRefPress,
                             const Scalar gocPc)
    {
        const int nsample = std::max(numPressurePts, 2);

        // Guard against a degenerate span, e.g. a region without any cells on
        // this process.
        const Scalar zbot = (zmax > zmin) ? zmax : zmin + 1.0;
        const std::array<Scalar, 2> span{zmin, zbot};

        reg.depthTab.resize(nsample + 1);
        const Scalar h = (zbot - zmin) / nsample;
        for (int i = 0; i <= nsample; ++i) {
            reg.depthTab[i] = zmin + i * h;
        }

        // Oil pressure: integrate up from the reference to zmin and down to zmax.
        const DensityOde oilOde(reg.oilComp, reg.tempVsDepth, oilPhaseIdx, eosType_, gravity);
        reg.oilPressTab = integratePressure(oilOde, oilRefDepth, oilRefPress, span, nsample);

        if (reg.twoPhase) {
            // Gas pressure above the contact using the constant vapour composition.
            const CompArray vapor = reg.vaporComp;
            const CompositionFunc gasComp = [vapor](Scalar /*depth*/) { return vapor; };
            const DensityOde gasOde(gasComp, reg.tempVsDepth, gasPhaseIdx, eosType_, gravity);
            const Scalar gasRefPress = oilRefPress + gocPc; // gas pressure at the contact
            reg.gasPressTab = integratePressure(gasOde, reg.gocDepth, gasRefPress, span, nsample);
        }
    }

    // Integrate dp/dz = f(depth, p) from the reference point to both ends of the
    // span, and evaluate the result on the equidistant depth table.
    template <class Ode>
    std::vector<Scalar> integratePressure(const Ode& ode,
                                          const Scalar refDepth,
                                          const Scalar refPress,
                                          const std::array<Scalar, 2>& span,
                                          const int nsample) const
    {
        const Scalar zmin = span[0];
        const Scalar zmax = span[1];
        const Scalar total = zmax - zmin;

        // Number of RK4 sub-steps for the upward (towards zmin) and downward
        // (towards zmax) integration, proportional to each segment's extent. A
        // segment may be empty when the reference sits at a span end.
        const bool hasUp = refDepth > zmin;
        const bool hasDown = refDepth < zmax;
        const int nUp = std::max(1, static_cast<int>(std::ceil(nsample * (refDepth - zmin) / total)));
        const int nDown = std::max(1, static_cast<int>(std::ceil(nsample * (zmax - refDepth) / total)));

        Details::RK4IVP<Scalar, Ode> up(ode, {refDepth, hasUp ? zmin : refDepth}, refPress, nUp);
        Details::RK4IVP<Scalar, Ode> down(ode, {refDepth, hasDown ? zmax : refDepth}, refPress, nDown);

        std::vector<Scalar> press(nsample + 1);
        const Scalar h = total / nsample;
        for (int i = 0; i <= nsample; ++i) {
            const Scalar depth = zmin + i * h;
            if (depth < refDepth && hasUp) {
                press[i] = up(depth);
            }
            else if (depth > refDepth && hasDown) {
                press[i] = down(depth);
            }
            else {
                press[i] = refPress;
            }
        }
        return press;
    }

    static Scalar interpTable(const std::vector<Scalar>& depthTab,
                              const std::vector<Scalar>& valTab,
                              const Scalar depth)
    {
        const std::size_t n = depthTab.size();
        if (n == 0) {
            return 0.0;
        }
        if (depth <= depthTab.front()) {
            return valTab.front();
        }
        if (depth >= depthTab.back()) {
            return valTab.back();
        }
        // Equidistant table: locate the bracketing interval directly.
        const Scalar h = depthTab[1] - depthTab[0];
        std::size_t i = static_cast<std::size_t>((depth - depthTab.front()) / h);
        if (i >= n - 1) {
            i = n - 2;
        }
        const Scalar t = (depth - depthTab[i]) / (depthTab[i + 1] - depthTab[i]);
        return valTab[i] * (1 - t) + valTab[i + 1] * t;
    }

    void assignCell(FluidState& fs, const RegionEquil& reg, const Scalar depth) const
    {
        const bool inGasZone = reg.twoPhase && (depth < reg.gocDepth);

        // Total composition and phase pressure.
        const CompArray z = inGasZone ? reg.vaporComp : reg.oilComp(depth);
        const Scalar press = inGasZone
            ? interpTable(reg.depthTab, reg.gasPressTab, depth)
            : interpTable(reg.depthTab, reg.oilPressTab, depth);
        const Scalar temp = reg.tempVsDepth.eval(depth, /*extrapolate=*/true);

        fs.setTemperature(temp);
        for (int p = 0; p < numPhases; ++p) {
            if (FluidSystem::phaseIsActive(p)) {
                fs.setPressure(p, press);
            }
        }
        for (int c = 0; c < numComponents; ++c) {
            fs.setMoleFraction(c, z[c]);
        }

        // Placeholder saturations. The actual phase split is recomputed by the flash
        // in the intensive quantities from the total composition; these values only
        // guard against reads of uninitialized state during setup.
        for (int p = 0; p < numPhases; ++p) {
            if (FluidSystem::phaseIsActive(p)) {
                fs.setSaturation(p, 0.0);
            }
        }
        fs.setSaturation(inGasZone ? gasPhaseIdx : oilPhaseIdx, 1.0);
    }

    TabulatedFunction makeTemperatureTable(const TableManager& tables, const int regionIdx) const
    {
        TabulatedFunction table;
        if (tables.hasTables("RTEMPVD")) {
            const auto& rtempvdTables = tables.getRtempvdTables();
            const auto& t = rtempvdTables.template getTable<RtempvdTable>(regionIdx);
            std::vector<Scalar> x(t.getDepthColumn().begin(), t.getDepthColumn().end());
            std::vector<Scalar> y(t.getTemperatureColumn().begin(), t.getTemperatureColumn().end());
            table.setXYContainers(x, y);
        } else {
            const Scalar rtemp = static_cast<Scalar>(tables.rtemp());
            std::vector<Scalar> x{0.0, 1.0};
            std::vector<Scalar> y{rtemp, rtemp};
            table.setXYContainers(x, y);
        }
        return table;
    }

    // Bubble-point pressure of a liquid of composition \c x at temperature \c T.
    // Also returns the equilibrium vapour composition \c yOut at the bubble point.
    Scalar bubblePointPressure(const CompArray& x, const Scalar T, CompArray& yOut) const
    {
        // Wilson correlation for the initial equilibrium ratios / pressure guess.
        auto wilsonK = [&](const Scalar p) {
            CompArray k;
            for (int c = 0; c < numComponents; ++c) {
                const Scalar pc = FluidSystem::criticalPressure(c);
                const Scalar tc = FluidSystem::criticalTemperature(c);
                const Scalar omega = FluidSystem::acentricFactor(c);
                k[c] = (pc / p) * std::exp(5.373 * (1.0 + omega) * (1.0 - tc / T));
            }
            return k;
        };

        // Full fugacity equality: iterate K = phi_liq / phi_vap at fixed pressure
        // (inner loop), then adjust the pressure so that sum(K*x) = 1 (outer loop).
        // The pressure is approached from below (starting well under the bubble
        // point) to avoid collapsing onto the trivial solution K == 1.
        Scalar p = 1.0e5; // 1 bar
        CompArray K = wilsonK(p);
        CompArray y = K;
        for (int outer = 0; outer < 500; ++outer) {
            for (int inner = 0; inner < 100; ++inner) {
                Scalar s = 0.0;
                for (int c = 0; c < numComponents; ++c) {
                    y[c] = K[c] * x[c];
                    s += y[c];
                }
                for (int c = 0; c < numComponents; ++c) {
                    y[c] /= s;
                }

                FluidState fs;
                fs.setTemperature(T);
                fs.setPressure(oilPhaseIdx, p);
                fs.setPressure(gasPhaseIdx, p);
                for (int c = 0; c < numComponents; ++c) {
                    fs.setMoleFraction(oilPhaseIdx, c, x[c]);
                    fs.setMoleFraction(gasPhaseIdx, c, y[c]);
                }
                typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType_);
                paramCache.updatePhase(fs, oilPhaseIdx);
                paramCache.updatePhase(fs, gasPhaseIdx);

                Scalar maxDiff = 0.0;
                for (int c = 0; c < numComponents; ++c) {
                    const Scalar phiL = FluidSystem::fugacityCoefficient(fs, paramCache, oilPhaseIdx, c);
                    const Scalar phiV = FluidSystem::fugacityCoefficient(fs, paramCache, gasPhaseIdx, c);
                    const Scalar kNew = phiL / phiV;
                    maxDiff = std::max(maxDiff, std::abs(kNew - K[c]));
                    K[c] = kNew;
                }
                if (maxDiff < 1e-12) {
                    break;
                }
            }

            Scalar s = 0.0;
            for (int c = 0; c < numComponents; ++c) {
                s += K[c] * x[c];
            }
            if (std::abs(s - 1.0) < 1e-11) {
                break;
            }
            p *= std::sqrt(s);
        }

        Scalar s = 0.0;
        for (int c = 0; c < numComponents; ++c) {
            yOut[c] = K[c] * x[c];
            s += yOut[c];
        }
        for (int c = 0; c < numComponents; ++c) {
            yOut[c] /= s;
        }
        return p;
    }

    EOSType eosType_;
    std::vector<ScalarFluidState> initialFluidStates_;
};

} // namespace EQUIL
} // namespace Opm

#endif // OPM_EQUIL_INITIALIZER_COMP_HPP

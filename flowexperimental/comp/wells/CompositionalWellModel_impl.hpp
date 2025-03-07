
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/utils/phaseUsageFromDeck.hpp>

namespace Opm
{

template <typename TypeTag>
CompositionalWellModel<TypeTag>::
CompositionalWellModel(Simulator& simulator)
    : WellConnectionModule(*this, simulator.gridView().comm())
    , simulator_(simulator)
    , schedule_(simulator.vanguard().schedule())
    , summary_state_(simulator.vanguard().summaryState())
    , ecl_state_(simulator.vanguard().eclState())
    , comm_(simulator.gridView().comm())
    , phase_usage_(phaseUsageFromDeck(ecl_state_))
    , comp_config_(ecl_state_.compositionalConfig())
    , comp_well_states_(phase_usage_, comp_config_)
{
    local_num_cells_ = simulator.gridView().size(0);
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
beginReportStep(unsigned report_step)
{
    // TODO: not considering the parallel running yet
    wells_ecl_ = schedule_.getWells(report_step);
    initWellConnectionData();
    initWellState();
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
beginTimeStep()
{
    createWellContainer();
    initWellContainer();
//    for (auto& well : well_container_) {
//        auto& well_state = comp_well_states_[well->name()];
//        well->calculateExplicitQuantities(simulator_, well_state);
//        well->updateSecondaryQuantities(simulator_);
//        well->solveWellEq(simulator_, well_state);
//    }
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
init()
{
    simulator_.model().addAuxiliaryModule(this);
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
createWellContainer()
{
    // const auto& schedule = simulator_.vanguard().schedule();
    const auto nw = wells_ecl_.size(); // not considering the parallel running yet
    well_container_.clear();
    for (auto w = 0 * nw; w < nw; ++w) {
        well_container_.emplace_back(std::make_shared<CompWell<TypeTag>>(wells_ecl_[w], w, well_connection_data_[w]));
    }
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellContainer()
{
    for (auto& well : well_container_) {
        well->init();
    }
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellConnectionData()
{
    // TODO: we need to consider the parallel running
    // we can refer to the BlackoilWellModelGeneric::initializeWellPerfData()
    well_connection_data_.resize(wells_ecl_.size());

    int well_index = 0;
    for (const auto& well : wells_ecl_) {
        int connection_index = 0;
        const auto& well_connections = well.getConnections();
        auto& well_connection_data = well_connection_data_[well_index];

        well_connection_data.reserve(well_connections.size());
        for (const auto& connection : well_connections) {
            const auto active_index =
                    this->compressedIndexForInterior(connection.global_index());

            const auto connIsOpen =
                    connection.state() == Connection::State::OPEN;

            if (connIsOpen && (active_index >= 0)) {
                auto& pd = well_connection_data_[well_index].emplace_back();

                pd.cell_index = active_index;
                pd.connection_transmissibility_factor = connection.CF();
                pd.satnum_id = connection.satTableId();
                pd.ecl_index = connection_index;
            }
            ++connection_index;
        }
        ++well_index;
    }

}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
initWellState()
{
    // TODO: not sure the following is correct
    const auto pressIx = [this]()
    {
        if (phase_usage_.phase_used[FluidSystem::oilPhaseIdx]) {
            return FluidSystem::oilPhaseIdx;
        }
        if (phase_usage_.phase_used[FluidSystem::gasPhaseIdx]) {
            return FluidSystem::gasPhaseIdx;
        }
        assert(false && "the usage of oil and gas phase is not correct");
        return FluidSystem::gasPhaseIdx;
    }();

    auto cell_pressure = std::vector<Scalar>(this->local_num_cells_, Scalar{0.});
    auto cell_mole_fractions = std::vector<std::vector<Scalar>>(this->local_num_cells_,
                                           std::vector<Scalar>(FluidSystem::numComponents, Scalar{0.}));

    auto elemCtx = ElementContext { this->simulator_ };
    const auto& gridView = this->simulator_.vanguard().gridView();

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

        const auto ix = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        const auto& fs = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0).fluidState();

        cell_pressure[ix] = fs.pressure(pressIx).value();
        for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
            cell_mole_fractions[ix][compIdx] = fs.moleFraction(compIdx).value();
        }
    }
    OPM_END_PARALLEL_TRY_CATCH("ComposotionalWellModel::initializeWellState() failed: ",
                           this->simulator_.vanguard().grid().comm());

    /* TODO: no prev well state for now */
    this->comp_well_states_.init(this->wells_ecl_,
                                 cell_pressure, cell_mole_fractions, this->well_connection_data_,
                                 this->summary_state_);
}


template <typename TypeTag>
std::size_t
CompositionalWellModel<TypeTag>::
compressedIndexForInterior(std::size_t cartesian_cell_idx) const {
    return simulator_.vanguard().compressedIndexForInterior(cartesian_cell_idx);
}

template <typename TypeTag>
std::vector<int>
CompositionalWellModel<TypeTag>::
getCellsForConnections(const Well& well) const
{
    std::vector<int> wellCells;
    // All possible connections of the well
    const auto& connectionSet = well.getConnections();
    wellCells.reserve(connectionSet.size());

    for (const auto& connection : connectionSet)
    {
        int compressed_idx = compressedIndexForInterior(connection.global_index());
        if (compressed_idx >= 0) { // Ignore connections in inactive/remote cells.
            wellCells.push_back(compressed_idx);
        }
    }

    return wellCells;

}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
beginIteration()
{
    // do we need to do every iteration here?
    const auto& grid = simulator_.vanguard().grid();
    const auto& gridView = grid.leafGridView();
    ElementContext elemCtx(simulator_);
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
    }

    assemble(simulator_.model().newtonMethod().numIterations(), simulator_.timeStepSize());
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
assemble(const int iterationIdx,
         const double dt)
{

    if (iterationIdx == 0) {
        this->calculateExplicitQuantities();
    }
    for (auto& well : well_container_) {
        auto& well_state = comp_well_states_[well->name()];
        well->iterateWellEq(simulator_, dt, well_state);
        // TODO: do we need this? can we use the one assembled previously?
        // well->assembleWellEq(simulator_, well_state, dt);
    }
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
calculateExplicitQuantities()
{
    for (auto& well : well_container_) {
        auto& well_state = comp_well_states_[well->name()];
        well->calculateExplicitQuantities(simulator_, well_state);
    }
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
computeTotalRatesForDof(RateVector& rate,
                        unsigned globalIdx) const {
    for (const auto& well: well_container_) {
        well->addCellRates(rate, globalIdx);
    }
}

template<typename TypeTag>
void
CompositionalWellModel<TypeTag>::
recoverWellSolutionAndUpdateWellState(const BVector& x)
{
    {
        for (const auto& well : well_container_) {
            const auto& cells = well->cells();
            x_local_.resize(cells.size());

            for (size_t i = 0; i < cells.size(); ++i) {
                x_local_[i] = x[cells[i]];
            }
            auto& ws = this->comp_well_states_[well->name()];
            well->recoverWellSolutionAndUpdateWellState(x_local_, ws);
        }
    }
}

template <typename TypeTag>
bool
CompositionalWellModel<TypeTag>::
getWellConvergence() const
{
    bool converged = true;
    for (const auto& well : this->well_container_) {
        converged = converged && well->getConvergence();
    }
    return converged;
}

template <typename TypeTag>
data::Wells
CompositionalWellModel<TypeTag>::
wellData() const
{
    return this->comp_well_states_.report();
}

//template <typename TypeTag>
//void
//CompositionalWellModel<TypeTag>::
//linearize(SparseMatrixAdapter&, GlobalEqVector&)
//{
//
//}


} // end of namespace Opm
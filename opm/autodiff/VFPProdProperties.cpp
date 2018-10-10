/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#include "config.h"

#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/autodiff/VFPHelpers.hpp>



namespace Opm {




VFPProdProperties::VFPProdProperties() {

}



VFPProdProperties::VFPProdProperties(const VFPProdTable* table){
    m_tables[table->getTableNum()] = table;
}




VFPProdProperties::VFPProdProperties(const std::map<int, std::shared_ptr<const VFPProdTable> >& tables) {
    for (const auto& table : tables) {
        m_tables[table.first] = table.second.get();
    }
}

VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
                                              const Wells& wells,
                                              const ADB& qs,
                                              const ADB& thp_arg,
                                              const ADB& alq) const {
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(wells.number_of_phases == 3);
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    return bhp(table_id, w, o, g, thp_arg, alq);
}


VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
                                              const ADB& aqua,
                                              const ADB& liquid,
                                              const ADB& vapour,
                                              const ADB& thp_arg,
                                              const ADB& alq) const {
    const int nw = thp_arg.size();

    std::vector<int> block_pattern = detail::commonBlockPattern(aqua, liquid, vapour, thp_arg, alq);

    assert(static_cast<int>(table_id.size()) == nw);
    assert(aqua.size()     == nw);
    assert(liquid.size()   == nw);
    assert(vapour.size()   == nw);
    assert(thp_arg.size()      == nw);
    assert(alq.size()      == nw);

    //Allocate data for bhp's and partial derivatives
    ADB::V value = ADB::V::Zero(nw);
    ADB::V dthp = ADB::V::Zero(nw);
    ADB::V dwfr = ADB::V::Zero(nw);
    ADB::V dgfr = ADB::V::Zero(nw);
    ADB::V dalq = ADB::V::Zero(nw);
    ADB::V dflo = ADB::V::Zero(nw);

    //Get the table for each well
    std::vector<const VFPProdTable*> well_tables(nw, nullptr);
    for (int i=0; i<nw; ++i) {
        if (table_id[i] >= 0) {
            well_tables[i] = detail::getTable(m_tables, table_id[i]);
        }
    }

    //Get the right FLO/GFR/WFR variable for each well as a single ADB
    const ADB flo = detail::combineADBVars<VFPProdTable::FLO_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB wfr = detail::combineADBVars<VFPProdTable::WFR_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB gfr = detail::combineADBVars<VFPProdTable::GFR_TYPE>(well_tables, aqua, liquid, vapour);

    //Compute the BHP for each well independently
    for (int i=0; i<nw; ++i) {
        const VFPProdTable* table = well_tables[i];
        if (table != nullptr) {
            //First, find the values to interpolate between
            //Value of FLO is negative in OPM for producers, but positive in VFP table
            auto flo_i = detail::findInterpData(-flo.value()[i], table->getFloAxis());
            auto thp_i = detail::findInterpData( thp_arg.value()[i], table->getTHPAxis());
            auto wfr_i = detail::findInterpData( wfr.value()[i], table->getWFRAxis());
            auto gfr_i = detail::findInterpData( gfr.value()[i], table->getGFRAxis());
            auto alq_i = detail::findInterpData( alq.value()[i], table->getALQAxis());

            detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

            value[i] = bhp_val.value;
            dthp[i] = bhp_val.dthp;
            dwfr[i] = bhp_val.dwfr;
            dgfr[i] = bhp_val.dgfr;
            dalq[i] = bhp_val.dalq;
            dflo[i] = bhp_val.dflo;
        }
        else {
            value[i] = -1e100; //Signal that this value has not been calculated properly, due to "missing" table
        }
    }

    //Create diagonal matrices from ADB::Vs
    ADB::M dthp_diag(dthp.matrix().asDiagonal());
    ADB::M dwfr_diag(dwfr.matrix().asDiagonal());
    ADB::M dgfr_diag(dgfr.matrix().asDiagonal());
    ADB::M dalq_diag(dalq.matrix().asDiagonal());
    ADB::M dflo_diag(dflo.matrix().asDiagonal());

    //Calculate the Jacobians
    const int num_blocks = block_pattern.size();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        //Could have used fastSparseProduct and temporary variables
        //but may not save too much on that.
        jacs[block] = ADB::M(nw, block_pattern[block]);

        if (!thp_arg.derivative().empty()) {
            jacs[block] += dthp_diag * thp_arg.derivative()[block];
        }
        if (!wfr.derivative().empty()) {
            jacs[block] += dwfr_diag * wfr.derivative()[block];
        }
        if (!gfr.derivative().empty()) {
            jacs[block] += dgfr_diag * gfr.derivative()[block];
        }
        if (!alq.derivative().empty()) {
            jacs[block] += dalq_diag * alq.derivative()[block];
        }
        if (!flo.derivative().empty()) {
            jacs[block] -= dflo_diag * flo.derivative()[block];
        }
    }

    ADB retval = ADB::function(std::move(value), std::move(jacs));
    return retval;
}



double VFPProdProperties::bhp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp_arg,
        const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg, alq);
    return retval.value;
}



double VFPProdProperties::bhpwithflo(const int table_id,
                      const double flo,
                      const double wfr,
                      const double gfr,
                      const double thp,
                      const double alq) const
{
    //Get the table
    const VFPProdTable* table = detail::getTable(m_tables, table_id);

    //First, find the values to interpolate between
    //Value of FLO is negative in OPM for producers, but positive in VFP table
    const auto flo_i = detail::findInterpData(-flo, table->getFloAxis());
    const auto thp_i = detail::findInterpData( thp, table->getTHPAxis()); // assume constant
    const auto wfr_i = detail::findInterpData( wfr, table->getWFRAxis());
    const auto gfr_i = detail::findInterpData( gfr, table->getGFRAxis());
    const auto alq_i = detail::findInterpData( alq, table->getALQAxis()); //assume constant

    detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

    return bhp_val.value;
}




double VFPProdProperties::
       calculateBhpWithTHPTarget(const std::vector<double>& rates1,
                                 const std::vector<double>& rates2,
                                 // bhp2 is the bhp limit, the bhp1 is the middle af bhp2 and cell pressure
                                 const double bhp1,
                                 const double bhp2,
                                 const double dp,
                                 const int table_id,
                                 const double thp,
                                 const double alq) const
{
    assert(table_id > 0);

    const VFPProdTable* table = detail::getTable(m_tables, table_id);

    const int Water = BlackoilPhases::Aqua;
    const int Oil = BlackoilPhases::Liquid;
    const int Gas = BlackoilPhases::Vapour;

    // FLO is the rate
    const double aqua1 = rates1[Water];
    const double liquid1 = rates1[Oil];
    const double vapour1 = rates1[Gas];
    const double flo_rate1 = detail::getFlo(aqua1, liquid1, vapour1, table->getFloType());

    const double aqua2 = rates2[Water];
    const double liquid2 = rates2[Oil];
    const double vapour2 = rates2[Gas];
    const double flo_rate2 = detail::getFlo(aqua2, liquid2, vapour2, table->getFloType());

    /* std::cout << "flo_rate1 is " << flo_rate1 << std::endl;
    std::cout << " bhp1 is " << bhp1 << std::endl;

    std::cout << "flo_rate2 is " << flo_rate2 << std::endl;
    std::cout << " bhp2 is " << bhp2 << std::endl; */

    const double wfr = detail::getWFR(aqua1, liquid1, vapour1, table);
    const double gfr = detail::getGFR(aqua1, liquid1, vapour1, table);

    const int sample_number = 1000;
    std::vector<double> bhp_samples(sample_number);
    std::vector<double> rate_samples(sample_number, 0.);

    // rate sampling interval
    const double rate_iterval = flo_rate2 / (sample_number - 1);
    for (int i = 0; i < sample_number; ++i) {
        rate_samples[i] =  i * rate_iterval;
    }

    // based on all the rate samples, let us calculate the bhp_samples
    for (int i = 0; i < sample_number; ++i) {
        bhp_samples[i] = bhpwithflo(table_id, rate_samples[i], wfr, gfr, thp, alq) - dp;
    }

    // std::cout << " the rate and bhp samples " << std::endl;
    // for (int i = 0; i < sample_number; ++i) {
    //     std::cout << rate_samples[i] << " " << bhp_samples[i] << std::endl;
    // }

    // interface needs to be re-designed
    double return_bhp = 0.;
    const bool found = detail::findIntersectionForBhp(rate_samples, bhp_samples, flo_rate1, flo_rate2, bhp1, bhp2, return_bhp);

    if (!found) {
        std::cout << " COULD NOT find an Intersection point, the well might need to be closed " << std::endl;
        std::cin.ignore();
    }

    // TODO: some sanity check for the obtained return_bhp value?
    return return_bhp;
}


void VFPProdProperties::
calculateBhpWithTHPTarget(const std::vector<double>& ipr_a,
                          const std::vector<double>& ipr_b,
                          const double bhp_limit,
                          const double thp_table_id,
                          const double thp_limit,
                          const double alq,
                          const double dp,
                          bool& obtain_solution_with_thp_limit,
                          bool& violate_bhp_limit_with_thp_limit) const
{
    // let us find a safe highest bhp to use based on IPR

    // begin with a big value
    double bhp_safe_limit = 1.e100;
    for (size_t i = 0; i < ipr_a.size(); ++i) {
        if (ipr_b[i] == 0.) continue;

        const double bhp = ipr_a[i] / ipr_b[i];
        if (bhp < bhp_safe_limit) {
            bhp_safe_limit = bhp;
        }
    }


    std::cout << " bhp_safe_limit " << bhp_safe_limit <<  std::endl;

    const double bhp_middle = (bhp_limit + bhp_safe_limit) / 2.0;

    // FLO is the rate
    // The two points correspond to the bhp values of bhp_limit and bhp_safe_limit

    // for producers, the rates are negative actually
    std::vector<double> rates_bhp_limit(ipr_a.size());
    std::vector<double> rates_bhp_middle(ipr_a.size());
    for (size_t i = 0; i < rates_bhp_limit.size(); ++i) {
        rates_bhp_limit[i] = bhp_limit * ipr_b[i] - ipr_a[i];
        rates_bhp_middle[i] = bhp_middle * ipr_b[i] - ipr_a[i];
    }

    // TODO: we need to be careful that there is nothings wrong related to the indices here
    const int Water = BlackoilPhases::Aqua;
    const int Oil = BlackoilPhases::Liquid;
    const int Gas = BlackoilPhases::Vapour;

    const VFPProdTable* table = detail::getTable(m_tables, thp_table_id);
    // FLO is the rate
    /* const double aqua_bhp_limit = rates_bhp_limit[Water];
    const double liquid_bhp_limit = rates_bhp_limit[Oil]; */
    const double aqua_bhp_limit = rates_bhp_limit[Oil];
    const double liquid_bhp_limit = rates_bhp_limit[Water];
    const double vapour_bhp_limit = rates_bhp_limit[Gas];
    const double flo_bhp_limit = detail::getFlo(aqua_bhp_limit, liquid_bhp_limit, vapour_bhp_limit, table->getFloType() );

    std::cout << " aqua_bhp_limit " << aqua_bhp_limit << " liquid_bhp_limit " << liquid_bhp_limit
              << " vapour_bhp_limit " << vapour_bhp_limit << " flo " << flo_bhp_limit << std::endl;

    /* const double aqua_bhp_middle = rates_bhp_middle[Water];
    const double liquid_bhp_middle = rates_bhp_middle[Oil]; */
    const double aqua_bhp_middle = rates_bhp_middle[Oil];
    const double liquid_bhp_middle = rates_bhp_middle[Water];
    const double vapour_bhp_middle = rates_bhp_middle[Gas];
    const double flo_bhp_middle = detail::getFlo(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table->getFloType() );

    std::cout << " aqua_bhp_middle " << aqua_bhp_middle << " liquid_bhp_middle " << liquid_bhp_middle
              << " vapour_bhp_middle " << vapour_bhp_middle << " flo " << flo_bhp_middle << std::endl;

    const double wfr = detail::getWFR(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table);
    const double gfr = detail::getGFR(aqua_bhp_middle, liquid_bhp_middle, vapour_bhp_middle, table);

    std::cout << " wfr " << wfr << " gfr " << gfr << std::endl;

    std::cout << " outputting the flo_sample points based on the points from the table " << std::endl;
    std::vector<double> flo_samples = table->getFloAxis();

    if (flo_samples[0] > 0.) {
        flo_samples.insert(flo_samples.begin(), 0.);
    }

    if (flo_samples.back() < std::abs(flo_bhp_limit)) {
        flo_samples.push_back(std::abs(flo_bhp_limit));
    }

    for (double& value : flo_samples) {
        value = -value;
    }

    std::vector<double> bhp_flo_samples(flo_samples.size());
    for (int i = 0; i < flo_samples.size(); ++i) {
        bhp_flo_samples[i] = bhpwithflo(thp_table_id, flo_samples[i], wfr, gfr, thp_limit, alq) - dp;
    }

    std::cout << " the rate and bhp samples based on the points in the table " << std::endl;
    for (int i = 0; i < flo_samples.size(); ++i) {
        std::cout << flo_samples[i] << " " << bhp_flo_samples[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << " the rate and bhp from the IPR line " << std::endl;
    std::cout << flo_bhp_middle << " " << bhp_middle << std::endl
               << flo_bhp_limit << " " << bhp_limit << std::endl;

    double return_bhp = 0.;
    obtain_solution_with_thp_limit = detail::findIntersectionForBhp(flo_samples, bhp_flo_samples, flo_bhp_middle, flo_bhp_limit,
                                                      bhp_middle, bhp_limit, return_bhp);

    if (obtain_solution_with_thp_limit) {
        violate_bhp_limit_with_thp_limit = (return_bhp < bhp_limit);
    } else {
        std::cout << " COULD NOT find an Intersection point, the well might need to be closed " << std::endl;
    }
}


double VFPProdProperties::thp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& bhp_arg,
        const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);
    const VFPProdTable::array_type& data = table->getTable();

    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = detail::getWFR(aqua, liquid, vapour, table);
    double gfr = detail::getGFR(aqua, liquid, vapour, table);
    // std::cout << " flo " << flo << " wfr " << wfr << " gfr " << gfr << std::endl;

    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     * Recall that flo is negative in Opm, so switch the sign
     */
    auto flo_i = detail::findInterpData(-flo, table->getFloAxis());
    auto wfr_i = detail::findInterpData( wfr, table->getWFRAxis());
    auto gfr_i = detail::findInterpData( gfr, table->getGFRAxis());
    auto alq_i = detail::findInterpData( alq, table->getALQAxis());
    // std::cout << " flo_i " << flo_i << " wfr_i " << wfr_i << " gfr_i " << gfr_i << " alq_i " << alq_i << std::endl;
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(data, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}






const VFPProdTable* VFPProdProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}







}

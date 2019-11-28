#ifdef ACL_DEVICE

#include <iostream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <iterator>
#include <string>

#include "algorithm.h"
#include "EncodingConstants.h"
#include "Condition.h"
#include "ConditionalTable.h"
#include "IndexTable.h"
#include "specific_types.h"

#include <Kernels/Interface.hpp>
#include <yateto.h>

using namespace seissol::initializers::binning;
using namespace seissol;


void check_key(conditional_table_t &table, ConditionalKey &key) {
    if (table.find(key) != table.end()) {
        // found a dublicated key. Hashing is not perfect
        throw std::string("ERROR::BINNING:: a table key conflict found. Problems with hashing.");
    }
}


void seissol::initializers::binning::test(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {

    kernels::LocalData::Loader local_data_loader;
    local_data_loader.load(handler, layer);
    LayerContainer &container = layer.getLayerContainer();
    auto &table = container.get_table_reference_to_init();

    // allocate counters for the main analysis
    size_t idofs_address_counter = 0;        // NOTE: reset at the very end
    std::unordered_map<index_t, index_t> idfos_address_registry{};
    size_t derivatives_address_counter = 0;  // NOTE: reset at the very end
    int distance = 0;

    // *************************************** time integral ***************************************

    // --------------------------------------- primary analysis::start ---------------------------------------
    // perform a primary analysis and find cells which have their own derivatives,
    // and cells which require to provide scratch memory for their derivatives computations
    std::vector<index_t> cells_without_derivatives{};
    std::vector<index_t> cells_with_derivatives{};

    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto data = local_data_loader.entry(cell);
        bool l_DerivativesProvided = ((data.cellInformation.ltsSetup >> 9) % 2) == 1;

        if (l_DerivativesProvided)
            cells_with_derivatives.push_back(cell);
        else
            cells_without_derivatives.push_back(cell);
    }

    // compute the check-sum
    index_t bin_check_sum = cells_with_derivatives.size() + cells_without_derivatives.size();
    assert(bin_check_sum == layer.getNumberOfCells() && "ERROR::BINNING::TIME_INTEGRATION::a bin check sum is different.");
    // if we pass the check above, it means that we prepare all cells for puting into appropriate bins
    // --------------------------------------- primary analysis::end ---------------------------------------

    // bin cells without their own derivatives buffers
    if (!cells_without_derivatives.empty()) {
        std::vector<index_t> dofs_indices{};
        std::vector<index_t> start_indices{};
        std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part
        std::vector<index_t> dQ_indices{}; //relative to the tmp buffer, namely: derivative part

      // debugging
      std::vector<index_t> elements_ids{};

        index_t base_cell_id = cells_without_derivatives[0];
        auto base_data = local_data_loader.entry(base_cell_id);

        for (auto cell: cells_without_derivatives) {
            auto current_data = local_data_loader.entry(cell);

            // dofs
            distance = std::distance(static_cast<real *>(base_data.dofs), static_cast<real *>(current_data.dofs));
            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            dofs_indices.push_back(distance);

            // idofs
            idofs_indices.push_back(idofs_address_counter);
            idfos_address_registry[cell] = idofs_address_counter;
            idofs_address_counter += tensor::I::size();

            // stars
            distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                                     static_cast<real *>(current_data.localIntegration.starMatrices[0]));

            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            start_indices.push_back(distance);

            // dQ; NOTE: the first derivatives is dofs itself; recording the address of the second derivative
            dQ_indices.push_back(derivatives_address_counter);

            // NOTE: the size is bigger than needed, namely: we take into account the first derivative
            derivatives_address_counter += yateto::computeFamilySize<tensor::dQ>();

            // debugging
            elements_ids.push_back(cell);
        }


       // compute a key for this particular condition
        ConditionalKey key(*KernelNames::time, *TimeComputationKind::without_derivatives);
        check_key(table, key);

        // record indices to the table
        table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
        table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
        table[key].variable_indices[*VariableID::derivatives] = new BasicIndices(dQ_indices);

        //debugging
        table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

        table[key].set_not_empty_flag();
    }


    //------------------------------- cells with derivatives ------------------------------------
    // bin cells with their own derivatives buffers
    if (!cells_with_derivatives.empty()) {

        std::vector<index_t> dofs_indices{};
        std::vector<index_t> start_indices{};
        std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part
        std::vector<index_t> dQ_indices{};

        // debugging
        std::vector<index_t> elements_ids{};

        index_t base_cell_id = cells_with_derivatives[0];
        auto base_data = local_data_loader.entry(base_cell_id);

        real ** derivatives = layer.var(handler.derivatives);
        real *base_derivative_address = layer.var(handler.derivatives)[base_cell_id];

        for (auto cell: cells_with_derivatives) {
            auto current_data = local_data_loader.entry(cell);

            // dofs
            distance = std::distance(static_cast<real *>(base_data.dofs), static_cast<real *>(current_data.dofs));
            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            dofs_indices.push_back(distance);

            // idofs
            idofs_indices.push_back(idofs_address_counter);
            idfos_address_registry[cell] = idofs_address_counter;
            idofs_address_counter += tensor::I::size();

            // stars
            distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                                     static_cast<real *>(current_data.localIntegration.starMatrices[0]));

            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            start_indices.push_back(distance);

            // dQ
            distance = std::distance(base_derivative_address, derivatives[cell]);
            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            dQ_indices.push_back(distance);

            // debugging
            elements_ids.push_back(cell);
        }

        ConditionalKey key(*KernelNames::time, *TimeComputationKind::with_derivatives);
        check_key(table, key);

        table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
        table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
        table[key].variable_indices[*VariableID::derivatives] = new RelativeIndices(dQ_indices, base_cell_id);

        //debugging
        table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

        table[key].set_not_empty_flag();
    }


    // *************************************** volume integral ***************************************
    // bin cells for the volume integral
    // NOTE: some clusters can have no elements inside
    if (layer.getNumberOfCells()) {
        std::vector<index_t> idofs_indices{};
        std::vector<index_t> dofs_indices{};
        std::vector<index_t> start_indices{};


      // debugging
        std::vector<index_t> elements_ids{};

        index_t base_cell_id = 0;
        auto base_data = local_data_loader.entry(base_cell_id);

        for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
            auto current_data = local_data_loader.entry(cell);

            // idofs
            idofs_indices.push_back(idfos_address_registry[cell]);

            // dofs
            distance = std::distance(static_cast<real *>(base_data.dofs),
                                     static_cast<real *>(current_data.dofs));
            assert((distance >= 0) && "ERROR::BINNING::VOLUME_INTEGRATION:: negative distance between elements detected");
            dofs_indices.push_back(distance);

            // stars
            distance = std::distance(static_cast<real *>(base_data.localIntegration.starMatrices[0]),
                                     static_cast<real *>(current_data.localIntegration.starMatrices[0]));

            assert((distance >= 0) && "ERROR::BINNING::VOLUME_INTEGRATION:: negative distance between elements detected");
            start_indices.push_back(distance);

            // debugging
            elements_ids.push_back(cell);
        }

        ConditionalKey key(*KernelNames::volume);
        check_key(table, key);
        table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
        table[key].variable_indices[*VariableID::start] = new RelativeIndices(start_indices, base_cell_id);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);

        //debugging
        table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

        table[key].set_not_empty_flag();
    }

    // *************************************** local flux integral ***************************************

    {
        for (unsigned face = 0; face < 4; ++face) {

            std::vector<index_t> idofs_indices{};
            std::vector<index_t> dofs_indices{};
            std::vector<index_t> AplusT_indices{};

            // debugging
            std::vector<index_t> elements_ids{};

            index_t base_cell_id = 0;
            auto base_data = local_data_loader.entry(base_cell_id);

            for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {

                auto current_data = local_data_loader.entry(cell);

                // no element local contribution in the case of dynamic rupture boundary conditions
                if(current_data.cellInformation.faceTypes[face] != dynamicRupture) {

                    // idofs
                    idofs_indices.push_back(idfos_address_registry[cell]);

                    // dofs
                    distance = std::distance(static_cast<real *>(base_data.dofs),
                                             static_cast<real *>(current_data.dofs));
                    assert((distance >= 0) && "ERROR::BINNING::LOCAL_FLUX:: negative distance between elements detected");
                    dofs_indices.push_back(distance);

                    // AplusT
                    distance = std::distance(static_cast<real *>(base_data.localIntegration.nApNm1[face]),
                                             static_cast<real *>(current_data.localIntegration.nApNm1[face]));
                    assert((distance >= 0) && "ERROR::BINNING::LOCAL_FLUX:: negative distance between elements detected");
                    AplusT_indices.push_back(distance);

                    elements_ids.push_back(cell);
                }
            }

            // NOTE: we can check any container, but we must check that a set is not empty!
            if (!dofs_indices.empty()) {
              ConditionalKey key(*KernelNames::local_flux, !FaceKinds::dynamicRupture, face);
              check_key(table, key);
              table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
              table[key].variable_indices[*VariableID::dofs] = new RelativeIndices(dofs_indices, base_cell_id);
              table[key].variable_indices[*VariableID::AplusT] = new RelativeIndices(AplusT_indices, base_cell_id);

              //debugging
              table[key].variable_indices[*VariableID::elements_ids] = new BasicIndices(elements_ids);

              table[key].set_not_empty_flag();
            }
        }
    }

    // *************************************** end of the local integral ***************************************
    // --------------------------------------- primary analysis::start ---------------------------------------

    // perform a primary analysis and find which cells have gts_buffers,
    // and which cells have lts buffers
    std::vector<index_t> cells_with_lts_buffers{};
    std::vector<index_t> cells_with_gts_buffers{};

    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        auto data = local_data_loader.entry(cell);
        bool l_BuffersProvided = ((data.cellInformation.ltsSetup >> 8) % 2) == 1;
        bool l_LtsBuffers = ((data.cellInformation.ltsSetup >> 10) % 2) == 1;

        if (l_BuffersProvided) {
            if (l_LtsBuffers)
                cells_with_lts_buffers.push_back(cell);
            else
                cells_with_gts_buffers.push_back(cell);
        }
    }
    // --------------------------------------- primary analysis::end ---------------------------------------

    //  bin cells with lts buffers: must be either overridden or accumulated, namely: run-time behaviour
    if (!cells_with_lts_buffers.empty()) {

        std::vector<index_t> buffers_indices{};
        std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part

        index_t base_cell_id = cells_with_lts_buffers[0];
        real ** buffers = layer.var(handler.buffers);
        real *base_buffer_address = layer.var(handler.buffers)[base_cell_id];

        for (auto cell: cells_with_lts_buffers) {
            auto current_data = local_data_loader.entry(cell);

            distance = std::distance(base_buffer_address, buffers[cell]);
            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            buffers_indices.push_back(distance);

            idofs_indices.push_back(idfos_address_registry[cell]);
        }

        ConditionalKey key(*KernelNames::time, *TimeComputationKind::with_lts_buffers);
        check_key(table, key);

        table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(buffers_indices, base_cell_id);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
        table[key].set_not_empty_flag();
    }

    // bin cells with gts buffers: must be only overridden: static behaviour, i.e. depends on a mesh
    if (!cells_with_gts_buffers.empty()) {

        std::vector<index_t> buffers_indices{};
        std::vector<index_t> idofs_indices{};  // relative to the tmp buffer, namely: idof part

        index_t base_cell_id = cells_with_gts_buffers[0];
        real ** buffers = layer.var(handler.buffers);
        real *base_buffer_address = layer.var(handler.buffers)[base_cell_id];

        for (auto cell: cells_with_gts_buffers) {
            auto current_data = local_data_loader.entry(cell);

            distance = std::distance(base_buffer_address, buffers[cell]);
            assert((distance >= 0) && "ERROR::BINNING::TIME_INTEGRATION:: negative distance between elements detected");
            buffers_indices.push_back(distance);

            idofs_indices.push_back(idfos_address_registry[cell]);
        }

        ConditionalKey key(*KernelNames::time, *TimeComputationKind::with_gts_buffers);
        check_key(table, key);

        table[key].variable_indices[*VariableID::buffers] = new RelativeIndices(buffers_indices, base_cell_id);
        table[key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_indices);
        table[key].set_not_empty_flag();
    }

    /*
    ConditionalKey cells_wd_key(*KernelNames::time,
                                *TimeComputationKind::with_derivatives);
    check_key(table, key);

     table[cells_wd_key].variable_indices[*VariableID::dofs] = new BasicIndices(dofs_wd_indices);

    table[cells_wd_key].variable_indices[*VariableID::start] = new BasicIndices(start_wd_indices);
    table[cells_wd_key].variable_indices[*VariableID::idofs] = new BasicIndices(idofs_wd_indices);
    table[cells_wd_key].variable_indices[*VariableID::derivatives] = new BasicIndices(dQ_wd_indices);
    */

    /*
    ConditionalKey cells_wod_key(*KernelNames::time,
                                 *TimeComputationKind::without_derivatives);
    check_key(table, key);
    table[cells_wod_key].variable_indices[*VariableID::dofs] = dofs_wod_indices;
    table[cells_wod_key].variable_indices[*VariableID::start] = start_wod_indices;
    table[cells_wod_key].variable_indices[*VariableID::idofs] = idofs_wod_indices;
    table[cells_wod_key].variable_indices[*VariableID::derivatives] = dQ_wod_indices;
    */

    /*
    kernels::NeighborData::Loader neighbor_data_loader;
    neighbor_data_loader.load(handler, layer);

    std::vector<std::pair<index_t, index_t>> regular_periodic_bins[*FaceRelations::Count];
    std::vector<std::pair<index_t, index_t>> dynamicRupture_bin;
    std::vector<std::pair<index_t, index_t>> freeSurface_bin;


    // Bin boundary conditions
    for (unsigned l_cell = 0; l_cell < layer.getNumberOfCells(); ++l_cell) {
        auto data = neighbor_data_loader.entry(l_cell);

        for(unsigned int l_face = 0; l_face < 4; l_face++) {

            auto item = std::make_pair(l_cell, l_face);
            switch (data.cellInformation.faceTypes[l_face]) {
                case regular:
                    // Fallthrough intended
                case periodic:
                {
                    //std::cout << l_cell << "," << l_face << ") regular and periodic" << std::endl;

                    // compute face type relation
                    unsigned bin_id = data.cellInformation.faceRelations[l_face][1]
                                      + 3 * data.cellInformation.faceRelations[l_face][0]
                                      + 12 * l_face;

                    regular_periodic_bins[bin_id].push_back(item);
                    break;
                }
                case freeSurface:
                {
                    //std::cout << l_cell << "," << l_face << ") freeSurface" << std::endl;
                    freeSurface_bin.push_back(item);
                    break;
                }
                case dynamicRupture:
                {
                    //std::cout << l_cell << "," << l_face << ") dynamicRupture" << std::endl;
                    dynamicRupture_bin.push_back(item);
                    break;
                }
                default:
                {
                    //std::cout << l_cell << "," << l_face << ") outflow" << std::endl;
                    break;
                }
            }
        }
    }

    // declare a conditional table for Neighbout Kernels
    std::unordered_map< ConditionalKey, IndexTable, ConditionalHash<ConditionalKey> > table;
    ConditionalKey dynamicRupture_key(*KernelNames::neighbor,
                                      *FaceKinds::dynamicRupture,
                                      *FaceRelations::any);
    check_key(table, key);
    table[dynamicRupture_key].element_face_pairs = dynamicRupture_bin;



    // fill the table with free surface
    ConditionalKey freeSurface_key(*KernelNames::neighbor,
                                   *FaceKinds::freeSurface,
                                   *FaceRelations::any);
    check_key(table, key);
    table[freeSurface_key].element_face_pairs = freeSurface_bin;


    // fill the table with regular and perioic
    for (unsigned relation = 0; relation < (*FaceRelations::Count); ++relation) {
        ConditionalKey key(*KernelNames::neighbor,
                           FaceKinds::regular || FaceKinds::periodic,
                           relation);
        check_key(table, key);
        table[key].element_face_pairs = regular_periodic_bins[relation];
    }

    */
    /*
    // DEBUGGING
    std::cout << "Free surface" << std::endl;
    for (auto pair: table[freeSurface_key].element_face_pairs) {
        std::cout << pair.first << "|" << pair.second << std::endl;
    }
    */

}
#endif
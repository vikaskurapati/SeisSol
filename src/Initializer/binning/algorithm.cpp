#include <iostream>
#include <utility>
#include <vector>
#include <unordered_map>

#include "algorithm.h"
#include "EncodingConstants.h"
#include "Condition.h"
#include "ConditionalTable.h"
#include "IndexTable.h"
#include "specific_types.h"

#include <Kernels/Interface.hpp>

using namespace seissol::initializers::binning;
using namespace seissol;

void seissol::initializers::binning::test(seissol::initializers::LTS &handler, seissol::initializers::Layer &layer) {

    kernels::NeighborData::Loader loader;
    loader.load(handler, layer);

    std::vector<std::pair<index_t, index_t>> regular_periodic_bins[*FaceRelations::Count];
    std::vector<std::pair<index_t, index_t>> dynamicRupture_bin;
    std::vector<std::pair<index_t, index_t>> freeSurface_bin;


    // Bin boundary conditions
    for (unsigned l_cell = 0; l_cell < layer.getNumberOfCells(); ++l_cell) {
        auto data = loader.entry(l_cell);

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

    table[dynamicRupture_key].element_face_pairs = dynamicRupture_bin;



    // fill the table with free surface
    ConditionalKey freeSurface_key(*KernelNames::neighbor,
                                   *FaceKinds::freeSurface,
                                   *FaceRelations::any);

    table[freeSurface_key].element_face_pairs = freeSurface_bin;


    // fill the table with regular and perioic
    for (unsigned relation = 0; relation < (*FaceRelations::Count); ++relation) {
        ConditionalKey key(*KernelNames::neighbor,
                           FaceKinds::regular || FaceKinds::periodic,
                           relation);

        table[key].element_face_pairs = regular_periodic_bins[relation];
    }

    layer.getLayerContainer().set_conditional_table(table);
    layer.getLayerContainer().get_conditional_table();
    /*
    // DEBUGGING
    std::cout << "Free surface" << std::endl;
    for (auto pair: table[freeSurface_key].element_face_pairs) {
        std::cout << pair.first << "|" << pair.second << std::endl;
    }
    */
}

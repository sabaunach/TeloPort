/* sequenceQuality.hpp
 * Contains declarations functions for parsing sequence reads and sequence data for sequence quality.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#ifndef SEQUENCEQUALITY_HPP
#define SEQUENCEQUALITY_HPP

#include "core/telomere_core.hpp"

namespace core = telomere_core;

namespace sequenceQuality {

/* Returns the number  of N's in the read.
 */
size_t countNs(core::Read read);

}

#endif // SEQUENCEQUALITY_HPP

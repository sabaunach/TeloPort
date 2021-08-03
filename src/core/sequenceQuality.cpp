/* sequenceQuality.cpp
 * Contains functions for parsing sequence reads and sequence data for sequence quality.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#include "core/telomere_core.hpp"
#include "core/sequenceQuality.hpp"

namespace core = telomere_core;

namespace sequenceQuality {

/* Returns the number  of N's in the read.
 */
size_t countNs(core::Read read) {
	size_t cnt = 0;
	for (char c: read.seq)
		if (c == 'N')
			cnt++;
	return cnt;
}

}  // end namespace sequenceQuality

/* telomereFinder.hpp
 * Contains declarations functions for mining input files for telomeric reads and outputting telomeric reads.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#ifndef TELOMEREFINDER_HPP
#define TELOMEREFINDER_HPP

#include "core/telomere_core.hpp"

namespace core = telomere_core;

namespace telomereFinder {

static unordered_map<bool, core::TelRepeatInfo> telRepeats = {{true, core::TelRepeatInfo("CCCTAA", true)},
	{false, core::TelRepeatInfo("TTAGGG", false)}};

/* Compare actual sequence to the sequence we would expect based on telRepeat and off.
 *   0 (fastest) : see if actual sequence exactly matches generated sequence
 *   1 (fast)    : precompute a match table s.t. for each substitution of a base with another base
 *     in position i, there is a cost c associated with that substitution
 *   2 (slow)    : do a global alignment of the expected sequence with the actual sequence
 *     (gaps should be high cost, and cost of substitution should be adjusted based on the
 *     frequency of that base within the string). This may or may not be better at recognizing
 *     telomeric sequences than mode 1; it allows for gaps to exist, but mode 1 also accounts
 *     better for the actual telomeric repeat sequence.
 *
 */
double scoreSubsequence(const string & actual, size_t start, size_t stop,
                        const core::TelRepeatInfo & telRepeatInfo, int off,
                        int mode);

/* Check sequence to see if it is a telomeric read.
 * First seek a telomere repeat. For forward repeats, this is the last repeat within the interval [start, stop].
 *  For reverse complement repeats, this is the first repeat within the interval [0, sequence.length() - start].
 *  This is because for reverse complement repeats, sequence quality drops drastically toward the end of the read,
 *      so we'd prefer to match earlier in the read.
 *  Also for reverse complement reads, we will try more matches of the telRepeat if they occur, to check for potential
 *      telomeres. This is to account for potentially having a random centrally-occurring telomere repeat in the read.
 *      For now this is limited to
 *
 * Then the read is scored based on the score() function, which compares actual sequence to what should be expected.
 * Returns: score telomeric
 */
double scoreTelomericRead(const core::Read & read, const core::TelRepeatInfo & telRepeatInfo,
                          int mode, int start);

} // end namespace telomereFinder

#endif // TELOMEREFINDER_HPP

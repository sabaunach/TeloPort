/* junctionFinder.hpp
 * Declarations for functions for finding the telomere junction location in a telomeric read.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#ifndef JUNCTIONFINDER_HPP
#define JUNCTIONFINDER_HPP

#include "core/telomere_core.hpp"

namespace core = telomere_core;

namespace junctionFinder {

static unordered_map<bool, core::TelRepeatInfo> telRepeats = {{true, core::TelRepeatInfo("CCCTAA", true)},
	{false, core::TelRepeatInfo("TTAGGG", false)}};
static int windowLen; // window length for first pass
static int stepLen; // step length for first pass
static double cutoff; // cutoff for first pass
static int shortWindowLen; // cutoff for second pass; to match the bad part.
static double shortCutoff; // cutoff for second pass; to match the bad part.
static int startTries; // amount of times to try skipping when first window fails
static int startSkipLen; // amount of bp to skip and try looking for a window if first window fails completely

/* Returns sum of insertions and deletions in first n positions of optimal WF edit path.
 *  Deletions are +1, insertions are -1.
 *  Uses traceback matrix
 */
int indelSum(vector<vector<pair<int, int> > > tb_matrix, int n);

/* Uses a sliding window approach to determine where the telomere junction ends.
 * Does two passes; one with a larger window size, one with a smaller size.
 *  0. Determine offset of initial sequence by matching first 2 repeats, aligning, picking best offset
 *  1. Run W-F algorithm on the actual sequence and expected sequence based on offset
 *		If penalty exceeds cutoff, go to 2
 *  	Otherwise, adjust offset (based on stepLen and also if there was an insertion or deletion in the step)
 *		Adjust window and repeat
 *	2. Take a window of 2 telomere repeats and read until quality drops off
 *		When quality drops off, look for next place where offset is 0 and return as junction
 *
 *	Returns 0 if no telomeric sequence was matched
 *  Returns INF ( > seq.size()) if the whole sequence was telomeric
 */
size_t findTelomereJunction(string seq, core::TelRepeatInfo telRepeatInfo,
                            size_t windowLen = windowLen, size_t stepLen = stepLen,
                            double cutoff = cutoff, double shortWindowLen = shortWindowLen,
                            double shortCutoff = shortCutoff, int startTries = startTries,
                            int startSkipLen = startSkipLen);

// includes information for last window matched for first and second pass
size_t findTelomereJunction_i(string seq, core::TelRepeatInfo telRepeatInfo,
                              size_t & firstWindowStart, size_t & secondWindowStart,
                              size_t windowLen = windowLen, size_t stepLen = stepLen,
                              double cutoff = cutoff, double shortWindowLen = shortWindowLen,
                              double shortCutoff = shortCutoff, int startTries = startTries,
                              int startSkipLen = startSkipLen);

} // end namespace junctionFinder

#endif // JUNCTIONFINDER_HPP

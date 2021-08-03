/* telomereFinder.cpp
 * Functions for mining input files for telomeric reads and outputting telomeric reads.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#include "telomereFinder/telomereFinder.hpp"

namespace telomereFinder {

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
                        int mode) {

	const string & telRepeat = telRepeatInfo.telRepeat;
	double penaltySum = 0;
	if (mode == 0 || mode == 1) {
		for (size_t i = start, j = off; i <= stop; ++i, ++j %= telRepeat.length()) {
			if (mode == 0) {
				if (actual[i] != telRepeat[j]) // check for any mismatch
					return numeric_limits<double>::max();
			}
			else if (mode == 1) {
				if (actual[i] != telRepeat[j])
					penaltySum += telRepeatInfo.cost_quick.at(j).at(actual[i]);
			}
		}
	} else if (mode == 2) {
		penaltySum = core::alignment_cost(actual, telRepeatInfo, off);
	}
	return penaltySum;
}

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
                          int mode, int start) {
	const string & sequence = read.seq;
	const string & telRepeat = telRepeatInfo.telRepeat;
	size_t loc;
	double penaltySum;
	if (telRepeatInfo.isForward) {
		vector<size_t> locs = {-telRepeat.length()}; // will wrap with + telRepeat.length()
		loc = -1; // will wrap with + 1
		while (true) {
			loc = sequence.find(telRepeat, loc + 1);
			if (loc == string::npos || loc > (sequence.length() - (size_t)start))
				break;
			locs.push_back(loc);
		}
		if (locs.size() == 1) // no repeats found
			return numeric_limits<double>::max();

		// score each sequence from bp after each repeat to the end of next repeat
		// do it from the beginning, adding each score in. once it looks telomeric, return
		penaltySum = 0;
		double minPenaltySum = numeric_limits<double>::max();
		for (size_t i = 1; i < locs.size(); i++) {
			size_t begin = locs[i - 1] + telRepeat.length();
			size_t end = locs[i] + telRepeat.length() - 1;
			// offset based on latter repeat. + telRepeat.length() * sequence.length() to ensure positive before mod
			size_t off = ((int)locs[i] - (int)begin + telRepeat.length() * sequence.length()) % telRepeat.length();
			penaltySum += scoreSubsequence(sequence, begin, end, telRepeatInfo, off, mode);
			if (end >= (size_t)start
			    && penaltySum / ((double)(end + 1) / (double)(telRepeat.length())) < minPenaltySum)
				// minPenaltySum is relative to # of tel repeats in subsequence
				minPenaltySum = penaltySum / ((double)(end + 1) / (double)(telRepeat.length()));
		}
		penaltySum = minPenaltySum;
	} else {
		vector<size_t> locs;
		loc = -1; // loc is size_t, this will wrap back around with + 1
		while (true) {
			loc = sequence.find(telRepeat, loc + 1);
			if (loc == string::npos || loc > (sequence.length() - (size_t)start))
				break;
			locs.push_back(loc);
		}
		if (locs.empty())
			return numeric_limits<double>::max();
		locs.push_back(sequence.length());
		// score each sequence from beginning of each repeat to the bp before next repeat
		// do it from the end, adding each score in. once it looks telomeric, return
		penaltySum = 0;
		double minPenaltySum = numeric_limits<double>::max();
		for (size_t i = locs.size() - 1; i > 0; i--) {
			size_t begin = locs[i - 1];
			size_t end = locs[i] - 1;

			penaltySum += scoreSubsequence(sequence, begin, end, telRepeatInfo, 0, mode);
			if (sequence.length() - begin >= (size_t)start
			    && penaltySum / ((double)(sequence.length() - begin + 1) / (double)(telRepeat.length())) < minPenaltySum)
				// minPenaltySum is relative to # of tel repeats in subsequence
				minPenaltySum = penaltySum / ((double)(sequence.length() - begin + 1) / (double)(telRepeat.length()));
		}
		penaltySum = minPenaltySum;
	}
	return penaltySum;
}

} // end namespace telomereFinder

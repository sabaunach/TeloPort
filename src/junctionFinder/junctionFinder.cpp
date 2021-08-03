/* junctionFinder.cpp
 * Functions for finding the telomere junction location in a telomeric read.
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#include "junctionFinder/junctionFinder.hpp"

namespace junctionFinder {

/* Returns sum of insertions and deletions in first n positions of optimal WF edit path.
 *  Deletions are +1, insertions are -1.
 *  Uses traceback matrix
 */
int indelSum(vector<vector<pair<int, int> > > tb_matrix, int n) {
	int i = tb_matrix.size() - 1, j = tb_matrix[i].size() - 1;
	vector<int> indels;
	while (true) {
		if (tb_matrix[i][j] == pair<int, int>(-1,  -1))
			break;
		else if (tb_matrix[i][j] == pair<int, int>(i - 1, j - 1))
			indels.push_back(0);
		else if (tb_matrix[i][j] == pair<int, int>(i - 1, j))
			indels.push_back(1);
		else if (tb_matrix[i][j] == pair<int, int>(i, j - 1))
			indels.push_back(-1);
		int itemp = i;
		i = tb_matrix[i][j].first;
		j = tb_matrix[itemp][j].second;
	}
	while (i-- > 0)
		indels.push_back(1);
	while (j-- > 0)
		indels.push_back(-1);

	int sum = 0;
	for (int k = 0; k < n; k++) {
		sum += indels[indels.size() - k - 1];
	}
	return sum;
}

size_t findTelomereJunction(string seq, core::TelRepeatInfo telRepeatInfo,
                            size_t windowLen, size_t stepLen, double cutoff,
                            double shortWindowLen, double shortCutoff,
                            int startTries, int startSkipLen) {
	size_t firstWindowStart, secondWindowStart;
	return findTelomereJunction_i(seq, telRepeatInfo, firstWindowStart, secondWindowStart,
	                              windowLen, stepLen, cutoff, shortWindowLen, shortCutoff,
	                              startTries, startSkipLen);
}

/* Uses a sliding window approach to determine where the telomere junction ends.
 * Does two passes; one with a larger window size, one with a smaller size.
 *  0. Determine offset of initial sequence by matching first 2 repeats, aligning, picking best offset.
 *      Do this for start windows beginning out from the sequence and ending at the beginning.
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
size_t findTelomereJunction_i(string seq, core::TelRepeatInfo telRepeatInfo,
                              size_t & firstWindowStart, size_t & secondWindowStart,
                              size_t windowLen, size_t stepLen, double cutoff,
                              double shortWindowLen, double shortCutoff,
                              int startTries, int startSkipLen) {
	if (!telRepeatInfo.isForward) {
		// reverse sequence (so that logic remains similar)
		reverse(seq.begin(), seq.end());
		string revTelRepeat = telRepeatInfo.telRepeat;
		reverse(revTelRepeat.begin(), revTelRepeat.end());
		core::TelRepeatInfo rev = core::TelRepeatInfo(revTelRepeat, false);
		telRepeatInfo = rev;
	}

	string telRepeat = telRepeatInfo.telRepeat;

	// 0
	double minPenalty = numeric_limits<double>::max();
	size_t start = startSkipLen * (startTries - 1);
	size_t off = 0;
	string actual;
	
	// in case we're starting beyond the telomere length
	while (start + telRepeat.length() * 2 - 1 >= seq.size()) {
		start -= startSkipLen;
		startTries--;
	}
	
	for (; minPenalty > .1 && startTries; startTries--, start -= startSkipLen) {
		minPenalty = numeric_limits<double>::max();
		actual = seq.substr(start, telRepeat.length() * 2);
		for (size_t j = 0; j < telRepeat.length(); j++) {
			double curPenalty = core::alignment_cost(actual, telRepeatInfo, j) / 2.0;
			if (curPenalty < minPenalty) {
				minPenalty = curPenalty;
				off = j;
			}
		}
	}
	if (start > seq.length()) // if start went below 0 and wrapped due to size_t
		return 0;
	// 1
	double penalty = -1;
	// penalty is relative to # of tel repeats
	double denom = (double)((double)(windowLen) / (double)telRepeat.length());
	// traceback matrix for tracing path of optimal alignment
	vector<vector<pair<int, int> > > tb_matrix;
	// slide until either quality drops, end of sequence, or no intial repeat found
	while (start + windowLen - 1 < seq.size()) {
		actual = seq.substr(start, windowLen);
		penalty = core::alignment_cost(actual, telRepeatInfo, off, tb_matrix) / denom;
		if (penalty > cutoff)
			break;

		int indelOff = indelSum(tb_matrix, stepLen);
		off = (off + stepLen + indelOff + telRepeat.length()) % telRepeat.length();
		start += stepLen;
	}
	// broke due to hitting end of sequence
	if (penalty <= cutoff)
		return numeric_limits<size_t>::max();

	firstWindowStart = start;

	// broke because penalty was exceeded: slide with short window
	// slide until qualtiy drops
	windowLen = shortWindowLen;
	stepLen = 1;
	denom = (double)((double)(windowLen) / (double)telRepeat.length());
	cutoff = shortCutoff;
	while (start + windowLen - 1 < seq.size()) {
		actual = seq.substr(start, windowLen);
		penalty = core::alignment_cost(actual, telRepeatInfo, off, tb_matrix) / denom;

		if (penalty > cutoff) { // quality dropped off
			break;
		}

		int indelOff = indelSum(tb_matrix, stepLen);
		off = (off + stepLen + indelOff + telRepeat.length()) % telRepeat.length();
		start += stepLen;
	}
	// broke on first read (no telomeric sequence matched)
	if (start == 0)
		return 0;
	secondWindowStart = start;

	// final check; slide until no longer exactly matching telomere repeat
	windowLen = telRepeat.length();
	cutoff = 0.0;
	while (start + windowLen - 1 < seq.size()) {
		actual = seq.substr(start, windowLen);
		penalty = core::alignment_cost(actual, telRepeatInfo, off);

		if (penalty != 0.0) { // does not exactly match
			break;
		}

		start += stepLen;
		off = (off + stepLen + telRepeat.length()) % telRepeat.length();
	}

	// slide until offset = 0
	while (start + windowLen - 1 < seq.size()) {
		actual = seq.substr(start, windowLen);
		core::alignment_cost(actual, telRepeatInfo, off, tb_matrix);

		int indelOff = indelSum(tb_matrix, stepLen);
		off = (off + stepLen + indelOff + telRepeat.length()) % telRepeat.length();
		start += stepLen;

		if (off == 0)
			break;
	}
	// broke due to hitting end of sequence; this really is unlikely to be hit.
	if (start + windowLen - 1 >= seq.size())
		return numeric_limits<size_t>::max();

	if (telRepeatInfo.isForward)
		return start;
	firstWindowStart = seq.size() - firstWindowStart;
	secondWindowStart = seq.size() - secondWindowStart;
	return seq.size() - start;
}

} // end namespace jf

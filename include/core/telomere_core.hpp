/* telomere_core.hpp
 * Declarations of core functions and datastructures for reading sequences from files and working with telomere patterns.
 * Author: Seth Baunach
 * Date: 6/30/2020
 */

#ifndef TELOMERE_CORE_HPP
#define TELOMERE_CORE_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using std::string;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using std::vector;
using std::unordered_map;
using std::numeric_limits;
using std::cout;
using std::cerr;
using std::endl;
using std::to_string;
using std::pair;
using std::exception;
using std::min;
using std::max;

namespace telomere_core {

const double INSBASE_DEF = .5; // default insert base cost
const double INSFACT_DEF = .5; // default insert frequency factor
const double DELBASE_DEF = .8; // default delete base cost
const double DELFACT_DEF = 0; // default delete frequency factor
const double SUBBASE_DEF = .25; // default substitute base cost
const double SUBFACT_DEF = .5; // default substitute frequency factor
const double NONBASE_DEF = 1.5; // default cost for non-occurring characters ('N' included)

/* Datastructure to hold information about a repeat sequence (forward or reverse).
 *  bool isFoward : whether the repeat is in forward or reverse complement orientation
 *  string telRepeat : repeat sequence
 *  vector<unordered_map<char, double> > cost_quick : cost table for quick comparison
 *      of a query sequence to the repeat sequence, based on the location of the
 *      characters and the frequency. Does not account for indels.
 *          cost_quick[loc][char] := cost
 *  unordered_map<string, unordered_map<char, double> > cost_alignment: cost table for
 *      edit distance costs for comapring a query sequence to the repeat sequence, based
 *      on frequency of the characters.
 *          cost_alignment["insert/delete/substitute"][char] := cost
 *  For both tables, the cost is based on the relative frequency of characters, not absolute frequency.
 *      Meaning you do not need to adjust for the length of the telomere repeat.
 */
struct TelRepeatInfo {
	bool isForward;
	string telRepeat;
	vector<unordered_map<char, double> > cost_quick;
	unordered_map<string, unordered_map<char, double> > cost_alignment;

	// Initialize a cost table for alignment based on a base cost and factor of frequency for each type of error.
	void init_cost_alignment(double insBase = INSBASE_DEF, double insFact = INSFACT_DEF,
	                         double delBase = DELBASE_DEF, double delFact = DELFACT_DEF,
	                         double subBase = SUBBASE_DEF, double subFact = SUBFACT_DEF,
	                         double nonBase = NONBASE_DEF);
	// Initialize a cost table for quick string comparison based on the frequency of repeats with default values.
	// Parameters really shouldn't be changed.
	void init_cost_quick();

	// Default constructor; nothing initialized.
	TelRepeatInfo();
	// Constructor with telRepeat sequence. Also initializes cost tables.
	TelRepeatInfo(const string & telRepeat, bool isForward);
};

/* Datastructure for holding basic information about a raw read.
 */
struct Read {
	string id = "";
	string seq = "";
	string phred = "";

	Read();
	Read(string id, string seq);
	Read(string id, string seq, string phred);
};

/* Split read at loc. [0 : loc - 1] <-> [loc : ]
 */
pair<Read, Read> splitRead(Read read, size_t loc);

/* Return a read which is the reverse of read.
 */
Read rev(Read read);

/* Return a read which is the complement of read.
 */
Read comp(Read read);

/* Return a read which is the reverse complement of read.
 */
Read revc(Read read);

/* compute Wagner-Fischer between x and y, return edit distance (bottom right cell)
 * result is not relative to length of sequence but is the absolute edit distance
 * based on costs defined in telRepeatInfo.cost_alignment
 * cost to get from expected to actual (x to y)
 */
double alignment_cost(const string & actual, const TelRepeatInfo & telRepeatInfo, int off,
                      vector<vector<pair<int, int> > > & t);

// overload when t is discarded
double alignment_cost(const string & actual, const TelRepeatInfo & telRepeatInfo, int off);

/* read a sequence in fasta format from f (buffer should be at beginning of seqID line)
 */
bool read_sequence_fasta(ifstream & f, Read & read);

/* write a sequence in fasta format to f (buffer should be at beginning of seqID line)
 */
bool write_sequence_fasta(ofstream & f, Read & read);

/* read a sequence in fastq format from f (buffer should be at beginning of seqID line)
 */
bool read_sequence_fastq(ifstream & f, Read & read);

/* write a sequence in fastq format to f (buffer should be at beginning of seqID line)
 */
bool write_sequence_fastq(ofstream & f, Read & read);

/* Write a sequence read to file f with format fmt.
 */
bool write_sequence(ofstream & f, Read & read, string fmt);

/* Read a sequence from file f with format fmt.
 */
bool read_sequence(ifstream & f, Read & read, string fmt);

/* Read from a file containing sequence information.
 * Store sequences in a vector.
 */
void read_seqFile(ifstream & f, string fmt, vector<Read> & reads);

/* Read from a file containing sequence information.
 * Store sequences in an unordered_map<id, read>
 */
void read_seqFile(ifstream & f, string fmt, unordered_map<string, Read> & reads);

/* Write to out file sequences occurring in both files using sequence information from second file (seq file).
 *  useFirstOrdering: if true, output sequences are ordered based on occurrence in ref file
 *                    if false, based on seq info file.
 *					  default: true
 */
bool write_ids_both_ref_seq(ifstream & ref, ifstream & seq, ofstream & out,
                            string refFmt = "fastq", string seqFmt = "fastq",
                            string outFmt = "fastq", bool useRefOrdering = true);

} // end namespace telomere_core

#endif // TELOMERE_CORE_HPP

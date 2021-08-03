/* telomere_core.cpp
 * Collection of core functions and datastructures for reading sequences from files and working with telomere patterns.
 * Author: Seth Baunach
 * Date: 6/30/2020
 * TODO: add logging functions to core
 */

#include "core/telomere_core.hpp"

namespace telomere_core {

TelRepeatInfo::TelRepeatInfo() {
}

TelRepeatInfo::TelRepeatInfo(const string & telRepeat, bool isForward) {
	this->telRepeat = telRepeat;
	this->isForward = isForward;
	init_cost_quick();
	init_cost_alignment();
}

Read::Read() {
}

Read::Read(string id, string seq) {
	this->id = id;
	this->seq = seq;
}

Read::Read(string id, string seq, string phred) {
	this->id = id;
	this->seq = seq;
	this->phred = phred;
}

/* Return a read which is the reverse of read.
 */
Read rev(Read read) {
	Read newRead = Read(read.id, read.seq, read.phred);
	reverse(newRead.seq.begin(), newRead.seq.end());
	reverse(newRead.phred.begin(), newRead.phred.end());
	return newRead;
}

/* Return a read which is the complement of read.
 */
Read comp(Read read) {
	Read newRead = Read(read.id, read.seq, read.phred);
	for (char & c: newRead.seq) {
		if (c == 'A')
			c = 'T';
		else if (c == 'T')
			c = 'A';
		else if (c == 'C')
			c = 'G';
		else if (c == 'G')
			c = 'C';
	}
	return newRead;
}

/* Return a read which is the reverse complement of read.
 */
Read revc(Read read) {
	return rev(comp(read));
}

/* Split read at loc. [0 : loc - 1] <-> [loc : ]
 */
pair<Read, Read> splitRead(Read read, size_t loc) {
	return {Read(read.id, read.seq.substr(0, loc), read.phred.substr(0, loc)),
	        Read(read.id, read.seq.substr(loc, string::npos), read.phred.substr(loc, string::npos))};
}

// Initialize a cost table for alignment based on a base cost and factor of frequency for each type of error.
void TelRepeatInfo::init_cost_alignment(double insBase, double insFact,
                                        double delBase, double delFact,
                                        double subBase, double subFact,
                                        double nonBase) {

	unordered_map<char, double> freqTable{{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};

	for (const char & c: telRepeat)
		freqTable[c]++;
	for (auto & kv: freqTable)
		kv.second = 1 - (kv.second/telRepeat.length());

	for (auto kv: freqTable) {
		char c = kv.first;
		cost_alignment["insert"][c] = kv.second == 1 ? nonBase : insBase + insFact*kv.second;
		cost_alignment["delete"][c] = delBase + delFact*kv.second;
		cost_alignment["substitute"][c] = kv.second == 1 ? nonBase : subBase + subFact*kv.second;
	}
	for (auto & m: cost_alignment) {
		m.second['N'] = nonBase;
	}
}

// Initialize a cost table for quick string comparison based on the frequency of repeats with default values.
// Parameters really shouldn't be changed.
void TelRepeatInfo::init_cost_quick() {

	unordered_map<char, double> freqTable{{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};

	for (const char & c: telRepeat)
		freqTable[c]++;
	for (auto & kv: freqTable)
		kv.second = 1 - (kv.second/telRepeat.length());

	cost_quick = vector<unordered_map<char, double> >(telRepeat.length(),
	                                                  unordered_map<char, double>(freqTable));

	// lazy algorithm, we don't really expect telRepeat to be huge so this is fine
	for (int i = 0; (size_t)i < telRepeat.length(); i++) {
		for (auto & kv: cost_quick[i]) {
			char base = kv.first;
			double * cost = &(kv.second);
			if (base == telRepeat[i]) {                                     // no substitution penalty, bases are same
				*cost = 0;
				continue;
			}
			if (freqTable[base] == 1.0) {                                     // base DNE in repeat sequence
				*cost = 2;
				continue;
			}
			for (int j = 1; (size_t)j <= (telRepeat.length())/2; j++) {
				if (telRepeat[(i + j) % telRepeat.length()] == base                                                 // wrap right to check for occurrence
				    || telRepeat[((i - j) + telRepeat.length()) % telRepeat.length()] == base) {                                                 // wrap left
					*cost += (double)j / telRepeat.length();                                                             // increase cost once occurrence is found
					break;
				}
			}
		}
	}
	for (auto & m: cost_quick) {
		m['N'] = 1;
	}
}

/* compute Wagner-Fischer between x and y, return edit distance (bottom right cell)
 * result is not relative to length of sequence but is the absolute edit distance
 * based on costs defined in telRepeatInfo.cost_alignment
 * cost to get from expected to actual (x to y)
 */
double alignment_cost(const string & actual, const TelRepeatInfo & telRepeatInfo, int off,
                      vector<vector<pair<int, int> > > & t) {

	string y = actual;
	string x(y.length(), '\0');
	// generate expected
	for (size_t i = 0, j = off; i < x.length(); i++, ++j %= telRepeatInfo.telRepeat.length()) {
		x[i] = telRepeatInfo.telRepeat[j];
	}

	/* WAGNER-FISCHER */

	size_t n = x.size();
	size_t m = y.size();

	// wagner-fischer matrix
	vector<vector<double> > p(n+1, vector<double>(m+1, 0.0));

	// traceback matrix
	t = vector<vector<pair<int, int> > >(n+1, vector<pair<int, int> >(m+1, {-1, -1}));

	// if we were to delete every character of x to get from x to blank string
	for (size_t i = 1; i <= n; i++)
		p[i][0] = p[i - 1][0] + telRepeatInfo.cost_alignment.at("delete").at(x[i - 1]);

	// if we were to insert every character of y to get from blank string to y
	for (size_t j = 1; j <= m; j++)
		p[0][j] = p[0][j - 1] + telRepeatInfo.cost_alignment.at("insert").at(y[j - 1]);

	// compute Wagner-Fischer matrix using x and y
	for (size_t i = 1; i <= n; i++) {
		for (size_t j = 1; j <= m; j++) {
			double cost_sub, cost_ins, cost_del;
			if (x[i - 1] == y[j - 1])
				cost_sub = p[i - 1][j - 1];
			else
				cost_sub = p[i - 1][j - 1] + telRepeatInfo.cost_alignment.at("substitute").at(y[i - 1]);
			cost_ins = p[i][j - 1] + telRepeatInfo.cost_alignment.at("insert").at(y[j - 1]);
			cost_del = p[i - 1][j] + telRepeatInfo.cost_alignment.at("delete").at(x[i - 1]);
			if (cost_sub <= cost_ins && cost_sub <= cost_del) {
				p[i][j] = cost_sub;
				t[i][j] = {i - 1, j - 1};
			} else if (cost_ins <= cost_del) {
				p[i][j] = cost_ins;
				t[i][j] = {i, j - 1};
			} else {
				p[i][j] = cost_del;
				t[i][j] = {i - 1, j};
			}
		}
	}
	/* output stuff
	   cout << "x (expected) : " << x << endl;
	   cout << "y (actual)   : " << y << endl;
	   cout << "wf-matrix: " << endl;
	   for (vector<double> row: p) {
	    for (double i: row) {
	        cout << std::setprecision(3);
	        cout << std::setw(5);
	        cout << i << " ";
	    }
	    cout << endl;
	   }
	   cout << "tb-matrix: " << endl;
	   for (size_t i = 0; i < t.size(); i++) {
	    for (size_t j = 0; j < t[i].size(); j++) {
	        if (t[i][j] == pair<int, int>({-1, -1}))
	            cout << "  ";
	        else if (t[i][j] == pair<int, int>({i - 1, j - 1}))
	            cout << "` ";
	        else if (t[i][j] == pair<int, int>({i, j - 1}))
	            cout << "< ";
	        else if (t[i][j] == pair<int, int>({i - 1, j}))
	            cout << "^ ";
	    }
	    cout << endl;
	   }*/
	return p[n][m];
}

double alignment_cost(const string & actual, const TelRepeatInfo & telRepeatInfo, int off) {
	vector<vector<pair<int, int> > > t;
	return alignment_cost(actual, telRepeatInfo, off, t);
}

bool read_sequence_fasta(ifstream & f, Read & read) {
	return read_sequence(f, read, "fasta");
}

bool write_sequence_fasta(ofstream & f, Read & read) {
	return write_sequence(f, read, "fasta");
}

bool read_sequence_fastq(ifstream & f, Read & read) {
	return read_sequence(f, read, "fastq");
}

bool write_sequence_fastq(ofstream & f, Read & read) {
	return write_sequence(f, read, "fastq");
}

/* Write a sequence read to file f with format fmt.
 */
bool write_sequence(ofstream & f, Read & read, string fmt) {
	if (fmt == "fasta") {
		f << '>' + read.id << "\n"
		  << read.seq << endl;
	} else if (fmt == "fastq") {
		f << '@' + read.id << "\n"
		  << read.seq
		  << "\n+\n"
		  << read.phred << endl;
	} else {
		cerr << "Error: Unrecognized write format: " << fmt;
		return false;
	}
	if (f.fail())
		return false;
	return true;
}

/* Read a sequence from file f with format fmt.
 */
bool read_sequence(ifstream & f, Read & read, string fmt) {
	if (fmt != "fastq" && fmt != "fasta") {
		cerr << "Error: Unrecognized read format: " << fmt;
		return false;
	}
	if (!getline(f, read.id)) {
		return false;
	}
	read.id = read.id.substr(1, read.id.length() - 1); // remove id prepending character
	if (!getline(f, read.seq)) {
		cerr << "Error while reading sequence from input file: error occurred (premature EOF) while reading sequence information from sequence with seqID:\n" << read.id << "\n";
		return false;
	}
	if (fmt == "fastq") {
		if (!getline(f, read.phred)) {
			cerr << "Error while reading sequence from input file: error occurred (premature EOF) while reading break line from sequence with seqID:\n" << read.id << "\n";
			return false;
		}
		if (!getline(f, read.phred)) {
			cerr << "Error while reading sequence from input file: error occurred (premature EOF) while reading PHRED information from sequence with seqID:\n" << read.id << "\n";
			return false;
		}
		if (read.seq.length() != read.phred.length()) {
			cerr << "Error while reading sequence from input file: error occurred (sequence length and PHRED length differ!) with sequence with seqID:\n" << read.id << "\n";
			return false;
		}
	}
	return true;
}

/* Read from a file containing sequence information.
 * Store sequences in a vector.
 */
void read_seqFile(ifstream & f, string fmt, vector<Read> & reads) {
	while (true) {
		Read read;
		if (!read_sequence(f, read, fmt))
			break;
		reads.push_back(read);
	}
}

/* Read from a file containing sequence information.
 * Store sequences in an unordered_map<id, read>
 */
void read_seqFile(ifstream & f, string fmt, unordered_map<string, Read> & reads) {
	while (true) {
		Read read;
		if (!read_sequence(f, read, fmt))
			break;
		reads.insert({read.id, read});
	}
}

/* Write to out file sequences occurring in both files using sequence information from second file (seq file).
 *  useFirstOrdering: if true, output sequences are ordered based on occurrence in ref file
 *                    if false, based on seq info file.
 *					  default: true
 */
bool write_ids_both_ref_seq(ifstream & ref, ifstream & seq, ofstream & out,
                            string refFmt, string seqFmt, string outFmt,
                            bool useRefOrdering) {
	vector<Read> order_vec;
	unordered_map<string, Read> other_map;

	if (useRefOrdering) {
		read_seqFile(ref, refFmt, order_vec);
		read_seqFile(seq, seqFmt, other_map);
	} else {
		read_seqFile(seq, seqFmt, order_vec);
		read_seqFile(ref, refFmt, other_map);
	}

	for (Read & read: order_vec) {
		auto match = other_map.find(read.id);
		if (match != other_map.end()) {
			if (useRefOrdering) { // seq is stored in other_map
				if (write_sequence(out, match->second, outFmt))
					return false; // error occurred
			} else { // seq is stored in order_vec
				if (write_sequence(out, read, outFmt))
					return false; // error occurred
			}
		}
	}
	return true;
}

} // end namespace

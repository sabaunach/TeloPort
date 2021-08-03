/* junctionFinder_driver.cpp
 * Functions for finding the telomere junction location in a telomeric read.
 * Author: Seth Baunach
 * Date: 6/30/2020
 */

#include "junctionFinder/junctionFinder.hpp"
#include "core/sequenceQuality.hpp"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace gg = boost::gregorian;
namespace pt = boost::posix_time;

namespace jf = junctionFinder;
namespace sq = sequenceQuality;

ofstream log_stats;
ofstream log_details;

bool o_opt_revc = false,
     o_opt_sort = true,
     o_opt_splitDir = true,
     o_opt_splitJunc = true,
     o_opt_pretty = false,
     o_opt_fullMatches = false;
int repeatsThrow = 4; // throw out reads with number of repeats in subTelSeq >= repeatsThrow

struct JunctionRead {
	core::Read read;
	bool isForward;
	size_t junc_loc;
	size_t firstWindowStart;
	size_t secondWindowStart;

	JunctionRead();

	JunctionRead(core::Read read, bool isForward, size_t junc_loc) {
		this->read = read;
		this->isForward = isForward;
		this->junc_loc = junc_loc;
	}

	JunctionRead(core::Read read, bool isForward, size_t junc_loc,
	             size_t firstWindowStart, size_t secondWindowStart) {
		this->read = read;
		this->isForward = isForward;
		this->junc_loc = junc_loc;
		this->firstWindowStart = firstWindowStart;
		this->secondWindowStart = secondWindowStart;
	}
};

void writePretty(JunctionRead r, ostream & out) {
	size_t i = 0;
	if (r.isForward) {
		out << "\033[4;31m";
		for (; i < r.firstWindowStart; i++)
			out << r.read.seq[i];
		out << "\033[4;35m";
		for (; i < r.secondWindowStart; i++)
			out << r.read.seq[i];
		out << "\033[4;36m";
		for (; i < r.junc_loc; i++)
			out << r.read.seq[i];
		out << "\033[24;36m";
		for (; i < r.secondWindowStart + jf::shortWindowLen; i++)
			out << r.read.seq[i];
		out << "\033[24;35m";
		for (; i < r.firstWindowStart + jf::windowLen; i++)
			out << r.read.seq[i];
		out << "\033[0m";
		for (; i < r.read.seq.length(); i++)
			out << r.read.seq[i];
		out << endl;
	} else {
		for (; (int)i < (int)(r.firstWindowStart - jf::windowLen); i++)
			out << r.read.seq[i];
		out << "\033[24;35m";
		for (; (int)i < (int)(r.secondWindowStart - jf::shortWindowLen); i++)
			out << r.read.seq[i];
		out << "\033[24;36m";
		for (; i < r.junc_loc; i++)
			out << r.read.seq[i];
		out << "\033[4;36m";
		for (; i < r.secondWindowStart; i++)
			out << r.read.seq[i];
		out << "\033[4;35m";
		for (; i < r.firstWindowStart; i++)
			out << r.read.seq[i];
		out << "\033[4;31m";
		for (; i < r.read.seq.length(); i++)
			out << r.read.seq[i];
		out << "\033[0m" << endl;
	}
}

/* Simple function to return the number of exact telomere repeats in seq.
 */
int countTelRepeats(string seq, core::TelRepeatInfo telRepeatInfo) {
	int repFind = 0;
	size_t loc = seq.find(telRepeatInfo.telRepeat, 0);
	for (; loc != string::npos; repFind++) {
		loc = seq.find(telRepeatInfo.telRepeat, loc + 1);
	}
	return repFind;
}

bool read_jf(ifstream & f_in, ofstream * telSeq_out_f, ofstream * subTelSeq_out_f,
             ofstream * telSeq_out_r, ofstream * subTelSeq_out_r) {
	int seqsread = 0, nonMatches = 0, fullMatches = 0, throwOut = 0;
	core::Read read;
	vector<JunctionRead> juncReads;
	while(true) {
		seqsread++;
		core::read_sequence_fastq(f_in, read);

		if (f_in.eof())
			break;

		size_t forward_firstWindowStart, forward_secondWindowStart,
		       reverse_firstWindowStart, reverse_secondWindowStart;
		size_t forward_res = jf::findTelomereJunction_i(read.seq, jf::telRepeats[true], forward_firstWindowStart, forward_secondWindowStart);
		size_t reverse_res = jf::findTelomereJunction_i(read.seq, jf::telRepeats[false], reverse_firstWindowStart, reverse_secondWindowStart);
		if (forward_res == numeric_limits<size_t>::max()
		    || reverse_res == numeric_limits<size_t>::max()) {
			// one matched whole telomere
			fullMatches++;
			if (o_opt_fullMatches) {
				if (forward_res == numeric_limits<size_t>::max()) {
					juncReads.push_back(JunctionRead(read, true, read.seq.length(), read.seq.length(), read.seq.length()));
				} else {
					juncReads.push_back(JunctionRead(read, false, 0, 0, 0));
				}
			}

		} else if (reverse_res == 0 && forward_res == 0) {
			// no match
			nonMatches++;
			continue;

		} else if (reverse_res == 0 || forward_res > read.seq.length() - reverse_res) {
			// forward match
			juncReads.push_back(JunctionRead(read, true, forward_res, forward_firstWindowStart, forward_secondWindowStart));

		} else {
			// reverse match
			if (o_opt_revc)
				juncReads.push_back(JunctionRead(core::revc(read), true, read.seq.length() - reverse_res,
				                                 reverse_firstWindowStart, reverse_secondWindowStart));
			else
				juncReads.push_back(JunctionRead(read, false, reverse_res, reverse_firstWindowStart, reverse_secondWindowStart));

		}
		seqsread++;
		writePretty(juncReads.back(), cout);
	}

	// write output
	if (o_opt_sort) {
		sort(juncReads.begin(), juncReads.end(),
		     [](const JunctionRead & a, const JunctionRead & b) -> bool {
			return a.junc_loc > b.junc_loc;
		});
	}
	if (o_opt_splitDir) {
		stable_partition(juncReads.begin(), juncReads.end(),
		                 [](const JunctionRead & a) -> bool {
			return a.isForward;
		});
	}
	for (JunctionRead r: juncReads) {

		pair<core::Read, core::Read> splitReads;
		splitReads = core::splitRead(r.read, r.junc_loc);

		if (r.isForward) {
			int repFind = countTelRepeats(splitReads.second.seq, jf::telRepeats[true]);
			if (repeatsThrow == 0 || repFind < repeatsThrow) {
				if (!o_opt_splitJunc) {
					core::write_sequence_fastq(*telSeq_out_f, r.read);
				} else {
					core::write_sequence_fastq(*telSeq_out_f, splitReads.first);
					core::write_sequence_fastq(*subTelSeq_out_f, splitReads.second);
				}
			} else {
				throwOut++;
				continue;
			}
		} else {
			int repFind = countTelRepeats(splitReads.first.seq, jf::telRepeats[false]);
			if (repeatsThrow == 0 || repFind < repeatsThrow) {
				if (!o_opt_splitJunc) {
					core::write_sequence_fastq(*telSeq_out_r, r.read);
				} else {
					core::write_sequence_fastq(*telSeq_out_r, splitReads.second);
					core::write_sequence_fastq(*subTelSeq_out_r, splitReads.first);
				}
			} else {
				throwOut++;
				continue;
			}
		}
		if (o_opt_pretty) {
			writePretty(r, cout);
		}

	}
	return true;
}

bool process_options(int argc, char** argv,
                     ifstream & f_in,
                     ofstream & telSeq_f_out, ofstream & subTelSeq_f_out,
                     ofstream & telSeq_r_out, ofstream & subTelSeq_r_out) {
	try {
		string i_filename;
		string out_dirname;

		po::options_description desc("Options");
		desc.add_options()
		    ("help,h", "produce this message")
		    ("in,i", po::value<string>(&i_filename), "specify input from a fastq file.")
		    ("out,o", po::value<string>(&out_dirname)->default_value("junctionFinder_out"), "specify output directory.")
		    ("startTries,T", po::value<int>(&jf::startTries)->default_value(10), "specify how many times to try a start window.")
		    ("startSkipLen,l", po::value<int>(&jf::startSkipLen)->default_value(12), "specify how many bp to skip if the first start fails.")
		    ("windowLen,W", po::value<int>(&jf::windowLen)->default_value(60), "specify window length on first pass.")
		    ("stepLen,L", po::value<int>(&jf::stepLen)->default_value(1), "specify step length on first pass.")
		    ("cutoff,C", po::value<double>(&jf::cutoff)->default_value(.5), "specify cutoff penalty for first pass.")
		    ("shortWindowLen,w", po::value<int>(&jf::shortWindowLen)->default_value(15), "specify window size for second pass.")
		    ("shortCutoff,c", po::value<double>(&jf::shortCutoff)->default_value(1.375), "specify cutoff penalty for second pass.")
		    ("repeatsThrow,r", po::value<int>(&repeatsThrow)->default_value(4), "ignore reads with at least u telomere repeats remaining after cut. (=0 to turn off)")
		    ("sort", po::value<bool>(&o_opt_sort)->default_value(true), "sort output based on telomere length")
		    ("splitDir", po::value<bool>(&o_opt_splitDir)->default_value(true), "split output based on telomere orientation (forward/reverse)")
		    ("splitJunc", po::value<bool>(&o_opt_splitJunc)->default_value(true), "split output based on junction location")
		    ("revc", po::value<bool>(&o_opt_revc)->default_value(false), "output reads in same orientation (by revc-ing revc reads)")
		    ("fullMatches", po::value<bool>(&o_opt_fullMatches)->default_value(false), "include reads that are entirely telomeric sequence")
		    ("pretty", po::value<bool>(&o_opt_pretty)->default_value(false), "print each cut with color to stdout")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			cout << "Explanation: \n"
			     << "The algorithm does two passes. The first pass should be used to determine roughly when the telomere repeat ends;\n"
			     << "            the second pass should be used to determine more precisely where it ends.\n\n"
			     << "There are two ways to modify the way the algorithm behaves when encountering errors:\n"
			     << "            Window length: A smaller window length will weight individual errors higher, because penalties are relative to window size.\n"
			     << "            Cutoff: A lower cutoff will make it so that errors are weighted higher.\n\n"
			     << "It is recommended that the first pass has a large window size (>=30) with a low cutoff (<=1). Ideally, the window should exceed the cutoff\n"
			     << "            while there is some telomeric sequence left at the beginning.\n"
			     << "After the second window exceeds the cutoff, the junction location will be determined by the NEXT position that would otherwise represent\n"
			     << "            the beginning of a telomere repeat. Therefore, it is ideal to have the second window exceed the cutoff when there is a little\n"
			     << "            bit of telomeric sequence in the beginning (less than a full repeat.) However, if it is too sensitive, it may cut off due to a\n"
			     << "            sequencing error instead.\n"
			     << "Therefore, the second window should have a small window size (<=18) with a higher cutoff (>=1).\n"
			     << "These values can be modified to achieve desired effects and sensitivity.\n\n\n\n"
			     << "Recommended settings:\n\n\n"
			     << "For accurately cutting very close to the last repeat sequence:\n\n"
			     << " -w 60 -c .5 --shortWindowLen 15 --shortCutoff 1.4\n\n"
			     << "This seems to work well in most cases. It can be slow with the default step length of 1, this can be increased with relatively minimal penalty.\n";
			return false;
		}

		// OUTPUT DIRECTORY
		fs::path p(out_dirname);
		if (fs::exists(p) && !fs::is_directory(p)) {
			cerr << "Error: On output: Output specified must be a directory!";
			return false;
		}
		if (!fs::exists(p)) {
			create_directory(p);
		}
		// OUTPUT FILES
		if (o_opt_splitDir && o_opt_splitJunc) {
			telSeq_f_out = ofstream(out_dirname + "/" + "telSeq_f.fastq");
			subTelSeq_f_out = ofstream(out_dirname + "/" + "telAdjSeq_f.fastq");
			telSeq_r_out = ofstream(out_dirname + "/" + "telSeq_r.fastq");
			subTelSeq_r_out = ofstream(out_dirname + "/" + "telAdjSeq_r.fastq");
		} else if (o_opt_splitJunc) {
			telSeq_f_out = ofstream(out_dirname + "/" + "telSeq.fastq");
			subTelSeq_f_out = ofstream(out_dirname + "/" + "telAdjSeq.fastq");
		} else if (o_opt_splitDir) {
			telSeq_f_out = ofstream(out_dirname + "/" + "telSeq_f.fastq");
			telSeq_r_out = ofstream(out_dirname + "/" + "telSeq_r.fastq");
		} else {
			telSeq_f_out = ofstream(out_dirname + "/" + "telSeq.fastq");
		}
		// OUTPUT LOGS
		pt::ptime curTime = pt::ptime(pt::second_clock::local_time());
		string isoe = pt::to_iso_extended_string(curTime);
		log_stats = ofstream(out_dirname + "/stats_" + isoe + ".log");
		log_details = ofstream(out_dirname + "/details_" + isoe + ".log");
		log_details << isoe << endl;

		// INPUT FILE
		f_in = ifstream(i_filename);
		if (!f_in.is_open()) {
			cerr << "Error: on input: error opening file " << i_filename;
			return false;
		}
	} catch (exception& e) {
		cerr << e.what() << "\n";
		return false;
	} catch (...) {
		cerr << "Unknown error" << "\n";
		return false;
	}
	return true;
}

int main(int argc, char** argv) {
	ifstream f_in;
	ofstream telSeq_f_out;
	ofstream subTelSeq_f_out;
	ofstream telSeq_r_out;
	ofstream subTelSeq_r_out;

	bool result = process_options(argc, argv, f_in, telSeq_f_out, subTelSeq_f_out,
	                              telSeq_r_out, subTelSeq_r_out);

	if (!result)
		return 1;

	if (o_opt_splitDir && o_opt_splitJunc) {
		read_jf(f_in, &telSeq_f_out, &subTelSeq_f_out, &telSeq_r_out, &subTelSeq_r_out);
	} else if (o_opt_splitJunc) {
		read_jf(f_in, &telSeq_f_out, &subTelSeq_f_out, &telSeq_f_out, &subTelSeq_f_out);
	} else if (o_opt_splitDir) {
		read_jf(f_in, &telSeq_f_out, &telSeq_f_out, &telSeq_r_out, &telSeq_r_out);
	} else {
		read_jf(f_in, &telSeq_f_out, &telSeq_f_out, &telSeq_f_out, &telSeq_f_out);
	}

	f_in.close();

	telSeq_f_out.close();
	subTelSeq_f_out.close();
	telSeq_r_out.close();
	subTelSeq_r_out.close();

	log_stats.close();
	log_details.close();
}

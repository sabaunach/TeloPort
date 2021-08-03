/* junctionFinder.cpp
 * Functions for dealing with sequence quality and trimming sequences.
 * Author: Seth Baunach
 * Date: 6/30/2020
 */

#include "core/telomere_core.hpp"
#include "core/sequenceQuality.hpp"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace gg = boost::gregorian;
namespace pt = boost::posix_time;

namespace sq = sequenceQuality;

ofstream log_stats;
ofstream log_details;

int opt_cut, opt_lengththrow;
string opt_ifmt, opt_ofmt;
bool opt_log = false;

bool read_sq(ifstream & f_in, ofstream & f_out) {
	int seqsread = 0;
	core::Read read;

	while(true) {
		seqsread++;
		if (opt_ifmt == "fastq")
			core::read_sequence_fastq(f_in, read);
		else
			core::read_sequence_fasta(f_in, read);

		if (f_in.eof())
			break;

		// cut
		if (opt_cut < 0) {
			read.seq = read.seq.substr(max((int)read.seq.length() - abs(opt_cut), 0), abs(opt_cut));
			read.phred = read.phred.substr(max((int)read.phred.length() - abs(opt_cut), 0), abs(opt_cut));
		} else if (opt_cut > 0) {
			read.seq = read.seq.substr(0, opt_cut);
			read.phred = read.phred.substr(0, opt_cut);
		}
		// length throw (skip sequences that aren't right length)
		if (opt_lengththrow < 0 && read.seq.length() > (size_t)abs(opt_lengththrow))
			continue;
		else if (opt_lengththrow > 0 && read.seq.length() < (size_t)opt_lengththrow)
			continue;

		if (opt_ofmt == "fastq")
			core::write_sequence_fastq(f_out, read);
		else
			core::write_sequence_fasta(f_out, read);

		seqsread++;
	}

	// write output
	/*if (o_opt_sort) {
	    sort(juncReads.begin(), juncReads.end(),
	         [](const JunctionRead & a, const JunctionRead & b) -> bool {
	        return a.junc_loc > b.junc_loc;
	    });
	   }
	   if (o_opt_split) {
	    stable_partition(juncReads.begin(), juncReads.end(),
	                     [](const JunctionRead & a) -> bool {
	        return a.isForward;
	    });
	   }*/
	return true;
}

bool process_options(int argc, char** argv,
                     ifstream & f_in, ofstream & f_out) {
	try {
		string i_filename;
		string o_filename;

		po::options_description desc("Options");
		desc.add_options()
		    ("help,h", "produce this message")
		    ("in,i", po::value<string>(&i_filename), "specify input file.")
		    ("out,o", po::value<string>(&o_filename)->default_value(""), "specify output file. (blank for default)")
		    ("ifmt,f", po::value<string>(&opt_ifmt)->default_value("fastq"), "specify input format. (fastq/fasta)")
		    ("ofmt,g", po::value<string>(&opt_ofmt)->default_value("fastq"), "specify output format. (fastq/fasta)")
		    ("cut,c", po::value<int>(&opt_cut)->default_value(0), "cut sequence to specified length. 0=no cut. negative=cut from end.")
		    ("lengththrow,l", po::value<int>(&opt_lengththrow)->default_value(0), "throw out sequences less than specified length. 0=no restriction. negative=throw out greater")
		    ("log", "enable logging")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return false;
		}
		if (vm.count("log")) {
			opt_log = true;
		}

		// FORMATS
		if ((opt_ifmt != "fastq" && opt_ifmt != "fasta") || (opt_ofmt != "fastq" && opt_ofmt != "fasta")) {
			cerr << "Error: file format specified must be fastq or fasta!";
			return false;
		}

		// OUTPUT FILE
		f_out = ofstream(o_filename);

		// OUTPUT LOGS
		if (opt_log) {
			pt::ptime curTime = pt::ptime(pt::second_clock::local_time());
			string isoe = pt::to_iso_extended_string(curTime);
			log_stats = ofstream("stats_" + isoe + ".log");
			log_details = ofstream("details_" + isoe + ".log");
			log_details << isoe << endl;
		}

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
	ofstream f_out;

	bool result = process_options(argc, argv, f_in, f_out);

	if (!result)
		return 1;

	read_sq(f_in, f_out);

	f_in.close();
	f_out.close();
	log_stats.close();
	log_details.close();
}

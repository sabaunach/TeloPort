/* telomericFinder_driver.cpp
 * Options:
 *      -separated [r1.fastq] [r2.fastq]        - Take input from separated paired-end read fastq files.
 *      -interleaved [interleaved.fastq]        - Take input from an interleaved paired-end fastq file.
 *      -out [output directory]                 - Specify output directory. Default: ./telomereFinder
 * Mines the input files for telomeric reads and outputs telomeric reads into two fastq files, tel_reads.fastq and subtel_reads.fastq, in directory specified by -out
 * Author: Seth Baunach
 * Date: 6/30/2020
 *
 */

#include "telomereFinder/telomereFinder.hpp"
#include "core/sequenceQuality.hpp"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace gg = boost::gregorian;
namespace pt = boost::posix_time;

namespace tf = telomereFinder;
namespace sq = sequenceQuality;

ofstream log_stats;
ofstream log_details;
int seqReadOutputBreak = 10000;
double nRatioThrow = .05; // throw out reads that have more than nRatioThrow N's
int cutoff = 1;
int mode = 1;
int start = 18;

/* Writes statistical output to the log file
 */
void write_output(int telR1forandR2rev, int telR1revandR2for, int telR1fornotR2rev,
                  int telR1revnotR2for, int telR2fornotR1rev, int telR2revnotR1for,
                  double sumPenalty, int totaltelomeres) {
	log_stats << telR1forandR2rev << "\tForward telomere in R1 AND reverse in R2\n"
	          << telR1revandR2for << "\tReverse telomere in R1 AND forward in R2\n"
	          << telR1fornotR2rev << "\tForward telomere in R1 AND NOT reverse in R2\n"
	          << telR1revnotR2for << "\tReverse telomere in R1 AND NOT forward in R2\n"
	          << telR2fornotR1rev << "\tForward telomere in R2 AND NOT reverse in R1\n"
	          << telR2revnotR1for << "\tReverse telomere in R2 AND NOT forward in R1\n"
	          << endl
	          << telR1forandR2rev + telR1fornotR2rev << "\tForward telomeres from R1\n"
	          << telR1revandR2for + telR1revnotR2for << "\tReverse telomeres from R1\n"
	          << telR1revandR2for + telR2fornotR1rev << "\tForward telomeres from R2\n"
	          << telR1forandR2rev + telR2revnotR1for << "\tReverse telomeres from R2\n"
	          << endl
	          << telR1forandR2rev + telR1fornotR2rev + telR1revandR2for + telR1revnotR2for
	          << "\tTotal telomeres from R1\n"
	          << telR1revandR2for + telR2fornotR1rev + telR1forandR2rev + telR2revnotR1for
	          << "\tTotal telomeres from R2\n"
	          << endl
	          << telR1fornotR2rev + telR1revnotR2for << "\tTelomeres from R1 and not R2\n"
	          << telR2fornotR1rev + telR2revnotR1for << "\tTelomeres from R2 and not R1\n"
	          << endl
	          << telR1forandR2rev + telR1revandR2for << "\tTelomeres from both R1 and R2\n"
	          << endl;
	log_stats << "Average penalty for a telomere: " << (double)((double)sumPenalty / (double)totaltelomeres)
	          << endl << endl;
}

bool read_tf(ifstream * r1_in, ifstream * r2_in, ofstream & telReads_out, ofstream & pairReads_out) {
	core::Read r1, r2;
	log_details << "ID\tRead\tR1_Forward\tR1_Reverse\tR2_Forward\tR2_Reverse\tR1_Penalty\tR2_Penalty\n";
	// LOGGING
	int seqsread = 0, totalTelomeres = 0;
	double sumPenalty = 0;
	int telR1forandR2rev = 0,     // forward telomeric read in R1 and reverse in R2
	    telR1revandR2for = 0,     // reverse telomeric read in R1 and forward in R2
	    telR1fornotR2rev = 0,     // forward telomeric read in R1 and no reverse in R2
	    telR1revnotR2for = 0,     // reverse telomeric read in R1 and no forward in R2
	    telR2fornotR1rev = 0,     // forward telomeric read in R2 and no reverse in R1
	    telR2revnotR1for = 0;     // reverse telomeric read in R2 and no forward in R1
	while (true) {
		if (seqsread % seqReadOutputBreak == 0) {
			cout << seqsread << " sequences read: " << totalTelomeres << " telomeres found\n";
		}
		core::read_sequence_fastq(*r1_in, r1);
		core::read_sequence_fastq(*r2_in, r2);
		if (r1_in->eof() ^ r2_in->eof()) {
			cerr << "Error while reading from input files: one file has less sequences than the other!\n";
			return false;
		}
		if (r1_in->eof() && r2_in->eof()) {
			break;
		}

		/* CHECK READS FOR TELOMERIC REPEATS AND WRITE TO FILES */
		double r1_forward_res = tf::scoreTelomericRead(r1, tf::telRepeats[true], mode, start),
		       r1_reverse_res = tf::scoreTelomericRead(r1, tf::telRepeats[false], mode, start),
		       r2_forward_res = tf::scoreTelomericRead(r2, tf::telRepeats[true], mode, start),
		       r2_reverse_res = tf::scoreTelomericRead(r2, tf::telRepeats[false], mode, start);
		bool r1_forward_isTel = r1_forward_res < cutoff,
		     r1_reverse_isTel = r1_reverse_res < cutoff,
		     r2_forward_isTel = r2_forward_res < cutoff,
		     r2_reverse_isTel = r2_reverse_res < cutoff;

		if (r1_forward_isTel  || r1_reverse_isTel || r2_forward_isTel || r2_reverse_isTel) {

			double penaltyR1 = -1.0,
			       penaltyR2 = -1.0;

			// For now: output any read which matched a telomeric sequence in telReads_out
			//      and output its mate-pair in pairReads_out
			if (r1_forward_isTel || r1_reverse_isTel) {
				if ((double)((double)sq::countNs(r1) / (double)r1.seq.length()) > nRatioThrow)
					continue;
				totalTelomeres++;
				core::write_sequence_fastq(telReads_out, r1);
				core::write_sequence_fastq(pairReads_out, r2);
				log_details << r1.id << "\tR1\t";
			}
			if (r2_forward_isTel || r2_reverse_isTel) {
				if ((double)((double)sq::countNs(r2) / (double)r2.seq.length()) > nRatioThrow)
					continue;
				totalTelomeres++;
				core::write_sequence_fastq(telReads_out, r2);
				core::write_sequence_fastq(pairReads_out, r1);
				if (!(r1_forward_isTel || r1_reverse_isTel))
					log_details << r2.id << "\tR2\t";
			}
			if (r1_forward_isTel) {
				penaltyR1 = r1_forward_res;
				log_details << "+\t-\t";
				if (r2_reverse_isTel)
					telR1forandR2rev++;
				else
					telR1fornotR2rev++;
			} else if (r1_reverse_isTel) {
				penaltyR1 = r1_reverse_res;
				log_details << "-\t+\t";
				if (r2_forward_isTel)
					telR1revandR2for++;
				else
					telR1revnotR2for++;
			} else {
				log_details << "-\t-\t";
			}
			if (r2_forward_isTel) {
				log_details << "+\t-\t";
				penaltyR2 = r2_forward_res;
				if (!r1_reverse_isTel)
					telR2fornotR1rev++;
			} else if (r2_reverse_isTel) {
				log_details << "-\t+\t";
				penaltyR2 = r2_reverse_res;
				if (!r1_forward_isTel)
					telR2revnotR1for++;
			} else {
				log_details << "-\t-\t";
			}
			log_details << (penaltyR1 == -1.0 ? "-" : to_string(penaltyR1)) << "\t"
			            << (penaltyR2 == -1.0 ? "-" : to_string(penaltyR2)) << endl;
			penaltyR1 = penaltyR1 == -1.0 ? 0 : penaltyR1;
			penaltyR2 = penaltyR2 == -1.0 ? 0 : penaltyR2;
			if (penaltyR1 != -1.0 && penaltyR2 != -1.0)
				sumPenalty += (penaltyR1 + penaltyR2)/2.0;
			else
				sumPenalty += penaltyR1 != -1.0 ? penaltyR1 : penaltyR2;
		}
		seqsread++;
	}
	write_output(telR1forandR2rev, telR1revandR2for, telR1fornotR2rev,
	             telR1revnotR2for, telR2fornotR1rev, telR2revnotR1for,
	             sumPenalty, totalTelomeres);
	return true;
}

bool process_options(int argc, char** argv,
                     ifstream & f1_in,
                     ifstream & f2_in,
                     bool & interleaved,
                     ofstream & telReads_out,
                     ofstream & pairReads_out) {
	try {
		vector<string> s_filenames;
		string i_filename;
		string out_dirname;

		po::options_description desc("Allowed options");
		desc.add_options()
		    ("help,h", "produce this message")
		    ("separated,s", po::value<vector<string> >()->multitoken(), "specify input from separated paired-end read fastq files (r1 then r2).")
		    ("interleaved,i", po::value<string>(&i_filename), "specify input from an interleaved paired-end read fastq file.")
		    ("out,o", po::value<string>(&out_dirname)->default_value("telomereFinder_out"), "specify output directory. default: telomereFinder_out/")
		    ("nRatio,n", po::value<double>(&nRatioThrow)->default_value(.05), "throw out reads with a ratio of N's greater than specified. default: 0.05")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return false;
		}

		if (!vm.count("separated") && !vm.count("interleaved")) {
			cerr << "Error: On input: At least one option -i or -s must be specified!\n";
			cerr << desc << "\n";
			return false;
		}

		// OUTPUT FILES
		fs::path p(out_dirname);
		if (fs::exists(p) && !fs::is_directory(p)) {
			cerr << "Error: On output: Output specified must be a directory!";
			return false;
		}
		if (!fs::exists(p)) {
			create_directory(p);
		}
		telReads_out = ofstream(out_dirname + "/" + "telReads.fastq");
		pairReads_out = ofstream(out_dirname + "/" + "pairReads.fastq");
		// OUTPUT LOGS
		pt::ptime curTime = pt::ptime(pt::second_clock::local_time());
		string isoe = pt::to_iso_extended_string(curTime);
		log_stats = ofstream(out_dirname + "/stats_" + isoe + ".log");
		log_details = ofstream(out_dirname + "/details_" + isoe + ".log");
		log_details << isoe << endl;

		// INPUT FILES
		// separate (r1.fastq and r2.fastq)
		if (vm.count("separated")) {
			s_filenames = vm["separated"].as<vector<string> >();
			if (s_filenames.size() != 2) {
				cerr << "Error: With option -s: Exactly two files must be specified! (r1 and r2)\n";
				return false;
			}
			interleaved = false;
			f1_in = ifstream(s_filenames[0]);
			f2_in = ifstream(s_filenames[1]);
			if (!f1_in.is_open() || !f2_in.is_open()) {
				cerr << "Error: on input: error opening file " << s_filenames[0] << " or " << s_filenames[1];
				return false;
			}
		} else {
			// interleaved
			interleaved = true;
			f1_in = ifstream(i_filename);
			if (!f1_in.is_open()) {
				cerr << "Error: on input: error opening file " << i_filename;
				return false;
			}
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
	ifstream f1_in;
	ifstream f2_in;
	bool interleaved;
	ofstream telReads_out;
	ofstream pairReads_out;

	bool result = process_options(argc, argv, f1_in, f2_in, interleaved, telReads_out, pairReads_out);

	if (!result)
		return 1;

	if (interleaved) {
		cout << "Reading from interleaved\n";
		read_tf(&f1_in, &f1_in, telReads_out, pairReads_out);
		f1_in.close();
	} else {
		cout << "Reading from separated\n";
		read_tf(&f1_in, &f2_in, telReads_out, pairReads_out);
		f1_in.close();
		f2_in.close();
	}

	telReads_out.close();
	pairReads_out.close();

	log_stats.close();
	log_details.close();
}

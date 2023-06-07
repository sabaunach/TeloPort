/* wcdInterrogate_driver.cpp
 * Functions for interrogating cluster output from wcd.
 * Author: Seth Baunach
 * Date: 6/30/2020
 */

#include "core/telomere_core.hpp"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace gg = boost::gregorian;
namespace pt = boost::posix_time;

namespace core = telomere_core;

ofstream log_stats;
ofstream log_details;

bool opt_log = false,
     opt_ids = false,
     opt_indices = false,
     opt_aux = false,
     opt_sort = false;
int opt_size;
string o_mfadirname;

bool read_wcd(ifstream & f_wcd, ifstream & f_wcdin, ifstream & f_seqin, string seqin_fmt,
              ifstream & f_cons, ofstream & f_out) {

	vector<core::Read> wcdin_reads; // reads from wcdin (indices should line up with wcd output)
	core::read_seqFile(f_wcdin, "fasta", wcdin_reads);

	unordered_map<string, core::Read> seqin_map; // sequence info from seqin
	core::read_seqFile(f_seqin, seqin_fmt, seqin_map);


	vector<vector<size_t> > clusters; // clusters represented by indices

	// Read from wcd output file, add cluster information to clusters
	string line;
	while (getline(f_wcd, line)) {
		if (line.length() == 0)
			continue;
		if (line[line.length() - 1] != '.') { // not a cluster information
			if (opt_aux) {
				cout << line << endl;
			}
			continue;
		}
		// add read indices to cluster
		clusters.push_back(vector<size_t>());
		istringstream iss(line.substr(0, line.length() - 1));
		string index;
		while (iss >> index) {
			clusters.back().push_back((size_t)stoi(index));
		}
	}

	// output
	if (opt_sort) {
		sort(clusters.begin(), clusters.end(),
		     [](const vector<size_t> & a, const vector<size_t> & b) -> bool {
			return a.size() > b.size();
		});
	}

	for (size_t i = 0; i < clusters.size(); i++) {
		vector<size_t> cluster = clusters[i];

		if (cluster.size() < (size_t)opt_size)
			continue;

		ofstream clusterFile;
		if (o_mfadirname != "")
			clusterFile = ofstream(o_mfadirname + "/cluster" + to_string(i) + ".fasta");

		f_out << "cluster " << i << "; size=" << cluster.size() << endl;
		sort(cluster.begin(), cluster.end());
		for (size_t index: cluster) {

			if (index >= wcdin_reads.size()) {
				f_out << "MISSING: index out of range in wcd in: " << index << endl;
				continue;
			}

			string id = wcdin_reads[index].id;
			auto readkv = seqin_map.find(id);

			if (readkv == seqin_map.end()) {
				f_out << "MISSING: id not found in seq info: " << id << endl;
				continue;
			}

			core::Read read = readkv->second;

			if (clusterFile.is_open()) {
				core::write_sequence(clusterFile, read, "fasta");
			}
			if (opt_indices) {
				f_out << index << ":" << endl;
			}
			if (opt_ids) {
				f_out << read.id << endl;
			}
			f_out << read.seq << endl;
		}
		f_out << endl;
		if (clusterFile.is_open())
			clusterFile.close();
	}

	return true;
}

bool process_options(int argc, char** argv,
                     ifstream & f_wcd, ifstream & f_wcdin,
                     ifstream & f_seqin, string & seqin_fmt,
                     ifstream & f_cons, ofstream & f_out) {
	try {
		string wcd_filename,
		       wcdin_filename,
		       seqin_filename,
		       cons_filename,
		       o_filename;

		po::options_description desc("Options");
		desc.add_options()
		    ("help,h", "produce this message")
		    ("wcd,w", po::value<string>(&wcd_filename), "specify wcd file.")
		    ("wcdin,i", po::value<string>(&wcdin_filename), "specify fasta file used for wcd input.")
		    ("out,o", po::value<string>(&o_filename), "specify output file for displaying clusters")
		    ("outmfadir,r", po::value<string>(&o_mfadirname)->default_value(""), "specify output directory for cluster multifasta files (blank for no do)")
		    ("seqin,s", po::value<string>(&seqin_filename)->default_value(""), "specify file with sequence info (blank for wcdin)")
		    ("seqfmt,f", po::value<string>(&seqin_fmt)->default_value("fasta"), "specify format of seqin")
		    ("cons,c", po::value<string>(&cons_filename)->default_value(""), "specify file with consensus sequence info and display consensus sequence info with clusters (not required)")
		    ("size", po::value<int>(&opt_size)->default_value(2), "only use clusters of size >= size")
		    ("ids", "display IDs of sequences in clusters")
		    ("indices", "display indices of sequences in clusters")
		    ("aux", "display aux information from wcd output in stdout")
		    ("sort", "sort clusters by cluster size")
		    ("log", "enable logging")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return false;
		}
		// TO-DO: this could be cleaned up with some datastructures, maybe add this stuff to core/po?
		if (vm.count("sort")) {
			opt_sort = true;
		}
		if (vm.count("aux")) {
			opt_aux = true;
		}
		if (vm.count("indices")) {
			opt_indices = true;
		}
		if (vm.count("ids")) {
			opt_ids = true;
		}
		if (vm.count("log")) {
			opt_log = true;
		}

		// OUTPUT FILES
		if (o_mfadirname != "") {
			fs::path p(o_mfadirname);
			if (fs::exists(p) && !fs::is_directory(p)) {
				cerr << "Error: On output (mfa dir): Output specified must be a directory!";
				return false;
			}
			if (!fs::exists(p)) {
				create_directories(p);
			}
		}
		// FORMATS
		if (seqin_fmt != "fastq" && seqin_fmt != "fasta") {
			cerr << "Error: file format specified must be fastq or fasta!";
			return false;
		}
		// OUTPUT LOGS
		if (opt_log) {
			pt::ptime curTime = pt::ptime(pt::second_clock::local_time());
			string isoe = pt::to_iso_extended_string(curTime);
			log_stats = ofstream("stats_" + isoe + ".log");
			log_details = ofstream("details_" + isoe + ".log");
			log_details << isoe << endl;
		}

		// OUTPUT FILE
		f_out = ofstream(o_filename);
		if (!f_out.is_open()) {
			cerr << "Error: on output: error opening file " << o_filename;
			return false;
		}

		// INPUT FILES
		f_wcd = ifstream(wcd_filename);
		if (!f_wcd.is_open()) {
			cerr << "Error: on input: error opening file " << wcd_filename;
			return false;
		}
		f_wcdin = ifstream(wcdin_filename);
		if (!f_wcdin.is_open()) {
			cerr << "Error: on input: error opening file " << wcdin_filename;
			return false;
		}
		if (seqin_filename != "") {
			f_seqin = ifstream(seqin_filename);
			if (!f_seqin.is_open()) {
				cerr << "Error: on input: error opening file " << seqin_filename;
				return false;
			}
		} else {
			f_seqin = ifstream(wcdin_filename);
		}
		if (cons_filename != "") {
			f_cons = ifstream(cons_filename);
			if (!f_cons.is_open()) {
				cerr << "Error: on input: error opening file " << cons_filename;
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
	ifstream f_wcd, f_wcdin, f_cons, f_seqin;
	ofstream f_out;
	string seqin_fmt;

	bool result = process_options(argc, argv, f_wcd, f_wcdin, f_seqin, seqin_fmt,
	                              f_cons, f_out);

	if (!result)
		return 1;

	read_wcd(f_wcd, f_wcdin, f_seqin, seqin_fmt, f_cons, f_out);

	f_wcd.close();
	f_wcdin.close();
	f_seqin.close();
	f_cons.close();
	f_out.close();
	log_stats.close();
	log_details.close();
}

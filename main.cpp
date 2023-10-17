#include "header.h"

#define EMBED 1

using namespace std;
unsigned g_kmer_len = 6;
unsigned g_kmer_len_spacer = 4;
unsigned g_step = 2;
uint8_t g_code[256];
uint32_t g_mask;
Embedding g_e;

Motif spacer, fwd_primer, rev_primer;
short spacer_cov_threshold, STANDARD_NB_MOTIFS;
char mode;
vector<short> motif_src_ids;

unsigned ncpus = 1;
int BATCH_SIZE = 120;
string ofile_name = "";
ofstream ofile;
IType g_itype = IType::Guppy;

atomic_long embed_motif_time = 0, index_motif_time = 0, embed_flank_time = 0,
    index_read_time = 0, index_read_flank_time = 0,
    seed_flank_time = 0, normalize_5pos_time = 0, cluster_pos_time = 0,
    locate_flank_pos_time = 0, shift_pos_time = 0, seed_payload_time = 0,
    edit_index_time = 0, embed_candidate_time = 0, build_output_time = 0, write_output_time = 0,
    input_time = 0, clear_time = 0, load_motif_time = 0;
atomic_long seed_key_time = 0, seed_sort_time = 0, seed_threshold_time = 0;
atomic_long best_match_time = 0;


void make_code(void) {
  for (int i = 0; i < 256; i++)
    g_code[i] = 4;

  g_code['A'] = g_code['a'] = 0;
  g_code['N'] = g_code['n'] = 0;
  g_code['C'] = g_code['c'] = 1;
  g_code['G'] = g_code['g'] = 2;
  g_code['T'] = g_code['t'] = 3;
}

std::istream &safeGetline(std::istream &is, std::string &t) {
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n': return is;
      case '\r':
        if (sb->sgetc() == '\n')
          sb->sbumpc();
        return is;
      case std::streambuf::traits_type::eof():
        // Also handle the case when the last line has no line ending
        if (t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default: t += (char) c;
    }
  }
}

void print_stats() {
//#if DBGPRINT
  cout << "Breakdown:\n" <<
       "embed motifs: " << embed_motif_time / 1000000.0 << "\n" <<
       "index motifs: " << index_motif_time / 1000000.0 << "\n" <<
       "embed 3', 5' flank seq: " << embed_flank_time / 1000000.0 << "\n" <<
       "load motif time: " << load_motif_time / 1000000.0 << "\n" <<
       "load read time: " << input_time / 1000000.0 << "\n" <<
       "clear read vector time: " << clear_time / 1000000.0 << "\n" <<
       "index read kmer5: " << index_read_time / 1000000.0 << "\n" <<
       "index read kmer4 to locate flank seq: " << index_read_flank_time / 1000000.0 << "\n" <<
       "locate each motif's start pos: "
       << (seed_flank_time + normalize_5pos_time + cluster_pos_time + locate_flank_pos_time + shift_pos_time)
           / 1000000.0 << "\n" <<
       "\t kmer4 overlap seed of flank seq: " << seed_flank_time / 1000000.0 << "\n" <<
       "\t normalize the 5' pos: " << normalize_5pos_time / 1000000.0 << "\n" <<
       "\t cluster the pos: " << cluster_pos_time / 1000000.0 << "\n" <<
       "\t locate each motif's start pos: " << locate_flank_pos_time / 1000000.0 << "\n" <<
       "\t shift each motif's start pos (embed): " << shift_pos_time / 1000000.0 << "\n" <<
       "align each motif: " << (seed_payload_time + best_match_time) / 1000000.0 << "\n" <<
       "\t kmer5 overlap seed of the payload seq: " << seed_payload_time / 1000000.0 << "\n" <<
       "\t\t seed_key_time : " << seed_key_time / 1000000000.0 << "\n" <<
       "\t\t seed_sort_time : " << seed_sort_time / 1000000000.0 << "\n" <<
       "\t\t seed_threshold_time : " << seed_threshold_time / 1000000000.0 << "\n" <<
       "\t find the best candidate: " << best_match_time / 1000000.0 << "\n" <<
       "\t\t select best candidate by edit distance : " << edit_index_time / 1000000000.0 << "\n" <<
       "\t\t select best candidate by embed distance: " << embed_candidate_time / 1000000000.0 << "\n" <<
       "build output time: " << build_output_time / 1000000.0 << "\n" <<
       "write output time: " << write_output_time / 1000000.0 << endl;
//#endif
}

/*
 * index motif to get : (kmer, vector(motif_id))
 * only index the payload sequence, exclude the 3' 5' flanking seq
 */
void Motif_dict::index_motif() {
  for (unsigned j = 0; j < motifs.size(); j++) {
    Motif m = motifs[j];

    unsigned len = m.seq.size();
    uint32_t k = 0;
//    int start = flank3.seq.length();
//    int end = len - flank5.seq.length();
    size_t start = 0;
    size_t end = len;
    for (unsigned i = start; i < start + g_kmer_len - 1; i++) {
      k = (k << 2) + *(g_code + m.seq[i]);
    }

    for (unsigned i = start + g_kmer_len - 1; i < end; i++) {
      k = (k << 2) + *(g_code + m.seq[i]);
      unsigned key = k & g_mask;
      index[key].push_back(j);
    }
  }
}

class Tbb_decoder {
  vector<Read> &reads;
  vector<Motif_dict> &motif_dict_list;
  ostream &out;

 public:
  Tbb_decoder(vector<Read> &_reads, vector<Motif_dict> &_motif_dict_list, ostream &_out) :
      reads(_reads), motif_dict_list(_motif_dict_list), out(_out) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      reads[i].decode_read(motif_dict_list, out);
    }
  }
};

ReadStats decode_reads(vector<Motif_dict> &motif_dict_list, const string &rfname, ostream &out) {
  ifstream in(rfname);
  if (!in.is_open()) {
    cerr << "ERROR: Invalid input file " << rfname << endl;
    return ReadStats{0, 0};
  }
  cerr << "Processing read file " << rfname << endl;

  int ntotal = 0, nmapped = 0;
  vector<Read> reads;
  reads.reserve(BATCH_SIZE);

  tbb::task_scheduler_init init(ncpus);

  do {
    auto start = std::chrono::high_resolution_clock::now();

    ntotal++;
    Read r{};
    in >> r;
    if (g_itype == IType::Other) {
      r.name = to_string(ntotal);
    }

    reads.push_back(r);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    input_time += elapsed.count();

    if (ntotal % BATCH_SIZE == 0) {

      tbb::parallel_for(tbb::blocked_range<size_t>(0, BATCH_SIZE), Tbb_decoder(reads, motif_dict_list, out));

      start = std::chrono::high_resolution_clock::now();

      reads.clear();

      end = std::chrono::high_resolution_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      clear_time += elapsed.count();
    }
  } while (in);

  if (ntotal % BATCH_SIZE != 0) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, ntotal % BATCH_SIZE),
                      Tbb_decoder(reads, motif_dict_list, out));
  }

  return ReadStats{ntotal, nmapped};
}

int main(int argc, char *argv[]) {
  auto start_all = std::chrono::high_resolution_clock::now();

  namespace po = boost::program_options;
  bool help{};
  po::options_description description{"motif-search [options]"};
  description.add_options()
      ("help,h", po::bool_switch(&help), "Display help")
      ("ncpus,t", po::value<unsigned>(), "Number of threads to use")
      ("motifs,m", po::value<vector<string>>()->multitoken(), "File containing motifs to be mapped")
      ("read,r", po::value<vector<string>>()->multitoken(), "File(s) containing reads")
      ("kmerlen,l", po::value<unsigned>(), "Length of kmer (5)")
      ("output,o", po::value<string>(), "Ouptut file(default stdout)")
      ("fp", po::value<string>(), "Forward primer")
      ("rp", po::value<string>(), "Reverse primer")
      ("s", po::value<string>(), "Spacer")
      ("mode", po::value<char>(), "Mode L(long) or S (short)")
      ("n", po::value<unsigned>(), "Number of motifs in each oligo")
      ("id", po::value<vector<unsigned>>()->multitoken(), "The motif source id list for each motif (0-index). No need to specify if all motifs come from sane motif source file.")
      ("itype,i", po::value<string>(), "Input format (bonito or other), default is guppy");
  po::command_line_parser parser{argc, argv};
  parser.options(description);
  auto parsed_result = parser.run();
  po::variables_map vm;
  po::store(parsed_result, vm);
  po::notify(vm);

  if (help) {
    cerr << description << endl;
    return 0;
  }

  if (!vm["n"].empty()) {
    STANDARD_NB_MOTIFS = vm["n"].as<unsigned>();
    cerr << "Number of motifs in each oligo: " << STANDARD_NB_MOTIFS << endl;
  }else{
    cerr << "Please specify number of motifs in each oligo" << endl;
    return 1;
  }

  if (!vm["id"].empty()) {
    for (auto id: vm["id"].as<vector<unsigned>>())
      motif_src_ids.push_back(id);
  } else {
    for (int i = 0; i < STANDARD_NB_MOTIFS; ++i)
      motif_src_ids.push_back(0);
  }

  if (!vm["mode"].empty()) {
    mode = vm["mode"].as<char>();
    cerr << "Mode: " << mode << endl;
  }else{
    cerr << "Please specify the mode L(long) or S (short)" << endl;
    return 1;
  }

  if (!vm["fp"].empty()) {
    fwd_primer = {"fwd_primer", vm["fp"].as<string>()};
    cerr << "Forward primer: " << vm["fp"].as<string>() << endl;
  }else{
    fwd_primer = {"fwd_primer", ""};
    cerr << "No forward primer" << endl;
  }

  if (!vm["rp"].empty()) {
    rev_primer = {"rev_primer", vm["rp"].as<string>()};
    cerr << "Reverse primer: " << vm["rp"].as<string>() << endl;
  }else{
    rev_primer = {"rev_primer", ""};
    cerr << "No reverse primer" << endl;
  }

  if (!vm["s"].empty()) {
    spacer = {"spacer", vm["s"].as<string>()};
    spacer_cov_threshold = spacer.seq.size() / g_kmer_len_spacer; //at least have 1/4 cov
    cerr << "Spacer: " << vm["s"].as<string>() << endl;
  }else{
    cout << "Must specify spacer" << endl;
    return 1;
  }

  if (vm["motifs"].empty() || vm["read"].empty()) {
    cout << "Must specify motif and at least one read file as input" << endl;
    return 1;
  }

  if (!vm["ncpus"].empty()) {
    ncpus = vm["ncpus"].as<unsigned>();
    cerr << "Using " << ncpus << " cpus" << endl;
  }

  if (!vm["kmerlen"].empty()) {
    g_kmer_len = vm["kmerlen"].as<unsigned>();
    cerr << "Using kmer of size " << g_kmer_len << endl;
  }
  assert(g_kmer_len < 16);
  g_mask = (1U << (g_kmer_len * 2)) - 1;

  if (!vm["itype"].empty()) {
    if (boost::algorithm::to_lower_copy(vm["itype"].as<string>()) == "bonito") {
      g_itype = IType::Bonito;
      cerr << "Setting input read format to bonito\n";
    } else if (boost::algorithm::to_lower_copy(vm["itype"].as<string>()) == "other") {
      g_itype = IType::Other;
      cerr << "Setting input read format to other \n";
    } else
      cerr << "Setting input read format to guppy \n";
  }

  if (!vm["output"].empty()) {
    ofile_name = vm["output"].as<string>();
    ofile.open(ofile_name);
  }

  make_code();

  vector<Motif_dict> motif_dict_list;
  motif_dict_list.reserve(vm["motifs"].as<vector<string>>().size());

  for (const string &motif_name: vm["motifs"].as<vector<string>>()) {
    /*
    * Load motif dictionary, e.g.
    * >00000000000000000
    * GATTACAACCTTATGAGGACGAATCTCCCGCTTATACTACAACGCA
    */
    auto start = std::chrono::high_resolution_clock::now();
    cerr << "Loading motif dictionary from " << motif_name << endl;
    ifstream mfile(motif_name);
    Motif_dict motif_dict;
    Motif m;
    while (mfile >> m)
      motif_dict.motifs.push_back(m);

    if (!motif_dict.motifs.size()) {
      cerr << "ERR: there is no motif in motif dictionary" << endl;
      return 1;
    }
    motif_dict.motif_len = motif_dict.motifs[0].seq.size();
    cerr << ">>> Found " << motif_dict.motifs.size() << " motifs, length is " << motif_dict.motif_len << endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    load_motif_time += elapsed.count();

    /*
   * Embed motifs
   */
    start = std::chrono::high_resolution_clock::now();
    for (Motif &m: motif_dict.motifs) {
      //spacer (rp+fp) + payload
      string seg = spacer.seq + m.seq;
      g_e.embed_string(seg, m.eseq);
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    embed_motif_time += elapsed.count();

    /*
     * index payload of motifs
     * get (kmer hash, vector(motif_id))
     */
    start = std::chrono::high_resolution_clock::now();
    motif_dict.index_motif();
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    index_motif_time += elapsed.count();

    motif_dict_list.push_back(motif_dict);
  }

  /*
   * embed spacer
   */
  auto start = std::chrono::high_resolution_clock::now();
  g_e.embed_string(spacer.seq, spacer.eseq);
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_flank_time += elapsed.count();

  /*
   * Decode reads
   */
  auto decode_start = std::chrono::high_resolution_clock::now();
  for (const string &rf : vm["read"].as<vector<string>>()) {
    decode_reads(motif_dict_list, rf, ofile.is_open() ? ofile : cout);
  }
  auto decode_end = std::chrono::high_resolution_clock::now();
  auto decode_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(decode_end - decode_start);
  cerr << "Saved decoded sequence in " << ofile_name << endl;
  cerr << "Decode read time: " << decode_elapsed.count() / 1000000.0 << endl;

  if (ofile.is_open())
    ofile.close();

  print_stats();
  auto end_all = std::chrono::high_resolution_clock::now();
  auto elapsed_all = std::chrono::duration_cast<std::chrono::microseconds>(end_all - start_all);

  cerr << "Total time: " << elapsed_all.count() / 1000000.0 << "\n Completed.";

  return 0;
}

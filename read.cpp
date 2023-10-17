//
// Created by yan on 21/12/2021.
//

#include "header.h"

using namespace std;
std::mutex sam_mutex;
char rcsymbol[6] = "TGCAN";

extern short spacer_cov_threshold;
extern short STANDARD_NB_MOTIFS;
extern char mode;
extern uint8_t g_code[256];
extern IType g_itype;
extern unsigned g_kmer_len_spacer;
extern unsigned g_kmer_len;
extern unsigned g_step;
extern Motif spacer;
extern Motif fwd_primer;
extern Motif rev_primer;
extern vector<short> motif_src_ids;
extern Embedding g_e;
extern atomic_long seed_key_time;
extern atomic_long seed_threshold_time;
extern atomic_long edit_index_time;
extern atomic_long embed_candidate_time;
extern atomic_long best_match_time;
extern atomic_long index_read_time;
extern atomic_long index_read_flank_time;
extern atomic_long seed_flank_time;
extern atomic_long normalize_5pos_time;
extern atomic_long cluster_pos_time;
extern atomic_long locate_flank_pos_time;
extern atomic_long shift_pos_time;
extern atomic_long seed_payload_time;
extern atomic_long build_output_time;
extern atomic_long write_output_time;

std::istream &operator>>(std::istream &in, Read &r) {
  if (g_itype == IType::Guppy) {
    std::string tmp;
    if (!getline(in, r.name)) return in;
    if (!getline(in, r.fwd)) return in;
    if (!getline(in, tmp)) return in;
    if (!getline(in, tmp)) return in;
  } else if (g_itype == IType::Bonito) {
    if (!safeGetline(in, r.name)) return in;
    std::string seq;
    r.fwd = "";
    do {
      safeGetline(in, seq);
      r.fwd += seq;
    } while (in && in.peek() != '>');
  } else {
    if (!getline(in, r.fwd)) return in;
  }

  r.rlen = r.fwd.size();
  return in;
}

void Read::get_rev_seq() {
  stringstream ss;

  for (int i = rlen - 1; i >= 0; i--) {
    uint8_t c = *(g_code + fwd[i]);
    ss << rcsymbol[c];
  }
  rev = ss.str();
}

void Read::index_read(string &seq, unsigned kmer_len, bool is_flank) {
  uint32_t k = 0;
  for (unsigned i = 0; i < kmer_len - 1; i++) {
    k = (k << 2) + *(g_code + seq[i]);
  }

  for (unsigned i = kmer_len - 1; i < rlen; i++) {
    k = (k << 2) + *(g_code + seq[i]);
    uint32_t mask = (1U << (kmer_len * 2)) - 1;
    unsigned key = k & mask;

    if (is_flank) {
      if (seq == fwd)
        index_spacer_kmer[key].push_back(i - (kmer_len - 1));
      else
        index_rev_spacer_kmer[key].push_back(i - (kmer_len - 1));
    } else {
      if (seq == fwd)
        index[key].push_back(i - (kmer_len - 1));
      else
        index_rev[key].push_back(i - (kmer_len - 1));
    }

  }
}

/*
 * find all substring's pos in the read which has same hash
 * - overlap seed on the flank sequence with flank kmer size
 */
void Read::find_all(const string &flank_seq, vector<Segment> &candidates,
                    map<uint32_t, std::vector<unsigned>> &spacer_map) {
  int threshold = 1;
  vector<int> all_candidates;
  for (unsigned i = 0; i + g_kmer_len_spacer <= flank_seq.size(); i++) {
    uint32_t k = 0;
    for (unsigned j = i; j < i + g_kmer_len_spacer; j++)
      k = (k << 2) + *(g_code + flank_seq[j]);

    uint32_t g_mask_flank = (1U << (g_kmer_len_spacer * 2)) - 1;
    assert(k == (k & g_mask_flank));
    for (auto c = spacer_map[k].begin(); c != spacer_map[k].end(); ++c) {
      all_candidates.push_back(*c > i ? *c - i : 0);
    }
  }

  sort(all_candidates.begin(), all_candidates.end());
  unsigned ncandidates = all_candidates.size();
  for (unsigned i = 0; i < ncandidates;) {
    unsigned j = i + 1;
    while (j < ncandidates && all_candidates[i] == all_candidates[j])
      ++j;
    if (j - i >= threshold) {
      Segment pos;
      pos.pos = all_candidates[i];
      pos.cnt = j - i;
      candidates.push_back(pos);
    }
    i = j;
  }
}

/*
 * shift the start pos, to find the pos of minimal embed dist of the flank seq
 */
void Read::shift_start_pos(string &seq, vector<Segment> &start_pos) {
  int len = start_pos.size();

  for (int i = 1; i < len; i++) {
    int min_shift = 0;
    int min_edist = numeric_limits<int>::max();

    for (int j = -flank_max_shift; j <= flank_max_shift; j++) {
//      string_view cview3(seq.data() + start_pos[i].pos + j, flank3.seq.length());
//      unsigned edist3 = g_e.embed_compare(cview3, flank3.eseq, min_edist);
//
//      string_view cview5(seq.data() + start_pos[i].pos + j - flank5.seq.length(), flank5.seq.length());
//      unsigned edist5 = g_e.embed_compare(cview5, flank5.eseq, min_edist);

      string_view cview(seq.data() + start_pos[i].pos + j, spacer.seq.length());
      short edist = g_e.embed_compare(cview, spacer.eseq, min_edist);

      if (edist < min_edist) {
        min_shift = j;
        min_edist = edist;
      }
    }

    start_pos[i].pos += min_shift;
  }
}

void Read::pigeonhole_query(int start, int end, string &seq, map<uint32_t, vector<unsigned>> &motif_index,
                            vector<unsigned> &candidates, unsigned threshold) {
  auto start_t = std::chrono::high_resolution_clock::now();

  int payload_start = start;
  int payload_end = end;
  assert(payload_start + g_kmer_len <= seq.length());
  assert(payload_end <= seq.length()); //it's [start, end)

  unsigned nkmers = (payload_end - payload_start - g_kmer_len) / g_step + 1;
  size_t ntotal_hits = 0;
  unsigned kmer_idx = 0;
  vector<unsigned>::iterator b[nkmers], e[nkmers];

  // Take overlapping seeds and find all hits
  for (unsigned i = payload_start; i + g_kmer_len <= payload_end; i += g_step) {
    uint32_t k = 0;
    for (unsigned j = i; j < i + g_kmer_len; j++)
      k = (k << 2) + *(g_code + seq[j]);

    uint32_t g_mask = (1U << (g_kmer_len * 2)) - 1;
    assert(k == (k & g_mask));
    b[kmer_idx] = motif_index[k].begin();
    e[kmer_idx] = motif_index[k].end();
    ntotal_hits += motif_index[k].size();
    kmer_idx++;
  }
  assert(kmer_idx == nkmers);
  auto end_t = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t - start_t);
  seed_key_time += elapsed.count();


  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  unsigned MAX_POS = numeric_limits<unsigned>::max(), last_pos = MAX_POS;
  size_t nprocessed = 0;
  int last_cov = 0;

  unsigned top_pos[nkmers];
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] != e[i])
      top_pos[i] = *(b[i]);
    else
      top_pos[i] = MAX_POS;
  }

  start_t = std::chrono::high_resolution_clock::now();
  while (nprocessed < ntotal_hits) {
    //find min
    unsigned *min_item = min_element(top_pos, top_pos + nkmers);
    unsigned min_pos = *min_item;
    int min_kmer = min_item - top_pos;

    if (min_pos == last_pos) {
      last_cov++;
    } else {
      if (last_cov >= threshold)
        candidates.push_back(last_pos);
      last_cov = 1;
      last_pos = min_pos;
    }

    // add next element
    b[min_kmer]++;
    unsigned next_pos = b[min_kmer] != e[min_kmer] ? *(b[min_kmer]) : MAX_POS;
    *min_item = next_pos;

    ++nprocessed;
  }

  if (last_cov >= threshold && last_pos != MAX_POS) {
    candidates.push_back(last_pos);
  }
  end_t = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t - start_t);
  seed_threshold_time += elapsed.count();
}

void ksw_align(const char *tseq, int tlen, const char *qseq, int qlen,
               int sc_mch, int sc_mis, int gapo, int gape, ksw_extz_t &ez) {
  int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  const uint8_t *ts = reinterpret_cast<const uint8_t *>(tseq);
  const uint8_t *qs = reinterpret_cast<const uint8_t *>(qseq);
  memset(&ez, 0, sizeof(ksw_extz_t));
  ksw_extz2_sse(0, qlen, qs, tlen, ts, 5, mat, gapo, gape, -1, -1, 0, 0, &ez);
}

int Read::get_as_ksw(string &motif_seq, string &rseq, int &start, int &end, string &cigar) {

  const char *ptr_ref = motif_seq.c_str();
  const char *ptr_read = rseq.c_str() + start;
  int ref_len = motif_seq.length();
  int rlen = end - start;

  ksw_extz_t ez;

  // penalty
  ksw_align(ptr_ref, ref_len, ptr_read, rlen, SC_MCH, SC_MIS, GAPO, GAPE, ez);
  free(ez.cigar);

  return ez.score;
}

//int get_dist_wfa(Read &R, vector<Motif> &motifs, int pos, Flank_start_pos &start,
//				 Flank_start_pos &end, string &cigar) {
//	int ref_len = motifs[pos].seq.length();
//	const char *pattern = motifs[pos].seq.c_str(); //ref
//
//	int rlen = end.pos - start.pos;
//	char *text = (char *) (R.seq.c_str() + start.pos);
//
//	// Allocate MM
//	mm_allocator_t *const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
//	// Set penalties
//	affine_penalties_t affine_penalties = {
//		.match = -1,
//		.mismatch = 1,
//		.gap_opening = 1,
//		.gap_extension = 1,
//	};
//
//	// Init Affine-WFA
//	affine_wavefronts_t *affine_wavefronts = affine_wavefronts_new_complete(
//		ref_len, rlen, &affine_penalties, NULL, mm_allocator);
//	// Align
//	affine_wavefronts_align(affine_wavefronts, pattern, ref_len, text, rlen);
//
//	// Display alignment
//	edit_cigar_t *edit_cigar = &affine_wavefronts->edit_cigar;
////	std::stringstream cigar_ss;
////	char last_op = edit_cigar->operations[edit_cigar->begin_offset];
////
////	int last_op_length = 1;
////	int i;
////	for (i = edit_cigar->begin_offset + 1; i < edit_cigar->end_offset; ++i) {
////		//convert to M if X
////		last_op = last_op == 'X' ? 'M' : last_op;
////		if (edit_cigar->operations[i] == last_op || (edit_cigar->operations[i] == 'X' && last_op == 'M')) {
////			++last_op_length;
////		} else {
////			cigar_ss << last_op_length << last_op;
////			last_op = edit_cigar->operations[i];
////			last_op_length = 1;
////		}
////	}
////
////	cigar = cigar_ss.str();
//	int dist = edit_cigar_score_gap_affine(edit_cigar, &affine_penalties);
//
//	// Free
//	affine_wavefronts_delete(affine_wavefronts);
//	mm_allocator_delete(mm_allocator);
//
//	return -dist;
//}

/**
 * get the best candidate with minimal distance
 * - use KSW for the first/second segments, because they are the read-id
 * - use embed dist for the later segments
 * @param i
 * @param start
 * @param end
 * @param candidates
 * @param motifs
 * @return
 */
void Read::decode_ith_segment(int i, vector<Segment> &start_pos, string &seq, vector<unsigned> &candidates,
                              vector<Motif> &motifs, Segment &d_seg, vector<Motif_dict> &motif_dict_list) {
  auto start_all = std::chrono::high_resolution_clock::now();

  // first: fp + payload + spacer
  // later: spacer + payload
  // last: spacer + payload + rp
  int start = start_pos[i].pos, end = 0;
  if (i == start_pos.size() - 1)
    end = min(int(start + spacer.seq.size() + motif_dict_list[motif_src_ids.back()].motif_len + rev_primer.seq.size()),
              rlen);
  else if (i == 0)
    end = start_pos[i + 1].pos + spacer.seq.size();
  else
    end = start_pos[i + 1].pos;

  int min_idx = 0, max_as = numeric_limits<int>::min();
  string cigar = "";

  for (unsigned idx : candidates) {
    int as;

    //if (i == 0 || i == 1) {
      auto start_t = std::chrono::high_resolution_clock::now();

      string ref_seq;
      if (i == 0) {
        ref_seq = fwd_primer.seq + motifs[idx].seq + spacer.seq;
      } else if (i == candidates.size()-1) {
        ref_seq = spacer.seq + motifs[idx].seq + rev_primer.seq;
      } else {
      	ref_seq = spacer.seq + motifs[idx].seq;
      }
      as = get_as_ksw(ref_seq, seq, start, end, cigar);
      //dist = get_dist_wfa(r, motifs, pos, start, end, cigar);

      auto end_t = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t - start_t);
      edit_index_time += elapsed.count();
    //} else {
    //  auto start_t = std::chrono::high_resolution_clock::now();

    //  string_view cview(seq.data() + start, end - start);
    //  int edist = g_e.embed_compare(cview, motifs[idx].eseq, -max_as);
    //  as = -edist;

    //  auto end_t = std::chrono::high_resolution_clock::now();
    //  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t - start_t);
    //  embed_candidate_time += elapsed.count();
    //}

    if (as > max_as) {
      min_idx = idx;
      max_as = as;
    }
  }

  auto end_all = std::chrono::high_resolution_clock::now();
  auto elapsed_all = std::chrono::duration_cast<std::chrono::microseconds>(end_all - start_all);
  best_match_time += elapsed_all.count();

  d_seg.m_name = motifs[min_idx].name;
  d_seg.m_seq = motifs[min_idx].seq;
  d_seg.pos = start;
  d_seg.as = max_as;
  d_seg.cigar = cigar;
}

void cluster_start_pos(vector<Segment> &candidates, vector<Segment> &clustered_pos) {
  int cluster_dist = 5;
  int i = 0;
  while (i < candidates.size()) {
    int j = i + 1, sum_cnt = candidates[i].cnt;
    int sum = candidates[i].pos, nb_val = 1;
    //sum = candidates[i].pos * candidates[i].cnt,
    while (j < candidates.size() && abs(candidates[j].pos - candidates[i].pos) < cluster_dist) {
      sum += candidates[j].pos;
      //      sum += candidates[j].pos * candidates[j].cnt;
      nb_val += 1;
      sum_cnt += candidates[j].cnt;
      ++j;
    }

    Segment pos;
    pos.cnt = sum_cnt;
    pos.pos = sum / nb_val;

    clustered_pos.push_back(pos);
    i = j;
  }
}

void Read::locate_clustered_start_pos_multi(vector<Segment> &clustered_pos,
                                            vector<vector<Segment>> &start_pos_list,
                                            vector<Motif_dict> &motif_dict_list) {
  int m = 0, m_bk = 0, min_motif_len;
  int nb_pos = clustered_pos.size();

  while (m < nb_pos) {
    vector<Segment> start_pos;

    //pre-run: determine the first (fp-m1) and second (spacer-m2) motif's pos
    int i = 0;
    while (m < nb_pos) {
      if (clustered_pos[m].pos > motif_dict_list[motif_src_ids[i]].motif_len
          && clustered_pos[m].cnt >= spacer_cov_threshold) {
        int pos = clustered_pos[m].pos - motif_dict_list[motif_src_ids[i]].motif_len - fwd_primer.seq.size();
        start_pos.push_back(Segment{pos});
        ++i;
        start_pos.push_back(Segment{clustered_pos[m].pos});
        m_bk = m;
        ++i;
        ++m;
        break;
      }
      ++m;
    }

    for (; i < STANDARD_NB_MOTIFS && m < nb_pos; i++) {
      short motif_len = motif_dict_list[motif_src_ids[i]].motif_len;
      int seg_len = i == 0 ? motif_len : motif_len + spacer.seq.size();
      min_motif_len = seg_len * (1 - MAX_INDEL);

      int last = start_pos.back().pos;
      int dist = numeric_limits<int>::max();

      int min_dist = numeric_limits<int>::max(), min_index = m, max_cnt = 0;
      while (m < nb_pos && clustered_pos[m].pos < last + 2 * min_motif_len) {
        if (clustered_pos[m].cnt < max_cnt || clustered_pos[m].cnt <= spacer_cov_threshold) {
          m++;
          continue;
        }

        dist = abs((int) clustered_pos[m].pos - (int) (last + seg_len));

        if (dist < MAX_INDEL * seg_len &&
            (clustered_pos[m].cnt > max_cnt ||
                (clustered_pos[m].cnt == max_cnt && dist < min_dist))) {
          max_cnt = clustered_pos[m].cnt;
          min_dist = dist;
          min_index = m;
        }
        m++;
      }
//    m = min_index + 1; //the next motif go the next pos

      if (min_dist != numeric_limits<int>::max()) {
        assert(clustered_pos[min_index].pos > 0);
        if (rlen - clustered_pos[min_index].pos > min_motif_len)
          start_pos.push_back(Segment{clustered_pos[min_index].pos});
        else
          break;
      }
    }
//    // check m from the next pos
//    if (m_bk + 1 < nb_pos)
//      m = m_bk + 1;

    //short oligo must find all motifs
    if (start_pos.size() == STANDARD_NB_MOTIFS)
      start_pos_list.push_back(start_pos);
  }

}

void Read::locate_clustered_start_pos(vector<Segment> &clustered_pos,
                                      vector<Segment> &start_pos, vector<Motif_dict> &motif_dict_list) {
  int m = 0, min_motif_len;

  // set the first segment start pos
  start_pos.push_back(Segment{(int) fwd_primer.seq.size(), 0});

  for (int i = 1; i < STANDARD_NB_MOTIFS && m < clustered_pos.size(); i++) {
    short motif_len = motif_dict_list[motif_src_ids[i - 1]].motif_len;

    int seg_len = i == 1 ? motif_len : motif_len + spacer.seq.size();
    min_motif_len = seg_len * (1 - MAX_INDEL);

    int last = start_pos.back().pos;
    int dist = numeric_limits<int>::max();

    int min_dist = numeric_limits<int>::max(), min_index = m, max_cnt = 0;
    while (m < clustered_pos.size() && clustered_pos[m].pos < last + 2 * min_motif_len) {
      if (clustered_pos[m].cnt < max_cnt) {
        m++;
        continue;
      }

      dist = abs((int) clustered_pos[m].pos - (int) (last + seg_len));

      if (dist < MAX_INDEL * seg_len &&
          (clustered_pos[m].cnt > max_cnt ||
              (clustered_pos[m].cnt == max_cnt && dist < min_dist))) {
        max_cnt = clustered_pos[m].cnt;
        min_dist = dist;
        min_index = m;
      }
      m++;
    }
//    m = min_index + 1; //the next motif go the next pos

    if (min_dist != numeric_limits<int>::max()) {
      assert(clustered_pos[min_index].pos > 0);
      if (rlen - clustered_pos[min_index].pos > min_motif_len)
        start_pos.push_back(Segment{clustered_pos[min_index].pos});
      else
        break;
    } else {
      // not find a pos
      if (last + 2 * min_motif_len < rlen) {
        int pp = last + seg_len;
        assert(pp > 0);
        if (rlen - pp > min_motif_len)  //the motif length could be at min_motif_len
          start_pos.push_back(Segment{pp});
        else
          break;
      } else
        break; //already reach the end
    }
  }

  // in case not find seed at the end
  int last = start_pos.size() == 0 ? fwd_primer.seq.size() - motif_dict_list[motif_src_ids[0]].motif_len
                                   : start_pos.back().pos;
  while (start_pos.size() < STANDARD_NB_MOTIFS && last + 2 * min_motif_len < rlen) {
    short motif_len = motif_dict_list[motif_src_ids[start_pos.size() - 1]].motif_len;
    min_motif_len = (motif_len + spacer.seq.size()) * (1 - MAX_INDEL);
    int seg_len = start_pos.size() == 1 ? motif_len : motif_len + spacer.seq.size();
    last += seg_len;
    assert(last > 0);
    if (rlen - last > min_motif_len)
      start_pos.push_back(Segment{last});
    else
      break;
  }

}

void Read::locate_motif_start_pos(string &seq, vector<Motif_dict> &motif_dict_list,
                                  vector<vector<Segment>> &start_pos) {
  /*
   * index read:
   * - by kmer_len
   * - by a small kmer_len, be more accurate to locate the flank seq
   */
  auto start = std::chrono::high_resolution_clock::now();
  index_read(seq, g_kmer_len, false);
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  index_read_time += elapsed.count();

  start = std::chrono::high_resolution_clock::now();
  index_read(seq, g_kmer_len_spacer, true);
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  index_read_flank_time += elapsed.count();


  /*
   * locate the 3', 5' flanking sequence by overlapping seeding
   */
  start = std::chrono::high_resolution_clock::now();
  vector<Segment> flank3_candidates;
  if (seq == fwd)
    find_all(spacer.seq, flank3_candidates, index_spacer_kmer);
  else
    find_all(spacer.seq, flank3_candidates, index_rev_spacer_kmer);

  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  seed_flank_time += elapsed.count();
#if DBGPRINT
  cout << "3' candidates: " << endl;
  for (auto a: flank3_candidates)
      cout << a.pos << ": " << a.cnt << ", ";
  cout << endl;
#endif

  /*
   * cluster the position within dist (e.g. 5), and represent by the central point
   */
  start = std::chrono::high_resolution_clock::now();
  vector<Segment> clustered_pos;
  cluster_start_pos(flank3_candidates, clustered_pos);
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  cluster_pos_time += elapsed.count();
#if DBGPRINT
  cout << " clustered candidates: " << endl;
  for (auto a: clustered_pos)
      cout << a.pos << ": " << a.cnt << ", ";
  cout << endl;
#endif

  /*
   * chain the start pos
   */
  start = std::chrono::high_resolution_clock::now();
  if (mode == 'L') {
    vector<Segment> start_pos_long;
    locate_clustered_start_pos(clustered_pos, start_pos_long, motif_dict_list);
    start_pos.push_back(start_pos_long);
  } else if (mode == 'S') {
    locate_clustered_start_pos_multi(clustered_pos, start_pos, motif_dict_list);
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  locate_flank_pos_time += elapsed.count();
#if DBGPRINT
  cout << "== selected start pos : " << endl;
  for (auto a: start_pos){
      for (auto aa: a)
        cout << aa.pos << ", ";
      cout << endl;
  }
  cout << endl;
#endif

  /*
   * shift the start pos, to get the minimal embed dist of the flank seq
   */
  start = std::chrono::high_resolution_clock::now();
  for (vector<Segment> &s_pos: start_pos)
    shift_start_pos(seq, s_pos);
  end = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  shift_pos_time += elapsed.count();
#if DBGPRINT
  cout << "== selected start pos after shiftig: " << endl;
  for (auto a: start_pos){
      for (auto aa: a)
        cout << aa.pos << ", ";
      cout << endl;
  }
  cout << endl;
#endif

  //print the segmented subreads
#if DBGPRINT
  for (auto a: start_pos){
    cout << "====" << endl;
    for (int i = 0; i < a.size(); i++){
      int start = a[i].pos;
      int end = i == a.size()-1 ? a[i].pos + motif_dict_list[motif_src_ids[i]].motif_len + spacer.seq.size(): a[i+1].pos;
      string_view subread(seq.data() + start, end-start);
//		cout << "@" << r.name << "-" << i << endl;
      cout << start << ": " << end << endl;
      cout << subread << endl;
//		cout << "+" << endl;
//		cout << subread << endl;
    }
  }
#endif
}

/**
 * get the sequence of best match motifs for a read
 * for each segments decided the start_pos
 * - get candidate motifs by seeding
 * @param motifs
 * @param motif_index
 * @param start_pos
 * @param decode_motifs
 */
void Read::map_motifs(string &seq, vector<Motif_dict> &motif_dict_list,
                      vector<vector<Segment>> &start_pos_multi,
                      vector<vector<Segment>> &decoded_segs_multi) {
  for (vector<Segment> &start_pos: start_pos_multi) {
    vector<Segment> decoded_segs;
    for (int i = 0; i < start_pos.size(); i++) {
      Motif_dict &motif_dict = motif_dict_list[motif_src_ids[i]];
      int payload_start = i == 0 ? start_pos[i].pos : start_pos[i].pos + spacer.seq.size();
      int payload_end = i == start_pos.size() - 1 ?
                        min(int(rlen - rev_primer.seq.size()), int(payload_start + motif_dict.motif_len))
                                                  : start_pos[i + 1].pos;

      if (payload_end - payload_start < (int) g_kmer_len)
        continue;

      auto start_t = std::chrono::high_resolution_clock::now();

      vector<unsigned> candidates;
      pigeonhole_query(payload_start, payload_end, seq, motif_dict.index, candidates, 2);
      if (candidates.empty())
        pigeonhole_query(payload_start, payload_end, seq, motif_dict.index, candidates, 1);

      auto end_t = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
      seed_payload_time += elapsed.count();

      if (candidates.empty())
        break;

      Segment d_seg;
      decode_ith_segment(i, start_pos, seq, candidates, motif_dict.motifs, d_seg, motif_dict_list);
      decoded_segs.push_back(d_seg);
    }

    if (decoded_segs.size() == STANDARD_NB_MOTIFS)
      decoded_segs_multi.push_back(decoded_segs);
  }
}

void Read::print_decoded_read(vector<Segment> &decoded_segs, ostream &out) {

  //rname, flag(+/-), pos, dist(e.g. 3-4, concat by '-'), motif
  auto start_t = std::chrono::high_resolution_clock::now();

  // short oligo mode: AS must be bigger than 0
  if (mode == 'S'){
    for (Segment &d_seg: decoded_segs)
      if (d_seg.as < 0)
        return;
  }
  
  stringstream ss;
  ss << "@@" << name << "\t" << strand << "\t" << decoded_segs[0].pos << "\t";

  for (Segment &d_seg: decoded_segs)
    ss << d_seg.as << ">";
  ss << "\t";

  for (int i = 0; i < decoded_segs.size(); i++) {
    ss << decoded_segs[i].m_name;
    if (i != decoded_segs.size() - 1)
      ss << "-";
  }
  ss << endl;

  auto end_t = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
  build_output_time += elapsed.count();

  {
    start_t = std::chrono::high_resolution_clock::now();

    std::lock_guard<std::mutex> guard(sam_mutex);
    out << ss.str();

    end_t = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_t - start_t);
    write_output_time += elapsed.count();
  }
}

void Read::check_overlap(vector<Motif_dict> &motif_dict_list, bool *flag, vector<vector<Segment>> &decoded_segs) {

  int i = 0;
  int min_last_seg_len = (spacer.seq.size() + motif_dict_list[motif_src_ids.back()].motif_len + rev_primer.seq.size())
      * (1 - MAX_INDEL);

  while (i < (int) decoded_segs.size() - 1) {

    int end = decoded_segs[i].back().pos + min_last_seg_len;

    int j = i + 1, min_index = i, max_as = 0;
    for (int k = 0; k < decoded_segs[i].size(); ++k) {
      max_as += decoded_segs[i][k].as;
    }

    while (j < decoded_segs.size() - 1) {
      int start = decoded_segs[j].begin()->pos;
      if (start >= end)
        break;

      int as = 0;
      for (int k = 0; k < decoded_segs[i].size(); ++k) {
        as += decoded_segs[j][k].as;
      }
      if (as > max_as) { //max for AS
        min_index = j;
        max_as = as;
      }
      ++j;
    }

    for (int k = i; k < j; ++k) {
      if (k == min_index)
        flag[k] = 1; //pick it
      else
        flag[k] = 0; // not pick
    }
    i = i == min_index ? j : min_index;
  }
}

void Read::decode_strand(vector<Motif_dict> &motif_dict_list, bool **flag, string &seq,
                         vector<vector<Segment>> &decoded_segs) {
  //locate each seg
  vector<vector<Segment>> start_pos;
  locate_motif_start_pos(seq, motif_dict_list, start_pos);

  // find best position for map motifs
  map_motifs(seq, motif_dict_list, start_pos, decoded_segs);

  // check whether overlap between the chains
  if (mode == 'S' && decoded_segs.size()) {
    bool *_flag = new bool[decoded_segs.size()]();
    memset(_flag, 1, sizeof(bool) * decoded_segs.size());
    check_overlap(motif_dict_list, _flag, decoded_segs);
    *flag = _flag;
  }
}

int avg_as(vector<vector<Segment>> &decoded_segs) {
  if (decoded_segs.empty())
    return INT_MIN;

  int dist = 0;

  for (auto decoded_seg: decoded_segs) {
    for (auto seg: decoded_seg) {
      dist += seg.as;
    }
  }

  dist /= decoded_segs.size();

  return dist;
}

void Read::decode_read(vector<Motif_dict> &motif_dict_list, ostream &out) {
  //ignore too short/long read
  int motifs_len_sum = 0;
  for (auto motif_src_id: motif_src_ids)
    motifs_len_sum += motif_dict_list[motif_src_id].motif_len;

  int oligo_len = motifs_len_sum + spacer.seq.size() * (STANDARD_NB_MOTIFS - 1) +
      fwd_primer.seq.size() + rev_primer.seq.size();

//    cerr << "Warning: Read " <<  name << " is not valid, it's too short or too long" << endl;
//  cerr << name << endl;
  if (mode == 'L') { //one oligo per read
    if (rlen < oligo_len * (1 - MAX_INDEL) || fwd.size() > oligo_len * (1 + MAX_INDEL))
      return;
  } else if (mode == 'S') { //may have multiple oligos per read
    if (rlen < oligo_len * (1 - MAX_INDEL))
      return;
  }

  get_rev_seq();
  vector<vector<Segment>> decoded_segs_fwd;
  vector<vector<Segment>> decoded_segs_rev;

  bool *flag_fwd = nullptr;
  bool *flag_rev = nullptr;
  decode_strand(motif_dict_list, &flag_fwd, fwd, decoded_segs_fwd);
  decode_strand(motif_dict_list, &flag_rev, rev, decoded_segs_rev);

  vector<vector<Segment>> *decoded_segs;
  bool *flag = nullptr;
  // use the average of AS of each chain
  if (avg_as(decoded_segs_fwd) > avg_as(decoded_segs_rev)) {
    decoded_segs = &decoded_segs_fwd;
    strand = '+';
    flag = flag_fwd;
  } else {
    decoded_segs = &decoded_segs_rev;
    strand = '-';
    flag = flag_rev;

    // adjust the pos for rev
    for (auto &segs: decoded_segs_rev) {
      for (auto &seg: segs) {
        seg.pos = rlen - oligo_len - seg.pos;
      }
    }
  }

  for (int i = 0; i < decoded_segs->size(); ++i) {
    if ((*decoded_segs)[i].size() && (mode == 'L' || flag[i])) {
      print_decoded_read((*decoded_segs)[i], out);
    }
  }

  if (flag_fwd)
    delete[] flag_fwd;
  if (flag_rev)
    delete[] flag_rev;
}

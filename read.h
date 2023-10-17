//
// Created by yan on 21/12/2021.
//

#ifndef MOTIF_SEARCH__READ_H_
#define MOTIF_SEARCH__READ_H_

using namespace std;

struct Segment {
  int pos, cnt;
  std::string m_name, m_seq, cigar;  //decoded motif name, seq
  short as; //seg start pos, alignment score/use neg embedding distance if embed

  bool operator()(const Segment &X, const Segment &Y) const {
    return X.pos < Y.pos;
  }
};


struct Read {
 private:
  const float MAX_INDEL = 0.3;
// max indel length count from 3', 5' direction, 5' is longer because from end to locate start is longer
  int flank_max_shift = 5;


 public:
  int rlen;
  char strand;
  std::string name, fwd, rev;
  std::map<uint32_t, std::vector<unsigned>> index, index_spacer_kmer, index_rev, index_rev_spacer_kmer;

  void get_rev_seq();
  void index_read(string &seq, unsigned kmer_len, bool is_flank);
  void find_all(const std::string &flank_seq, std::vector<Segment> &candidates,
                map<uint32_t, std::vector<unsigned>> &spacer_map);
  void shift_start_pos(string &seq, std::vector<Segment> &start_pos);
  void pigeonhole_query(int start, int end, string &seq, std::map<uint32_t, std::vector<unsigned>> &motif_index,
                        std::vector<unsigned> &candidates, unsigned threshold);
  int get_as_ksw(string &motif_seq, string &rseq, int &start, int &end, string &cigar);
  void decode_ith_segment(int i, vector<Segment> &start_pos, string &seq, vector<unsigned> &candidates,
                          vector<Motif> &motifs, Segment &d_seg, vector<Motif_dict> &motif_dict_list);
  void locate_motif_start_pos(string &seq, vector<Motif_dict> &motif_dict_list,
                              vector<vector<Segment>> &start_pos);
  void locate_clustered_start_pos(vector<Segment> &clustered_pos,
                                  vector<Segment> &start_pos, vector<Motif_dict> &motif_dict_list);
  void locate_clustered_start_pos_multi(vector<Segment> &clustered_pos,
                                        vector<vector<Segment>> &start_pos_list, vector<Motif_dict> &motif_dict_list);
  void map_motifs(string &seq, vector<Motif_dict> &motif_dict_list, vector<vector<Segment>> &start_pos_multi,
                  vector<vector<Segment>> &decoded_segs_multi);
  void print_decoded_read(vector<Segment> &decoded_segs, ostream &out);

  void check_overlap(vector<Motif_dict> &motif_dict_list, bool *flag, vector<vector<Segment>> &decoded_segs);

  void decode_read(vector<Motif_dict> &motif_dict_list, ostream &out);
  void decode_strand(vector<Motif_dict> &motif_dict_list, bool **flag, string &seq,
                     vector<vector<Segment>> &decoded_segs);

  friend std::istream &operator>>(std::istream &in, Read &r);
};

#endif //MOTIF_SEARCH__READ_H_

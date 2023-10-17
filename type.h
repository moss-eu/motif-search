//
// Created by yan on 20/12/2021.
//

#ifndef MOTIF_SEARCH__TYPE_H_
#define MOTIF_SEARCH__TYPE_H_

enum class IType {
  Bonito,
  Guppy,
  Other
};

extern IType g_itype;

std::istream &safeGetline(std::istream &is, std::string &t);


struct Motif {
  std::string name, seq;
  std::vector<std::string> eseq; //embedded motif seq for each random value

  friend std::istream &operator>>(std::istream &in, Motif &m) {
    if (!safeGetline(in, m.name)) return in; // avoid the \r at the end of string
    return safeGetline(in >> std::ws, m.seq); //remove the whitespace in the beginning of line
  }
};

struct Motif_dict{
  short motif_len;
  std::vector<Motif> motifs;
  std::map<uint32_t, std::vector<unsigned>> index;

  void index_motif();
};

struct ReadStats {
  int ntotal, nmapped;
};

#endif //MOTIF_SEARCH__TYPE_H_

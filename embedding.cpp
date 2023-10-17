#include "header.h"
//#define CGK2_EMBED 1

using namespace std;
extern uint8_t g_code[256];

#ifdef CGK2_EMBED
	constexpr static int EFACTOR = 2;
#else
	constexpr static int EFACTOR = 3;
#endif

Embedding::Embedding() {
//	time_t seed = time(NULL);
//	srand(seed);
//	cerr << "Embedding using random seed " << seed << endl;

	srand(1559063236);
	cerr << "Embedding using random seed 1559063236 " << endl;

	for (int n = 0; n < N_RNDSTR; n++)
		for (int i = 0; i < N_ECHAR; i++)
			for (int j = 0; j < MAX_ELEN; j++)
				rnd_str[BITPOS(n, j, i)] = rand() % 2;
}

void Embedding::embed_string(const string_view &str, vector<string> &evec) {
	for (int i = 0; i < N_RNDSTR; i++)
		#ifdef CGK2_EMBED
			evec.push_back(cgk2_embed_string(str, i));
		#else
			evec.push_back(embed_string(str, i));
		#endif
}

unsigned Embedding::embed_compare(const string_view &to_embed,
								  const vector<string> &evec, unsigned threshold) {
	for (int i = 0; i < N_RNDSTR; i++) {
		#ifdef CGK2_EMBED
			unsigned nmismatch = cgk2_embed_nmismatch(to_embed, evec[i], i, threshold);
		#else
			unsigned nmismatch = embed_nmismatch(to_embed, evec[i], i, threshold);
		#endif
		if (nmismatch < threshold)
			threshold = nmismatch;
	}

	return threshold;
}

//unsigned Embedding::embed_compare_reverse(const string_view &to_embed,
//								  const vector<string> &evec, unsigned threshold) {
//	for (int i = 0; i < N_RNDSTR; i++) {
//#ifdef CGK2_EMBED
//		unsigned nmismatch = cgk2_embed_nmismatch_reverse(to_embed, evec[i], i, threshold);
//#else
//		//TODO: not implemented yet
//		unsigned nmismatch = embed_nmismatch_reverse(to_embed, evec[i], i, threshold);
//#endif
//		if (nmismatch < threshold)
//			threshold = nmismatch;
//	}
//
//	return threshold;
//}

string Embedding::embed_string(const string_view &str, int str_id) {
	unsigned len = str.size();
	unsigned elen = len * EFACTOR;
	string estr(elen, 0);
	unsigned i = 0;
	for (unsigned j = 0; j < elen; j++) {
		uint8_t s = i < len ? str[i] : EPAD;
		estr[j] = s;
		i = i + rnd_str[BITPOS(str_id, j, *(g_code + s))];
	}

	return estr;
}

unsigned Embedding::embed_nmismatch(const string_view &to_embed,
								  const string &ref, int str_id, unsigned threshold) {
	int len = to_embed.length();
	int elen = len * EFACTOR;
	int ref_len = ref.length();
//	assert(elen == ref.size());

	int i = 0, nmismatch = 0;
	for (int j = 0; j < min(elen, ref_len); j++) {
		uint8_t s = i < len ? to_embed[i] : EPAD;
		nmismatch += (ref[j] == s ? 0 : 1);
		if (nmismatch > threshold) {
			break;
		}
		i = i + rnd_str[BITPOS(str_id, j, *(g_code + s))];
	}

	nmismatch += abs(elen - ref_len);

	return nmismatch;
}

string Embedding::cgk2_embed_string(const string_view &str, int str_id) {
	int j = 0;
	unsigned rlen = str.size();
	int elen = EFACTOR * rlen;
	string estr(elen, 0);
	for (unsigned i = 0; i < rlen; i++) {
		uint8_t s = str[i];
		char bit = rnd_str[BITPOS(str_id, j, *(g_code + s))];
		if (!bit) {
			estr[j] = s;
			j++;
		} else {
			// both jth and j+1th value are s.
			estr[j + 1] = estr[j] = s;
			j += 2;
		}
	}

	//append the rest with EMBED_PAD
	for (; j < elen; j++) {
		estr[j] = EPAD;
	}
	assert(j <= elen);

	return estr;
}

int Embedding::cgk2_embed_nmismatch(const string_view &to_embed,
									const string &ref, int str_id, unsigned threshold) {
	int nmismatch = 0, j = 0;
	int rlen = to_embed.length();
	int elen = EFACTOR * rlen;
	int ref_len = ref.length();
	int cmp_len = min(elen, ref_len);

	for (unsigned i = 0; i < rlen && j + 1 < cmp_len; i++) {
		uint8_t s = to_embed[i];
		char bit = rnd_str[BITPOS(str_id, j, *(g_code + s))];
		if (!bit) {
			nmismatch += (ref[j] == s ? 0 : 1);
			if (nmismatch > threshold)
				goto end;
			j++;
		} else {
			// here, jth and j+1th value are both s.
			nmismatch += (ref[j] == s ? 0 : 1);
			nmismatch += (ref[j + 1] == s ? 0 : 1);
			if (nmismatch > threshold)
				goto end;
			j += 2;
		}
	}

	for (; j < cmp_len; j++)
		nmismatch += (ref[j] == EPAD ? 0 : 1);

	nmismatch += abs(elen - ref_len);

	end:
	return nmismatch;
}


int Embedding::cgk2_embed_nmismatch_reverse(const string_view &to_embed,
									const string &ref, int str_id, unsigned threshold) {
	int nmismatch = 0, j = 0;
	int rlen = to_embed.length();
	int elen = EFACTOR * rlen;
	int ref_len = ref.length();
	int cmp_len = min(elen, ref_len);

	for (unsigned i = 0; i < rlen && j + 1 < cmp_len; i++) {
		uint8_t s = to_embed[rlen - i];
		char bit = rnd_str[BITPOS(str_id, j, *(g_code + s))];
		if (!bit) {
			nmismatch += (ref[ref_len - j] == s ? 0 : 1);
			if (nmismatch > threshold)
				goto end;
			j++;
		} else {
			// here, jth and j+1th value are both s.
			nmismatch += (ref[ref_len - j] == s ? 0 : 1);
			nmismatch += (ref[ref_len - j - 1] == s ? 0 : 1);
			if (nmismatch > threshold)
				goto end;
			j += 2;
		}
	}

	nmismatch += max(elen, ref_len) - j;

	end:
	return nmismatch;
}

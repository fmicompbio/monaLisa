// Biostrings C API:
// https://github.com/Bioconductor/Biostrings/blob/master/inst/include/Biostrings_interface.h
//
// outdated:
// https://stat.ethz.ch/pipermail/bioc-devel/2012-December/003941.html
// https://stat.ethz.ch/pipermail/bioc-devel/2016-January/008517.html
//
// useful:
// base encoding:
//    > Biostrings:::.DNA_CODES
//      A  C  G  T  M  R  W  S  Y  K  V  H  D  B  N  -  +
//      1  2  4  8  3  5  9  6 10 12  7 11 13 14 15 16 32
// converted by:
// char DNAencode(char c);
// char DNAdecode(char code);
//


// R includes
#include <Rcpp.h>

// system includes
#include <math.h>
#include <string>

// package includes
#include "Biostrings_interface.h"
#include "_Biostrings_stubs.c"

#include "XVector_interface.h"
// #include "_XVector_stubs.c"

// recursive function to enumerate k-mers
void populate(const int k, int current_depth, std::string base, std::string* kmers, int* position) {
    if(current_depth == k) {
        kmers[*position].assign(base);
        (*position)++;
    } else {
        static char bases[] = { 'A', 'C', 'G', 'T' };
        for(int i = 0; i < 4; i++)
            populate(k, current_depth + 1, base + bases[i], kmers, position);
    }
}

// mapping table to map encoded bases in a DNAString to a numerical encoding
// A = 1, C = 2, G = 3, T = 4
int baseval[16] = {-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1};

// calculate k-mer index for k-mer starting at char*
int kmer_index_at(const char* seq_char, int k, int* pow4) {
    int p, index = 0, val = 0;
    const char *base;
    for (p = k - 1, base = seq_char; p >= 0; p--, base++) {
        val = baseval[(int)(*base)];
        if (val < 0) { // non-standard-base (A,C,G,T) character
            index = -1;
            break;
        }
        index += (pow4[p] * val);
    }
    return index;
}

// Code for walking along each sequence of an XStringSet object
// and read 1 character at a time.

//' Count k-mer pairs in sequences
//'
//' Over a set of sequences, for each k-mer a, count occurences of all k-mers b
//' within a defined distance downstream of the start of a.
//'
//' @param x DNAStringSet with sequences.
//' @param k An integer scalar with the length of sequences.
//' @param n An integer scalar defining the maximum downstream distance of
//'     second k-mers, relative to the start position of the first k-mer.
//' @param zoops (not implemented yet)
//'
//' @examples
//' count_kmer_pairs(DNAStringSet(c("AACCGGTT")), k = 2, n = 1)
//'
//' @return A numeric matrix with observed k-mer pairs counts.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix count_kmer_pairs(SEXP x,
                                     int k = 6,
                                     int n = 5,
                                     bool zoops = false) {
    // pre-flight checks
    if (! ::Rf_inherits(x, "DNAStringSet"))
        ::Rf_error("'x' must be a DNAStringSet");
    if (k <= 1)
        ::Rf_error("'k' must be greater than 1");
    if (k > 10)
        ::Rf_warning("'k' (%d) is large - this might take long an use a lot of memory", k);
    if (n < 1)
        ::Rf_error("'n' must be greater than 0");

    // prepare
    int i = 0, j = 0, l = 0, idx1 = 0, idx2 = 0;
    int *pow4 = new int[k + 1];
    for (i = 0; i <= k; i++)
        pow4[i] = pow(4, i);
    int nk = pow4[k];
    std::string* strkmers = new std::string[nk];
    int position = 0;
    populate(k, 0, "", strkmers, &position);
    Rcpp::CharacterVector kmers(nk);
    for (i = 0; i < nk; i++)
        kmers(i) = strkmers[i];

    Rcpp::NumericMatrix m(nk, nk);
    rownames(m) = kmers;
    colnames(m) = kmers;

    // loop through sequences
    const XStringSet_holder X = hold_XStringSet(x);
    const int x_len = get_XStringSet_length(x);

    Chars_holder seq;
    const char *seq_char1, *seq_char2;
    for (i = 0; i < x_len; i++) {
        seq = get_elt_from_XStringSet_holder(&X, i);
        if (seq.length < k)
            continue;
        //Rprintf("seq %d (%d bp): ", i+1, seq.length);
        for (j = 0, seq_char1 = seq.ptr; j < seq.length - k; j++, seq_char1++) {
            //Rprintf("%x ", *seq_char1);
            idx1 = kmer_index_at(seq_char1, k, pow4);
            if (idx1 < 0)
                continue; // ignore k-mers with non-ACGT characters
            //Rprintf("%d(%s) ", idx, strkmers[idx].c_str());
            for (l = j + 1, seq_char2 = seq_char1 + 1; l <= j+n && l <= seq.length - k; l++, seq_char2++) {
                idx2 = kmer_index_at(seq_char2, k, pow4);
                if (idx2 < 0)
                    continue; // ignore k-mers with non-ACGT characters
                m(idx1,idx2) += 1;
            }
        }
        //Rprintf("\n");
    }

    // clean up
    delete [] strkmers;
    delete [] pow4;

    return m;
}

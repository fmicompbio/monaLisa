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



#include <Rcpp.h>

#include <math.h>
#include <string>

#include "Biostrings_interface.h"

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

/*
cachedCharSeq X = cache_XRaw(x);
int x_len = X.length;  // get the length of the sequence in 'x'
const char *x_seq = X.seq;  // get a pointer to the sequence in 'x'
*/

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
//[[Rcpp::export]]
Rcpp::NumericMatrix count_kmer_pairs(SEXP x,
                                     int k = 6,
                                     int n = 5,
                                     bool zoops = false) {
    // pre-flight checks
    if (k <= 1) {
        ::Rf_error("'k' must be greater than 1");
    }
    if (k > 10) {
        ::Rf_warning("'k' (%d) is large - this might take long an use a lot of memory", k);
    }
    if (n < 1) {
        ::Rf_error("'n' must be greater than 0");
    }

    // prepare
    int nk = pow(4, k);
    std::string* strkmers = new std::string[nk];
    int position = 0;
    populate(k, 0, "", strkmers, &position);
    Rcpp::CharacterVector kmers(nk);
    for (int i = 0; i < nk; i++)
        kmers(i) = strkmers[i];
    delete [] strkmers;

    Rcpp::NumericMatrix m(nk, nk);
    rownames(m) = kmers;
    colnames(m) = kmers;

    // loop through sequences
    /*
cachedXStringSet X = cache_XStringSet(x);
int x_len = get_cachedXStringSet_length(&X);
for (int i = 0; i < x_len; i++) {
    cachedCharSeq X_elt = get_cachedXStringSet_elt(&X, i);
    int x_elt_len = X_elt.length;
    const char *x_elt_seq = X_elt.seq;
    for (int j = 0; j < x_elt_len; j++) {
        char c = x_elt_seq[j];
    }
     */
    return m;
}

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

// calculate k-mer indices for a given XString sequence
// and store them in kidx (grow size if necessary)
void calc_kmer_indices_for_seq(std::vector<int> &kidx,
                               Chars_holder &seq, int k, int* pow4) {
    int i;
    const char *seq_char;

    if (kidx.size() < (size_t)seq.length)
        kidx.resize(seq.length);

    for (i = 0, seq_char = seq.ptr; i < seq.length - k + 1; i++, seq_char++) {
        kidx[i] = kmer_index_at(seq_char, k, pow4);
    }
}


//' Count k-mer pairs in sequences
//'
//' Over a set of sequences, for each k-mer a, count occurences of all k-mers b
//' within a defined distance downstream of the start of a.
//'
//' @param x DNAStringSet with sequences.
//' @param k An integer scalar with the length of sequences.
//' @param n An integer scalar defining the maximum downstream distance of
//'     second k-mers, relative to the start position of the first k-mer.
//' @param zoops A logical scalar. If TRUE, count each observed k-mer pair
//'     only once per sequence.
//'
//' @examples
//' countKmerPairs(Biostrings::DNAStringSet(c("AACCGGTT")), k = 2, n = 1)
//'
//' @return A numeric matrix with observed k-mer pairs counts.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix countKmerPairs(SEXP x,
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
    int i = 0, j = 0, l = 0, a = 0, b = 0, idx1 = 0, idx2 = 0;
    std::vector<int> kidx;
    
    int *pow4 = new int[k + 1];
    for (i = 0; i <= k; i++)
        pow4[i] = pow(4, i);
    int nk = pow4[k];
    
    std::string* strkmers = new std::string[nk];
    int position = 0;
    populate(k, 0, "", strkmers, &position);
    
    Rcpp::CharacterVector kmers(nk);
    for (i = 0; i < nk; i++) {
        kmers(i) = strkmers[i];
    }
    
    Rcpp::NumericMatrix m(nk, nk);
    rownames(m) = kmers;
    colnames(m) = kmers;
    
    // loop through sequences
    const XStringSet_holder X = hold_XStringSet(x);
    const int x_len = get_XStringSet_length(x);
    
    Chars_holder seq;

    if (zoops) { // zoops == true
        
        bool **seen = new bool*[nk];
        for (a = 0; a < nk; a++)
            seen[a] = new bool[nk];
        
        for (i = 0; i < x_len; i++) {
            
            seq = get_elt_from_XStringSet_holder(&X, i);
            if (seq.length < k)
                continue;
            
            // pre-calculate k-mer indices
            calc_kmer_indices_for_seq(kidx, seq, k, pow4);
            
            // initialize seen matrix
            for (a = 0; a < nk; a++)
                for (b = 0; b < nk; b++)
                    seen[a][b] = false; 
            
            //Rprintf("seq %d (%d bp)\n", i+1, seq.length);
            for (j = 0; j < seq.length - k; j++) {
                idx1 = kidx[j];
                if (idx1 < 0)
                    continue; // ignore k-mers with non-ACGT characters
                for (l = j + 1; l <= j+n && l <= seq.length - k; l++) {
                    idx2 = kidx[l];
                    if (idx2 < 0 || seen[idx1][idx2])
                        continue; // ignore seen k-mer pairs (zoops = true) and k-mers with non-ACGT characters
                    seen[idx1][idx2] = true; 
                    m(idx1, idx2) += 1;
                }
            }
        }
        
        for (a = 0; a < nk; a++)
            delete [] seen[a];
        delete [] seen;
        
    } else {     // zoops == false
        for (i = 0; i < x_len; i++) {
            
            seq = get_elt_from_XStringSet_holder(&X, i);
            if (seq.length < k)
                continue;
            
            // pre-calculate k-mer indices
            calc_kmer_indices_for_seq(kidx, seq, k, pow4);

            //Rprintf("seq %d (%d bp)\n", i+1, seq.length);
            for (j = 0; j < seq.length - k; j++) {
                idx1 = kidx[j];
                if (idx1 < 0)
                    continue; // ignore k-mers with non-ACGT characters
                for (l = j + 1; l <= j+n && l <= seq.length - k; l++) {
                    idx2 = kidx[l];
                    if (idx2 < 0)
                        continue; // ignore k-mers with non-ACGT characters
                    m(idx1, idx2) += 1;
                }
            }
        }
    }
    
    // clean up
    delete [] strkmers;
    delete [] pow4;
    
    return m;
}


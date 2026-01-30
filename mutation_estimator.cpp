#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <cmath>
#include <vector>

using namespace std;

void print_usage() {
    cout << "Mutation Rate Estimator\n"
         << "=======================\n\n"
         << "Usage: mutation_estimator [options]\n\n"
         << "Required arguments:\n"
         << "  --mode <mode>      Estimation mode: pp, pc, or cc\n"
         << "  -m <length>        Sequence length (number of base pairs)\n"
         << "  -k <kmer_size>     K-mer size\n\n"
         << "Mode-specific arguments:\n"
         << "  --pp mode (presence-to-presence):\n"
         << "    -n <count>       Number of novel distinct k-mers in t's spectrum\n\n"
         << "  --pc mode (presence-to-count):\n"
         << "    -n <count>       Number of novel k-mers with multiplicity in t\n\n"
         << "  --cc mode (count-to-count):\n"
         << "    -n <count>       Number of novel k-mers with multiplicity in t\n"
         << "    -s <file>        FASTA file containing original sequence s\n\n"
         << "Optional arguments:\n"
         << "  -t <theta>         Sketch rate (default: 1.0, no sketching)\n"
         << "                     Use this if your novel k-mer count comes from\n"
         << "                     FracMinHash\n\n"
         << "Examples:\n"
         << "  # Using distinct k-mers of t(pp mode)\n"
         << "  mutation_estimator --mode pp -m 1000000 -k 21 -n 15000\n\n"
         << "  # Using k-mer counts of t (pc mode)\n"
         << "  mutation_estimator --mode pc -m 1000000 -k 21 -n 18000\n\n"
         << "  # Using k-mer counts of s and t (cc mode)\n"
         << "  mutation_estimator --mode cc -m 1000000 -k 21 -n 18000 -s seq.fasta\n\n"
         << "  # With FracMinHash sketch (theta = 0.01)\n"
         << "  mutation_estimator --mode pc -m 1000000 -k 21 -n 180 -t 0.01\n";
}

// Structure to hold input parameters
struct EstimatorInput {
    int L;              // Sequence length (will be converted to L-k+1)
    int k;              // k-mer size
    double theta;       // Sketch rate (1.0 for no sketching)
    
    // For --pp mode
    int novel_spectrum; // Number of novel kmers in t's spectrum (distinct)
    
    // For --pc and --cc modes
    int novel_multi;    // Number of novel kmers with multiplicity
    
    // For --cc mode only
    vector<string> s_sequences;  // Original sequences from FASTA
};

// Check if a k-mer contains N
bool contains_N(const string& kmer) {
    return kmer.find('N') != string::npos;
}

// Compute h1 values for all kmers efficiently using masked hash table
unordered_map<string, int> compute_h1_values(const set<string>& kmer_set, int k) {
    unordered_map<string, int> h1_map;
    
    // Build masked hash table for fast lookup
    unordered_map<string, set<string>> mask_to_kmers;
    
    for (const auto& kmer : kmer_set) {
        for (int pos = 0; pos < k; pos++) {
            string masked = kmer;
            masked[pos] = 'N';
            mask_to_kmers[masked].insert(kmer);
        }
    }
    
    // For each kmer, count how many hamming-1 neighbors exist in kmer_set
    for (const auto& kmer : kmer_set) {
        int h1_count = 0;
        
        for (int pos = 0; pos < k; pos++) {
            string masked = kmer;
            masked[pos] = 'N';
            
            if (mask_to_kmers.find(masked) != mask_to_kmers.end()) {
                h1_count += mask_to_kmers[masked].size() - 1;
            }
        }
        
        h1_map[kmer] = h1_count;
    }
    
    return h1_map;
}

// Get k-spectrum from sequence, excluding k-mers with N
set<string> kspectrum(const string& s, int k) {
    set<string> kmers;
    for (size_t i = 0; i <= s.length() - k; i++) {
        string kmer = s.substr(i, k);
        if (!contains_N(kmer)) {
            kmers.insert(kmer);
        }
    }
    return kmers;
}

// Count kmer occurrences, excluding k-mers with N
unordered_map<string, int> kmer_counter(const string& s, int k) {
    unordered_map<string, int> kmer_count;
    for (size_t i = 0; i <= s.length() - k; i++) {
        string kmer = s.substr(i, k);
        if (!contains_N(kmer)) {
            kmer_count[kmer]++;
        }
    }
    return kmer_count;
}

// Compute sum of occurrence * h1 for all kmers across all sequences
double compute_sum_occ_h1(const vector<string>& sequences, int k) {
    // Collect all k-mers from all sequences
    set<string> kmer_set;
    unordered_map<string, int> kmer_count;
    
    for (const auto& s : sequences) {
        set<string> seq_kmers = kspectrum(s, k);
        kmer_set.insert(seq_kmers.begin(), seq_kmers.end());
        
        unordered_map<string, int> seq_counts = kmer_counter(s, k);
        for (const auto& [kmer, count] : seq_counts) {
            kmer_count[kmer] += count;
        }
    }
    
    unordered_map<string, int> h1_map = compute_h1_values(kmer_set, k);
    
    double sum_occ_h1 = 0.0;
    for (const auto& [kmer, h1] : h1_map) {
        sum_occ_h1 += kmer_count[kmer] * h1;
    }
    
    return sum_occ_h1;
}

// Estimate mutation rate using r_ss formula (--pp mode)
// q_ss = novel_spectrum / L
// r_ss = 1 - (1 - q_ss)^(1/k)
double estimate_r_pp(const EstimatorInput& input) {
    int L_effective = input.L;  // Already L-k+1
    double novel = input.novel_spectrum / input.theta;
    
    // Cap novel at L_effective
    if (novel > L_effective) {
        novel = L_effective;
    }
    
    double q_ss = novel / L_effective;
    
    // Cap q at 1.0
    if (q_ss >= 1.0) {
        return 1.0;
    }
    
    double r_ss = 1.0 - pow(1.0 - q_ss, 1.0 / input.k);
    return r_ss;
}

// Estimate mutation rate using r_sm formula (--pc mode)
// q_sm = novel_multi / L
// r_sm = 1 - (1 - q_sm)^(1/k)
double estimate_r_pc(const EstimatorInput& input) {
    int L_effective = input.L;  // Already L-k+1
    double novel = input.novel_multi / input.theta;
    
    // Cap novel at L_effective
    if (novel > L_effective) {
        novel = L_effective;
    }
    
    double q_sm = novel / L_effective;
    
    // Cap q at 1.0
    if (q_sm >= 1.0) {
        return 1.0;
    }
    
    double r_sm = 1.0 - pow(1.0 - q_sm, 1.0 / input.k);
    return r_sm;
}

// Estimate mutation rate using r_strong formula (--cc mode)
// q_strong = (novel_multi + sum_occ_h1 * (1-r_sm)^(k-1) * r_sm / 3) / L
// r_strong = 1 - (1 - q_strong)^(1/k)
double estimate_r_cc(const EstimatorInput& input) {
    int L_effective = input.L;  // Already L-k+1
    double novel = input.novel_multi / input.theta;
    
    // Cap novel at L_effective
    if (novel > L_effective) {
        novel = L_effective;
    }
    
    // First compute r_sm as initial estimate
    double q_sm = novel / L_effective;
    if (q_sm >= 1.0) {
        return 1.0;
    }
    double r_sm = 1.0 - pow(1.0 - q_sm, 1.0 / input.k);
    
    // Compute sum_occ_h1 from all sequences
    double sum_occ_h1 = compute_sum_occ_h1(input.s_sequences, input.k);
    
    // Compute q_strong
    double correction = sum_occ_h1 * pow(1.0 - r_sm, input.k - 1) * r_sm / 3.0;
    double q_strong = (novel + correction) / L_effective;
    
    // If q_strong > 1.0, fall back to simple estimate
    if (q_strong >= 1.0) {
        return 1.0 - pow(novel / L_effective, 1.0 / input.k);
    }
    
    double r_strong = 1.0 - pow(1.0 - q_strong, 1.0 / input.k);
    return r_strong;
}

// Read all sequences from FASTA file
vector<string> read_fasta(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    
    vector<string> sequences;
    string line;
    string current_sequence;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Save previous sequence if exists
            if (!current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
        } else {
            // Convert to uppercase and append
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            current_sequence += line;
        }
    }
    
    // Don't forget the last sequence
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }
    
    file.close();
    return sequences;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    string mode = "";
    EstimatorInput input;
    input.L = -1;
    input.k = -1;
    input.theta = 1.0;
    input.novel_spectrum = -1;
    input.novel_multi = -1;
    // s_sequences is a vector, no need to initialize
    string s_file = "";
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage();
            return 0;
        }
        else if (arg == "--mode") {
            if (i + 1 < argc) {
                mode = argv[++i];
            }
        }
        else if (arg == "-m") {
            if (i + 1 < argc) {
                input.L = stoi(argv[++i]);
            }
        }
        else if (arg == "-k") {
            if (i + 1 < argc) {
                input.k = stoi(argv[++i]);
            }
        }
        else if (arg == "-n") {
            if (i + 1 < argc) {
                int n = stoi(argv[++i]);
                input.novel_spectrum = n;
                input.novel_multi = n;
            }
        }
        else if (arg == "-t") {
            if (i + 1 < argc) {
                input.theta = stod(argv[++i]);
            }
        }
        else if (arg == "-s") {
            if (i + 1 < argc) {
                s_file = argv[++i];
            }
        }
    }
    
    // Validate inputs
    if (mode != "pp" && mode != "pc" && mode != "cc") {
        cerr << "Error: Mode must be one of: pp, pc, cc\n";
        print_usage();
        return 1;
    }
    
    if (input.L <= 0) {
        cerr << "Error: Sequence length -m must be positive\n";
        return 1;
    }
    
    if (input.k <= 0) {
        cerr << "Error: K-mer size -k must be positive\n";
        return 1;
    }
    
    if (input.k > input.L) {
        cerr << "Error: K-mer size cannot be larger than sequence length\n";
        return 1;
    }
    
    if (mode == "pp" && input.novel_spectrum < 0) {
        cerr << "Error: --pp mode requires -n (novel distinct k-mer count)\n";
        return 1;
    }
    
    if (mode == "pc" && input.novel_multi < 0) {
        cerr << "Error: --pc mode requires -n (novel k-mer count with multiplicity)\n";
        return 1;
    }
    
    if (mode == "cc") {
        if (input.novel_multi < 0) {
            cerr << "Error: --cc mode requires -n (novel k-mer count with multiplicity)\n";
            return 1;
        }
        if (s_file.empty()) {
            cerr << "Error: --cc mode requires -s (original sequence file)\n";
            return 1;
        }
    }
    
    if (input.theta <= 0.0 || input.theta > 1.0) {
        cerr << "Error: Theta must be in range (0, 1]\n";
        return 1;
    }
    
    // Convert L to effective length (L - k + 1)
    input.L = input.L - input.k + 1;
    
    // Load sequences for cc mode
    if (mode == "cc") {
        input.s_sequences = read_fasta(s_file);
        if (input.s_sequences.empty()) {
            cerr << "Error: No sequences found in " << s_file << "\n";
            return 1;
        }
        
        // Check if any sequence is too short
        bool has_short_seq = false;
        for (const auto& seq : input.s_sequences) {
            if (seq.length() < static_cast<size_t>(input.k)) {
                has_short_seq = true;
                break;
            }
        }
        
        if (has_short_seq) {
            cerr << "Warning: Some sequences in " << s_file 
                 << " are shorter than k-mer size and will be skipped for those positions\n";
        }
        
        cout << "Loaded " << input.s_sequences.size() << " sequence(s) from " << s_file << "\n";
    }
    
    // Compute mutation rate estimate
    double r_estimate = 0.0;
    
    if (mode == "pp") {
        r_estimate = estimate_r_pp(input);
        cout << "Mode: presence-to-presence (--pp)\n";
    }
    else if (mode == "pc") {
        r_estimate = estimate_r_pc(input);
        cout << "Mode: presence-to-count (--pc)\n";
    }
    else if (mode == "cc") {
        r_estimate = estimate_r_cc(input);
        cout << "Mode: count-to-count (--cc)\n";
    }
    
    // Output results
    cout << "Parameters:\n";
    cout << "  Sequence length -k +1: " << input.L << "\n";
    cout << "  K-mer size: " << input.k << "\n";
    cout << "  Novel k-mers: " << (mode == "pp" ? input.novel_spectrum : input.novel_multi) << "\n";
    if (input.theta < 1.0) {
        cout << "  Sketch rate (theta): " << input.theta << "\n";
    }
    cout << "\nEstimated mutation rate: " << fixed << setprecision(6) << r_estimate << "\n";
    
    return 0;
}
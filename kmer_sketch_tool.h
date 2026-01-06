#ifndef KMER_SKETCH_TOOL_H
#define KMER_SKETCH_TOOL_H

#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <cstdint>
#include <fstream>
#include <memory>

using namespace std;

// ============================================================================
// 2-bit encoding for DNA sequences
// ============================================================================

class DNAEncoder {
public:
    // Encode DNA character to 2-bit
    // Returns 0xFF for invalid characters (N, etc.)
    static inline uint8_t encode(char c) {
        switch(c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0xFF; // Invalid character
        }
    }
    
    // Check if character is valid DNA base
    static inline bool is_valid_base(char c) {
        return (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
                c == 'a' || c == 'c' || c == 'g' || c == 't');
    }
    
    // Check if k-mer string is valid (no N or other non-ACGT)
    static bool is_valid_kmer(const string& kmer) {
        for (char c : kmer) {
            if (!is_valid_base(c)) return false;
        }
        return true;
    }
    
    // Decode 2-bit to DNA character
    static inline char decode(uint8_t code) {
        const char bases[] = "ACGT";
        return bases[code & 0x03];
    }
    
    // Complement base
    static inline uint8_t complement(uint8_t code) {
        return (~code) & 0x03; // A<->T (0<->3), C<->G (1<->2)
    }
    
    // Encode k-mer to 64-bit integer (supports k up to 32)
    // Returns UINT64_MAX if k-mer contains invalid characters
    static uint64_t encode_kmer(const string& kmer) {
        uint64_t encoded = 0;
        for (char c : kmer) {
            uint8_t code = encode(c);
            if (code == 0xFF) return UINT64_MAX; // Invalid k-mer
            encoded = (encoded << 2) | code;
        }
        return encoded;
    }
    
    // Decode 64-bit to k-mer string
    static string decode_kmer(uint64_t encoded, int k) {
        string kmer(k, 'A');
        for (int i = k - 1; i >= 0; i--) {
            kmer[i] = decode(encoded & 0x03);
            encoded >>= 2;
        }
        return kmer;
    }
    
    // Get reverse complement of encoded k-mer
    static uint64_t reverse_complement(uint64_t kmer, int k) {
        uint64_t rc = 0;
        for (int i = 0; i < k; i++) {
            rc = (rc << 2) | complement(kmer & 0x03);
            kmer >>= 2;
        }
        return rc;
    }
    
    // Get canonical k-mer (minimum of forward and reverse complement)
    static uint64_t canonical(uint64_t kmer, int k) {
        uint64_t rc = reverse_complement(kmer, k);
        return (kmer < rc) ? kmer : rc;
    }
};

// ============================================================================
// KmerSketch: Stores sketched k-mers and optional h1 statistics
// ============================================================================

struct KmerSketch {
    string id;                           // Sample/file ID
    int k;                              // k-mer size
    double theta;                       // Sketch fraction
    size_t total_kmers;                 // Total k-mers in sequence
    
    // Sketched canonical k-mers (encoded as 64-bit)
    unordered_map<uint64_t, uint16_t> kmer_counts;  // kmer -> count
    
    // Optional: h1 value for strong estimator
    bool has_h1;
    double h1_avg;                      // Average h1 value: sum(occ*h1) / L
    
    KmerSketch() : k(0), theta(1.0), total_kmers(0), has_h1(false), h1_avg(0.0) {}
    
    // Serialize to binary file
    void save(const string& filename) const;
    
    // Deserialize from binary file
    static KmerSketch load(const string& filename);
    
    // Get memory usage in bytes
    size_t memory_usage() const {
        return sizeof(KmerSketch) + 
               kmer_counts.size() * (sizeof(uint64_t) + sizeof(uint16_t)) +
               id.capacity();
    }
};

// ============================================================================
// SketchDatabase: Collection of sketches
// ============================================================================

class SketchDatabase {
private:
    vector<KmerSketch> sketches;
    string db_path;
    
public:
    SketchDatabase() = default;
    explicit SketchDatabase(const string& path) : db_path(path) {}
    
    void add_sketch(const KmerSketch& sketch) {
        sketches.push_back(sketch);
    }
    
    const vector<KmerSketch>& get_sketches() const {
        return sketches;
    }
    
    size_t size() const {
        return sketches.size();
    }
    
    // Save database to directory
    void save(const string& dir_path) const;
    
    // Load database from directory
    static SketchDatabase load(const string& dir_path);
    
    // Get total memory usage
    size_t total_memory() const {
        size_t total = 0;
        for (const auto& sketch : sketches) {
            total += sketch.memory_usage();
        }
        return total;
    }
};

// ============================================================================
// SketchBuilder: Build sketches from FASTA files
// ============================================================================

class SketchBuilder {
private:
    int k;
    double theta;
    bool compute_h1;
    uint32_t seed;
    int num_threads;
    
    // Hash function for sketching
    bool should_keep_kmer(uint64_t kmer_encoded) const;
    
    // Compute h1 values efficiently using bit operations
    unordered_map<uint64_t, int> compute_h1_values(
        const unordered_map<uint64_t, uint16_t>& kmer_counts
    ) const;
    
public:
    SketchBuilder(int k_size, double theta_val, bool calc_h1, 
                  uint32_t seed_val = 42, int threads = 1)
        : k(k_size), theta(theta_val), compute_h1(calc_h1), 
          seed(seed_val), num_threads(threads) {}
    
    // Build sketch from a single FASTA file
    KmerSketch build_sketch(const string& fasta_file, const string& sketch_id) const;
    
    // Build sketches from multiple FASTA files (parallel)
    SketchDatabase build_database(const vector<string>& fasta_files) const;
    
    // Build sketch from string (for testing)
    KmerSketch build_from_string(const string& seq, const string& id) const;
};

// ============================================================================
// QueryEngine: Query sequences against sketch database
// ============================================================================

struct QueryResult {
    string query_id;
    string target_id;
    size_t shared_kmers;      // Intersection size
    size_t novel_kmers;       // D_obs (with multiplicity)
    size_t query_total;       // Total k-mers in query
    double r_sm;              // Naive estimator
    double r_strong;          // Strong estimator (if available)
    bool has_strong;
    
    QueryResult() : shared_kmers(0), novel_kmers(0), query_total(0),
                   r_sm(0.0), r_strong(0.0), has_strong(false) {}
};

enum class QueryMode {
    PER_FILE,      // Concatenate all sequences in FASTA, one result per file
    PER_SEQUENCE   // Separate result for each sequence in FASTA
};

class QueryEngine {
private:
    const SketchDatabase& db;
    int k;
    double theta;
    uint32_t seed;
    int num_threads;
    
    // Hash function for sketching
    bool should_keep_kmer(uint64_t kmer_encoded) const;
    
    // Compute mutation rate estimates
    double compute_r_sm(size_t D_obs, size_t L, int k) const;
    double compute_r_strong(size_t D_obs, size_t L, int k, double h1_avg) const;
    
public:
    QueryEngine(const SketchDatabase& database, int k_size, 
                double theta_val, uint32_t seed_val = 42, int threads = 1)
        : db(database), k(k_size), theta(theta_val), 
          seed(seed_val), num_threads(threads) {}
    
    // Query a single sequence against all sketches in database
    vector<QueryResult> query(const string& sequence, const string& query_id) const;
    
    // Query from FASTA file - TWO MODES
    // Mode 1 (PER_FILE): Concatenate all sequences, one result per file
    // Mode 2 (PER_SEQUENCE): Separate result for each sequence
    vector<QueryResult> query_fasta(const string& fasta_file, 
                                     QueryMode mode = QueryMode::PER_FILE) const;
    
    // Batch query multiple sequences (parallel)
    vector<vector<QueryResult>> batch_query(
        const vector<string>& sequences, 
        const vector<string>& query_ids
    ) const;
};

// ============================================================================
// Utility functions
// ============================================================================

// Read FASTA file and return sequences
vector<pair<string, string>> read_fasta(const string& filename);

// Read single sequence from FASTA (first sequence only)
string read_fasta_single(const string& filename);

// Format memory size
string format_memory(size_t bytes);

// MurmurHash3 for k-mer sketching (from your existing code)
uint64_t murmurhash3_64(uint64_t kmer, uint32_t seed);

#endif // KMER_SKETCH_TOOL_H
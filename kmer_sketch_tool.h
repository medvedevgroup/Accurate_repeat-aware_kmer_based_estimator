#ifndef KMER_SKETCH_TOOL_H
#define KMER_SKETCH_TOOL_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <cmath>

using namespace std;

#ifdef __SIZEOF_INT128__
    typedef __uint128_t kmer_t;
    #define KMER_MAX_K 64
#else
    typedef uint64_t kmer_t;
    #define KMER_MAX_K 32
    #warning "128-bit integers not supported, k limited to 32"
#endif

namespace KmerUtils {
    inline string kmer_to_string(kmer_t val) {
        if (val == 0) return "0";
        
        char buffer[40];
        char* p = buffer + 39;
        *p = '\0';
        
        kmer_t tmp = val;
        while (tmp > 0) {
            *--p = '0' + (tmp % 10);
            tmp /= 10;
        }
        return string(p);
    }
    
    struct KmerHash {
        size_t operator()(const kmer_t& k) const {
            uint64_t low = static_cast<uint64_t>(k);
            uint64_t high = static_cast<uint64_t>(k >> 64);
            return hash<uint64_t>()(low ^ high);
        }
    };
}

uint64_t murmurhash3_64(kmer_t kmer, uint32_t seed);
vector<pair<string, string>> read_fasta(const string& filename);
string read_fasta_single(const string& filename);
string format_memory(size_t bytes);

class DNAEncoder {
public:
    static kmer_t encode_kmer(const string& kmer) {
        if (kmer.length() > KMER_MAX_K) {
            return ~kmer_t(0);
        }
        
        kmer_t encoded = 0;
        for (char c : kmer) {
            encoded <<= 2;
            switch(c) {
                case 'A': case 'a': encoded |= 0; break;
                case 'C': case 'c': encoded |= 1; break;
                case 'G': case 'g': encoded |= 2; break;
                case 'T': case 't': encoded |= 3; break;
                default: return ~kmer_t(0);
            }
        }
        return encoded;
    }
    
    static string decode_kmer(kmer_t encoded, int k) {
        string kmer(k, 'A');
        for (int i = k - 1; i >= 0; i--) {
            int base = encoded & 0x03;
            kmer[i] = "ACGT"[base];
            encoded >>= 2;
        }
        return kmer;
    }
    
    static kmer_t reverse_complement(kmer_t kmer, int k) {
        kmer_t rc = 0;
        for (int i = 0; i < k; i++) {
            rc <<= 2;
            rc |= (3 - (kmer & 0x03));
            kmer >>= 2;
        }
        return rc;
    }
    
    static kmer_t canonical(kmer_t kmer, int k) {
        kmer_t rc = reverse_complement(kmer, k);
        return (kmer < rc) ? kmer : rc;
    }
    
    static bool is_valid_kmer(const string& kmer) {
        for (char c : kmer) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                return false;
            }
        }
        return true;
    }
};

struct KmerSketch {
    string id;
    int k;
    double theta;
    uint32_t seed;
    size_t total_kmers;
    unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash> kmer_counts;
    
    bool has_h1;
    double h1_avg;
    
    KmerSketch() : k(0), theta(1.0), seed(42), total_kmers(0), has_h1(false), h1_avg(0.0) {}
    
    void save(const string& filename) const;
    static KmerSketch load(const string& filename);
    
    size_t memory_usage() const {
        return sizeof(*this) + 
               kmer_counts.size() * (sizeof(kmer_t) + sizeof(uint16_t)) +
               id.capacity();
    }
};

struct DatabaseMetadata {
    string db_path;
    vector<pair<string, string>> sketch_entries;
    int k;
    double theta;
    uint32_t seed;
    bool has_h1;
    
    DatabaseMetadata() : k(0), theta(1.0), seed(42), has_h1(false) {}
    
    static DatabaseMetadata load(const string& dir_path);
    size_t size() const { return sketch_entries.size(); }
};

class SketchDatabase {
private:
    vector<KmerSketch> sketches;
    string db_path;
    
public:
    SketchDatabase() {}
    SketchDatabase(const string& path) : db_path(path) {}
    
    void add_sketch(const KmerSketch& sketch) {
        sketches.push_back(sketch);
    }
    
    const vector<KmerSketch>& get_sketches() const {
        return sketches;
    }
    
    size_t size() const {
        return sketches.size();
    }
    
    size_t total_memory() const {
        size_t total = 0;
        for (const auto& sketch : sketches) {
            total += sketch.memory_usage();
        }
        return total;
    }
    
    void save(const string& dir_path) const;
    static SketchDatabase load(const string& dir_path);
};

enum class SketchMode {
    INDIVIDUAL,
    COMBINED
};

class SketchBuilder {
private:
    int k;
    double theta;
    kmer_t hash_threshold;
    bool compute_h1;
    uint32_t seed;
    int num_threads;
    
    bool should_keep_kmer(kmer_t kmer_encoded) const {
        if (theta >= 1.0) return true;
        uint64_t hash_val = murmurhash3_64(kmer_encoded, seed);
        return hash_val <= hash_threshold;
    }
    
    unordered_map<kmer_t, int, KmerUtils::KmerHash> compute_h1_values(
        const unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash>& kmer_counts
    ) const;
    
public:
    SketchBuilder(int kmer_size, double theta_val, bool h1_flag = false,
                  uint32_t hash_seed = 42, int threads = 1)
        : k(kmer_size), theta(theta_val), compute_h1(h1_flag), 
          seed(hash_seed), num_threads(threads) {
        
        if (theta >= 1.0) {
            hash_threshold = UINT64_MAX;
        } else {
            hash_threshold = static_cast<uint64_t>(theta * static_cast<double>(UINT64_MAX));
        }
    }
    
    KmerSketch build_from_string(const string& seq, const string& id) const;
    
    KmerSketch build_sketch(const string& fasta_file, const string& sketch_id = "") const;
    
    void build_database_incremental(const vector<string>& fasta_files,
                                     const string& output_dir,
                                     SketchMode mode = SketchMode::INDIVIDUAL) const;
};

struct QueryResult {
    string query_id;
    string target_id;
    size_t shared_kmers;
    size_t novel_kmers;
    size_t query_total;
    double r_sm;
    double r_strong;
    bool has_strong;
    
    QueryResult() : shared_kmers(0), novel_kmers(0), query_total(0),
                    r_sm(0.0), r_strong(0.0), has_strong(false) {}
};

enum class QueryMode {
    PER_FILE,
    PER_SEQUENCE,
    BATCH
};

struct QueryKmers {
    string query_id;
    unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash> kmer_counts;
    size_t total_kmers;
    size_t sketched_kmer_count;
    
    QueryKmers() : total_kmers(0), sketched_kmer_count(0) {}
};

class QueryEngine {
private:
    DatabaseMetadata metadata;
    int k;
    double theta;
    uint32_t seed;
    kmer_t hash_threshold;
    int num_threads;
    
    bool should_keep_kmer(kmer_t kmer_encoded) const {
        if (theta >= 1.0) return true;
        uint64_t hash_val = murmurhash3_64(kmer_encoded, seed);
        return hash_val <= hash_threshold;
    }
    
    double compute_r_sm(size_t D_obs, size_t L, int k) const;
    double compute_r_strong(size_t D_obs, size_t L, int k, double h1_avg) const;
    
    QueryKmers extract_query_kmers(const string& sequence, const string& query_id) const;
    
    QueryResult compare_with_sketch(const QueryKmers& query, const KmerSketch& target) const;
    
public:
    QueryEngine(const DatabaseMetadata& meta, int threads = 1);
    
    double get_theta() const { return theta; }
    int get_k() const { return k; }
    uint32_t get_seed() const { return seed; }
    
    void query_streaming(const vector<string>& fasta_files,
                        QueryMode mode,
                        const string& output_file) const;
};

#endif // KMER_SKETCH_TOOL_H
#include "kmer_sketch_tool.h"
#include "MurmurHash3.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <thread>
#include <mutex>
#include <sstream>
#include <iostream>

namespace fs = std::filesystem;

// ============================================================================
// MurmurHash3 wrapper for 64-bit k-mers
// ============================================================================

uint64_t murmurhash3_64(uint64_t kmer, uint32_t seed) {
    uint64_t hash[2];
    MurmurHash3_x64_128(&kmer, sizeof(uint64_t), seed, hash);
    return hash[0];
}

// ============================================================================
// Utility functions
// ============================================================================

vector<pair<string, string>> read_fasta(const string& filename) {
    vector<pair<string, string>> sequences;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return sequences;
    }
    
    string line, current_id, current_seq;
    
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Save previous sequence
            if (!current_id.empty() && !current_seq.empty()) {
                sequences.emplace_back(current_id, current_seq);
            }
            
            // Start new sequence
            current_id = line.substr(1);
            // Remove everything after first space
            size_t space_pos = current_id.find(' ');
            if (space_pos != string::npos) {
                current_id = current_id.substr(0, space_pos);
            }
            current_seq.clear();
        } else {
            // Append sequence (convert to uppercase)
            for (char c : line) {
                current_seq += toupper(c);
            }
        }
    }
    
    // Save last sequence
    if (!current_id.empty() && !current_seq.empty()) {
        sequences.emplace_back(current_id, current_seq);
    }
    
    file.close();
    return sequences;
}

string read_fasta_single(const string& filename) {
    auto sequences = read_fasta(filename);
    return sequences.empty() ? "" : sequences[0].second;
}

string format_memory(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_idx = 0;
    double size = static_cast<double>(bytes);
    
    while (size >= 1024.0 && unit_idx < 4) {
        size /= 1024.0;
        unit_idx++;
    }
    
    ostringstream oss;
    oss << fixed << setprecision(2) << size << " " << units[unit_idx];
    return oss.str();
}

// ============================================================================
// KmerSketch serialization
// ============================================================================

void KmerSketch::save(const string& filename) const {
    ofstream out(filename, ios::binary);
    
    if (!out.is_open()) {
        cerr << "Error: Cannot write to " << filename << endl;
        return;
    }
    
    // Write header
    size_t id_len = id.length();
    out.write(reinterpret_cast<const char*>(&id_len), sizeof(id_len));
    out.write(id.c_str(), id_len);
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&theta), sizeof(theta));
    out.write(reinterpret_cast<const char*>(&total_kmers), sizeof(total_kmers));
    out.write(reinterpret_cast<const char*>(&has_h1), sizeof(has_h1));
    out.write(reinterpret_cast<const char*>(&h1_avg), sizeof(h1_avg));
    
    // Write k-mer counts
    size_t num_kmers = kmer_counts.size();
    out.write(reinterpret_cast<const char*>(&num_kmers), sizeof(num_kmers));
    
    for (const auto& [kmer, count] : kmer_counts) {
        out.write(reinterpret_cast<const char*>(&kmer), sizeof(kmer));
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    }
    
    out.close();
}

KmerSketch KmerSketch::load(const string& filename) {
    KmerSketch sketch;
    ifstream in(filename, ios::binary);
    
    if (!in.is_open()) {
        cerr << "Error: Cannot read from " << filename << endl;
        return sketch;
    }
    
    // Read header
    size_t id_len;
    in.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
    sketch.id.resize(id_len);
    in.read(&sketch.id[0], id_len);
    in.read(reinterpret_cast<char*>(&sketch.k), sizeof(sketch.k));
    in.read(reinterpret_cast<char*>(&sketch.theta), sizeof(sketch.theta));
    in.read(reinterpret_cast<char*>(&sketch.total_kmers), sizeof(sketch.total_kmers));
    in.read(reinterpret_cast<char*>(&sketch.has_h1), sizeof(sketch.has_h1));
    in.read(reinterpret_cast<char*>(&sketch.h1_avg), sizeof(sketch.h1_avg));
    
    // Read k-mer counts
    size_t num_kmers;
    in.read(reinterpret_cast<char*>(&num_kmers), sizeof(num_kmers));
    
    for (size_t i = 0; i < num_kmers; i++) {
        uint64_t kmer;
        uint16_t count;
        in.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
        in.read(reinterpret_cast<char*>(&count), sizeof(count));
        sketch.kmer_counts[kmer] = count;
    }
    
    in.close();
    return sketch;
}

// ============================================================================
// SketchDatabase
// ============================================================================

void SketchDatabase::save(const string& dir_path) const {
    // Create directory if not exists
    fs::create_directories(dir_path);
    
    // Save index file
    ofstream index(dir_path + "/index.txt");
    for (size_t i = 0; i < sketches.size(); i++) {
        string sketch_file = "sketch_" + to_string(i) + ".bin";
        index << sketches[i].id << "\t" << sketch_file << "\n";
        
        // Save individual sketch
        sketches[i].save(dir_path + "/" + sketch_file);
    }
    index.close();
    
    cout << "Saved " << sketches.size() << " sketches to " << dir_path << endl;
    cout << "Total memory: " << format_memory(total_memory()) << endl;
}

SketchDatabase SketchDatabase::load(const string& dir_path) {
    SketchDatabase db(dir_path);
    
    // Read index file
    ifstream index(dir_path + "/index.txt");
    if (!index.is_open()) {
        cerr << "Error: Cannot find index.txt in " << dir_path << endl;
        return db;
    }
    
    string line;
    while (getline(index, line)) {
        istringstream iss(line);
        string id, sketch_file;
        iss >> id >> sketch_file;
        
        KmerSketch sketch = KmerSketch::load(dir_path + "/" + sketch_file);
        db.add_sketch(sketch);
    }
    
    index.close();
    
    cout << "Loaded " << db.size() << " sketches from " << dir_path << endl;
    cout << "Total memory: " << format_memory(db.total_memory()) << endl;
    
    return db;
}

// ============================================================================
// SketchBuilder - Part 1
// ============================================================================

bool SketchBuilder::should_keep_kmer(uint64_t kmer_encoded) const {
    if (theta >= 1.0) return true;  // Keep all k-mers
    
    uint64_t hash_val = murmurhash3_64(kmer_encoded, seed);
    double hash_frac = static_cast<double>(hash_val) / static_cast<double>(UINT64_MAX);
    return hash_frac <= theta;
}

unordered_map<uint64_t, int> SketchBuilder::compute_h1_values(
    const unordered_map<uint64_t, uint16_t>& kmer_counts
) const {
    unordered_map<uint64_t, int> h1_map;
    
    if (k > 32) {
        cerr << "Warning: k > 32, h1 computation may be slow" << endl;
        // Fallback to slower method for large k
        return h1_map;
    }
    
    // Build masked hash table
    // For 2-bit encoded k-mers, we can efficiently generate all 1-mutation neighbors
    unordered_map<uint64_t, vector<uint64_t>> mask_to_kmers;
    
    for (const auto& [kmer, _] : kmer_counts) {
        // For each position, create masked version
        uint64_t mask = 0x03;  // 2-bit mask
        for (int pos = 0; pos < k; pos++) {
            uint64_t masked = kmer & ~(mask << (2 * pos));
            mask_to_kmers[masked | ((uint64_t)4 << (2 * pos))].push_back(kmer);
            mask <<= 2;
        }
    }
    
    // Count h1 neighbors for each k-mer
    for (const auto& [kmer, _] : kmer_counts) {
        int h1_count = 0;
        
        uint64_t mask = 0x03;
        for (int pos = 0; pos < k; pos++) {
            uint64_t masked = kmer & ~(mask << (2 * pos));
            uint64_t key = masked | ((uint64_t)4 << (2 * pos));
            
            if (mask_to_kmers.find(key) != mask_to_kmers.end()) {
                // Count neighbors, excluding self
                h1_count += mask_to_kmers[key].size();
                // Check if kmer itself is in the list
                for (uint64_t neighbor : mask_to_kmers[key]) {
                    if (neighbor == kmer) {
                        h1_count--;
                        break;
                    }
                }
            }
            mask <<= 2;
        }
        
        h1_map[kmer] = h1_count;
    }
    
    return h1_map;
}
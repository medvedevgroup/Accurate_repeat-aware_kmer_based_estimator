#include "kmer_sketch_tool.h"
#include "MurmurHash3.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <thread>
#include <mutex>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace fs = std::filesystem;

uint64_t murmurhash3_64(kmer_t kmer, uint32_t seed) {
    uint64_t hash[2];
    
    uint64_t kmer_low = static_cast<uint64_t>(kmer);
    
#ifdef __SIZEOF_INT128__
    uint64_t kmer_high = static_cast<uint64_t>(kmer >> 64);
    if (kmer_high != 0) {
        kmer_low ^= kmer_high;
    }
#endif
    
    MurmurHash3_x64_128(&kmer_low, sizeof(uint64_t), seed, hash);
    return hash[0];
}

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
            if (!current_id.empty() && !current_seq.empty()) {
                sequences.emplace_back(current_id, current_seq);
            }
            
            current_id = line.substr(1);
            size_t space_pos = current_id.find(' ');
            if (space_pos != string::npos) {
                current_id = current_id.substr(0, space_pos);
            }
            current_seq.clear();
        } else {
            for (char c : line) {
                current_seq += toupper(c);
            }
        }
    }
    
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

void KmerSketch::save(const string& filename) const {
    ofstream out(filename, ios::binary);
    
    if (!out.is_open()) {
        cerr << "Error: Cannot write to " << filename << endl;
        return;
    }
    
    size_t id_len = id.length();
    out.write(reinterpret_cast<const char*>(&id_len), sizeof(id_len));
    out.write(id.c_str(), id_len);
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&theta), sizeof(theta));
    out.write(reinterpret_cast<const char*>(&seed), sizeof(seed));
    out.write(reinterpret_cast<const char*>(&total_kmers), sizeof(total_kmers));
    out.write(reinterpret_cast<const char*>(&has_h1), sizeof(has_h1));
    out.write(reinterpret_cast<const char*>(&h1_avg), sizeof(h1_avg));
    
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
    
    size_t id_len;
    in.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
    sketch.id.resize(id_len);
    in.read(&sketch.id[0], id_len);
    in.read(reinterpret_cast<char*>(&sketch.k), sizeof(sketch.k));
    in.read(reinterpret_cast<char*>(&sketch.theta), sizeof(sketch.theta));
    in.read(reinterpret_cast<char*>(&sketch.seed), sizeof(sketch.seed));
    in.read(reinterpret_cast<char*>(&sketch.total_kmers), sizeof(sketch.total_kmers));
    in.read(reinterpret_cast<char*>(&sketch.has_h1), sizeof(sketch.has_h1));
    in.read(reinterpret_cast<char*>(&sketch.h1_avg), sizeof(sketch.h1_avg));
    
    size_t num_kmers;
    in.read(reinterpret_cast<char*>(&num_kmers), sizeof(num_kmers));
    
    for (size_t i = 0; i < num_kmers; i++) {
        kmer_t kmer;
        uint16_t count;
        in.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
        in.read(reinterpret_cast<char*>(&count), sizeof(count));
        sketch.kmer_counts[kmer] = count;
    }
    
    in.close();
    return sketch;
}

void SketchDatabase::save(const string& dir_path) const {
    fs::create_directories(dir_path);
    
    ofstream index(dir_path + "/index.txt");
    for (size_t i = 0; i < sketches.size(); i++) {
        string sketch_file = "sketch_" + to_string(i) + ".bin";
        index << sketches[i].id << "\t" << sketch_file << "\n";
        
        sketches[i].save(dir_path + "/" + sketch_file);
    }
    index.close();
    
    cout << "Saved " << sketches.size() << " sketches to " << dir_path << endl;
    cout << "Total memory: " << format_memory(total_memory()) << endl;
}

SketchDatabase SketchDatabase::load(const string& dir_path) {
    SketchDatabase db(dir_path);
    
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

DatabaseMetadata DatabaseMetadata::load(const string& dir_path) {
    DatabaseMetadata meta;
    meta.db_path = dir_path;
    
    ifstream index(dir_path + "/index.txt");
    if (!index.is_open()) {
        cerr << "Error: Cannot find index.txt in " << dir_path << endl;
        return meta;
    }
    
    string line;
    while (getline(index, line)) {
        istringstream iss(line);
        string id, sketch_file;
        iss >> id >> sketch_file;
        meta.sketch_entries.emplace_back(id, sketch_file);
    }
    index.close();
    
    // Load first sketch to get k, theta, seed, has_h1
    if (!meta.sketch_entries.empty()) {
        string first_sketch = dir_path + "/" + meta.sketch_entries[0].second;
        KmerSketch sample = KmerSketch::load(first_sketch);
        meta.k = sample.k;
        meta.theta = sample.theta;
        meta.seed = sample.seed;
        meta.has_h1 = sample.has_h1;
    }
    
    cout << "Loaded metadata: " << meta.sketch_entries.size() 
         << " sketches, k=" << meta.k << ", theta=" << meta.theta 
         << ", seed=" << meta.seed << endl;
    
    return meta;
}

unordered_map<kmer_t, int, KmerUtils::KmerHash> SketchBuilder::compute_h1_values(
    const unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash>& kmer_counts
) const {
    unordered_map<kmer_t, int, KmerUtils::KmerHash> h1_map;
    
    if (k > KMER_MAX_K) {
        cerr << "Warning: k > " << KMER_MAX_K << ", h1 computation may be slow" << endl;
        return h1_map;
    }
    
    // For each k-mer, count how many Hamming-1 neighbors exist in the set
    for (const auto& [kmer, _] : kmer_counts) {
        int h1_count = 0;
        
        // Try mutating each position
        for (int pos = 0; pos < k; pos++) {
            // Extract current base at this position (2 bits)
            // Position 0 is the rightmost (least significant) 2 bits
            kmer_t shift_amount = 2 * pos;
            kmer_t base_mask = kmer_t(0x03) << shift_amount;
            kmer_t current_base = (kmer & base_mask) >> shift_amount;
            
            // Try all 4 possible bases (A=0, C=1, G=2, T=3)
            for (kmer_t new_base = 0; new_base < 4; new_base++) {
                // Skip if it's the same as current base
                if (new_base == current_base) {
                    continue;
                }
                
                // Create mutated k-mer by replacing the base at this position
                kmer_t mutated = (kmer & ~base_mask) | (new_base << shift_amount);
                
                // Check if this mutated k-mer exists in the set
                if (kmer_counts.find(mutated) != kmer_counts.end()) {
                    h1_count++;
                }
            }
        }
        
        h1_map[kmer] = h1_count;
    }
    
    return h1_map;
}
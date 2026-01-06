// Continuation of kmer_sketch_tool implementation
// Part 2: SketchBuilder and QueryEngine

#include "kmer_sketch_tool.h"
#include <cmath>
#include <thread>
#include <mutex>
#include <iostream>

// ============================================================================
// SketchBuilder - Part 2
// ============================================================================

KmerSketch SketchBuilder::build_from_string(const string& seq, const string& id) const {
    KmerSketch sketch;
    sketch.id = id;
    sketch.k = k;
    sketch.theta = theta;
    sketch.has_h1 = compute_h1;
    
    if (seq.length() < static_cast<size_t>(k)) {
        cerr << "Warning: Sequence too short for k=" << k << endl;
        return sketch;
    }
    
    // Count total valid k-mer positions (excluding those with N)
    size_t valid_kmers = 0;
    size_t skipped_kmers = 0;
    
    // Extract and encode k-mers
    unordered_map<uint64_t, uint16_t> all_kmers;  // Before sketching
    
    for (size_t i = 0; i <= seq.length() - k; i++) {
        string kmer_str = seq.substr(i, k);
        
        // Skip k-mers with N or other non-ACGT characters
        if (!DNAEncoder::is_valid_kmer(kmer_str)) {
            skipped_kmers++;
            continue;
        }
        
        valid_kmers++;
        
        // Encode k-mer
        uint64_t kmer_encoded = DNAEncoder::encode_kmer(kmer_str);
        
        // Double check (should not happen after is_valid_kmer check)
        if (kmer_encoded == UINT64_MAX) {
            skipped_kmers++;
            continue;
        }
        
        // Get canonical k-mer
        uint64_t canonical_kmer = DNAEncoder::canonical(kmer_encoded, k);
        
        // Store in all_kmers (for h1 computation if needed)
        if (compute_h1) {
            all_kmers[canonical_kmer]++;
        }
        
        // Check if should keep in sketch
        if (should_keep_kmer(canonical_kmer)) {
            sketch.kmer_counts[canonical_kmer]++;
        }
    }
    
    sketch.total_kmers = valid_kmers;
    
    if (skipped_kmers > 0) {
        cerr << "Info: Skipped " << skipped_kmers << " k-mers with N/non-ACGT in " 
             << id << " (" << valid_kmers << " valid k-mers)\n";
    }
    
    // Compute h1 if requested
    if (compute_h1 && !all_kmers.empty()) {
        // Compute h1 on ALL k-mers (not just sketched ones)
        auto h1_map = compute_h1_values(all_kmers);
        
        // Compute average: sum(occ * h1) / L
        double sum_occ_h1 = 0.0;
        for (const auto& [kmer, occ] : all_kmers) {
            if (h1_map.find(kmer) != h1_map.end()) {
                sum_occ_h1 += static_cast<double>(occ) * static_cast<double>(h1_map[kmer]);
            }
        }
        
        sketch.h1_avg = sum_occ_h1 / static_cast<double>(sketch.total_kmers);
    }
    
    return sketch;
}

KmerSketch SketchBuilder::build_sketch(const string& fasta_file, const string& sketch_id) const {
    auto sequences = read_fasta(fasta_file);
    
    if (sequences.empty()) {
        cerr << "Error: No sequences found in " << fasta_file << endl;
        return KmerSketch();
    }
    
    // Use first sequence or concatenate all?
    // For now, concatenate all sequences
    string combined_seq;
    for (const auto& [id, seq] : sequences) {
        combined_seq += seq;
    }
    
    string id = sketch_id.empty() ? sequences[0].first : sketch_id;
    return build_from_string(combined_seq, id);
}

SketchDatabase SketchBuilder::build_database(const vector<string>& fasta_files) const {
    SketchDatabase db;
    
    if (num_threads <= 1) {
        // Single-threaded
        for (const auto& file : fasta_files) {
            cout << "Processing " << file << "..." << endl;
            KmerSketch sketch = build_sketch(file, "");
            db.add_sketch(sketch);
        }
    } else {
        // Multi-threaded
        vector<KmerSketch> sketches(fasta_files.size());
        vector<thread> threads;
        mutex cout_mutex;
        
        auto worker = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; i++) {
                {
                    lock_guard<mutex> lock(cout_mutex);
                    cout << "Thread processing " << fasta_files[i] << "..." << endl;
                }
                sketches[i] = build_sketch(fasta_files[i], "");
            }
        };
        
        size_t chunk_size = (fasta_files.size() + num_threads - 1) / num_threads;
        for (int t = 0; t < num_threads; t++) {
            size_t start = t * chunk_size;
            size_t end = min(start + chunk_size, fasta_files.size());
            if (start < end) {
                threads.emplace_back(worker, start, end);
            }
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
        
        for (auto& sketch : sketches) {
            db.add_sketch(sketch);
        }
    }
    
    return db;
}

// ============================================================================
// QueryEngine
// ============================================================================

double QueryEngine::compute_r_sm(size_t D_obs, size_t L, int k) const {
    if (L == 0) return 0.0;
    
    double q_sm = static_cast<double>(D_obs) / static_cast<double>(L);
    
    // Clamp q to [0, 1)
    if (q_sm >= 1.0) q_sm = 0.999999;
    if (q_sm < 0.0) q_sm = 0.0;
    
    double r_sm = 1.0 - pow(1.0 - q_sm, 1.0 / k);
    return r_sm;
}

double QueryEngine::compute_r_strong(size_t D_obs, size_t L, int k, double h1_avg) const {
    if (L == 0) return 0.0;
    
    // First get r_sm as initial estimate
    double r_sm = compute_r_sm(D_obs, L, k);
    
    // Bias correction term
    double bias_correction = h1_avg * pow(1.0 - r_sm, k - 1) * r_sm / 3.0;
    
    // Corrected q
    double q_strong = (static_cast<double>(D_obs) / static_cast<double>(L)) + bias_correction;
    
    // Clamp q to [0, 1)
    if (q_strong >= 1.0) q_strong = 0.999999;
    if (q_strong < 0.0) q_strong = 0.0;
    
    double r_strong = 1.0 - pow(1.0 - q_strong, 1.0 / k);
    return r_strong;
}

vector<QueryResult> QueryEngine::query(const string& sequence, const string& query_id) const {
    vector<QueryResult> results;
    
    if (sequence.length() < static_cast<size_t>(k)) {
        cerr << "Warning: Query sequence too short for k=" << k << endl;
        return results;
    }
    
    // Extract k-mers from query with multiplicity
    unordered_map<uint64_t, uint16_t> query_kmers;
    size_t valid_kmers = 0;
    size_t skipped_kmers = 0;
    
    for (size_t i = 0; i <= sequence.length() - k; i++) {
        string kmer_str = sequence.substr(i, k);
        
        // Skip k-mers with N or other non-ACGT characters
        if (!DNAEncoder::is_valid_kmer(kmer_str)) {
            skipped_kmers++;
            continue;
        }
        
        valid_kmers++;
        
        uint64_t kmer_encoded = DNAEncoder::encode_kmer(kmer_str);
        
        // Double check
        if (kmer_encoded == UINT64_MAX) {
            skipped_kmers++;
            continue;
        }
        
        uint64_t canonical_kmer = DNAEncoder::canonical(kmer_encoded, k);
        
        if (should_keep_kmer(canonical_kmer)) {
            query_kmers[canonical_kmer]++;
        }
    }
    
    size_t total_query_kmers = valid_kmers;
    
    if (skipped_kmers > 0) {
        cerr << "Info: Skipped " << skipped_kmers << " k-mers with N/non-ACGT in query " 
             << query_id << " (" << valid_kmers << " valid k-mers)\n";
    }
    
    // Scale up if using sketch (theta < 1)
    size_t L = total_query_kmers;
    if (theta < 1.0) {
        L = static_cast<size_t>(total_query_kmers / theta);
    }
    
    // Query against each sketch in database
    for (const auto& target_sketch : db.get_sketches()) {
        QueryResult result;
        result.query_id = query_id;
        result.target_id = target_sketch.id;
        result.query_total = total_query_kmers;
        
        // Compute intersection and difference
        size_t shared = 0;
        size_t novel = 0;
        
        for (const auto& [kmer, count] : query_kmers) {
            auto it = target_sketch.kmer_counts.find(kmer);
            if (it != target_sketch.kmer_counts.end()) {
                // K-mer is shared
                shared += min(count, it->second);
            } else {
                // K-mer is novel
                novel += count;
            }
        }
        
        // Scale up if using sketch
        if (theta < 1.0) {
            novel = static_cast<size_t>(novel / theta);
            shared = static_cast<size_t>(shared / theta);
        }
        
        result.shared_kmers = shared;
        result.novel_kmers = novel;
        
        // Compute r_sm
        result.r_sm = compute_r_sm(novel, L, k);
        
        // Compute r_strong if available
        if (target_sketch.has_h1) {
            result.r_strong = compute_r_strong(novel, L, k, target_sketch.h1_avg);
            result.has_strong = true;
        } else {
            result.r_strong = 0.0;
            result.has_strong = false;
        }
        
        results.push_back(result);
    }
    
    return results;
}

bool QueryEngine::should_keep_kmer(uint64_t kmer_encoded) const {
    if (theta >= 1.0) return true;
    
    uint64_t hash_val = murmurhash3_64(kmer_encoded, seed);
    double hash_frac = static_cast<double>(hash_val) / static_cast<double>(UINT64_MAX);
    return hash_frac <= theta;
}

vector<QueryResult> QueryEngine::query_fasta(const string& fasta_file, QueryMode mode) const {
    auto sequences = read_fasta(fasta_file);
    vector<QueryResult> all_results;
    
    if (mode == QueryMode::PER_FILE) {
        // Mode 1: Concatenate all sequences in the file
        string combined_seq;
        string combined_id = fasta_file;
        
        // Extract just the filename without path
        size_t last_slash = fasta_file.find_last_of("/\\");
        if (last_slash != string::npos) {
            combined_id = fasta_file.substr(last_slash + 1);
        }
        
        // Concatenate all sequences
        for (const auto& [id, seq] : sequences) {
            combined_seq += seq;
        }
        
        // Query once with combined sequence
        auto results = query(combined_seq, combined_id);
        all_results.insert(all_results.end(), results.begin(), results.end());
        
    } else {
        // Mode 2: Query each sequence separately
        for (const auto& [id, seq] : sequences) {
            auto results = query(seq, id);
            all_results.insert(all_results.end(), results.begin(), results.end());
        }
    }
    
    return all_results;
}

vector<vector<QueryResult>> QueryEngine::batch_query(
    const vector<string>& sequences,
    const vector<string>& query_ids
) const {
    vector<vector<QueryResult>> all_results(sequences.size());
    
    if (num_threads <= 1) {
        // Single-threaded
        for (size_t i = 0; i < sequences.size(); i++) {
            all_results[i] = query(sequences[i], query_ids[i]);
        }
    } else {
        // Multi-threaded
        vector<thread> threads;
        mutex cout_mutex;
        
        auto worker = [&](size_t start, size_t end) {
            for (size_t i = start; i < end; i++) {
                all_results[i] = query(sequences[i], query_ids[i]);
                {
                    lock_guard<mutex> lock(cout_mutex);
                    cout << "Processed query " << query_ids[i] << endl;
                }
            }
        };
        
        size_t chunk_size = (sequences.size() + num_threads - 1) / num_threads;
        for (int t = 0; t < num_threads; t++) {
            size_t start = t * chunk_size;
            size_t end = min(start + chunk_size, sequences.size());
            if (start < end) {
                threads.emplace_back(worker, start, end);
            }
        }
        
        for (auto& thread : threads) {
            thread.join();
        }
    }
    
    return all_results;
}
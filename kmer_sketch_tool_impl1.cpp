#include "kmer_sketch_tool.h"
#include <cmath>
#include <thread>
#include <mutex>
#include <iostream>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

QueryEngine::QueryEngine(const DatabaseMetadata& meta, int threads, bool pp_flag)
    : metadata(meta), num_threads(threads), pp_mode(pp_flag) {
    
    // Read all parameters from metadata
    k = metadata.k;
    theta = metadata.theta;
    seed = metadata.seed;
    
    if (theta >= 1.0) {
        hash_threshold = UINT64_MAX;
    } else {
        hash_threshold = static_cast<uint64_t>(theta * static_cast<double>(UINT64_MAX));
    }
    
    cout << "Information initialized from database:\n";
    cout << "  k = " << k << "\n";
    cout << "  theta = " << theta << "\n";
    cout << "  seed = " << seed << "\n";
    cout << "  has_h1 = " << (metadata.has_h1 ? "yes" : "no") << "\n";
    cout << "  pp_mode = " << (pp_mode ? "yes" : "no") << "\n";
}

KmerSketch SketchBuilder::build_from_string(const string& seq, const string& id) const {
    KmerSketch sketch;
    sketch.id = id;
    sketch.k = k;
    sketch.theta = theta;
    sketch.seed = seed;
    sketch.has_h1 = compute_h1;
    
    if (seq.length() < static_cast<size_t>(k)) {
        cerr << "Warning: Sequence too short for k=" << k << endl;
        return sketch;
    }
    
    size_t total_positions = seq.length() - k + 1;
    
    if (num_threads <= 1 || total_positions < 10000) {
        size_t valid_kmers = 0;
        size_t skipped_kmers = 0;
        unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash> all_kmers;
        
        for (size_t i = 0; i < total_positions; i++) {
            string kmer_str = seq.substr(i, k);
            
            if (!DNAEncoder::is_valid_kmer(kmer_str)) {
                skipped_kmers++;
                continue;
            }
            
            valid_kmers++;
            kmer_t kmer_encoded = DNAEncoder::encode_kmer(kmer_str);
            if (kmer_encoded == ~kmer_t(0)) {
                skipped_kmers++;
                continue;
            }
            
            kmer_t canonical_kmer = DNAEncoder::canonical(kmer_encoded, k);
            
            if (compute_h1) {
                all_kmers[canonical_kmer]++;
            }
            
            if (should_keep_kmer(canonical_kmer)) {
                sketch.kmer_counts[canonical_kmer]++;
            }
        }
        
        sketch.total_kmers = valid_kmers;
        
        if (skipped_kmers > 0) {
            cerr << "Info: Skipped " << skipped_kmers << " k-mers with N/non-ACGT in " 
                 << id << " (" << valid_kmers << " valid k-mers)\n";
        }
        
        if (compute_h1 && !all_kmers.empty()) {
            auto h1_map = compute_h1_values(all_kmers);
            double sum_occ_h1 = 0.0;
            for (const auto& [kmer, occ] : all_kmers) {
                if (h1_map.find(kmer) != h1_map.end()) {
                    sum_occ_h1 += static_cast<double>(occ) * static_cast<double>(h1_map[kmer]);
                }
            }
            sketch.h1_avg = sum_occ_h1 / static_cast<double>(sketch.total_kmers);
        }
        
    } else {
        // cout << "  Using " << num_threads << " threads for k-mer extraction...\n";
        
        vector<unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash>> thread_sketch_kmers(num_threads);
        vector<unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash>> thread_all_kmers(num_threads);
        vector<size_t> thread_valid_kmers(num_threads, 0);
        vector<size_t> thread_skipped_kmers(num_threads, 0);
        
        vector<thread> threads;
        size_t chunk_size = (total_positions + num_threads - 1) / num_threads;
        
        auto worker = [&](int thread_id) {
            size_t start = thread_id * chunk_size;
            size_t end = std::min(start + chunk_size, total_positions);
            
            for (size_t i = start; i < end; i++) {
                string kmer_str = seq.substr(i, k);
                
                if (!DNAEncoder::is_valid_kmer(kmer_str)) {
                    thread_skipped_kmers[thread_id]++;
                    continue;
                }
                
                thread_valid_kmers[thread_id]++;
                kmer_t kmer_encoded = DNAEncoder::encode_kmer(kmer_str);
                if (kmer_encoded == ~kmer_t(0)) {
                    thread_skipped_kmers[thread_id]++;
                    continue;
                }
                
                kmer_t canonical_kmer = DNAEncoder::canonical(kmer_encoded, k);
                
                if (compute_h1) {
                    thread_all_kmers[thread_id][canonical_kmer]++;
                }
                
                if (should_keep_kmer(canonical_kmer)) {
                    thread_sketch_kmers[thread_id][canonical_kmer]++;
                }
            }
        };
        
        for (int t = 0; t < num_threads; t++) {
            threads.emplace_back(worker, t);
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        // cout << "  Merging results from " << num_threads << " threads...\n";
        
        size_t total_valid = 0;
        size_t total_skipped = 0;
        unordered_map<kmer_t, uint16_t, KmerUtils::KmerHash> all_kmers;
        
        for (int t = 0; t < num_threads; t++) {
            total_valid += thread_valid_kmers[t];
            total_skipped += thread_skipped_kmers[t];
            
            for (const auto& [kmer, count] : thread_sketch_kmers[t]) {
                sketch.kmer_counts[kmer] += count;
            }
            
            if (compute_h1) {
                for (const auto& [kmer, count] : thread_all_kmers[t]) {
                    all_kmers[kmer] += count;
                }
            }
        }
        
        sketch.total_kmers = total_valid;
        
        if (total_skipped > 0) {
            cerr << "Info: Skipped " << total_skipped << " k-mers with N/non-ACGT in " 
                 << id << " (" << total_valid << " valid k-mers)\n";
        }
        
        if (compute_h1 && !all_kmers.empty()) {
            cout << "  Computing h1 statistics...\n";
            
            auto h1_map = compute_h1_values(all_kmers);
            double sum_occ_h1 = 0.0;
            for (const auto& [kmer, occ] : all_kmers) {
                if (h1_map.find(kmer) != h1_map.end()) {
                    sum_occ_h1 += static_cast<double>(occ) * static_cast<double>(h1_map[kmer]);
                }
            }
            sketch.h1_avg = sum_occ_h1 / static_cast<double>(sketch.total_kmers);
        }
    }
    
    return sketch;
}

KmerSketch SketchBuilder::build_sketch(const string& fasta_file, const string& sketch_id) const {
    auto sequences = read_fasta(fasta_file);
    
    if (sequences.empty()) {
        cerr << "Error: No sequences found in " << fasta_file << endl;
        return KmerSketch();
    }
    
    string combined_seq;
    for (const auto& [id, seq] : sequences) {
        combined_seq += seq;
    }
    
    string id;
    if (!sketch_id.empty()) {
        id = sketch_id;
    } else {
        size_t last_slash = fasta_file.find_last_of("/\\");
        size_t start = (last_slash == string::npos) ? 0 : last_slash + 1;
        id = fasta_file.substr(start);
        
        size_t last_dot = id.find_last_of('.');
        if (last_dot != string::npos) {
            id = id.substr(0, last_dot);
        }
    }
    
    return build_from_string(combined_seq, id);
}

void SketchBuilder::build_database_incremental(const vector<string>& fasta_files,
                                                const string& output_dir,
                                                SketchMode mode) const {
    if (mode == SketchMode::COMBINED) {
        cout << "Building combined sketch from " << fasta_files.size() << " files...\n";
        
        string combined_seq;
        string combined_id = "combined";
        
        for (const auto& file : fasta_files) {
            cout << "  Reading " << file << "...\n";
            auto sequences = read_fasta(file);
            for (const auto& [id, seq] : sequences) {
                combined_seq += seq;
            }
        }
        
        cout << "Building sketch from combined sequence (" 
             << combined_seq.length() << " bp)...\n";
        KmerSketch sketch = build_from_string(combined_seq, combined_id);
        
        fs::create_directories(output_dir);
        sketch.save(output_dir + "/sketch_0.bin");
        
        ofstream index(output_dir + "/index.txt");
        index << sketch.id << "\tsketch_0.bin\n";
        index.close();
        
        cout << "\nSaved 1 sketch to " << output_dir << endl;
        cout << "Total memory: " << format_memory(sketch.memory_usage()) << endl;
        
    } else {
        cout << "Building individual sketches from " << fasta_files.size() << " files...\n";
        
        fs::create_directories(output_dir);
        ofstream index(output_dir + "/index.txt");
        
        if (num_threads <= 1) {
            size_t total_memory = 0;
            
            for (size_t i = 0; i < fasta_files.size(); i++) {
                cout << "Processing [" << (i+1) << "/" << fasta_files.size() << "] " 
                     << fasta_files[i] << "..." << endl;
                
                KmerSketch sketch = build_sketch(fasta_files[i], "");
                
                if (sketch.k > 0) {
                    string sketch_file = "sketch_" + to_string(i) + ".bin";
                    sketch.save(output_dir + "/" + sketch_file);
                    index << sketch.id << "\t" << sketch_file << "\n";
                    
                    total_memory += sketch.memory_usage();
                    
                    cout << "  Written sketch '" << sketch.id << "' (" 
                         << format_memory(sketch.memory_usage()) << ")\n";
                }
            }
            
            index.close();
            cout << "\nSaved " << fasta_files.size() << " sketches to " << output_dir << endl;
            cout << "Total size on disk: " << format_memory(total_memory) << endl;
            
        } else {
            // cout << "  Using " << num_threads << " threads with batch processing\n";
            
            size_t batch_size = num_threads * 2;
            size_t sketch_counter = 0;
            size_t total_memory = 0;
            mutex index_mutex;
            
            for (size_t batch_start = 0; batch_start < fasta_files.size(); batch_start += batch_size) {
                size_t batch_end = min(batch_start + batch_size, fasta_files.size());
                size_t current_batch_size = batch_end - batch_start;
                
                // cout << "  Processing batch " << (batch_start / batch_size + 1) 
                //      << " (files " << batch_start << "-" << (batch_end-1) << ")\n";
                
                vector<KmerSketch> batch_sketches(current_batch_size);
                vector<thread> worker_threads;
                mutex cout_mutex;
                
                auto worker = [&](size_t local_start, size_t local_end) {
                    for (size_t i = local_start; i < local_end; i++) {
                        size_t file_idx = batch_start + i;
                        {
                            lock_guard<mutex> lock(cout_mutex);
                            cout << "    Thread processing [" << (file_idx+1) << "/" 
                                 << fasta_files.size() << "] " 
                                 << fasta_files[file_idx] << "..." << endl;
                        }
                        batch_sketches[i] = build_sketch(fasta_files[file_idx], "");
                    }
                };
                
                size_t chunk_size = (current_batch_size + num_threads - 1) / num_threads;
                for (int t = 0; t < num_threads; t++) {
                    size_t local_start = t * chunk_size;
                    size_t local_end = min(local_start + chunk_size, current_batch_size);
                    if (local_start < local_end) {
                        worker_threads.emplace_back(worker, local_start, local_end);
                    }
                }
                
                for (auto& t : worker_threads) {
                    t.join();
                }
                
                // cout << "  Writing batch to disk...\n";
                for (const auto& sketch : batch_sketches) {
                    if (sketch.k > 0) {
                        string sketch_file = "sketch_" + to_string(sketch_counter) + ".bin";
                        sketch.save(output_dir + "/" + sketch_file);
                        
                        lock_guard<mutex> lock(index_mutex);
                        index << sketch.id << "\t" << sketch_file << "\n";
                        
                        total_memory += sketch.memory_usage();
                        sketch_counter++;
                    }
                }
                
                // cout << "  Batch written, total sketches so far: " << sketch_counter << "\n";
            }
            
            index.close();
            cout << "\nSaved " << sketch_counter << " sketches to " << output_dir << endl;
            cout << "Total size on disk: " << format_memory(total_memory) << endl;
        }
    }
}

double QueryEngine::compute_r_pc(size_t D_obs, size_t L, int k) const {
    if (L == 0) return 0.0;
    
    double q_pc = static_cast<double>(D_obs) / static_cast<double>(L);
    
    if (q_pc >= 1.0) q_pc = 1.0;
    if (q_pc < 0.0) q_pc = 0.0;
    
    double r_pc = 1.0 - pow(1.0 - q_pc, 1.0 / k);
    return r_pc;
}

double QueryEngine::compute_r_cc(size_t D_obs, size_t L, int k, double h1_avg) const {
    if (L == 0) return 0.0;
    
    double r_pc = compute_r_pc(D_obs, L, k);
    
    double bias_correction = h1_avg * pow(1.0 - r_pc, k - 1) * r_pc / 3.0;
    
    double q_cc = (static_cast<double>(D_obs) / static_cast<double>(L)) + bias_correction;
    if (q_cc >= 1.0) q_cc = (static_cast<double>(D_obs) / static_cast<double>(L));
    if (q_cc < 0.0) q_cc = 0.0;
    
    double r_cc = 1.0 - pow(1.0 - q_cc, 1.0 / k);
    return r_cc;
}

double QueryEngine::compute_r_pp(size_t D_obs_sketched, size_t L, double theta) const {
    if (L == 0 || theta <= 0.0) return 0.0;
    
    // q = novel_sketched_kmers / (L * theta)
    double q_pp = static_cast<double>(D_obs_sketched) / (static_cast<double>(L) * theta);
    
    if (q_pp >= 1.0) q_pp = 1.0;
    if (q_pp < 0.0) q_pp = 0.0;
    
    double r_pp = 1.0 - pow(1.0 - q_pp, 1.0 / k);
    return r_pp;
}

QueryKmers QueryEngine::extract_query_kmers(const string& sequence, const string& query_id) const {
    QueryKmers result;
    result.query_id = query_id;
    
    if (sequence.length() < static_cast<size_t>(k)) {
        cerr << "Warning: Query sequence too short for k=" << k << endl;
        return result;
    }
    
    size_t valid_kmers = 0;
    size_t skipped_kmers = 0;
    
    for (size_t i = 0; i <= sequence.length() - k; i++) {
        string kmer_str = sequence.substr(i, k);
        
        if (!DNAEncoder::is_valid_kmer(kmer_str)) {
            skipped_kmers++;
            continue;
        }
        
        valid_kmers++;
        
        kmer_t kmer_encoded = DNAEncoder::encode_kmer(kmer_str);
        if (kmer_encoded == ~kmer_t(0)) {
            skipped_kmers++;
            continue;
        }
        
        kmer_t canonical_kmer = DNAEncoder::canonical(kmer_encoded, k);
        
        if (should_keep_kmer(canonical_kmer)) {
            result.kmer_counts[canonical_kmer]++;
        }
    }
    
    result.total_kmers = valid_kmers;
    
    for (const auto& [kmer, count] : result.kmer_counts) {
        result.sketched_kmer_count += count;
    }
    
    if (skipped_kmers > 0) {
        cerr << "Info: Skipped " << skipped_kmers << " k-mers with N/non-ACGT in query " 
             << query_id << " (" << valid_kmers << " valid k-mers)\n";
    }
    
    return result;
}

QueryResult QueryEngine::compare_with_sketch(const QueryKmers& query, 
                                              const KmerSketch& target) const {
    QueryResult result;
    result.query_id = query.query_id;
    result.target_id = target.id;
    result.query_total = query.total_kmers;
    
    size_t shared_sketched = 0;
    size_t novel_sketched = 0;
    
    for (const auto& [kmer, count] : query.kmer_counts) {
        auto it = target.kmer_counts.find(kmer);
        if (it != target.kmer_counts.end()) {
            shared_sketched += count;
        } else {
            novel_sketched += count;
        }
    }
    
    // Scale back to estimate full k-mer counts
    size_t shared = shared_sketched;
    size_t novel = novel_sketched;
    
    if (theta < 1.0) {
        shared = static_cast<size_t>(static_cast<double>(shared_sketched) / theta);
        novel = static_cast<size_t>(static_cast<double>(novel_sketched) / theta);
    }
    
    // Ensure shared + novel = total (handle rounding errors)
    size_t estimated_total = shared + novel;
    if (estimated_total > query.total_kmers) {
        double scale = static_cast<double>(query.total_kmers) / estimated_total;
        shared = static_cast<size_t>(shared * scale);
        novel = query.total_kmers - shared;
    } else if (estimated_total < query.total_kmers) {
        novel = query.total_kmers - shared;
    }
    
    result.shared_kmers = shared;
    result.novel_kmers = novel;
    if (novel >= query.total_kmers){
        novel = query.total_kmers - shared;
    }
    
    result.r_pc = compute_r_pc(novel, query.total_kmers, k);
    
    if (target.has_h1) {
        result.r_cc = compute_r_cc(novel, query.total_kmers, k, target.h1_avg);
        result.has_cc = true;
    } else {
        result.r_cc = 0.0;
        result.has_cc = false;
    }
    
    // NEW: Compute r_pp if pp_mode is enabled
    if (pp_mode) {
        result.r_pp = compute_r_pp(novel_sketched, query.total_kmers, theta);
        result.has_pp = true;
    } else {
        result.r_pp = 0.0;
        result.has_pp = false;
    }
    
    return result;
}

void QueryEngine::query_streaming(const vector<string>& fasta_files,
                                   QueryMode mode,
                                   const string& output_file) const {
    cout << "Step 1: Extracting k-mers from query sequences...\n";
    vector<QueryKmers> all_queries;
    
    for (size_t i = 0; i < fasta_files.size(); i++) {
        cout << "  Processing query file [" << (i+1) << "/" << fasta_files.size() 
             << "]: " << fasta_files[i] << "\n";
        
        auto sequences = read_fasta(fasta_files[i]);
        
        if (mode == QueryMode::PER_FILE || mode == QueryMode::BATCH) {
            string combined_seq;
            string combined_id;
            
            size_t last_slash = fasta_files[i].find_last_of("/\\");
            size_t start = (last_slash == string::npos) ? 0 : last_slash + 1;
            combined_id = fasta_files[i].substr(start);
            
            size_t last_dot = combined_id.find_last_of('.');
            if (last_dot != string::npos) {
                combined_id = combined_id.substr(0, last_dot);
            }
            
            for (const auto& [id, seq] : sequences) {
                combined_seq += seq;
            }
            
            all_queries.push_back(extract_query_kmers(combined_seq, combined_id));
            
        } else {
            for (const auto& [id, seq] : sequences) {
                all_queries.push_back(extract_query_kmers(seq, id));
            }
        }
    }
    
    cout << "  Extracted k-mers from " << all_queries.size() << " queries\n";
    
    size_t query_memory = 0;
    for (const auto& q : all_queries) {
        query_memory += q.kmer_counts.size() * (sizeof(kmer_t) + sizeof(uint16_t));
    }
    cout << "  Query k-mers memory: " << format_memory(query_memory) << "\n\n";
    
    ofstream out;
    bool use_stdout = output_file.empty();
    
    if (!use_stdout) {
        out.open(output_file);
        if (!out.is_open()) {
            cerr << "Error: Cannot write to " << output_file << "\n";
            return;
        }
    }
    
    ostream& output_stream = use_stdout ? cout : out;
    
    output_stream << "query_id\ttarget_id\tshared_kmers\tnovel_kmers\tquery_total\t"
                  << "r_pc\tr_cc\thas_cc";
    
    if (pp_mode) {
        output_stream << "\tr_pp\thas_pp\n";
    } else {
        output_stream << "\n";
    }
    
    // cout << "Step 2: Streaming through database (" << metadata.size() << " sketches)...\n";
    
    size_t total_comparisons = 0;
    
    for (size_t sketch_idx = 0; sketch_idx < metadata.sketch_entries.size(); sketch_idx++) {
        const auto& [sketch_id, sketch_file] = metadata.sketch_entries[sketch_idx];
        
        cout << "  Loading sketch [" << (sketch_idx+1) << "/" << metadata.size() 
             << "]: " << sketch_id << "\n";
        
        string sketch_path = metadata.db_path + "/" + sketch_file;
        KmerSketch target = KmerSketch::load(sketch_path);
        
        if (target.k == 0) {
            cerr << "    Warning: Failed to load sketch " << sketch_file << "\n";
            continue;
        }
        
        for (const auto& query : all_queries) {
            QueryResult result = compare_with_sketch(query, target);
            
            output_stream << result.query_id << "\t"
                         << result.target_id << "\t"
                         << result.shared_kmers << "\t"
                         << result.novel_kmers << "\t"
                         << result.query_total << "\t"
                         << scientific << setprecision(6) << result.r_pc << "\t";
            
            if (result.has_cc) {
                output_stream << scientific << setprecision(6) << result.r_cc << "\tyes";
            } else {
                output_stream << "NA\tno";
            }
            
            if (pp_mode) {
                output_stream << "\t";
                if (result.has_pp) {
                    output_stream << scientific << setprecision(6) << result.r_pp << "\tyes\n";
                } else {
                    output_stream << "NA\tno\n";
                }
            } else {
                output_stream << "\n";
            }
            
            total_comparisons++;
        }
    }
    
    if (!use_stdout) {
        out.close();
        cout << "\nResults written to " << output_file << "\n";
    }
    
    cout << "Completed " << total_comparisons << " query-target comparisons\n";
}
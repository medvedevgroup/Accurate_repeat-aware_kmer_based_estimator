#include "kmer_sketch_tool.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

void print_usage() {
    cout << "Repeat-Robust Estimator - Fast mutation rate estimation\n\n";
    cout << "Usage:\n";
    cout << "  Build database:\n";
    cout << "    repeat_robust_estimator sketch -i <input> -o <output_dir> [options]\n\n";
    cout << "  Query database:\n";
    cout << "    repeat_robust_estimator query -d <database_dir> -q <query_fasta> [options]\n\n";
    cout << "Options:\n";
    cout << "  Common:\n";
    cout << "    -p INT        Number of threads (default: 1)\n";
    cout << "    -h, --help    Print this help message\n\n";
    cout << "  Sketch mode:\n";
    cout << "    -k INT        K-mer size (default: 21, max: " << KMER_MAX_K << ")\n";
    cout << "    -s INT        Random seed (default: 42)\n";
    cout << "    -i PATH       Input can be:\n";
    cout << "                    - Directory containing FASTA files\n";
    cout << "                    - Single FASTA file (.fa, .fasta, .fna, etc.)\n";
    cout << "                    - Text file with list of FASTA paths (one per line)\n";
    cout << "                    - Multiple files (can repeat -i or use wildcards)\n";
    cout << "                  Note: Each FASTA file becomes one sketch (all sequences\n";
    cout << "                        within a file are concatenated). Sketch ID is the\n";
    cout << "                        filename without path and extension.\n";
    cout << "                  Memory efficient: sketches are written incrementally.\n";
    cout << "    -o DIR        Output database directory\n";
    cout << "    -t FLOAT      Sketch fraction theta (default: 1.0, no sketching)\n";
    cout << "    --h1          Compute h1 statistics (for r_cc estimator, default: off)\n";
    cout << "    --sketch-mode MODE  Sketch mode (default: individual)\n";
    cout << "                    individual: Each file -> separate sketch\n";
    cout << "                    combined:   All files -> one sketch\n\n";
    cout << "  Query mode:\n";
    cout << "    -d DIR        Database directory\n";
    cout << "    -q FILE(S)    Query FASTA file(s). Can specify multiple -q options\n";
    cout << "                  or use shell wildcards (e.g., -q *.fasta)\n";
    cout << "    -o FILE       Output results file (default: stdout)\n";
    cout << "    --pp          Enable presence-presence mode (compute r_pp)\n";
    cout << "                    r_pp: presence-presence (query uses sets, DB uses sets)\n";
    cout << "                    r_pc: presence-count (query uses sets, DB uses counts)\n";
    cout << "                    r_cc: count-count (both use counts with h1 correction)\n";
    cout << "    --mode MODE   Query mode (default: file)\n";
    cout << "                    file:     Concatenate all sequences in each query file\n";
    cout << "                              Query ID is filename without path/extension\n";
    cout << "                    sequence: Query each sequence separately\n";
    cout << "                              Query ID is the sequence header\n";
    cout << "                    batch:    Process multiple files (same as file mode)\n";
    cout << "    --top INT     Show only top N results per query (default: all)\n";
    cout << "    --no-streaming  Load entire database into memory (faster but uses more RAM)\n";
    cout << "\n";
    cout << "NOTE: Query mode automatically reads k, theta, and seed from the database.\n";
    cout << "      No need to specify these parameters when querying!\n";
    cout << "\n";
    cout << "Examples:\n";
    cout << "  # Build from directory (each file -> one sketch, incremental writing)\n";
    cout << "  repeat_robust_estimator sketch -i ./genomes/ -o db/ -k 21 -t 0.001 --h1 -s 42 -p 8\n";
    cout << "    Result: genome_A.fasta -> sketch ID 'genome_A'\n";
    cout << "            genome_B.fasta -> sketch ID 'genome_B'\n";
    cout << "            Each sketch is written to disk immediately after creation\n";
    cout << "            Parameters k=21, theta=0.001, seed=42 are saved in sketches\n\n";
    cout << "  # Build from wildcards\n";
    cout << "  repeat_robust_estimator sketch -i /path/*.fasta -o db/ -k 21 --h1\n\n";
    cout << "  # Build combined sketch (all files merged into one)\n";
    cout << "  repeat_robust_estimator sketch -i ./genomes/ -o db/ -k 21 --sketch-mode combined\n";
    cout << "    Result: All files -> single sketch ID 'combined'\n\n";
    cout << "  # Query single file (parameters auto-read from database)\n";
    cout << "  repeat_robust_estimator query -d db/ -q assembly.fasta -o results.tsv -p 4\n";
    cout << "    Result: Query ID is 'assembly' (filename without extension)\n";
    cout << "            Uses k=21, theta=0.001, seed=42 from database automatically\n";
    cout << "            Database sketches loaded one at a time (streaming, memory efficient)\n\n";
    cout << "  # Query with --pp mode (enable r_pp calculation)\n";
    cout << "  repeat_robust_estimator query -d db/ -q assembly.fasta -o results.tsv --pp\n";
    cout << "    Result: Outputs r_pc, r_cc, and r_pp\n\n";
    cout << "  # Query multiple files (batch mode)\n";
    cout << "  repeat_robust_estimator query -d db/ -q *.fasta --mode batch -o results.tsv\n\n";
    cout << "  # Query with full database in memory (faster but more RAM)\n";
    cout << "  repeat_robust_estimator query -d db/ -q query.fasta --no-streaming -o results.tsv\n\n";
    cout << "  # Query each sequence separately (sequence mode)\n";
    cout << "  repeat_robust_estimator query -d db/ -q contigs.fasta --mode sequence -o results.tsv\n";
    cout << "    Result: Each sequence header becomes a separate query ID\n\n";
    cout << "Multi-sequence FASTA files:\n";
    cout << "  During sketch building:\n";
    cout << "    - All sequences in a file are concatenated into one sketch\n";
    cout << "    - Sketch ID = filename (without path and extension)\n";
    cout << "    - Example: chr1.fasta containing chr1, chr2, chr3 -> ID 'chr1'\n";
    cout << "  During query:\n";
    cout << "    - file/batch mode: concatenate all sequences, use filename as ID\n";
    cout << "    - sequence mode: query each sequence separately, use headers as IDs\n\n";
}

void print_sketch_summary(const SketchDatabase& db) {
    cout << "\n" << string(80, '=') << "\n";
    cout << "Sketch Database Summary\n";
    cout << string(80, '=') << "\n";
    cout << "Number of sketches: " << db.size() << "\n";
    cout << "Total memory: " << format_memory(db.total_memory()) << "\n";
    
    if (!db.get_sketches().empty()) {
        const auto& first = db.get_sketches()[0];
        cout << "K-mer size: " << first.k << "\n";
        cout << "Theta: " << first.theta << "\n";
        cout << "Seed: " << first.seed << "\n";
        cout << "H1 stats: " << (first.has_h1 ? "Yes" : "No") << "\n";
    }
    cout << string(80, '=') << "\n\n";
}

void write_result(ostream& out, const QueryResult& result, bool pp_mode) {
    out << result.query_id << "\t"
        << result.target_id << "\t"
        << result.shared_kmers << "\t"
        << result.novel_kmers << "\t"
        << result.query_total << "\t"
        << scientific << setprecision(6) << result.r_pc << "\t";
    
    if (result.has_cc) {
        out << scientific << setprecision(6) << result.r_cc << "\tyes";
    } else {
        out << "NA\tno";
    }
    
    if (pp_mode) {
        out << "\t";
        if (result.has_pp) {
            out << scientific << setprecision(6) << result.r_pp << "\tyes\n";
        } else {
            out << "NA\tno\n";
        }
    } else {
        out << "\n";
    }
}

void print_query_results(const vector<QueryResult>& results, ostream& out, bool pp_mode, int top_n = -1) {
    out << "query_id\ttarget_id\tshared_kmers\tnovel_kmers\tquery_total\t"
        << "r_pc\tr_cc\thas_cc";
    
    if (pp_mode) {
        out << "\tr_pp\thas_pp\n";
    } else {
        out << "\n";
    }
    
    vector<QueryResult> sorted_results = results;
    std::sort(sorted_results.begin(), sorted_results.end(),
         [](const QueryResult& a, const QueryResult& b) {
             return a.r_pc < b.r_pc;
         });
    
    int count = 0;
    for (const auto& result : sorted_results) {
        if (top_n > 0 && count >= top_n) break;
        write_result(out, result, pp_mode);
        count++;
    }
}

bool is_fasta_file(const string& filename) {
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == string::npos) return false;
    
    string ext = filename.substr(dot_pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    return (ext == "fa" || ext == "fasta" || ext == "fna" || 
            ext == "ffn" || ext == "faa" || ext == "frn");
}

vector<string> get_fasta_files_from_dir(const string& dir_path) {
    vector<string> fasta_files;
    
    try {
        for (const auto& entry : fs::directory_iterator(dir_path)) {
            if (entry.is_regular_file()) {
                string filepath = entry.path().string();
                if (is_fasta_file(filepath)) {
                    fasta_files.push_back(filepath);
                }
            }
        }
    } catch (const exception& e) {
        cerr << "Error reading directory " << dir_path << ": " << e.what() << "\n";
    }
    
    return fasta_files;
}

int sketch_mode(int argc, char* argv[]) {
    vector<string> input_paths;
    string output_dir;
    int k = 21;
    double theta = 1.0;
    int num_threads = 1;
    uint32_t seed = 42;
    bool compute_h1 = false;
    string sketch_mode_str = "individual";
    
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                input_paths.push_back(argv[++i]);
            }
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            theta = atof(argv[++i]);
        } else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--h1") == 0) {
            compute_h1 = true;
        } else if (strcmp(argv[i], "--sketch-mode") == 0 && i + 1 < argc) {
            sketch_mode_str = argv[++i];
        }
    }
    
    if (input_paths.empty() || output_dir.empty()) {
        cerr << "Error: -i and -o are required for sketch mode\n";
        print_usage();
        return 1;
    }
    
    if (k < 1 || k > KMER_MAX_K) {
        cerr << "Error: k must be between 1 and " << KMER_MAX_K << "\n";
        return 1;
    }
    
    if (theta <= 0.0 || theta > 1.0) {
        cerr << "Error: theta must be in (0, 1]\n";
        return 1;
    }
    
    vector<string> fasta_files;
    
    for (const string& input_path : input_paths) {
        if (fs::is_directory(input_path)) {
            cout << "Scanning directory: " << input_path << "\n";
            auto dir_files = get_fasta_files_from_dir(input_path);
            fasta_files.insert(fasta_files.end(), dir_files.begin(), dir_files.end());
            
        } else if (fs::is_regular_file(input_path)) {
            if (is_fasta_file(input_path)) {
                fasta_files.push_back(input_path);
            } else {
                cout << "Reading file list: " << input_path << "\n";
                ifstream list_file(input_path);
                if (!list_file.is_open()) {
                    cerr << "Error: Cannot open file: " << input_path << "\n";
                    return 1;
                }
                
                string line;
                while (getline(list_file, line)) {
                    line.erase(0, line.find_first_not_of(" \t\r\n"));
                    line.erase(line.find_last_not_of(" \t\r\n") + 1);
                    
                    if (!line.empty() && line[0] != '#') {
                        fasta_files.push_back(line);
                    }
                }
                list_file.close();
            }
        } else {
            cerr << "Warning: Input path does not exist: " << input_path << "\n";
        }
    }
    
    if (fasta_files.empty()) {
        cerr << "Error: No FASTA files found\n";
        return 1;
    }
    
    cout << "Found " << fasta_files.size() << " FASTA file(s) total\n";
    
    SketchMode smode;
    if (sketch_mode_str == "individual" || sketch_mode_str == "separate") {
        smode = SketchMode::INDIVIDUAL;
    } else if (sketch_mode_str == "combined" || sketch_mode_str == "merge") {
        smode = SketchMode::COMBINED;
    } else {
        cerr << "Error: Invalid sketch mode '" << sketch_mode_str 
             << "'. Use 'individual' or 'combined'\n";
        return 1;
    }
    
    cout << "Building sketch database...\n";
    cout << "  Sketch mode: " << (smode == SketchMode::INDIVIDUAL ? "individual" : "combined") << "\n";
    cout << "  K-mer size: " << k << "\n";
    cout << "  Theta: " << theta << "\n";
    cout << "  Seed: " << seed << "\n";
    cout << "  Compute h1: " << (compute_h1 ? "Yes" : "No") << "\n";
    cout << "  Threads: " << num_threads << "\n";
    cout << "  Input files: " << fasta_files.size() << "\n";
    if (smode == SketchMode::INDIVIDUAL) {
        cout << "  Note: Using incremental disk writing (memory efficient).\n";
        cout << "        Each sketch is written immediately after creation.\n";
        cout << "        Parameters (k, theta, seed) are saved in each sketch.\n";
    }
    cout << "\n";
    
    SketchBuilder builder(k, theta, compute_h1, seed, num_threads);
    builder.build_database_incremental(fasta_files, output_dir, smode);
    
    SketchDatabase db = SketchDatabase::load(output_dir);
    print_sketch_summary(db);
    
    return 0;
}

int query_mode(int argc, char* argv[]) {
    string db_dir, output_file;
    vector<string> query_files;
    int num_threads = 1;
    int top_n = -1;
    string mode_str = "file";
    bool use_streaming = true;
    bool pp_mode = false;  // NEW
    
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            db_dir = argv[++i];
        } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                query_files.push_back(argv[++i]);
            }
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
            top_n = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            mode_str = argv[++i];
        } else if (strcmp(argv[i], "--no-streaming") == 0) {
            use_streaming = false;
        } else if (strcmp(argv[i], "--pp") == 0) {
            pp_mode = true;  // NEW
        }
    }
    
    if (db_dir.empty() || query_files.empty()) {
        cerr << "Error: -d and -q are required for query mode\n";
        print_usage();
        return 1;
    }
    
    QueryMode qmode;
    if (mode_str == "file" || mode_str == "per-file") {
        qmode = QueryMode::PER_FILE;
    } else if (mode_str == "sequence" || mode_str == "per-sequence") {
        qmode = QueryMode::PER_SEQUENCE;
    } else if (mode_str == "batch") {
        qmode = QueryMode::BATCH;
    } else {
        cerr << "Error: Invalid mode '" << mode_str << "'. Use 'file', 'sequence', or 'batch'\n";
        return 1;
    }
    
    if (use_streaming) {
        cout << "Loading database metadata from " << db_dir << "...\n";
        DatabaseMetadata metadata = DatabaseMetadata::load(db_dir);
        
        if (metadata.size() == 0) {
            cerr << "Error: Database is empty or failed to load\n";
            return 1;
        }
        
        cout << "\n" << string(80, '=') << "\n";
        cout << "Database Metadata\n";
        cout << string(80, '=') << "\n";
        cout << "Number of sketches: " << metadata.size() << "\n";
        cout << "K-mer size: " << metadata.k << "\n";
        cout << "Theta: " << metadata.theta << "\n";
        cout << "Seed: " << metadata.seed << "\n";
        cout << "H1 stats: " << (metadata.has_h1 ? "Yes" : "No") << "\n";
        cout << string(80, '=') << "\n\n";
        
        cout << "Using STREAMING query mode (memory efficient)\n";
        cout << "  - Query k-mers will be extracted and kept in memory\n";
        cout << "  - Database sketches will be loaded one at a time\n";
        cout << "  - Results will be written incrementally\n";
        cout << "  - All parameters (k, theta, seed) read from database automatically\n";
        
        if (pp_mode) {
            cout << "  - PP mode ENABLED: Will compute r_pp (presence-presence)\n";
        }
        
        if (qmode == QueryMode::PER_FILE || qmode == QueryMode::BATCH) {
            cout << "  - Query mode: Concatenate all sequences per file\n";
            cout << "  - Query ID: Filename (without path/extension)\n";
        } else {
            cout << "  - Query mode: Each sequence separately\n";
            cout << "  - Query ID: Sequence header\n";
        }
        cout << "\n";
        
        QueryEngine engine(metadata, num_threads, pp_mode);  // Pass pp_mode
        
        cout << "Query will use:\n";
        cout << "  k = " << engine.get_k() << " (from database)\n";
        cout << "  theta = " << engine.get_theta() << " (from database)\n";
        cout << "  seed = " << engine.get_seed() << " (from database)\n";
        cout << "  pp_mode = " << (engine.get_pp_mode() ? "yes" : "no") << "\n";
        cout << "\n";
        
        engine.query_streaming(query_files, qmode, output_file);
        
    } else {
        cerr << "Error: Non-streaming mode not implemented in this version.\n";
        cerr << "Please use streaming mode (default) or remove --no-streaming flag.\n";
        return 1;
    }
    
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    string mode = argv[1];
    
    if (mode == "-h" || mode == "--help") {
        print_usage();
        return 0;
    } else if (mode == "sketch") {
        return sketch_mode(argc, argv);
    } else if (mode == "query") {
        return query_mode(argc, argv);
    } else {
        cerr << "Error: Unknown mode '" << mode << "'\n";
        print_usage();
        return 1;
    }
    
    return 0;
}
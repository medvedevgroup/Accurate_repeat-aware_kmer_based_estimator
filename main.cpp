// main.cpp - Command-line tool for k-mer sketching and querying
//
// Usage:
//   1. Build database: ./kmer_tool sketch -i input_files.txt -o database/ -k 21 -t 0.001 --h1
//   2. Query: ./kmer_tool query -d database/ -q query.fasta -k 21 -t 0.001 -o results.tsv

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
    cout << "K-mer Sketch Tool - Fast mutation rate estimation\n\n";
    cout << "Usage:\n";
    cout << "  Build database:\n";
    cout << "    kmer_tool sketch -i <input> -o <output_dir> [options]\n\n";
    cout << "  Query database:\n";
    cout << "    kmer_tool query -d <database_dir> -q <query_fasta> [options]\n\n";
    cout << "Options:\n";
    cout << "  Common:\n";
    cout << "    -k INT        K-mer size (default: 21, max: 32)\n";
    cout << "    -t FLOAT      Sketch fraction theta (default: 1.0, no sketching)\n";
    cout << "    -p INT        Number of threads (default: 1)\n";
    cout << "    -s INT        Random seed (default: 42)\n";
    cout << "    -h, --help    Print this help message\n\n";
    cout << "  Sketch mode:\n";
    cout << "    -i PATH       Input can be:\n";
    cout << "                    - Directory containing FASTA files\n";
    cout << "                    - Single FASTA file (.fa, .fasta, .fna, etc.)\n";
    cout << "                    - Text file with list of FASTA paths (one per line)\n";
    cout << "    -o DIR        Output database directory\n";
    cout << "    --h1          Compute h1 statistics (for strong estimator)\n";
    cout << "    --no-h1       Do not compute h1 (faster, only r_sm available)\n\n";
    cout << "  Query mode:\n";
    cout << "    -d DIR        Database directory\n";
    cout << "    -q FILE       Query FASTA file\n";
    cout << "    -o FILE       Output results file (default: stdout)\n";
    cout << "    --mode MODE   Query mode (default: file)\n";
    cout << "                    file (or per-file):     One result per FASTA file (concatenate sequences)\n";
    cout << "                    sequence (or per-sequence): One result per sequence in FASTA\n";
    cout << "    --top INT     Show only top N results (default: all)\n\n";
    cout << "Examples:\n";
    cout << "  # Build from directory\n";
    cout << "  kmer_tool sketch -i ./genomes/ -o db/ -k 21 -t 0.001 --h1 -p 8\n\n";
    cout << "  # Build from single FASTA\n";
    cout << "  kmer_tool sketch -i genome.fasta -o db/ -k 21 --h1\n\n";
    cout << "  # Build from file list\n";
    cout << "  kmer_tool sketch -i genomes.txt -o db/ -k 21 -t 1.0 --h1\n\n";
    cout << "  # Query (per-file mode - concatenate all sequences)\n";
    cout << "  kmer_tool query -d db/ -q query.fasta -k 21 --mode file -o results.tsv\n\n";
    cout << "  # Query (per-sequence mode - separate result for each sequence)\n";
    cout << "  kmer_tool query -d db/ -q query.fasta -k 21 --mode sequence -o results.tsv\n\n";
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
        cout << "H1 stats: " << (first.has_h1 ? "Yes" : "No") << "\n";
    }
    cout << string(80, '=') << "\n\n";
}

void print_query_results(const vector<QueryResult>& results, ostream& out, int top_n = -1) {
    // Header
    out << "query_id\ttarget_id\tshared_kmers\tnovel_kmers\tquery_total\t"
        << "r_sm\tr_strong\thas_strong\n";
    
    // Sort by r_sm (ascending - lower is more similar)
    vector<QueryResult> sorted_results = results;
    std::sort(sorted_results.begin(), sorted_results.end(),
         [](const QueryResult& a, const QueryResult& b) {
             return a.r_sm < b.r_sm;
         });
    
    // Output results
    int count = 0;
    for (const auto& result : sorted_results) {
        if (top_n > 0 && count >= top_n) break;
        
        out << result.query_id << "\t"
            << result.target_id << "\t"
            << result.shared_kmers << "\t"
            << result.novel_kmers << "\t"
            << result.query_total << "\t"
            << scientific << setprecision(6) << result.r_sm << "\t";
        
        if (result.has_strong) {
            out << scientific << setprecision(6) << result.r_strong << "\t"
                << "yes\n";
        } else {
            out << "NA\tno\n";
        }
        
        count++;
    }
}

// Helper function to check if file is FASTA
bool is_fasta_file(const string& filename) {
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == string::npos) return false;
    
    string ext = filename.substr(dot_pos + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    return (ext == "fa" || ext == "fasta" || ext == "fna" || 
            ext == "ffn" || ext == "faa" || ext == "frn");
}

// Helper function to get all FASTA files from directory
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
    string input_path, output_dir;
    int k = 21;
    double theta = 1.0;
    int num_threads = 1;
    uint32_t seed = 42;
    bool compute_h1 = false;
    
    // Parse arguments
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_path = argv[++i];
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
        } else if (strcmp(argv[i], "--no-h1") == 0) {
            compute_h1 = false;
        }
    }
    
    // Validate
    if (input_path.empty() || output_dir.empty()) {
        cerr << "Error: -i and -o are required for sketch mode\n";
        print_usage();
        return 1;
    }
    
    if (k < 1 || k > 32) {
        cerr << "Error: k must be between 1 and 32\n";
        return 1;
    }
    
    if (theta <= 0.0 || theta > 1.0) {
        cerr << "Error: theta must be in (0, 1]\n";
        return 1;
    }
    
    // Determine input type and collect FASTA files
    vector<string> fasta_files;
    
    if (fs::is_directory(input_path)) {
        // Input is a directory - get all FASTA files
        cout << "Input is directory, scanning for FASTA files...\n";
        fasta_files = get_fasta_files_from_dir(input_path);
        
        if (fasta_files.empty()) {
            cerr << "Error: No FASTA files found in directory: " << input_path << "\n";
            return 1;
        }
        
        cout << "Found " << fasta_files.size() << " FASTA files\n";
        
    } else if (fs::is_regular_file(input_path)) {
        // Check if it's a FASTA file or a list file
        if (is_fasta_file(input_path)) {
            // Single FASTA file
            cout << "Input is single FASTA file\n";
            fasta_files.push_back(input_path);
            
        } else {
            // Assume it's a text file with list of FASTA files
            cout << "Input is file list\n";
            ifstream list_file(input_path);
            if (!list_file.is_open()) {
                cerr << "Error: Cannot open file: " << input_path << "\n";
                return 1;
            }
            
            string line;
            while (getline(list_file, line)) {
                // Remove leading/trailing whitespace
                line.erase(0, line.find_first_not_of(" \t\r\n"));
                line.erase(line.find_last_not_of(" \t\r\n") + 1);
                
                if (!line.empty() && line[0] != '#') {
                    fasta_files.push_back(line);
                }
            }
            list_file.close();
            
            if (fasta_files.empty()) {
                cerr << "Error: No input files found in list: " << input_path << "\n";
                return 1;
            }
        }
        
    } else {
        cerr << "Error: Input path does not exist or is not accessible: " << input_path << "\n";
        return 1;
    }
    
    cout << "Building sketch database...\n";
    cout << "  K-mer size: " << k << "\n";
    cout << "  Theta: " << theta << "\n";
    cout << "  Compute h1: " << (compute_h1 ? "Yes" : "No") << "\n";
    cout << "  Threads: " << num_threads << "\n";
    cout << "  Input files: " << fasta_files.size() << "\n\n";
    
    // Build database
    SketchBuilder builder(k, theta, compute_h1, seed, num_threads);
    SketchDatabase db = builder.build_database(fasta_files);
    
    // Save database
    db.save(output_dir);
    
    print_sketch_summary(db);
    
    return 0;
}

int query_mode(int argc, char* argv[]) {
    string db_dir, query_file, output_file;
    int k = 21;
    double theta = 1.0;
    int num_threads = 1;
    uint32_t seed = 42;
    int top_n = -1;
    string mode_str = "file";  // Default: per-file mode
    
    // Parse arguments
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            db_dir = argv[++i];
        } else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            query_file = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            theta = atof(argv[++i]);
        } else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
            top_n = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            mode_str = argv[++i];
        }
    }
    
    // Validate
    if (db_dir.empty() || query_file.empty()) {
        cerr << "Error: -d and -q are required for query mode\n";
        print_usage();
        return 1;
    }
    
    // Parse mode
    QueryMode qmode;
    if (mode_str == "file" || mode_str == "per-file") {
        qmode = QueryMode::PER_FILE;
    } else if (mode_str == "sequence" || mode_str == "per-sequence") {
        qmode = QueryMode::PER_SEQUENCE;
    } else {
        cerr << "Error: Invalid mode '" << mode_str << "'. Use 'file' or 'sequence'\n";
        return 1;
    }
    
    // Load database
    cout << "Loading database from " << db_dir << "...\n";
    SketchDatabase db = SketchDatabase::load(db_dir);
    
    if (db.size() == 0) {
        cerr << "Error: Database is empty or failed to load\n";
        return 1;
    }
    
    print_sketch_summary(db);
    
    // Query
    cout << "Querying sequences from " << query_file << "...\n";
    cout << "Query mode: " << (qmode == QueryMode::PER_FILE ? "per-file" : "per-sequence") << "\n";
    
    QueryEngine engine(db, k, theta, seed, num_threads);
    vector<QueryResult> results = engine.query_fasta(query_file, qmode);
    
    cout << "Found " << results.size() << " query-target pairs\n\n";
    
    // Output results
    if (output_file.empty()) {
        print_query_results(results, cout, top_n);
    } else {
        ofstream out(output_file);
        if (!out.is_open()) {
            cerr << "Error: Cannot write to " << output_file << "\n";
            return 1;
        }
        print_query_results(results, out, top_n);
        out.close();
        cout << "Results written to " << output_file << "\n";
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
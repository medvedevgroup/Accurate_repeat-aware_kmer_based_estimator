# Accurate_repeat-aware_kmer_based_estimator

A tool for estimating mutation rates between genomic sequences using k-mer sketching. Designed to handle repetitive sequences robustly with multiple estimation modes.

## Features

- **Three estimation modes:**
  - **r_pc** (presence-count): Query uses k-mer sets, database uses k-mer counts
  - **r_cc** (count-count): Both use k-mer counts with h1 bias correction
  - **r_pp** (presence-presence): Both use k-mer sets

## Installation

### Requirements

- C++17 compatible compiler
- MPFR library (`libmpfr-dev`)
- GMP library (`libgmp-dev`)
- OpenMP

### Build

```bash
make
```

This will create the `repeat_robust_estimator` executable.

## Quick Start

We provide a toy example in the `toy_example/` directory to help you get started quickly.

### Directory Structure

```
toy_example/
├── ref_genomes/       # Reference genomes for building the database
│   ├── genome1.fna
│   └── genome2.fna
├── query_genomes/     # Query genomes for testing
│   └── query1.fna
└── toy_db/            # Output database (created after step 1)
```

### Step 1: Build a sketch database from reference genomes

```bash
./repeat_robust_estimator sketch \
    -i ./toy_example/ref_genomes/ \
    -o ./toy_example/toy_db/ \
    -k 21 \
    -t 0.1 \
    --h1 \
    -p 8
```

**Parameters:**
- `-i`: Input directory containing reference genomes
- `-o`: Output directory for the sketch database
- `-k 21`: Use 21-mers
- `-t 0.1`: FracMinHash sampling rate
- `--h1`: Compute h1 statistics (needed for r_cc estimator)
- `-p 8`: Use 8 threads

### Step 2: Query genomes against the database

```bash
./repeat_robust_estimator query \
    -d ./toy_example/toy_db/ \
    -q ./toy_example/query_genomes/*.fna \
    -o results.tsv
# Query with presence-presence mode
./repeat_robust_estimator query \
    -d ./toy_example/toy_db/ \
    -q ./toy_example/query_genomes/*.fna \
    -o results_with_pp.tsv \
    --pp
```


## Usage

### Sketch Mode

Build a k-mer sketch database from genomic sequences.

```bash
repeat_robust_estimator sketch [options]
```

**Required options:**

- `-i PATH`: Input FASTA file(s) or directory
- `-o DIR`: Output database directory

**Optional parameters:**

- `-k INT`: K-mer size (default: 21, max: 32 or 64 depending on build)

- `-t FLOAT`: Sketch fraction theta (default: 1.0, no sketching)

- `-s INT`: Random seed (default: 42)

- `-p INT`: Number of threads (default: 1)

- `--h1`: Compute h1 statistics (enables r_cc estimator)

- ```
  --sketch-mode MODE
  ```

  - `individual`: Each file → separate sketch (default)
  - `combined`: All files → one merged sketch

**Input formats:**

- Single FASTA file
- Directory containing FASTA files
- Text file listing FASTA paths (one per line)
- Multiple files via wildcards: `-i *.fasta`

**Examples:**

```bash
# Build from directory with h1 statistics
./repeat_robust_estimator sketch -i genomes/ -o db/ -k 21 -t 0.001 --h1 -p 8

# Build from wildcards
./repeat_robust_estimator sketch -i /path/*.fasta -o db/ -k 21 --h1

# Build combined sketch (merge all sequences)
./repeat_robust_estimator sketch -i genomes/ -o db/ -k 21 --sketch-mode combined
```

### Query Mode

Query sequences against a sketch database to estimate mutation rates.

```bash
repeat_robust_estimator query [options]
```

**Required options:**

- `-d DIR`: Database directory
- `-q FILE(S)`: Query FASTA file(s)

**Optional parameters:**

- `-o FILE`: Output results file (default: stdout)

- `-p INT`: Number of threads (default: 1)

- `--pp`: Enable presence-presence mode (compute r_pp)

- ```
  --mode MODE
  ```

  - `file`: Concatenate sequences per file, use filename as ID (default)
  - `sequence`: Query each sequence separately, use header as ID
  - `batch`: Same as file mode

- `--top INT`: Show only top N results per query

**Note:** Query mode automatically reads k, theta, and seed from the database.

**Examples:**

```bash
# Basic query
./repeat_robust_estimator query -d db/ -q assembly.fasta -o results.tsv -p 4

# Query with presence-presence mode
./repeat_robust_estimator query -d db/ -q assembly.fasta -o results.tsv --pp

# Query multiple files
./repeat_robust_estimator query -d db/ -q *.fasta --mode batch -o results.tsv

# Query each sequence separately
./repeat_robust_estimator query -d db/ -q contigs.fasta --mode sequence -o results.tsv
```

## Estimation Methods

This tool implements three k-mer based estimators for mutation rate estimation. Each estimator uses different levels of information from the query sequence (s) and the database sequence (t).

### Estimator Comparison Table

| Estimator      | Knowledge of k-mers in s (query) | Knowledge of k-mers in t (database) | Formula for q                                                |
| -------------- | -------------------------------- | ----------------------------------- | ------------------------------------------------------------ |
| $\hat{q}_{pp}$ | Presence/absence                 | Presence/absence                    | $\vert sp(t) \setminus sp(s) \vert $                         |
| $\hat{q}_{pc}$ | Presence/absence                 | Counts                              | $\sum_{\tau \in sp(t) \setminus sp(s) } occ(\tau,t) /L$      |
| $\hat{q}_{cc}$ | Counts                           | Counts                              | $q_pc + (1-\hat{r}\_{pc})^{k-1} · (\hat{r}\_{pc}/3L) \cdot \sum_{\tau \in sp(s)} occ(\tau,s) · h_1(\tau,s)$ |

**Notation:**
- $s$: Query sequence with $L = \vert s \vert - k + 1$
- $t$: Database sequence (considered as result of mutation process on s with substitution rate r)
- $sp(s)$: Set of k-mers present in sequence $s$ (k-spectrum)
- $occ(\tau, t)$: Number of occurrences of k-mer $\tau$ in sequence $t$
- $h_1(\tau, s)$: Number of k-mers in $sp(s)$ with Hamming distance 1 from $\tau$
- $q$: Shorthand for $1 - (1-r)^k$, where r is the substitution rate
- **Conversion**: Any estimator $\hat{r}$ is defined by $\hat{r} = 1 - (1 - \hat{q})^{1/k}$

## File Organization

When building a sketch database, the tool creates:

```
db/
├── index.txt          # Maps sketch IDs to files
├── sketch_0.bin       # Binary sketch file
├── sketch_1.bin
└── ...
```

Each binary sketch contains:

- Sketch ID
- Parameters (k, theta, seed)
- K-mer counts
- h1 statistics (if computed)

## Citation

If you use this tool in your research, please cite:

[TBD]


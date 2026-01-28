# Repeat-robust mutation rate estimators by $k$-mer sketchs

This tool is a prototype of three **repeat-robust**, **$k$-mer-based** estimators. We consider the following random **substitution process**, parameterized by a rate $0 \le r \le 1$. Given a string $s$, the character at each position mutates to one of the three other nucleotides with probability $r/3$ per nucleotide independently. We denote the mutated string as $t$. The set of all the distinct k-mers of the string $s$ is called a **spectrum** of $s$. When k-mers are associated with their **occurrence counts** in string $s$, we call it a **multiplicity** of $s$. The estimators are applicable in different settings, based on whether they need count information of the sequences $s$ and $t$. Note here the roles of $s$ and $t$ is **not symmetric**, especially for highly repetitive sequences. The three estimator corresponds to three modes in our tool:

- **--pp**: presence-to-presence (uses distinct k-mers of $s$ and $t$)
- **--pc**: presence-to-count (uses distinct k-mers of $s$  and k-mer counts of $t$)
- **--cc**: count-to-count (uses k-mer counts of $s$ and $t$). 

## Installation

```bash
# Compile the tool
make
```

Requirements: C++17 compatible compiler. 

## Usage

### Basic 

```bash
mutation_estimator --mode <pp|pc|cc> -L <length> -k <kmer_size> -n <novel_kmer_count> [options]
```

### Required arguments

- `--mode <mode>`: Estimation mode (pp, pc, or cc)
- `-L <length>`: Total sequence length in base pairs of $t$
- `-k <size>`: K-mer size
- `-n <count of novel kmers>`: Number of novel k-mers (see mode-specific details below)

### Optional arguments

- `-t <theta>`: FracMinHash sampling rate (default: 1.0)
- `-s <file>`: FASTA file with original sequence (required for --cc mode)

## Modes Explained

### Mode 1: --pp (presence-to-presence)

Use distinct k-mer sets of $s$ and $t$. Suitable when you have the k-mer spectra of $s$ and $t$.

**Required input**:

- `-n`: Number of distinct novel k-mers in $t$'s spectrum (k-mers present in $t$ but not in $s$)

**Example**:

```bash
# You have 10,000 distinct novel k-mers
mutation_estimator --mode pp -L 1000000 -k 21 -n 10000
```

### Mode 2: --pc (presence-to-count)

Use distinct k-mer set of $s$ and k-mer counts of $t$. More accurate than --pp when dealing with repetitive sequences.

**Required input**:

- `-n`: Total number of novel k-mer occurrences in t (with multiplicity)

**Example**:

```bash
# You have 10,000 total novel k-mer occurrences
mutation_estimator --mode pc -L 1000000 -k 21 -n 10000
```

### Mode 3: --cc (count-to-count)

Use k-mer counts of $s$ and $t$. Most accurate mode for repetitive sequences. We need $d_1$ value of each k-mer in $s$, where $d_1$ is the number of k-mers with hamming distance $1$ in the spectrum of $s$. 

**Required input**:

- `-n`: Total number of novel k-mer occurrences in t
- `-s`: FASTA file containing the original sequence s

**Example**:

```bash
# You have 10,000 total novel k-mer occurrences and original sequence s
mutation_estimator --mode cc -L 1000000 -k 21 -n 10000 -s original.fasta
```

## Using with FracMinHash

If you count your novel kmer based on FracMinHash, use the `-t` parameter to specify the sampling rate:

```bash
# If you sketched with theta=0.01
# and found 150 novel sketched k-mers
mutation_estimator --mode pp -L 1000000 -k 21 -n 150 -t 0.01
```

### How to get k-mer counts

You can obtain novel k-mer counts from various tools, like [KMC](https://github.com/refresh-bio/KMC), [Jellyfish](https://github.com/gmarcais/Jellyfish), or estimate by FracMinHash by using [sourmash](https://github.com/sourmash-bio/sourmash)

## Citation

If you use this tool, please cite:

```
Haonan Wu, Paul Medvedev. The gift of creation: repeat-robust estimators of substitution rates. 2026 under review. 
```

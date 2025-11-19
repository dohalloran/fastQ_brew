# fastQ_brew v2.1

[![GitHub license](https://img.shields.io/badge/license-GPL_2.0-orange.svg)](https://raw.githubusercontent.com/dohalloran/fastQ_brew/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/dohalloran/fastQ_brew.svg)](https://github.com/dohalloran/fastQ_brew/issues)

**fastQ_brew** is a lightweight and flexible FASTQ pre-processing toolkit, written in Perl and designed to run on macOS and Linux without heavy or non-portable dependencies.

It provides read filtering, trimming, adapter removal, format conversions, reverse-complementing, demultiplexing validation, and now **parallel chunk-based processing** for large FASTQ datasets.

---

## ⭐ Features

- 🚀 **FASTQ preprocessing**: length filtering, quality filtering, trimming  
- ✂️ **Adapter removal** (TruSeq, 1 mismatch allowed)  
- 🧬 **Convert FASTQ → FASTA**  
- 🔄 **Reverse complement reads**  
- 🔁 **DNA→RNA conversion** (T/t → U/u)  
- 🔍 **Demultiplexing validation** (`--plex`)  
- ⚡ **Parallel mode for large FASTQs** (no XS modules required)  
- 🧹 **Core-only Perl modules** (Moose + core utilities)

---

## 📦 1. Installation / Setup

### **Option A — Run directly (recommended)**  
No installation is required. Clone the repo:

```bash
git clone https://github.com/dohalloran/fastQ_brew.git
cd fastQ_brew
````

Make the main scripts executable:

```bash
chmod +x fastQ_brew.pl fastQ_brew_parallel.pl
```

Run:

```bash
./fastQ_brew.pl ...
./fastQ_brew_parallel.pl ...
```

---

### **Option B — Install as a Perl module (optional)**

```bash
perl Makefile.PL
make
make test
make install
```

---

## 📚 2. Requirements

* Perl **5.10+**
* Modules (minimal):

  * `Moose`
  * `Term::ANSIColor`
  * all other modules are Perl core

Install missing modules with:

```bash
cpanm Moose Term::ANSIColor
```

---

## 🚀 3. Quick Start Examples

### **3.1 Serial FASTQ processing**

```bash
perl fastQ_brew.pl \
  --i input.fastq \
  --o filtered.fastq \
  --lib sanger \
  --qf 50 \
  --lf 25 \
  --trim_l 5 \
  --trim_r 3 \
  --truseq
```

#### **Key Options**

| Option                  | Meaning                                 |
| ----------------------- | --------------------------------------- |
| `--i`                   | Input FASTQ (required)                  |
| `--o`                   | Output FASTQ (default: `filtered.fq`)   |
| `--lib`                 | Quality encoding (`sanger`, `illumina`) |
| `--qf`                  | Max percent error probability (1–100)   |
| `--lf`                  | Minimum post-trim read length           |
| `--trim_l` / `--trim_r` | Trim N bases from left / right          |
| `--truseq`              | Remove TruSeq adapters                  |

---

### **3.2 Parallel processing**

Use this for very large FASTQs (multi-GB).
The wrapper:

1. Splits the FASTQ into chunks
2. Runs `fastQ_brew.pl` on each chunk in parallel
3. Joins all output chunks into a final filtered FASTQ

```bash
perl fastQ_brew_parallel.pl \
  --i big.fastq \
  --o big.filtered.fastq \
  --lib sanger \
  --qf 50 \
  --lf 25 \
  --trim_l 5 \
  --trim_r 3 \
  --truseq \
  --jobs 8 \
  --reads_per_chunk 200000
```

#### Parallel Options

| Option              | Meaning                                     |
| ------------------- | ------------------------------------------- |
| `--jobs`            | Number of worker processes                  |
| `--reads_per_chunk` | FASTQ reads per chunk (each read = 4 lines) |

> Parallel mode skips `--plex` because demux validation requires both full FASTQ files.

---

## 🧪 4. Full Command-Line Reference

### **Filtering**

```bash
--qf 50        # Max % error probability (lower = stricter)
--lf 25        # Minimum read length after trimming
--trim_l 5     # Trim 5 bases from left
--trim_r 3     # Trim 3 bases from right
--truseq       # Remove TruSeq adapters (allows 1 mismatch)
```

---

### **Conversion / Utility Options**

```bash
--fasta        # Write FASTA instead of FASTQ
--dna_rna      # Convert T/t → U/u
--rev_comp     # Reverse-complement sequences
```

---

### **Demultiplexing Validation**

```bash
--plex \
  --i sample1.fastq \
  --x sample2.fastq
```

Reports matching and unique barcodes between two files.

**Not supported in parallel mode.**

---

### **General Options**

```bash
--i <file>     # Input FASTQ (required)
--o <file>     # Output FASTQ (default: filtered.fq)
--lib sanger   # Quality encoding (default)
--lib illumina # Alt quality encoding
```

---

## 🧬 5. Using fastQ_brew Programmatically (Perl API)

```perl
use strict;
use warnings;

use lib 'lib';
use fastQ_brew;

my $app = fastQ_brew->new(
    i       => 'input.fastq',
    o       => 'filtered.fastq',
    lib     => 'sanger',
    qf      => 50,
    lf      => 25,
    trim_l  => 5,
    trim_r  => 3,
    truseq  => 1,
);

$app->run_fastQ_brew();
```

---

## 📄 6. References

1. Cock PJ, et al. *The Sanger FASTQ file format…* NAR 2010.
2. Ewing B & Green P. *Phred I & II accuracy models.* Genome Res 1998.

---

## 🤝 7. Contributing

Pull requests and issues are welcome:

👉 [https://github.com/dohalloran/fastQ_brew/issues](https://github.com/dohalloran/fastQ_brew/issues)

---

## 🛟 8. Support

Open an issue on GitHub for help.

---

## 📜 9. License

GNU General Public License v2.0 (GPL-2.0)


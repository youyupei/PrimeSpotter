# PrimSpotter: Detecting Internal Priming Artifacts in RNA-seq Data

PrimSpotter is a tool to identify and characterize internal priming artifacts in RNA-seq data, especially those caused by A-rich genomic regions that lead to spurious priming events.

## Purpose

Reverse transcription in RNA-seq can occasionally initiate at internal A-rich regions of transcripts rather than true poly(A) tails, leading to false-positive transcript ends. This artifact can affect transcript assembly, quantification, and downstream analysis.

This tool identifies such artifacts by scanning the genome for A-rich regions and cross-referencing them with RNA-seq read alignments to flag likely internally primed reads.

## Algorithm Overview

### Step 1: Detect A-Enriched Genomic Regions (Both Strands)

1.1. **Sliding Window Search**  
Scan the genome using a sliding window of size 10 nucleotides.  
Select windows that contain **≥ 8 As** (on the plus strand) or **≥ 8 Ts** (on the minus strand).

1.2. **Merge Overlapping Windows**  
Merge overlapping or adjacent enriched windows into larger contiguous regions.

1.3. **Region Extension**  
Extend each merged region by **25 nucleotides** upstream and downstream to accommodate read alignment variability.

### Step 2: Flag Internally Primed Reads

2.1. **Read Categorization**  
Label a read as an **internally primed read** if:
- It **ends** within an A-rich region, or  
- It **starts** within a T-rich region

## Usage

To be added: script usage examples, input format, and dependencies.

## Visualization

![Workflow Diagram](PrimSpotter.svg)

## License

GPL3 License
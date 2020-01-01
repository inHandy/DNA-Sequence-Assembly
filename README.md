# Project: DNA-Sequence-Assembly

An interesting project in which data structures as well as a De bruijn Graph are used to create and assemble DNA sequences from large files efficiently.

## Getting Started

Some knowledge about the De Bruijn Graph and DNA sequencing is required.

## Overview

A hash map is used to store the k-mers (sequences of k DNA strands), using a hashing function specific to this use case.
By combining common sequences, it becomes more efficient to conduct operations on the data contained in large file, as with the depth first search.

## Future Objectives

- Refactor and add the existing Unit Tests to GitHub.
- Test different hashing functions. 

# Project: DNA-Sequence-Assembly

An interesting project in which data structures as well as a Debruijn Graph are used to create and assemble DNA sequences from large files efficiently.

## Getting Started

Some knownledge about the De Bruijn Graph and DNA sequencing is required.

## Overview

A hash map is used to store the k-mers (sequences of k strands), using a hash function specific to this use-case.
By combining common sequence, it becomes more efficient to conduct operations on the data contained in large file, as with the depth first search.

## Future Objectives

- Refactor and add the existing Unit Tests to GitHub.
- Test different hashing functions. 

# Project: DNA-Sequence-Assembly

An interesting project in which data structures as well as a De bruijn Graph are used to create and assemble DNA 
sequences from large files efficiently.

## Overview

A De Bruijn Graph stored inside a hash map is used to contain the given k-mers (sequences of k DNA strands), indexed 
by a hashing function specific to this use case.

Includes a Depth-First search method that generates contigs.
For longer contigs, it's recommended to choose a smaller k-mer size .

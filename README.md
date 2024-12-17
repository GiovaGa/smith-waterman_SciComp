# Smith-Waterman algorithm implementations for Introduction to Scientific Computing course

This repository contains code of some simple implementations and testing code for the Smith-Waterman algorithm. It was created as a final project for the Introduction to Scientific Computing course, taught in autumn 2024 at the [Galielian School of Higher Education](https://scuolagalileiana.unipd.it/)).  

## Introduction to the problem and the Smith-Waterman algorithm
Genome and protein alignment are fundamental bio-informatics techniques used to compare sequences of DNA, RNA, or proteins to identify regions of similarity.
These alignments reveal evolutionary relationships, functional similarities, or conserved regions across species or within different genes and proteins of the same organism.
Genomes can consist of billions of base pairs (the human genome has over 3 billion base pairs!). Aligning such large datasets requires significant computational resources.

The Smith-Waterman algorithm is a dynamic programming algorithm used specifically for local sequence alignment, that is to identify the most similar subsequences between two sequences.

## This repository

We have experimented mainly with three different approaches using OpenMP:
1. Using locks, a block division using locks to syncronize.
2. Using atomic counters
3. Using tasks, a Divide&Conquer recursive approach.

The relative implementations are in the `src/` folder, along with some testing code.

## Future work

Some nice things could be tried in the future:
A number of different approaches are possible:
1. **MPI** could be well suited to solve this problem: the communication between processes is much smaller than the computation done by each of them.
2. **GPUs** are widely used in the literature for this algorithm, though this is probably very tricky.

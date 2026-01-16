# 2025 Summer Research Project
This is the code produced from the research performed by myself (Alon) and supervising professor Dr. Robert Craigen from the Univserity of Manitoba.

### Background
A Hadamard matrix is a square grid of 1s and -1s in which the dot product of any two rows is 0. The dot product between two rows can be defined as the sum of the products of elements in the same column.

Two Hadamard matrices of the same dimensions are equivalent if there exists a set of row swaps, column swaps, row negations, and column negations such that you can permute one matrix into the other.

In this specific domain, the rank of a matrix can be defined as the number of rows which are distinct (with a row being the negation of another not counting as distinct)

### Problem Statement
Given a Hadamard matrix, partition it into 3x3 submatrices by inserting dividing lines. We wish to find a Hadamard matrix of every class that has the same rank within each of the 3x3 submatrices, if they even exist.

### Nauty:
Nauty is a very helpful tool for testing the automorphism of graphs quickly. There is a formulation that allows you to test the equivalence of Hadamard matrices by testing the automorphism of corresponding graphs, so Nauty was used for this.

### Folders:
`collections` contains lists of matrices of specific types, for easy retrieval of common matrices.

`tasks` is specific code to perform a specific task. Small variations on the code will easily modify the task.

`tools` is the collection of useful methods and algorithms used throughout the project.

### Aside:
The details and descriptions of the algorithms are described in a private overleaf project.
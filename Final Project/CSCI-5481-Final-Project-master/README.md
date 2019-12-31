# CSCI-5481-Final-Project

Project Option 2.

Download a chromosome from the human genome (http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/). Simulate whole-genome shotgun data by breaking it into random overlapping substrings of length k. Break it down further into a list of k-mers for some small k. Create the De Bruijn graph using k-mers as nodes, where a Hamiltonian circuit would be needed to find a valid path through the genome. Then convert this to a k-1-dimensional De Bruijn graph in which an Eulerian circuit would traverse all k-mers. Implement an algorithm to find an Eulerian circuit, and report the assembled chromosome. Compare it to the original chromosome by plotting the coordinates of the assembled location of each base against the known location of each base.

Try this first with a short section of 25-50 base pairs with several different values of k. Generate plots of the k-dimensional and k-1-dimensional graphs. Then try it on larger and larger portions of the chromosome. Report the performance (fraction correctly assembled and runtime). Then try simulating sequencing error in the data (1 random error per 100 bases) and report the effects.

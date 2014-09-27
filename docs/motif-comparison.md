# Motif comparison

CollecTF allows users to perform pair-wise comparisons of arbitrary searches on
the database. Query results can be compared using site- and motif-based
statistics relying on dynamic alignment of query results and contextualized
using permutation tests.

The first step in motif comparison is to perform two searches, results of which
will be compared. Then, two motifs are compared using motif- and site-based
similarity measures.

## Motif-based similarity

### Motif alignment

The alignment for two motifs are computed using ungapped [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) local
alignment algorithm. To score column similarity, [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
is used.

### Motif similarity

The statistics for similarity of two motifs are computed using several
similarity functions, defined for a pair of columns. To score motif similarity,
the alignment that gives to best score is used.

-   *Pearson correlation coefficient* gives a measure of correlation between base
    frequencies of two motif columns.
-   *Average log-likelihood ratio* between two motif columns is the sum of
    log-likelihood ratios.
-   *Kullback-Leibler divergence* (KLD) is a non-symmetric measure of the difference
    between two probability distributions. Symmetric form of KLD is used to
    measure the difference between two motif columns.
-   *Euclidean distance* is the distance between two points in n-dimensional
    space.

For details on adopting these functions for motif similarity, see [Mahony and
Benos, 2007](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933206/).

For each similarity function, permutation tests are applied to assess the
significance of the motif similarity. For each test, 100 replicates of each
motif are obtained by permuting columns of motifs. The histogram of 100 scores
are presented with the true similarity score of the two motifs.

### Unaligned/aligned sites

Motif comparison page also contains two binding site lists, the search results,
side by side for comparison.

## Site-based similarity

In addition to motif-based similarity/distance measures, two motifs can be
compared by measuring the similarity/distance between all pairs of sites from
the two motifs. CollecTF uses [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) (also known as "edit
distance") to measure the distance between two sites.

## Compared motifs

The tab "compared motifs" on the result page contains the list of TFs and
species of all binding motifs for the comparison. Individual motifs and detailed
regulation reports are available through given links.

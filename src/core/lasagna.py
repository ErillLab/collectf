"""Wrapper for LASAGNA binding site alignment."""

import lasagna_lib
from bioutils import reverse_complement


def lasagna(site_instances):
    """Aligns SiteInstance sequences using LASAGNA and returns alignment.

    Given a list of models.SiteInstance objects, run LASAGNA algorithm and
    return the aligned sequences. Here is the link to the paper:
    http://www.biomedcentral.com/1471-2105/14/108
    """
    aligned, idxAligned, strands = lasagna_lib.LASAGNA(
        [site_instance.sequence.lower() for site_instance in site_instances], 0)
    aligned = [seq.upper() for seq in aligned]

    assert (map(int, idxAligned) == range(len(idxAligned)),
            "Lasagna sites are not sorted")

    return [fill_gaps(site_instance, aligned_site, aligned_strand)
            for (site_instance, aligned_site, aligned_strand) in
            zip(site_instances, aligned, strands)]


def fill_gaps(site_instance, aligned_site, aligned_strand):
    """Fills the gaps in the LASAGNA alignment."""
    genome_seq = site_instance.get_genome_sequence()
    original_seq = (site_instance.seq if aligned_strand == '+'
                    else reverse_complement(site_instance.seq))
    num_left_gaps = 0
    num_right_gaps = 0
    # Count left and right gaps
    while aligned_site[num_left_gaps] == '-':
        num_left_gaps += 1
    while aligned_site[-num_right_gaps-1] == '-':
        num_right_gaps += 1
    # Check if correct
    assert ('-' * num_left_gaps + str(original_seq) +
            '-'*num_right_gaps == aligned_site)
    # Extend site to both sides.
    max_gaps = max(num_left_gaps, num_right_gaps)
    extended_site = extend_site(site_instance, genome_seq, max_gaps)
    if aligned_strand == '-':
        extended_site = reverse_complement(extended_site)
    if num_right_gaps == max_gaps:
        recovered_site = extended_site[max_gaps-num_left_gaps:]
    else:
        recovered_site = extended_site[max_gaps-num_left_gaps:
                                       num_right_gaps-max_gaps]
    # Sanity check
    assert (len(recovered_site) == len(aligned_site))
    assert recovered_site in extended_site
    return str(recovered_site)


def extend_site(site_instance, genome_seq, n=250):
    """Extends the SiteInstance sequence by n bases both sides."""
    seq = genome_seq[site_instance.start-n: site_instance.end+n+1]
    return seq if site_instance.strand == 1 else reverse_complement(seq)

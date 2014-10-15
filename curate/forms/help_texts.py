
# -*- coding: utf-8 -*-

"""This file contains all help texts that appear somewhere during paper
submission and curation processes. The purpose of keeping long help texts in a
separate file is to keep the code clean and make it easy to change those texts
later, if required."""

pubmed_publication_form = dict(
    pmid="Paste the PubMed ID obtained from the NCBI website.",
    reported_TF="""Type the name of the transcription factor(s) reported in the
    manuscript.""",
    reported_species="Type the name of the species reported in the manuscript.",
    contains_promoter_data="""The paper provides experimental data on the
structure and sequence of TF-regulated promoter.""",
    contains_expression_data="""The paper provides experimental support for
TF-mediated regulation of genes.""",
    submission_notes="""Include any additional
details about the submission. For instance, you might indicate the approximate
number of sites reported, whether high-throughput techniques are used or any
other factor that might help prioritize curation.""",
)

non_pubmed_publication_form = dict(
    reported_TF="""Type the name of the transcription factor(s) reported in the
    manuscript.""",
    reported_species="Type the name of the species reported in the manuscript.",
    contains_promoter_data="""The paper provides experimental data on the
    structure and sequence of TF-regulated promoter.""",
    contains_expression_data="""The paper provides experimental support for
    TF-mediated regulation of genes.""",
    submission_notes="""Include any additional details about the
submission. For instance, you might indicate the approximate number of sites
reported, whether high-throughput techniques are used or any other factor that
might help prioritize curation.""",
)

# Curation help texts
publication_form = dict(
    pub='',

    no_data= """Check this button if, after examining the paper, you find that
    the paper does not have data on binding sites. Checking this button will
    mark the paper as having no binding site data and set it to the 'curation
    complete' status. Also, the curation process will be ended as the paper has
    no data to be curated.""",
)

genome_form = dict( 

    TF="""Select the transcription factor you are curating on from list. If not
    in list, please add the corresponding TF/TF-family using Data submission
    menu.""",

    TF_type="""If specified in the manuscript, select the quaternary structure
    for the transcription factor when binding to the sites reported in this
    curation.""",

    TF_function="""If specified in the manuscript, select the mode of operation
    for the TF on the sites reported in this curation.""",

    genome_accession="""Paste the NCBI GenBank genome accession number for the
    species closest to the reported species/strain.
    (e.g. <code>NC_000913.2</code>) You can add more than one chromosome. """,

    TF_species_same="""
    Check if the reported strain and selected RefSeq strain are same.""",

    site_species_same="""
    Check if the reported strain and selected RefSeq strain are same for the TF.
    """,

    TF_accession="""Paste the NCBI TF protein accession number for the species
    closest to the reported species/strain.  (e.g. <code>NP_799324</code>) You
    can add more than one TF.""",

    TF_species="""If the work you are reporting uses a strain different from the
selected RefSeq genome, please type/paste the original strain
(e.g. <code>Pseudomonas sp. ADP</code>). This allows us to
keep track of the correspondence between reported and mapped strains.
""",

    site_species="""If the work you are reporting uses a strain different from
the selected RefSeq genome, please type/paste the original strain 
(e.g. <code>Pseudomonas putida plasmid pEST1226</code>). This allows
us to keep track of the correspondence between reported and mapped strains.
""",

    contains_promoter_data="""Check if the paper provides experimental data on
the structure and sequence of a TF-regulated promoter""",

    contains_expression_data="""Check if the paper provides experimental support
    for TF-mediated regulation of genes.  Please make sure that this field is
    checked if you plan to report differential gene expression associated with
    TF activity.""")

techniques_form = dict(
    techniques="""Select as many as apply to sites reported in this submission. Hover over any technique to see the
    description.""",

    experimental_process="""Write a concise, intuitive description of the experimental process to ascertain
    binding/induced expression.

<a class="" data-toggle="modal" data-target="#experiment_modal" href="#">
    [examples]
</a>

<!-- Modal -->
<div class="modal fade" id="experiment_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">Experimantal description examples</h4>
      </div>
      <div class="modal-body">
        <h4>Experimental description example 1</h4>
        <p>
          LuxT was identified as binding the smcR promoter through
          SDS-PAGE. EMSAs confirmed that it binds specifically and DNAse
          footprint identified a protected region with sequence
          AGTGCAATACGCTATTTACTATCACA. A luxT- mutant was shown to induce smcR
          expression through Western blot and luciferase report
        </p>
        <h4>Experimental description example 2</h4>
        <p>
          ChIP-chip analysis was performed on a csgD deletion strains
          transformed with an arabinose-inducible csgD-expressing plasmid
          (pBADcsgD) or control plasmid (pBAD18). DNAse footprinting was
          performed on six peaks (csgD-csgB, fliE-fliF, wrbA-ymdF, nlpA-yicS,
          yccU-yccT, and yhbT-yhbU) and identified several protected regions. A
          motif resembling the previously reported CsgD motif in S. enterica was
          derived by multiple sequence alignment of footprinted regions. Using
          LacZ fusions, promoters were tested for repression (fliE, yhbT
          promoter) or activation (yccT adrA). Individual binding to sites in
          the csgB promoter was confirmed by EMSA and site-directed mutagenesis.
        </p>
      </div>
    </div>
  </div>
</div>
""",

    external_db_type="""Select type of external database containing data
    (e.g. DNA-array data) reported in paper""",

    external_db_accession="""Type the accession number for external database
    referenced in paper.""",

    forms_complex="""Check if the manuscript reports characterization of the
    interaction of the TF with another protein""",

    complex_notes="""Provide brief description of the proteins involved in the
    complex and how it affects binding"""
)

site_entry_form = dict(
    is_motif_associated="""Check this option if sites are reported by the
    authors as being associated with a known TF-binding motif. Uncheck if sites
    are reported solely as DNA fragments shown to be bound by the TF, without
    reporting any specific binding sequences.""",

    site_type="""


<a class="" data-toggle="modal" data-target="#motif_type_modal" href="#">
    [motif examples]
</a>

<!-- Modal -->
<div class="modal fade" id="motif_type_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">Motif type examples</h4>
      </div>
      <div class="modal-body">
        <img src="/static/Motif_associated_example.jpg" width="70%"
             class="center-block">
        <hr/>
        <img src="/static/Motif_associated_example.jpg" width="70%"
        class="center-block">
        <hr/>
        <img src="/static/Variable_motif_associated_example.jpg" width="70%"
        class="center-block">
        <hr/>
        <img src="/static/Variable_motif_associated_example2.jpg" width="70%"
        class="center-block">
        
      </div>
    </div>
  </div>
</div>
""",

    sites="""
 Enter the list of sites in FASTA format, raw sequence or coordinate format (one
site per line). Sequence entries must be in unambiguous DNA code (A, C, G or T;
no degenerate IUPAC codes (e.g. W, Y, N...) or gap symbols (-, _, *). FASTA
format does not support quantitative data entry. Quantitative data (q-val) can
be added to raw sequence or coordinate entries. All fields (i.e. site & q-val or
coordinates & q-val) must be either space or tab separated.
    
<a class="" data-toggle="modal" data-target="#site_modal" href="#">
    [examples]
</a>    

<!-- Modal -->
<div class="modal fade" id="site_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">Site Entry Examples</h4>
      </div>
      <div class="modal-body">
    <h4>FASTA example</h4>
    <pre>
>Binding site 1
AGAACATTTGTTCC
>Binding site 2
AGAACTCATGTTCG
>Binding site 3
AGAACATTCGTTCT
>Binding site 4
AGAACGTATGTTTT
>Binding site 5
AGAACGTACATTCC
>Binding site 6
AGAACGTGCATTCG
    </pre>
    <h4>RAW example (sequence)</h4>
    <pre>
AGAACATTTGTTCC
AGAACTCATGTTCG
AGAACATTCGTTCT
AGAACGTATGTTTT
AGAACGTACATTCC
AGAACGTGCATTCG
    </pre>
    <h4>Raw sequence example (sites + quantitative information)</h4>
    <pre>
AGAACATTTGTTCC	3.2
AGAACTCATGTTCG	1.2
AGAACATTCGTTCT	2.2
AGAACGTATGTTTT	2.3
AGAACGTACATTCC	0.8
AGAACGTGCATTCG	1.9
    </pre>
    <h4>Raw sequence example (sites + quantitative information)</h4>
    <pre>
AGAACATTTGTTCC	3.2
AGAACTCATGTTCG	1.2
AGAACATTCGTTCT	2.2
AGAACGTATGTTTT	2.3
AGAACGTACATTCC	0.8
AGAACGTGCATTCG	1.9
    </pre>
    <h4>Coordinates example (sites + quantitative information)</h4>
    <pre>
155970	155983	3.2
608806	608819	1.2
655172	655159	2.2
719335	719348	2.3
1056255	1056268	0.8
1064783	1064796	1.9
    </pre>   
      </div>
    </div>
  </div>
</div>
""",
    
    quantitative_data_format="""If the manuscript reports quantitative values
    associated with sites, please enter the quantitative data format here. If
    not, you can leave this field empty.

<a class="" data-toggle="modal" data-target="#q_modal" href="#">[example]</a>

<!-- Modal -->
<div class="modal fade" id="q_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">Quantitative field example</h4>
      </div>
      <div class="modal-body">
        <p>
          ChIP-Seq fold enrichment: change in detection of a region comparing
          the sequences complexed to the induced AmrZ to the input DNA from that
          sample. Range: 3.00 to 35.53
        </p>
      </div>
    </div>
  </div>
</div>
""",

    peaks="""Enter the peak data (in either coordinate or sequence mode). If
    there is any quantitative data associated with the peak data, they will be
    automatically mapped to entered sites. Mapped peak intensity values will be
    displayed for review before curation submission.""",

    assay_conditions="""Describe the conditions of the high-throughput
    experiment that capture the specifics of the in-vivo setting for
    cross-linking. Were cells at exponetial-phase? Was the system induced? How
    were cells grown?

<a class="" data-toggle="modal" data-target="#assay_conditions_modal" href="#">
    [examples]
</a>

<!-- Modal -->
<div class="modal fade" id="assay_conditions_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close"
        data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">ChIP assay conditions examples</h4>
      </div>
      <div class="modal-body">
        <h4>ChIP assay conditions example 1</h4>
        <p>Strains were grown at 37uC on LANS (LBNS with 1.5% agar) or
        Pseudomonas Isolation Agar (Difco, Detroit, MI) agar plates. Cultures
        were induced with 0.5% arabinose at an OD600 of 0.1 and allowed to grow
        for two hours at 37uC in a roller.</p>
        <h4>ChIP assay conditions example 2</h4>
        <p>Wild-type Anabaena sp. PCC 7120 cells growing in bubbled cultures
        with ammonium as the N source were subjected to incubation in a combined
        N-depleted medium for 3 hours, after which the cultures were treated
        with formaldehyde to fix the proteins bound to DNA. After cell lysis and
        DNA fragmentation, the extracts were treated with an anti-NtcA antibody
        to specifically immunoprecipitate the NtcA-bound DNA. The
        immunoprecipitated material was then incubated at 65ºC to reverse the
        crosslinking, and the DNA was isolated. A sample of total DNA was also
        isolated prior to anti-NtcA treatment of the extracts to serve as the
        control input sample. Cells of Anabaena sp. (also known as Nostoc sp.)
        strain PCC 7120 growing exponentially (3-5 μg Chl/ml) in the light (75
        μE·m-2·s-1) at 30°C in BG110 medium supplemented with 10 mM NaHCO3
        (referred to as BG110C) containing 6 mM NH4Cl and 12 mM TES and bubbled
        with a mixture of air + 1% CO2 were collected, washed with BG110C,
        resuspended in BG110C, and incubated in the same conditions for 3
        h. Formaldehyde was then added to the cultures to a final concentration
        of 1%, and the cultures were incubated for 15 min (no aeration,
        occasional shaking). Glycine was added at 125 mM final concentration and
        the incubation was continued for 5 min to stop the fixing reaction. The
        cells were then filtered, washed with cold TBS (20 mM Tris–HCl, pH 7.4,
        140 mM NaCl) and collected in tubes (25 ml of culture per tube). The
        pellets were frozen in liquid nitrogen and stored at -20°C until used.
        </p>
      </div>
    </div>
  </div>
</div>
""",

    method_notes="""Describe (use copy-paste if appropriate) the
    high-throughput protocol. What antibodies were used? What chip/sequencer and
    using what parameters? Etc.

<a class="" data-toggle="modal" data-target="#chip_notes_modal" href="#">[example]</a>

<!-- Modal -->
<div class="modal fade" id="chip_notes_modal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close"
        data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="myModalLabel">ChIP notes example</h4>
      </div>
      <div class="modal-body">
        <p>
Protein-DNA complexes were cross-linked by addition of formaldehyde- to a final
concentration of 1.0% and incubated at room temperature for ten
minutes. Cross-linking was quenched by addition of glycine (final concentration
250 mM). The final OD600 was recorded and cells were collected from 1 OD600 of
culture via centrifugation and washed once in LBNS. The supernatant was removed
and pellets were stored for further processing at -80uC. Cell pellets were
resuspended in 1.0 mL of lysis buffer (20 mM HEPES, pH 7.9; 50 mM KCl; 0.5 mM
DTT; 500 mM NaCl; 10 mM imidazole; 1% BSA; 1 mg/mL leupeptin/pepstatin; and 400
mM PMSF) per 1 OD600 of culture. Samples were sonicated on Covaris with the
following conditions: Duty Cycle 20%, Intensity 8, Cycles per burst 200, with
frequency sweeping 20 min total shearing time (60 sec cycles, 20 cycles). Lysate
was cleared via centrifugation (20,0006 g, 30 minutes, 4uC) and the supernatant
was transferred to a fresh tube as the input sample. Magne-HIS beads (Promega
V8560) were blocked at room temperature in lysis buffer for 30 minutes, and then
500 ml of the input sample was added to the beads. After 30 minutes of binding
at room temperature with agitation, the supernatant was removed from the beads
via magnetic separation. Beads were washed five times in wash buffer (100 mM
HEPES, pH 7.5, 10 mM imidazole, 500 mM NaCl, and 1% BSA). Elution buffer (100 mM
HEPES, pH 7.5; and 500 mM imidazole) was added to the beads and incubated at
room temperature for 30 minutes. Supernatant was collected after magnetic
separation and combined with SDS (1.25% final concentration), then heated to
70uC for 30 minutes to reverse cross-links. DNA was purified via
phenol:chloroform extraction and ethanol precipitation. The chip DNA was
quantified with Qubit 2 flurometer (Life Technologies) using Qubit dsDNA BR
Assay. 10 ng of DNA was used to construct each Chip sequencing library,
following NEXTflex ChIP-Seq kit (Bioo Scientific) instruction. NEXTflex ChIP-Seq
Barcodes (Bioo Scientific) were used to index the library. The final DNA
libraries were validated with Agilent 2100 Bioanalyzer using Agilent High
Sensitivity DNA Kit. And the library concentrations were determined by Q-PCR
using KAPA SYBR Fast qPCR kit. The libraries were then run on Single End
flowcell on HiSeq200. HiSeq2000 sequencing was performed, resulting in
approxmately 255 million total single-end 52 bp reads from the six control and
eight treatment samples. Reads were aligned using bwa (0.5.10) to the
Pseudomonas aeruginosa PAO1 reference genome. Approximately 220 million reads
aligned uniquely to the reference (86.3%). A TDF file was created for each
sample for visualization in IGV, which was scaled to reads per 10 million data
using bedtools (2.17.0) and igvtools (2.3.3). ChIP-Seq analysis was performed
using HOMER (4.2). First, aligned data was transformed into a
platform-independent data structure for further HOMER analyses using the
makeTagDirectory function. Secondly, HOMER’s findPeaks-style factor was utilized
to identify peaks, or regions of the genome where more reads are present than
random. Lastly, HOMER’s findMotifsGenome.pl was used to analyze genomic
positions for de novo enriched motif regions of length 50 or 200 and identified
peaks were annotated with the motifs using the annotatePeaks.pl function.
        </p>
      </div>
    </div>
  </div>
</div>
""",

    peak_techniques="""Select all techniques that have been used to identify
    high-throughput data. Note that selected techniques are for peaks only. You will be able to
    select used experimental techniques for each binding site, individually.""",
)

site_exact_match_form = dict()

site_soft_match_form = dict()

site_annotation_form = dict()

gene_regulation_form = dict()

curation_review_form = dict(
    revision_reasons="""Select, if needed, the reason why this curation may
    require revision.  See detailed list of reasons in the curation guide.""",

    confidence="""Check this if experimental techniques and results meet the
    standards specified in the curation guide""",

    NCBI_submission_ready="""A curation is ready for submission if: (a) the
    identified genome sequence matches the reported one or (b) identified and
    reported genomes match at the species level and at least 90% of reported
    sites are located as exact matches.""",

    paper_complete="""Check this box if there are no more curations pending for
    this paper (additional sites, sites supported by different techniques, sites
    for other TFs, etc.""",

    notes="""Type in any additional notes on the curation process. For instance,
    if reported sites were left out for some reason, what prompted selection of
    a surrogate genome instead of another, general comments on the experimental
    process, etc.
    """,

    confirm="Check to submit when you click \"next step\"",
)

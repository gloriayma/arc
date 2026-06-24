1. DATA 
(see /Users/gloria/dev/research/arc/sc_nona/own_writing/4_big_picture_diploid.md)

first test:
single cell line
~50 scRNA tracks from that line
aligned onto the genotype (haploid version)



steps:
    1. download metadata
    - cell metadata from zenodo: cell_id map to donor
    - find out which donors have open genomes (how?)
    - download a single donor's VCF (what's VCF stand for) 
        VCF = Variant Call Format — text format listing variants (chrom, pos, REF, ALT, genotype). The phased version uses | (e.g. 0|1) instead of / (ie, it knows "haplotype 1 has these alleles" rather than "both these alleles are present, but haplotype 1/2 is unresolved").
    - download CRAM manifest so you get the links to cell CRAMs to download


    2. download reference files
    - reference genome hs37d5 (because CRAM was encoded against this reference genome)
    - GENCODE v19 gene annotations for intron/exon splice sites, so we can reference them when aligning our scRNA to genome

    3. build the donor's reference FASTA (what does FASTA stand for)
        FASTA = "FAST-All" — plain text sequence format with >header lines + sequence lines. The de facto format for reference genomes.
    
    - build using bcftools

    i want to build the haploid version. by haploid, i mean: 
        (c) Consensus / pseudo-haploid: homozygous ALT applied, het sites left as REF. Coordinates match hs37d5 (skip indels at hets).
    the reason for wanting things to align to hs37d5 is just for easy sanity checking at this stage. 

    3.5 build the genome STAR, for later alignment of scRNA onto it 
    4. download crams, for a single donor, all cells from that donor (or maybe, all cells from that donor of a certain day/maturity).
        - create some figure / table that tells me: how many donors (genotypes), how many cells per donor, and how many cells per day per donor. i envision like, a figure where each row is one donor; the x-axis is number of cells, adn each row is partitioned into 3 (or however many "age" categories there are)
    
    5. CRAM → FASTQ (per-cell thing of the scRNA) through unpacking using samtools fastq → STAR (does this do alignment on your generated genome) (pipe straight through with samtools) → bedGraph -> BigWig 
        STAR does alignment? → Yes, splice-aware short-read aligner. Maps reads (including ones that span introns) back to a genome.
        bedGraph is a plain text format: chrom, start, end, value, one interval per line. BigWig is the same information indexed and binary-compressed for fast random access by genomic region. For human-genome-scale per-base data, a bedGraph might be 1–5 GB uncompressed; the equivalent BigWig is 5–50 MB. Models that need to slice into windows quickly want BigWig.
        
        STAR actually emits four bedGraph files when run with --outWigType bedGraph:
        - Signal.Unique.str1.out.bg — uniquely-mapped reads, strand 1
        - Signal.Unique.str2.out.bg — uniquely-mapped reads, strand 2 (empty because smart-seq2 is unstranded, unlike smart-seq3)
        - Signal.UniqueMultiple.str1.out.bg — uniquely + multi-mapped, strand 1
        - Signal.UniqueMultiple.str2.out.bg — uniquely + multi-mapped, strand 2

        For Smart-seq2 (unstranded — see point 7) the str1 files contain all the coverage you want. Use UniqueMultiple if you want multi-mappers included with their EM-apportioned weights (recommended — single-cell data is sparse enough that you usually want to keep multi-mappers).


    9. Sanity check after STAR runs

        Before committing to the (L × 51) matrix and trying to overfit a model, you want to know the BigWigs are sensible. The cheapest check is to look at one of them at a known housekeeping gene where you know you should see strong, exon-structured coverage:

        - GAPDH at chr12:6,643,584–6,647,536 (GRCh37 coordinates, matches your donor consensus under option c).
        - ACTB at chr7:5,566,782–5,570,304.
        - HPRT1 at chrX:133,594,175–133,634,698.

        You'd expect to see:
        - High coverage over exonic regions (typically hundreds to thousands at peaks).
        - Near-zero coverage over introns.
        - Sharp drops at exon boundaries (this is the splicing signal point 5 in the original chat answer talked about).
        - Coverage over UTRs that may be lower than coding exons but still nonzero (Smart-seq2 is full-length so 5' and 3' UTRs are covered, unlike 10x).

        Quick check:
        import pyBigWig
        bw = pyBigWig.open("cell.bw")
        print(bw.stats("12", 6643584, 6647536, type="mean"))
        print(bw.values("12", 6643584, 6647536))

        If you see flat zeros, something is broken (wrong reference for FASTQ extraction, STAR index built against wrong genome, etc.). If you see ~uniform coverage across introns and exons, something is broken (STAR fell back to ungapped alignment somehow). The exon-structured profile is what tells you the splice-aware pipeline is working.

        You can also pull the matching GENCODE v19 record for the gene and overlay exon intervals on the coverage plot — if exon boundaries don't line up with coverage drops, that's the smoking gun.

----

2. ARCHITECTURE 
(/Users/gloria/dev/research/arc/sc_nona/own_writing/7_big_picture_arch.md)

3. OVERFIT RUN

For the overfit run, pick a single ~1 Mb window with a well-expressed gene cluster in your chosen donor. Good candidates:
- The HLA region (chr6:28-33Mb) — extremely gene-dense, lots of expression, but also highly variable and a notorious alignment hazard, so maybe not first.
- The β-globin locus (chr11:5.2-5.3Mb) — classic enhancer-promoter pairing example, but β-globin specifically may not be expressed in endoderm. Pick something endoderm-relevant.
- For day 3 endoderm: look at SOX17, FOXA2, GATA6, EOMES — these should all be strongly induced. Pick a 1 Mb window centered on one of them, ideally one with several neighboring expressed genes so the axial-attention "covariation between RNA expression of different genes" claim has something to work with.

Practical scaffolding: in your dataloader, pre-window the BigWigs once (extract the 1 Mb region from each of the 50 cells + the donor consensus FASTA, save as a .npz with arrays of shape (L,) for DNA tokens and (L, 50) for coverage). That avoids re-querying BigWig files every training step.
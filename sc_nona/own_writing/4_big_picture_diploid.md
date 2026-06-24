STAR maps the fastq file type? to align with genome

if diploid
- Using STAR
    - Concatenate the two haplotypes end-to-end. 
        - Some reads will map to two places on the genome. Allow that
        - This is messy because you don't teach it like those two places on the genome are actually the same place

the ideal version in my head is like, 
- You definitely don't concatenate the two genomes end-to-end, that is dumb
- Something more like, you build an alignment of the two genomes (allow the insertions/deletions as gap tokens)
- And then the RNA counts are 3 tracks. One track is if it can't be resolved between the two haplotypes. But if it's belonging to a certain haplotype then it goes into the respective track (2 or 3).


^^^ bookmarked for now. 

At minimum, I probably don't need to do diploid alignment. I can just do naive alignment to a single haplotype (would I still recover the haplotype in the same way as if i was recovering the two haplotypes?)
and whatever doesn't align just gets... thrown away or something. 

And then you just have your base DNA track (single haplotype), and many (~50) bigwig-style tracks above, one per cell (with a count per base.)


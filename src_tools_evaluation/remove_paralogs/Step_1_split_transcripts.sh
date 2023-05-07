grep -E "^NC_000001.11|^NC_000009.12" ../../Dataset/refseq_all_transcripts.gff > ../../Dataset/refseq_chr1_chr9_transcripts.gff
grep -Ev "^NC_000001.11|^NC_000009.12" ../../Dataset/refseq_all_transcripts.gff > ../../Dataset/refseq_chr1_chr9_rev_transcripts.gff

library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# # Mane only
# dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_MANE.shuffle.fa")
# dir.create("MANE")
# # Donor
# ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
# ggsave(
#   "MANE/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
# ggsave(
#   "MANE/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )


# # Refseq alternative only
# dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_ALTS.shuffle.fa")
# dir.create("ALTS")
# # Donor
# ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
# ggsave(
#   "ALTS/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
# ggsave(
#   "ALTS/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )


# # Refseq alternative only
# dna <- readDNAStringSet("../../src/INPUTS/800bp/input_neg_1.shuffle.fa")
# dir.create("neg_1")
# # Donor
# ggseqlogo(as.character(subseq(dna[1:50000], 580, 620)), method = 'bits' )
# ggsave(
#   "neg_1/acceptor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )
# # Acceptor
# ggseqlogo(as.character(subseq(dna[1:50000], 180, 220)), method = 'bits' )
# ggsave(
#   "neg_1/donor.png",
#   width = 15,
#   height = 3,
#   dpi = 300
# )


# Refseq alternative only
dna <- readDNAStringSet("../../../src/INPUTS/800bp/input_neg_random.shuffle.fa")
dir.create("neg_random")
# Donor
print("as.character(subseq(dna[1:50000], 580, 620)): ", as.character(subseq(dna[1:50000], 580, 620)))
ggseqlogo(as.character(subseq(dna[1:50000], 580, 620)), method = 'probability', seq_type='dna', col_scheme='nucleotide')
ggseqlogo(as.character(subseq(dna[1:50000], 580, 620)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
ggsave(
  "neg_random/acceptor.png",
  width = 15,
  height = 3,
  dpi = 300
)
# Acceptor
ggseqlogo(as.character(subseq(dna[1:50000], 180, 220)), method = 'probability', seq_type='dna', col_scheme='nucleotide')
ggseqlogo(as.character(subseq(dna[1:50000], 180, 220)), method = 'bits', seq_type='dna', col_scheme='nucleotide')
ggsave(
  "neg_random/donor.png",
  width = 15,
  height = 3,
  dpi = 300
)

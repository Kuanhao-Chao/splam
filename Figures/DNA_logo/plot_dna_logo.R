library(ggplot2)
library(ggseqlogo)

library(Biostrings)

# Mane only
dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_MANE.shuffle.fa")
dir.create("MANE")
# Donor
ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
ggsave(
  "MANE/acceptor.png",
  width = 25,
  height = 6,
  dpi = 300
)
# Acceptor
ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
ggsave(
  "MANE/donor.png",
  width = 25,
  height = 6,
  dpi = 300
)


# Refseq alternative only
dna <- readDNAStringSet("../../src/INPUTS/800bp/input_pos_ALTS.shuffle.fa")
dir.create("ALTS")
# Donor
ggseqlogo(as.character(subseq(dna, 580, 620)), method = 'bits' )
ggsave(
  "ALTS/acceptor.png",
  width = 25,
  height = 6,
  dpi = 300
)
# Acceptor
ggseqlogo(as.character(subseq(dna, 180, 220)), method = 'bits' )
ggsave(
  "ALTS/donor.png",
  width = 25,
  height = 6,
  dpi = 300
)
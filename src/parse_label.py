import pprint
from BCBio import GFF
# from BCBio.GFF import GFFExaminer

in_file = "../Dataset/GCF_000001405.40_GRCh38.p14_genomic.primary.gff"
# examiner = GFFExaminer()
# in_handle = open(in_file)
# pprint.pprint(examiner.parent_child_map(in_handle))
# in_handle.close()

print("Hello world!")
in_handle = open(in_file)
for rec in GFF.parse(in_handle):
    print(rec)
in_handle.close()
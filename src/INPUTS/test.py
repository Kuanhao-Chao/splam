
import random
import sys
from pyfaidx import Fasta

n = 10

pos_shuf = "./input_pos_shuf.fa"
fa = Fasta("./input_pos.fa")

# fw = open(pos_shuf, "w")

seq_num = len(fa.keys()) 
idices = [*range(seq_num)]

random.shuffle(idices)
# print("idices: ", len(idices))

fa = fa[idices]

# for idx in idices:
#     fw.write(">" + fa[idx].name + "\n")
#     fw.write(str(fa[idx]) + "\n")
# fw.close()
# print("idices: ", (ls_idx))

# print(fa[278936])
# counter = 0
# for record in fa:
#     counter += 1
#     # print(record.name)
#     # print(record)
# print("Counter: ", counter)
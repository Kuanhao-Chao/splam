import pandas as pd
import matplotlib.pyplot as plt

d_a_df = pd.read_csv("d_a.bed", sep="\t", header=None)
print(d_a_df)
seq_len_ls = d_a_df.iloc[:, 2] - d_a_df.iloc[:, 1]
seq_len_sum = seq_len_ls.sum()

splam_seq_len_intron = len(seq_len_ls) * 400
splam_seq_len_flanking = len(seq_len_ls) * 400
splam_seq_len = splam_seq_len_intron + splam_seq_len_flanking

spliceai_10k_noN_intron = seq_len_sum
spliceai_10k_noN_flanking = len(seq_len_ls) * 5400
spliceai_10k_noN = spliceai_10k_noN_intron + spliceai_10k_noN_flanking

spliceai_10k_N_intron = seq_len_sum
spliceai_10k_N_flanking = len(seq_len_ls) * 400
spliceai_10k_N = spliceai_10k_N_intron + spliceai_10k_N_flanking

intron_len = [splam_seq_len_intron, spliceai_10k_noN_intron, spliceai_10k_N_intron]
flanking_len = [splam_seq_len_flanking, spliceai_10k_noN_flanking, spliceai_10k_N_flanking]

print("splam_seq_len       : ", splam_seq_len)
print("spliceai_noN_seq_len: ", spliceai_10k_noN)
print("spliceai_N_seq_len  : ", spliceai_10k_N)
print("\n\n")
print("spliceai_noN_seq_len / splam_seq_len : ", spliceai_10k_noN / splam_seq_len)
print("spliceai_N_seq_len / splam_seq_len   : ", spliceai_10k_N / splam_seq_len)
ratio = [1, spliceai_10k_noN/splam_seq_len, spliceai_10k_N/splam_seq_len]
height = [splam_seq_len, spliceai_10k_noN, spliceai_10k_N]
tools = ["SPLAM", "SpliceAI-10k", "SpliceAI-10k-Ns"]
fig, ax = plt.subplots()
ax.bar(tools, intron_len, label="intronic_sequence", linewidth = 1, edgecolor = "black", width=0.7)
ax.bar(tools, flanking_len, bottom=intron_len, label="flanking_sequence", linewidth = 1, edgecolor = "black", width=0.7)
ax.legend()
ax.set_ylabel('Number of nucleotide')

totals = ratio
# Set an offset that is used to bump the label up a bit above the bar.
y_offset = 20000000

for i, total in enumerate(totals):
    print(i, total)
    ax.text(i, height[i] + y_offset, str(round(total, 4)), ha='center',  weight='bold')
    #ax.text(totals.index[i], total + y_offset, round(total), ha='center',  weight='bold')



    #x, y = p.get_xy()
    #ax.annotate(f'{height:.0%}', str(atio[1]), ha='center')
    #ax.annotate(f'{height:.0%}', (splam_seq_len), ha='center')
    #ax.annotate(f'{height:.0%}', (splam_seq_len, spliceai_10k_noN, spliceai_10k_N), ha='center')

plt.savefig('information.png', dpi=300)
plt.show()

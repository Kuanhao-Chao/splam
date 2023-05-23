mkdir /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/

for sample in R2826  R2835  R2839  R2845  R2855  R2857  R2869  R2874  R2894  R2895; do
    mkdir /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/$sample/
    echo scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA_STAR/${sample}.bamAligned.sortedByCoord.out/fasta/ /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/$sample
    scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA_STAR/${sample}.bamAligned.sortedByCoord.out/fasta/ /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/$sample
    echo scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA_STAR/${sample}.bamAligned.sortedByCoord.out/junction.bed /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/$sample
    scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA_STAR/${sample}.bamAligned.sortedByCoord.out/junction.bed /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA_STAR/$sample
done

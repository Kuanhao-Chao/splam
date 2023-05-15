for sample in R2826  R2835  R2839  R2845  R2855  R2857  R2869  R2874  R2894  R2895; do
    echo scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA/$sample/fasta/ /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA/$sample/
    scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA/$sample/fasta/ /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA/$sample/
    echo scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA/$sample/junction.bed /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA/$sample/
    scp -r s3:/ccb/salz3/kh.chao/SPLAM/results/polyA/$sample/junction.bed /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/polyA/$sample/
done

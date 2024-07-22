module load samtools

mkdir -p test_output

python3 ../PrimeSpotter/PrimeSpotter.py --bam_file SF3B1.bam \
                                        --gtf_file SF3B1.gtf \
                                        --output-summary test_output/test_summary.txt \
                                        --genome-ref chr2.fasta \
                                        --processes 5 | samtools view -S -b | samtools sort > test_output/test.sorted.bam

samtools index test_output/test.sorted.bam                                    
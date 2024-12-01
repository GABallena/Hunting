#!/bin/bash

# Question 2: Number of alignments with mate unmapped in the region
echo "Question 2: Number of alignments with mate unmapped in the region:"
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '*'  # Mate unmapped
echo ""

# Question 3: Number of alignments with a deletion (D) in the region
echo "Question 3: Number of alignments with a deletion (D) in the region:"
samtools view athal_wu_0_A.region.bam | cut -f6 | grep -c 'D'  # Deletion
echo ""

# Question 4: Number of alignments with mate mapped to the same chromosome in the region
echo "Question 4: Number of alignments with mate mapped to the same chromosome in the region:"
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '='  # Mate mapped to the same chromosome
echo ""

# Question 7: Number of alignments with mate unmapped in the region
echo "Question 7: Number of alignments with mate unmapped in the region:"
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '*'  # Mate unmapped
echo ""

# Question 8: Number of alignments with a deletion (D) in the region
echo "Question 8: Number of alignments with a deletion (D) in the region:"
samtools view athal_wu_0_A.region.bam | cut -f6 | grep -c 'D'  # Deletion
echo ""

# Question 9: Number of alignments with mate mapped to the same chromosome in the region
echo "Question 9: Number of alignments with mate mapped to the same chromosome in the region:"
samtools view athal_wu_0_A.region.bam | cut -f7 | grep -c '='  # Mate mapped to the same chromosome
echo ""

# Question 10: Number of spliced alignments in the region
echo "Question 10: Number of spliced alignments in the region:"
samtools view athal_wu_0_A.region.bam | grep -c 'N'  # Spliced alignments (indicating introns)
echo ""

# Question 11: What alignment tool was used?
echo "Question 11: What alignment tool was used?"
samtools view -H athal_wu_0_A.bam | grep -m 1 -i "program"  # Aligning tool (program information)
echo ""

# Question 12: What is the read identifier (name) for the first alignment?
echo "Question 12: What is the read identifier (name) for the first alignment?"
samtools view athal_wu_0_A.bam | head -n 1 | cut -f1  # Read name for the first alignment
echo ""

# Question 13: What is the start position of this read’s mate on the genome?
echo "Question 13: What is the start position of this read’s mate on the genome?"
samtools view athal_wu_0_A.bam | head -n 1 | cut -f8  # Mate position (column 8)
echo ""

# Question 14: How many overlaps (each overlap is reported on one line) are reported?
echo "Question 14: How many overlaps (each overlap is reported on one line) are reported?"
bedtools intersect -a athal_wu_0_A.bam -b athal_wu_0_A_annot.gtf | wc -l  # Overlaps between BAM and GTF
echo ""

# Question 16: How many overlaps are reported?
echo "Question 16: How many overlaps are reported?"
bedtools intersect -a athal_wu_0_A.bam -b athal_wu_0_A_annot.gtf | wc -l  # Total overlaps
echo ""

# Question 17: How many of these overlaps are 10 bases or longer?
echo "Question 17: How many of these overlaps are 10 bases or longer?"
bedtools intersect -a athal_wu_0_A.bam -b athal_wu_0_A_annot.gtf -f 0.10 | wc -l  # Overlaps >=10 bases
echo ""

# Question 18: How many alignments overlap the annotations?
echo "Question 18: How many alignments overlap the annotations?"
bedtools intersect -a athal_wu_0_A.bam -b athal_wu_0_A_annot.gtf -wa | wc -l  # Alignments that overlap
echo ""

# Question 19: Conversely, how many exons have reads mapped to them?
echo "Question 19: How many exons have reads mapped to them?"
bedtools intersect -a athal_wu_0_A.bam -b athal_wu_0_A_annot.gtf -wa | cut -f9 | grep -c "exon"  # Exons with mapped reads
echo ""

# Question 20: How many BED records would be generated from the transcript annotations?
echo "Question 20: How many BED records would be generated from the transcript annotations?"
awk '$3 == "transcript" {print $1"\t"$4-1"\t"$5"\t"$9}' athal_wu_0_A_annot.gtf | wc -l  # BED record count for transcripts
echo ""

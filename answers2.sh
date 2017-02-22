datasets="/vol2/home/hsparks8/data-sets"

#question 1. Use bed tools to identify the size of the largest overlap
#between CTCF and H3K4me3 locations

TFBS="$datasets/bed/encode.tfbs.chr22.bed.gz"
H3K="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"
CTCF="$datasets/bedtools/ctcf.hela.chr22.bg.gz"
hg19="$datasets/fasta/hg19.chr22.fa"
bed="$datasets/bed/hg19interval.bed"
TSS="$datasets/bed/tss.hg19.chr22.bed.gz"
genes="$datasets/bed/genes.hg19.bed.gz"
genome="$datasets/genome/hg19.genome"

answer_1=$(bedtools intersect -a $TFBS -b $H3K \
    | awk 'BEGIN {OFS="\t"} ($4=="CTCf") {print 40 $3-$2}' \
    | sort -k5nr \
    | cut -f 5 \
    | head -n 1)

echo "answer-1: $answer_1"

#question 2.use bedtools to calculate the GC content of nucleotides  to on chr22 of hg19 genome build. report the GC content
#as fraction

answer_2=$(echo -e "chr22\t19000000\t19000500\n" > new.bed \
    | bedtools nuc -fi $hg19 -bed new.bed \
    | grep -v '^#' \
    | cut -f 5)

echo "answer-2: $answer_2"

#question 3. identify the length of the CTCF chip-seq peak with the
#largest mean signal

answer_3=$(bedtools map -a $TFBS -b $CTCF -c 4 -o mean \
    | sort -k5nr \
    | awk '($4=="CTCF")' \
    | awk 'BEGIN {OFS="\t"} (NR==1) {print $0, $3-$2}' \
    | cut -f 6)
    
echo "answer-3: $answer_3"

#question 4. id the gene promoter with highest median signal in
#cfcf.hela.chr22.bg.gz, report gene name

answer_4=$(bedtools slop -i $TSS -g $genome -l 1000 -r 0 -s \
    | bedtools sort \
    | bedtools map -a - -b $CTCF -c 4 -o median \
    | sort -k7nr \
    | cut -f 4 \
    | head -n 1)

echo "answer-4: $answer_4"

#question 5. id the longest interval on chr22 that is not covered by
#genes.hg19.bed.gz report interval: chr1:100-500


answer_5=$(bedtools sort -i $genes \
    | bedtools complement -i - -g $genome \
    | awk 'BEGIN {OFS="\t"} ($1=="chr22") {print $1, $2, $3, $3-$2}' \
    | sort -k4nr \
    | head -1 \
    | awk '{print $1":"$2"-"$3}')

echo "answer-5: $answer_5"




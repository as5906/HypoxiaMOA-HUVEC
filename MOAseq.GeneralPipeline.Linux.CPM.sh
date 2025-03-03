#!/bin/bash
helpFunction()
{
   echo ""
   echo "Usage: $0 -a [Paired-End File1] -b [Paired-End File2] -c [Reference Genome FASTA File] -d [Genome BWA Index Path and Prefix]"
   echo ""
   echo -e "\tOptional Parameters:"
   echo -e "\t\t-e [Effective Genome Size] ; Option to enter your own ideal effective genome size"
   echo -e "\t\t-f [int MAPQ] ; Option will skip alignments with MAPQ scores lower than this user provided integer when converting SAM file to BAM file. Default is 0"
   echo -e "\t\t-g ; Option to skip flash and not merge paired-ends"
   exit 1 # Exit script after printing help
}
mkdir tempFiles
MAPQ=0
doFlash=true
while getopts "a:b:c:d:e::f::g" opt
do
   case "$opt" in
      a ) Pair1="$OPTARG" ;;
      b ) Pair2="$OPTARG" ;;
      c ) Genome="$OPTARG" ;;
      d ) IndexPath="$OPTARG" ;;
      e ) EffectiveGenome="$OPTARG" ;;
      f ) MAPQ="$OPTARG" ;;
      g ) doFlash=false ;;
      ? ) helpFunction ;;
   esac
done

samtools faidx $Genome -o tmp1.fai
awk '{OFS="\t" ; print $1,$2}' tmp1.fai > $Genome.fai
rm tmp1.fai

echo "Starting the Run. Removing adapters from MOA raw reads..."

	cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o Pair1.trim.gz -p Pair2.trim.gz $Pair1 $Pair2

#       SeqPurge \
#	-min_len 20 \
#	-threads 48 \
#	-qcut 0 \
#	-prefetch 50000 \
#	-in1 $Pair1 \
#	-in2 $Pair2 \
#	-out1 Pair1.trim.gz \
#	-out2 Pair2.trim.gz \
#	-summary Summary.trim.stats \
#	-progress 120000 \
#	-qc ReadQC.trim.qc

if [ $doFlash = true ];
then 
	echo -e "Merging Paired-Ends.."
	flash \
 	      	$Pair1 \
        	$Pair2 \
        	-t 40 \
        	--cap-mismatch-quals \
        	-M 100 \
        	-d Flash_Output \
		-q
fi

if [ $doFlash = true ];
then
	echo -e "BWA aligning to genome and then converting to sorted BAM.."
	bwa-mem2 mem \
        	-M -t 40 \
        	$IndexPath \
        	Flash_Output/out.extendedFrags.fastq  > output.sam
	samtools view \
		output.sam \
		-b \
		-q $MAPQ \
		-@ 2 | samtools sort \
		-@ 2 \
		-m 10G \
		-o output.sort.bam
else
	echo -e "BWA aligning to genome and then converting to sorted BAM.."
	bwa-mem2 mem \
                $IndexPath \
                Pair1.trim.gz \
                Pair2.trim.gz  | samtools view \
                -b \
                -q $MAPQ \
                -@ 2 \
                | samtools sort
                -@ 2 \
                -m 10G \
                -o output.sort.bam
fi

echo -e "Generating BAM statistics file..."
samtools stats output.sort.bam > output.sort.bam.stats.txt
frg_length=`gawk -v OFS='\t' '{if($2=="average" && $3=="first" && $4=="fragment"){print $6}}' output.sort.bam.stats.txt`

if [[ -z $EffectiveGenome ]]
then
	echo -e "Calculating Effective Genome Size with fragment length of  ${frg_length} base pairs"
	unique-kmers.py -q -k ${frg_length} -R effectiveGenome.txt ${Genome}
	EG=`awk '{if(NR==1){print $1}}' effectiveGenome.txt`

	echo -e "Effective genome size is ${EG}"
else
	EG=$EffectiveGenome
	echo -e "Effective genome size is ${EG}"
fi

echo -e "Converting BAM to BED..."
bedtools bamtobed -i output.sort.bam > output.bed

echo -e "Producing a bed file with fragment centers of 20 base pairs (shortenedReads)..."
gawk -v OFS='\t' -v FL=20 '{if(($2+$3) %2 == 0){a=int(($2+$3)/2+0.5) ; if(a-(FL/2)>0){$2=a-(FL/2)} else{$2=0} ; $3=a+(FL/2) ; print $0} else{a=int(($2+$3)/2+rand()) ; if(a-(FL/2)>0){$2=a-(FL/2)} else{$2=0}  ; $3=a+(FL/2) ; print $0}}' output.bed > output.shortenedReads.bed

echo -e "Sorting BED files..."
LC_COLLATE=C sort -k1,1 -k2,2n output.bed > output.sort.bed
rm output.bed
LC_COLLATE=C sort -k1,1 -k2,2n output.shortenedReads.bed > output.shortenedReads.sort.bed
rm output.shortenedReads.bed 

unireads=`wc -l output.sort.bed | awk '{print $1}'`

x=$(printf "scale =10; %d/%d \n" $EG $unireads | bc);
#x=$( printf "%s%d\n" 'scale = 10; 1000000/' ((( $allreads-$unireads )/2)+ $unireads) | bc);

echo -e "Genertaing BEDGRAPHS..."
genomeCoverageBed -scale $x -split -bg -i output.sort.bed -g $Genome.fai > output.sort.bg
genomeCoverageBed -scale $x -split -bg -i output.shortenedReads.sort.bed -g $Genome.fai > output.shortenedReads.sort.bg

echo -e "Generating BIGWIGS..."
bedGraphToBigWig output.sort.bg $Genome.fai output.sort.bw
bedGraphToBigWig output.shortenedReads.sort.bg $Genome.fai output.shortenedReads.sort.bw

mv output.sam output.sort.bam output.sort.bam.stats.txt tempFiles

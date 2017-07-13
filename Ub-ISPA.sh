command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed.  Aborting." >&2; exit 1; }
echo "cutadapt found"

command -v bowtie2 >/dev/null 2>&1 || { echo "I require bowtie2 but it's not installed.  Aborting." >&2; exit 1; }
echo "bowtie2 found"

if [ $# -eq 0 ]; then
	
echo "Please select option (1 or 2)"
echo "1- Single End Reads"
echo "2- Paired End Reads"
read Option

echo "Directory name to save files"
read dir1
mkdir $dir1

echo "Please select host Organism for alignment"
echo "h-human"
echo "o-other"
read host

if test "$host" = 'h'; then
	echo "hg19 for alignment? (y/n)"
	read res
	if test "$res" = 'y'; then
		indexfile=hg19
		bedfile=refGene_hg19_genebody.bed
		
	else
		echo "Please enter the indexed (Bowtie2) host genome file"
		read indexfile
		echo "Please enter the bed file"
		read bedfile	
	fi

elif test "$host" = 'o'; then
	echo "Please enter the indexed (Bowtie2) host genome file"
	read indexfile
	echo "Please enter the bed file"
	read bedfile

else 
	echo "Wrong selection"
	exit
fi
	

if test $Option -eq 1; then
	
	echo "Plesae enter the Fastq file"
	read File
	
	echo "Please enter the Forward Primer"
	read FP

	echo "Please enter the Reverse Primer"
	read RP

# to cut the adaptor of 5' end and select the trimmed reads only
 cutadapt --discard-untrimmed -g ${FP} -o ${dir1}/adapter_trimmed.fastq $File
if [ $? -ne 0 ]; then
    echo "5' primer trimming and selection failed"
exit
fi


#to cut the adaptor at 3' end
cutadapt --minimum-length 30 -a ${RP} -o ${dir1}/adapter_trimmed_3prime.fastq ${dir1}/adapter_trimmed.fastq

if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

#bowtie2 alignment
bowtie2 -q -k 2 -x ${indexfile} -U ${dir1}/adapter_trimmed_3prime.fastq -S ${dir1}/alignment.sam --al ${dir1}/concord.fastq

if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${dir1}/concord.fastq | uniq -c | grep -w -v 2 > ${dir1}/aligned_one_reads.txt
sed 's/1//' ${dir1}/aligned_one_reads.txt | sed 's/1:N:0:1//' > ${dir1}/aligned_one_reads_without_line_no.txt

#to convert fastq in single line per sequence
perl -pne 'if($.%4){s/\n/\t/;}' ${dir1}/adapter_trimmed_3prime.fastq > ${dir1}/adapter_trimmed_3prime_single_line.txt

awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${dir1}/aligned_one_reads_without_line_no.txt ${dir1}/adapter_trimmed_3prime_single_line.txt > ${dir1}/aligned_only_once.txt

sed 's/\t/\n/g' ${dir1}/aligned_only_once.txt > ${dir1}/once_aligned_reads.fastq

#rm once_aligned_reads.fastq
#for a in `cat aligned_one_reads_without_line_no.txt`; do 
	#grep -A3 $a adapter_trimmed_3prime.fastq >> once_aligned_reads.fastq
#done

echo "align read allowing one mismatch"
bowtie2 -q -k 2 -N 1 -x ${indexfile} -U ${dir1}/once_aligned_reads.fastq -S ${dir1}/alignment_1mismatch.sam --al ${dir1}/concord_1mismatch.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${dir1}/concord_1mismatch.fastq | uniq -c | grep -w -v 2 > ${dir1}/aligned_one_reads_1mismatch.txt
sed 's/1//' ${dir1}/aligned_one_reads_1mismatch.txt | sed 's/1:N:0:1//' > ${dir1}/aligned_one_reads_without_line_no_1mismatch.txt

awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${dir1}/aligned_one_reads_without_line_no_1mismatch.txt ${dir1}/adapter_trimmed_3prime_single_line.txt > ${dir1}/aligned_only_once_1mismatch.txt

sed 's/\t/\n/g' ${dir1}/aligned_only_once_1mismatch.txt > ${dir1}/once_aligned_reads_1mismatch.fastq

echo "align read allowing two mismatch"
bowtie2 -q -k 2 --mp 2,0 -x ${indexfile} -U ${dir1}/once_aligned_reads_1mismatch.fastq -S ${dir1}/alignment_2mismatch.sam --al ${dir1}/concord_2mismatch.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${dir1}/concord_2mismatch.fastq | sed 's/@//' | uniq -c | grep -w -v 2 > ${dir1}/aligned_one_reads_2mismatch.txt
sed 's/1//' ${dir1}/aligned_one_reads_2mismatch.txt | sed 's/1:N:0:1//' > ${dir1}/aligned_one_reads_without_line_no_2mismatch.txt

#to prepare final sam files
rm ${dir1}/final.sam
head -86 ${dir1}/alignment_2mismatch.sam > ${dir1}/final.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${dir1}/aligned_one_reads_without_line_no_2mismatch.txt ${dir1}/alignment_2mismatch.sam >> ${dir1}/final.sam

#to prepare bed files
samtools view -bS ${dir1}/final.sam > ${dir1}/final.bam
bedtools bamtobed -i ${dir1}/final.bam > ${dir1}/aligned.bed

elif test $Option -eq 2; then

	echo "Plesae enter the read1 Fastq file"
	read Read1

	echo "Please enter the Forward Primer"
	read FP1

	echo "Please enter the Reverse Primer"
	read RP1
	
	echo "Plesae enter the read2 Fastq file"
	read Read2
	
	echo "Please enter the Forward Primer"
	read FP2

	echo "Please enter the Reverse Primer"
	read RP2

echo "to remove distal adaptors of R1"
cutadapt --no-indels -a ${RP1} -o ${dir1}/reads1_DA_trimmed.fastq $Read1

if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

echo "to remove distal adaptors of R2"
cutadapt --no-indels -a ${RP2} -o ${dir1}/reads2_DA_trimmed.fastq $Read2

if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

echo "to cut the adaptor of 5' end and select the trimmed reads only"
cutadapt --minimum-length 30 --discard-untrimmed -g ${FP1} -G ${FP2} -o ${dir1}/trimmed_R1.fastq -p ${dir1}/trimmed_R2.fastq ${dir1}/reads1_DA_trimmed.fastq ${dir1}/reads2_DA_trimmed.fastq

if [ $? -ne 0 ]; then
    echo "5' primer trimming and selection failed"
exit
fi

echo "bowtie2 command"
bowtie2 -q --no-discordant -k 2 -x ${indexfile} -1 ${dir1}/trimmed_R1.fastq -2 ${dir1}/trimmed_R2.fastq -S ${dir1}/bowtie4.sam --al-conc ${dir1}/concord_1.fastq

if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

#to extract the reads aligned only once concordantly
sed -n '1~4p' ${dir1}/concord_1.1.fastq | sed 's/@//' | sed 's/1:N:0:1//' | uniq -c | grep -w -v 2 > ${dir1}/aligned_one_reads.txt
sed 's/1//' ${dir1}/aligned_one_reads.txt > ${dir1}/aligned_one_reads_without_line_no.txt

rm ${dir1}/concord_aligned_one.sam
head -86 ${dir1}/bowtie4.sam > ${dir1}/concord_aligned_one.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${dir1}/aligned_one_reads_without_line_no.txt ${dir1}/bowtie4.sam >> ${dir1}/concord_aligned_one.sam

samtools view -bS ${dir1}/concord_aligned_one.sam > ${dir1}/concord_aligned_one.bam

bedtools bamtobed -i ${dir1}/concord_aligned_one.bam > ${dir1}/concord_aligned_one.bed
grep /1 ${dir1}/concord_aligned_one.bed > ${dir1}/aligned.bed

else

	echo "Wrong Selection"
	exit 

fi
###################################################################
else	
	mkdir $2
	indexfile=$3
	bedfile=$4
	dir1=$2
if [ $1 -eq 1 ]; then
if [ $# -eq 7 ]; then
	

# to cut the adaptor of 5' end and select the trimmed reads only
cutadapt --discard-untrimmed -g $6 -o ${2}/adapter_trimmed.fastq $5

if [ $? -ne 0 ]; then
    echo "5' primer trimming and selection failed"
exit
fi


#to cut the adaptor at 3' end
cutadapt --minimum-length 30 -a $7 -o ${2}/adapter_trimmed_3prime.fastq ${2}/adapter_trimmed.fastq
if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

#bowtie2 alignment
bowtie2 -q -k 2 -x ${indexfile} -U ${2}/adapter_trimmed_3prime.fastq -S ${2}/alignment.sam --al ${2}/concord.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${2}/concord.fastq | uniq -c | grep -w -v 2 > ${2}/aligned_one_reads.txt
sed 's/1//' ${2}/aligned_one_reads.txt | sed 's/1:N:0:1//' > ${2}/aligned_one_reads_without_line_no.txt

#to convert fastq in single line per sequence
perl -pne 'if($.%4){s/\n/\t/;}' ${2}/adapter_trimmed_3prime.fastq > ${2}/adapter_trimmed_3prime_single_line.txt

awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${2}/aligned_one_reads_without_line_no.txt ${2}/adapter_trimmed_3prime_single_line.txt > ${2}/aligned_only_once.txt

sed 's/\t/\n/g' ${2}/aligned_only_once.txt > ${2}/once_aligned_reads.fastq

#rm once_aligned_reads.fastq
#for a in `cat aligned_one_reads_without_line_no.txt`; do 
	#grep -A3 $a adapter_trimmed_3prime.fastq >> once_aligned_reads.fastq
#done

echo "align read allowing one mismatch"
bowtie2 -q -k 2 -N 1 -x ${indexfile} -U ${2}/once_aligned_reads.fastq -S ${2}/alignment_1mismatch.sam --al ${2}/concord_1mismatch.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${2}/concord_1mismatch.fastq | uniq -c | grep -w -v 2 > ${2}/aligned_one_reads_1mismatch.txt
sed 's/1//' ${2}/aligned_one_reads_1mismatch.txt | sed 's/1:N:0:1//' > ${2}/aligned_one_reads_without_line_no_1mismatch.txt

awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${2}/aligned_one_reads_without_line_no_1mismatch.txt ${2}/adapter_trimmed_3prime_single_line.txt > ${dir1}/aligned_only_once_1mismatch.txt

sed 's/\t/\n/g' ${2}/aligned_only_once_1mismatch.txt > ${2}/once_aligned_reads_1mismatch.fastq

echo "align read allowing two mismatch"
bowtie2 -q -k 2 --mp 2,0 -x ${indexfile} -U ${2}/once_aligned_reads_1mismatch.fastq -S ${2}/alignment_2mismatch.sam --al ${2}/concord_2mismatch.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

# to extract reads aligned only once
sed -n '1~4p' ${2}/concord_2mismatch.fastq | sed 's/@//' | uniq -c | grep -w -v 2 > ${2}/aligned_one_reads_2mismatch.txt
sed 's/1//' ${2}/aligned_one_reads_2mismatch.txt | sed 's/1:N:0:1//' > ${2}/aligned_one_reads_without_line_no_2mismatch.txt

#to prepare final sam files
rm ${2}/final.sam
head -86 ${2}/alignment_2mismatch.sam > ${2}/final.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${2}/aligned_one_reads_without_line_no_2mismatch.txt ${2}/alignment_2mismatch.sam >> ${2}/final.sam

#to prepare bed files
samtools view -bS ${2}/final.sam > ${2}/final.bam
bedtools bamtobed -i ${2}/final.bam > ${2}/aligned.bed

else
echo "Incorrent number of commands given"
	exit 
fi

elif [ $1 -eq 2 ]; then
if [ $# -eq 10 ]; then
	Read1=$5
	FP1=$6
	RP1=$7
	Read2=$8
	FP2=$9
	RP2=$10
	

echo "to remove distal adaptors of R1"
cutadapt --no-indels -a $7 -o ${2}/reads1_DA_trimmed.fastq $5
if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

echo "to remove distal adaptors of R2"
cutadapt --no-indels -a $10 -o ${2}/reads2_DA_trimmed.fastq $8
if [ $? -ne 0 ]; then
    echo "3' primer not trimmed"
exit
fi

echo "to cut the adaptor of 5' end and select the trimmed reads only"
cutadapt --minimum-length 30 --discard-untrimmed -g $6 -G $9 -o ${2}/trimmed_R1.fastq -p ${2}/trimmed_R2.fastq ${2}/reads1_DA_trimmed.fastq ${2}/reads2_DA_trimmed.fastq
if [ $? -ne 0 ]; then
    echo "5' primer trimming and selection failed"
exit
fi

echo "bowtie2 command"
bowtie2 -q --no-discordant -k 2 -x ${indexfile} -1 ${2}/trimmed_R1.fastq -2 ${2}/trimmed_R2.fastq -S ${2}/bowtie4.sam --al-conc ${2}/concord_1.fastq
if [ $? -ne 0 ]; then
    echo "problem with bowtie2 alignment"
exit
fi

#to extract the reads aligned only once concordantly
sed -n '1~4p' ${2}/concord_1.1.fastq | sed 's/@//' | sed 's/1:N:0:1//' | uniq -c | grep -w -v 2 > ${2}/aligned_one_reads.txt
sed 's/1//' ${2}/aligned_one_reads.txt > ${2}/aligned_one_reads_without_line_no.txt

rm ${2}/concord_aligned_one.sam
head -86 ${2}/bowtie4.sam > ${2}/concord_aligned_one.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${2}/aligned_one_reads_without_line_no.txt ${2}/bowtie4.sam >> ${2}/concord_aligned_one.sam

samtools view -bS ${2}/concord_aligned_one.sam > ${2}/concord_aligned_one.bam

bedtools bamtobed -i ${2}/concord_aligned_one.bam > ${2}/concord_aligned_one.bed
grep /1 ${2}/concord_aligned_one.bed > ${2}/aligned.bed

else
echo "Incorrent number of commands given"
	exit 
fi

else

	echo "Wrong Selection"
	exit 

fi
fi


#######################################################

#to add integrations sites to bed file   
grep -w - ${dir1}/aligned.bed | cut -f 1,2,3,6 | uniq > ${dir1}/no_of_reads_aligned_atIS_-strand.txt
grep -w + ${dir1}/aligned.bed | cut -f 1,2,3,6 | uniq > ${dir1}/no_of_reads_aligned_atIS_+strand.txt

rm ${dir1}/A_pos.bed
rm ${dir1}/A_neg.bed
for a in `cat ${dir1}/no_of_reads_aligned_atIS_+strand.txt| cut -f2`; do
      b=`expr $a - 1`
	echo $b >> ${dir1}/A_pos.bed
done

paste ${dir1}/no_of_reads_aligned_atIS_+strand.txt ${dir1}/A_pos.bed > ${dir1}/final_+_bed_file.bed
  

for c in `cat ${dir1}/no_of_reads_aligned_atIS_-strand.txt| cut -f3`; do
      	echo $c >> ${dir1}/A_neg.bed
done
paste ${dir1}/no_of_reads_aligned_atIS_-strand.txt ${dir1}/A_neg.bed > ${dir1}/final_-_bed_file.bed

cat ${dir1}/final_+_bed_file.bed ${dir1}/final_-_bed_file.bed > ${dir1}/final_bed_file.bed

#to create file for TSS proximal, integenic and intragenic calling
cut -f 1,4,5 ${dir1}/final_bed_file.bed | sort | uniq | awk '{print $1 "\t" $3 "\t" $3 "\t" $2}' > ${dir1}/unique_IS_final_file.bed

#to create bed files
sort -k1,1 -k2,2n ${bedfile} > ${dir1}/sorted_genebody.bed
cut -f 1,2 ${dir1}/sorted_genebody.bed > ${dir1}/temp.bed
cut -f 2,4,5,6 ${dir1}/sorted_genebody.bed > ${dir1}/temp1.bed
paste ${dir1}/temp.bed ${dir1}/temp1.bed > ${dir1}/TSSbedfile.bed
sort -k1,1 -k2,2n ${dir1}/TSSbedfile.bed > ${dir1}/sorted_TSS_bedfile.bed
rm ${dir1}/temp.bed ${dir1}/temp1.bed

#to find the TSS proximity IS
sort -k1,1 -k2,2n ${dir1}/unique_IS_final_file.bed > ${dir1}/sorted_unique_IS_final_file.bed

bedtools closest -a ${dir1}/sorted_unique_IS_final_file.bed -b ${dir1}/sorted_TSS_bedfile.bed -t first -d > ${dir1}/TSS_location.bed

awk ' $11 <= 2500 ' ${dir1}/TSS_location.bed | cut -f 1,2,4,8,9 | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/TSS_proximal.txt

#to identify intragenic IS
awk ' $11 <= 2500 ' ${dir1}/TSS_location.bed | cut -f 1,2,4 | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/TSS_proximal1.txt

bedtools intersect -a ${dir1}/sorted_unique_IS_final_file.bed -b ${dir1}/sorted_genebody.bed -wa -wb | cut -f 1,2,4 | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/intragenic_TSS.txt

awk 'NR==FNR{a[$0];next} !($0 in a)' ${dir1}/TSS_proximal1.txt ${dir1}/intragenic_TSS.txt | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/intragenic.txt

#to identify intergenic IS
cat ${dir1}/intragenic.txt ${dir1}/TSS_proximal1.txt | sort -u -k 1,1 -k 2,2 -k 3,3 | cut -f 1,2,3 > ${dir1}/temp3.txt
cut -f 1,2,4 ${dir1}/sorted_unique_IS_final_file.bed | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/temp4.txt
awk 'NR==FNR{a[$0];next} !($0 in a)' ${dir1}/temp3.txt ${dir1}/temp4.txt | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/intergenic.txt

#to annotate intragenic file
awk '{print $1 "\t" $2 "\t" $2 "\t" $3}' ${dir1}/intragenic.txt > ${dir1}/intragenic.bed
sort -k1,1 -k2,2n ${dir1}/intragenic.bed > ${dir1}/sorted_intragenic.bed
bedtools intersect -a ${dir1}/sorted_intragenic.bed -b ${dir1}/sorted_genebody.bed -wa -wb | cut -f 1,2,4,8,9 | sort -u -k 1,1 -k 2,2 -k 3,3 > ${dir1}/final_intragenic.txt

#rm adapter_trimmed.fastq adapter_trimmed_3prime.fastq alignment.sam concord.fastq aligned_one_reads.txt aligned_one_reads_without_line_no.txt adapter_trimmed_3prime_single_line.txt aligned_only_once.txt once_aligned_reads.fastq alignment_1mismatch.sam concord_1mismatch.fastq aligned_one_reads_1mismatch.txt aligned_one_reads_without_line_no_1mismatch.txt aligned_only_once_1mismatch.txt aligned_only_once_1mismatch.txt alignment_2mismatch.sam concord_2mismatch.fastq aligned_one_reads_2mismatch.txt aligned_one_reads_without_line_no_2mismatch.txt final.sam  final.bam aligned.bed no_of_reads_aligned_atIS_-strand.txt no_of_reads_aligned_atIS_+strand.txt unique_IS_final_file.bed sorted_unique_IS_final_file.bed TSS_location.bed TSS_proximal1.txt intragenic_TSS.txt temp3.txt temp4.txt intragenic.bed  sorted_intragenic.bed A_pos.bed A_neg.bed final_+_bed_file.bed final_-_bed_file.bed final_bed_file.bed once_aligned_reads.fastq once_aligned_reads_1mismatch.fastq


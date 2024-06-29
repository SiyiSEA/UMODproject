# This pipeline use some code from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095


#! /bin/bash
#PBS -l nodes=1:ppn=4:centos6,cput=24:00:00,walltime=48:00:00
#PBS -N non_discovery_q10
#PBS -d /export/home/biostuds/2549343w/My_Project/scripts
#PBS -m abe
#PBS -M SYwangSEA@gmail.com
#PBS -q bioinf-stud

# decompress the raw reads files and rename the file
path="/export/biostuds/2549343w/My_Project/raw_data"
for sample in 888THPWtHS_S1_L001_R1_001 911THPKOHS_S3_L001_R1_001 720THPWtNS_S5_L001_R1_001 721THPWtNS_S6_L001_R1_001 722THPWtNS_S7_L001_R1_001 756THPKONS_S8_L001_R1_001 757THPKONS_S9_L001_R1_001 758THPKONS_S10_L001_R1_001 884THPWtHS_S11_L001_R1_001 887THPWtHS_S12_L001_R1_001 888THPWtHS_S1_L001_R1_001 910THPKOHS_S2_L001_R1_001 911THPKOHS_S3_L001_R1_001 912THPKOHS_S4_L001_R1_001
do
	gzip -d ${path}/${sample}.fastq.gz
	name=${sample:3}
	new_name=${name%????????????}
	mv ${path}/${sample}.fastq ${path}/$new_name.fastq
done

# build the hisat2 index
ref_genome = "/export/biostuds/2549343w/My_Project/hisat2_index/Mus_musculus.GRCm39.dna.primary_assembly.fa"
mv ${ref_genome} genome.fa
hisat2-build -p 16 genome.fa genome

#resource file
data_path="/export/biostuds/2549343w/My_Project/raw_data"
adapter="/export/biostuds/2549343w/My_Project/raw_data_test/adapter_2.fa"
hisat2_index="/export/biostuds/2549343w/My_Project/hisat2_index/genome"
GTF="/export/biostuds/2549343w/My_Project/Ensembl/Mus_musculus.GRCm39.106_no_gene.gtf"

# make few subdirs unless they exist
#fastqc_before="/export/biostuds/2549343w/My_Project/fastqc_before"
fastqc_after_10q="/export/biostuds/2549343w/My_Project/fastqc_after_10q"
new_trimmed="/export/biostuds/2549343w/My_Project/quality10/trimmed"
hisat2_result="/export/biostuds/2549343w/My_Project/quality10/nondiscovery/hisat2_result"
stringtie_result="/export/biostuds/2549343w/My_Project/quality10/nondiscovery/stringtie_result"
matrix_result="/export/biostuds/2549343w/My_Project/quality10/matrix_result"

mkdir -p ${new_trimmed}
#mkdir -p ${fastqc_before}
mkdir -p ${fastqc_after_10q}
mkdir -p ${hisat2_result}
mkdir -p ${stringtie_result}
mkdir -p ${matrix_result}

gtf_list="${matrix_result}/nondis_list.gtf.txt"
rm -f $gtf_list

for sample in THPWtHS_S1 THPKOHS_S2 THPKOHS_S3 THPKOHS_S4 THPWtNS_S5 THPWtNS_S6 THPWtNS_S7 THPKONS_S8 THPKONS_S9 THPKONS_S10  THPWtHS_S11 THPWtHS_S12
do
	# make few variables to store file names
	raw="${data_path}/${sample}.fastq"
	# check quality before trimming
	#fastqc -o ${fastqc_before} ${raw}
	trim1="${new_trimmed}/${sample}_tr1.fastq"
	trim2="${new_trimmed}/${sample}_tr2.fastq"
	# running scythe/sickle
	scythe -o ${trim1} -a ${adapter} -q sanger ${raw}
	sickle se -f ${trim1} -t sanger -o ${trim2} -q 10 -l 67
	fastqc -o ${fastqc_after_10q} ${trim2}
	rm $trim1
	
	# align reads to the reference genome
	hisat2 --rna-strandness R -p 10 --phred33 -x $hisat2_index -U $trim2 -S "${hisat2_result}/${sample}.sam" 
	# comvert sam file to bam file
	samtools view -bS -o "${hisat2_result}/${sample}.bam" "${hisat2_result}/${sample}.sam"
	# sort the bam file
	samtools sort "${hisat2_result}/${sample}.bam" -o "${hisat2_result}/${sample}_sort.bam"
	# index the sorted bam file
	samtools index "${hisat2_result}/${sample}_sort.bam"
	# delete no need files
	#rm "${hisat2_result}/${sample}.bam" "${hisat2_result}/${sample}.sam"
	# assemble reads to transcripts
	sub_stringtie="${stringtie_result}/${sample}"
	mkdir -p $sub_stringtie
	stringtie "${hisat2_result}/${sample}_sort.bam" -p 10 -B -e -G $GTF --rf -o "${sub_stringtie}/${sample}.gtf"
	# make a list of gtf file
	gtf_line="${sample} ${sub_stringtie}/${sample}.gtf"
	echo ${gtf_line} >> ${gtf_list}

done
gene="${matrix_result}/gene_count_matrix_nondiscovery_10q.csv"
trans="${matrix_result}/transcript_count_matrix_nondiscovery_10q.csv"
python2.7 /export/projects/polyomics/App/prepDE.py -i ${gtf_list} -g $gene -t $trans


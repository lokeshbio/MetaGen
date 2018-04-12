-   [Sequence processing](#sequence-processing)
    -   [Obtaining reads](#obtaining-reads)
    -   [Quality check](#quality-check)
    -   [Quality trim, filter and check (again)](#quality-trim-filter-and-check-again)
-   [Assembly](#assembly)
    -   [SortMeRNA](#sortmerna)
    -   [Assembly program](#assembly-program)
    -   [Assembly quality](#assembly-quality)
-   [Files for Binning](#files-for-binning)
    -   [Mapping](#mapping)
    -   [Coverage](#coverage)
    -   [Essential genes and their taxonomy](#essential-genes-and-their-taxonomy)
    -   [Quick check on the abundance of our favourite taxa groups](#quick-check-on-the-abundance-of-our-favourite-taxa-groups)
-   [Binning](#binning)
    -   [Coverage bi-plot](#coverage-bi-plot)
    -   [Visualization of TACK only](#visualization-of-tack-only)
-   [Shearing methods for library prep](#shearing-methods-for-library-prep)

Sequence processing
-------------------

This script will give the guidelines how the sequencing data is processed right from the time in obtaining the data from the sequencing facility up until obatining genomes from the metagenomes and so on!

### Obtaining reads

The reads come in 'bam' files from the sequencing facility and then we use samtools and bedtools to obtaing the PE sequncing data from the bam files.

``` bash
declare -A samples=([60858]=Sec-20-NSS [60859]=Sec-20-F1 [60860]=H2-2N [60861]=H2-9N [60862]=H2-17N [60863]=H2-2K [60864]=H2-9K [60865]=H2-17K [60866]=BF3-8N [60867]=BF3-8L [60868]=NV-NP [60870]=H2-2N-Enz [60871]=H2-2K-Enz [60872]=BF3-8N-Enz [60873]=BF3-8L-Enz)

for filename in ./*bam; do
    bin_name=`basename $filename|cut -f 1 -d '_'`
    #echo "${samples[$bin_name]}"
    samtools sort -n -@ 8 -o ${samples[$bin_name]}'.qsort.bam' $filename
    bedtools bamtofastq -i ${samples[$bin_name]}'.qsort.bam' -fq ${samples[$bin_name]}'.R1.fq' -fq2 ${samples[$bin_name]}'.R2.fq'
    #prokka $filename --outdir $bin_name --prefix $bin_name --kingdom Archaea --cpus 16 --fast --rnammer
    #nice -5 diamond blastp --query $filename --db /proj/Lokesh/Databases/arCOGs/arCOG_proteins --daa $bin_name'.arCOG.daa' --threads 10
done
```

here in this case we had about 15 samples and based on the sample numbers, I also renamed the sequence files to their corresponding biologically relevant sample names.

### Quality check

Then these reads are processed through 'fastqc' which is a way to look into the sequencing qulaity of the samples separately and then we can also use 'multiqc' to see the quality of the samples together! here we can alsio use multi-threading to increase the speed of

``` bash
fastqc -t 8 *fq
multiqc -s -i Logan_NV -o Logan_NV_multiqc Logan_NV
```

here in the multiqc commnad Logan\_NV at the end of the command denotes the folder containing all the fastqc results you want to combine together in one multiqc plot!

### Quality trim, filter and check (again)

Now we want to use these sequencing data to go through quality control steps! This involves using 'trimmomatic' for teimming the reads based on their quality and then to use 'prinseq' for filtering the reads based on quality and length. Then we check how they look using the 'fastqc' and 'multiqc' again.

``` bash
trimmomatic PE -threads 16 -phred33 -trimlog $bin_name'.trim.log'  $bin_name'.R1.fq.gz' $bin_name'.R2.fq.gz' $bin_name'.paired.R1.fq.gz' $bin_name'.unpaired.R1.fq.gz' $bin_name'.paired.R2.fq.gz' $bin_name'.unpaired.R2.fq.gz' ILLUMINACLIP:/apps/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50 HEADCROP:6
```

This will create 4 files as shown from the script above 2 paired and 2 unpaired files with reads! Now when I checked this with 'fastqc' and 'multiqc', I found that there were wtill some sequences with homopolymer errors just with 'G's in them! this could be due to a poor library preparation method! But these sequences only constitute a minor proportion of the entire dataset and weshould be able to work with them well anyway!

So, now we will use 'PrinSeq-Lite' to filter these sequences that have a lot of homopolymers in them. You can do this by the following statement.

``` bash
    prinseq-lite.pl -fastq temp_paired_r1 -fastq2 temp_paired_r2 -log prinseq_filter.log -min_len 50 -custom_params "A 80%;T 80%;G 80%;C 80%;G 15;C 15" -out_good $bin_name'.paired' -out_bad null
    prinseq-lite.pl -fastq temp_unpaired_r1 -log prinseq_filter.log -min_len 50 -custom_params "A 80%;T 80%;G 80%;C 80%;G 15;C 15" -out_good $bin_name'.unpaired.R1' -out_bad null
    prinseq-lite.pl -fastq temp_unpaired_r2 -log prinseq_filter.log -min_len 50 -custom_params "A 80%;T 80%;G 80%;C 80%;G 15;C 15" -out_good $bin_name'.unpaired.R2' -out_bad null
```

As you can see, we need three different 'prinseq' statements to be able to process all the four different files from the trimmomatic. Once this filteration process is done, we can run 'fastqc' and 'multiqc' again to make sure that the filtered reads are of good quality!

Assembly
--------

After this step the quality of the sequences looked pretty good!

### SortMeRNA

Just out of curiosity, we wanted to check quickly how the community of the microbes looked in the metagenomic dataset. So, we screened for rRNA signatures already in the read level on the metagenomic dataset as it is faster. So we ran this program called 'sortmerna'.

``` bash
    pigz -dc Sequences/Step2_filtered/*/$bin_name'.paired_1.fastq.gz' > temp_paired_r1
    pigz -dc Sequences/Step2_filtered/*/$bin_name'.paired_2.fastq.gz' > temp_paired_r2
    merge-paired-reads.sh temp_paired_r1 temp_paired_r2 temp_paired_int
    sortmerna --ref /apps/sortmerna/2.1/rRNA_databases/silva-arc-16s-id95.fasta,/apps/sortmerna/2.1/index/silva-arc-16s-db:/apps/sortmerna/2.1/rRNA_databases/silva-arc-23s-id98.fasta,/apps/sortmerna/2.1/index/silva-arc-23s-db:/apps/sortmerna/2.1/rRNA_databases/silva-bac-16s-id90.fasta,/apps/sortmerna/2.1/index/silva-bac-16s-db:/apps/sortmerna/2.1/rRNA_databases/silva-bac-23s-id98.fasta,/apps/sortmerna/2.1/index/silva-bac-23s-db  --reads temp_paired_int --sam --num_alignments 1 --fastx --aligned SortmeRNA_files/$bin_name'SortMeRNA__sam_rRNA' -a 16 --log 
    rm temp*
```

As the program will only take fastq files without 'gzipped' we ran the above command for all the metagenomic samples we had over a bash script. The reference database for the rRNA dataset indexed dor sortmerna was also found in the 'CUBE' server in the location mentioned above in the script.

You get sam files by running this script and then you get the 'Reference ids' in the sam file and I made a list of the sortmeRNA reference ids and their respective taxonomic path in a file. Then we can actually calculate the taxonomy abundance by the following command.

``` bash
tail -n +3 Sec-20-NSSSortMeRNA__sam_rRNA.sam|cut -f 3|sort|uniq -c|perl -pe 's/;/\t/g' > Sec-20-NSS.sam.abund.txt
multiple_replace.py Ref_sortmerna.tab Sec-20-NSS.sam.abund.txt Sec-20-NSS.sam.taxa.abund.txt
perl -pe 's/ +/ /g' Sec-20-NSS.taxa.abund.txt|perl -pe 's/ //g' |perl -pe 's/ /\t/g'|perl -pe 's/;/\t/g' > Sec-20-NSS.sam.abund.txt1
mv Sec-20-NSS.sam.abund.txt1 Sec-20-NSS.sam.abund.txt
ktImportText *txt -o ../Piran_SortMeRNA.html
```

We can compare the the taxa abundances of the different samples in the 'Krona' plot in the 'html' file that is created in the last command.

### Assembly program

Now, we want to move on to doing the actual assembly! Here we will use the assembly program 'metaspades' and 'megahit' to start with. Spades especially is known to work well and perform better than the other assemblers in most of the recent publications in comparing the different assemblers!!

To do the assembly, we take all the sequences that is available. So, for each sample (library) there should be 4 files. Two files that still has R1 and R2 paired end reads. Then two or four single end reads (only R1 or R2) that lost one of the pair due to lower qulaity from trimming and/or the filtering step.

So, first we run the Spades assembly by using the metaspades option, since it is the version of spades that is to work better for Metganomic sequencing. Here we use all the read that is available!

``` bash
metaspades.py -k 21,33,55,77 -m 800 -o Piran_spades -t 32 --pe1-1 $dir'Sec-20-F1.paired_1.fastq.gz' --pe1-1 $dir'Sec-20-NSS.paired_1.fastq.gz' --pe1-2 $dir'Sec-20-F1.paired_2.fastq.gz' --pe1-2 $dir'Sec-20-NSS.paired_2.fastq.gz' --pe1-s $dir'Sec-20-F1.SE_1.fastq.gz' --pe1-s $dir'Sec-20-F1.SE_2.fastq.gz' --pe1-s $dir'Sec-20-NSS.SE_1.fastq.gz' --pe1-s $dir'Sec-20-NSS.SE_2.fastq.gz'
```

We also run 'megahit' assembly just so that it is better to combine different assemblies togther when you want to bin Genomes in the end!

``` bash
megahit --k-min 27 --k-max 77 --k-step 10 -m 100000000000 --mem-flag 2 -t 16 -o BF3_MegaH_out --out-prefix BF3_asm_MH -1 $dir'BF3-8L.paired_1.fastq.gz' -1 $dir'BF3-8N.paired_1.fastq.gz' -2 $dir'BF3-8L.paired_2.fastq.gz' -2 $dir'BF3-8N.paired_2.fastq.gz' -r $dir'BF3-8L.SE_1.fastq.gz' -r $dir'BF3-8L.SE_2.fastq.gz' -r $dir'BF3-8N.SE_1.fastq.gz' -r $dir'BF3-8N.SE_2.fastq.gz'
```

### Assembly quality

We can check the assembly quality by running aprogram called 'quast' which gives all the statistics about the assembly that has been made. We can also compare many different assemblies togther, so that we can compare the analysis. This is run by the following code.

``` bash
quast.py -o Quast_Alex_Spades_a_MH -s -L -t 8 /proj/Lokesh/Piran_a_Alex_samples/H5GHCBGX5_merged/Assembly/Alex/*/scaffolds.fasta /scratch/Lokesh_scratch/Assembly_MH/Alex_MH/*/*fa
```

Then from these assemblies, I would only use the scaffolds that were more 1000bp for further analyses which is he binning and so on. To do this I first make sure that the fasta files are one-line fasta files and then I run my own code to extract the sequences that are &gt;= 1000bp long

``` bash
cd /proj/Lokesh/Piran_a_Alex_samples/H5GHCBGX5_merged/Assembly/Alex/H2-2_spades
myown_oneLine.pl scaffolds.fasta
sel_seq_length.pl 1000 scaffolds.fasta.oneline.fasta H2-9_spades_1k.fasta
```

Files for Binning
-----------------

Now that we have the assemblies ready, it is time for binning genomes from metagenomes. So, in all the samples that we have in this dataset, we have done differential coverage treatment. This basically means that from each sample two different DNA extarctions were made and the assemblies were correspondingly by combining the reads per sample. You can check this at the assembly section.

So, the idea here is that we map each of the two different Sequence library to the assembly and we should get different coverage for the different scaffolds that could help us in binning our genomes of interest. You can read more about in the "Nature Biotechnology paper by Mads Albertsson on Differential Coverage Binning".

### Mapping

So, lets map the reads first to the corresponding assembly! Below is an example of one such sample. here 'bbmap' which is amapping tool that is written and it is faster and we can calculate the average coverage of each scaffold from it very easily.

``` bash
mkdir MH_map
mkdir Spades_map
module load bbmap
Alex_reads="/proj/Lokesh/Piran_a_Alex_samples/H5GHCBGX5_merged/Sequences/Step2_filtered/Alex_samples/"
Piran_reads="/proj/Lokesh/Piran_a_Alex_samples/H5GHCBGX5_merged/Sequences/Step2_filtered/Piran_samples/"
Alex_MH="/scratch/Lokesh_scratch/Assembly_MH/Alex_MH/"
Alex_SP="/proj/Lokesh/Piran_a_Alex_samples/H5GHCBGX5_merged/Assembly/Alex/"

bbmap.sh in1=$Alex_reads'H2-9K.paired_1.fastq.gz' in2=$Alex_reads'H2-9K.paired_2.fastq.gz' ref=$Alex_MH'H2-9_MegaH_out/H2-9_MH_1k.fasta' out=MH_map/H2-9K_MH.mapped.bam t=16 idfilter=90
bbmap.sh in1=$Alex_reads'H2-9N.paired_1.fastq.gz' in2=$Alex_reads'H2-9N.paired_2.fastq.gz' out=MH_map/H2-9N_MH.mapped.bam t=16 idfilter=90

bbmap.sh in1=$Alex_reads'H2-9K.paired_1.fastq.gz' in2=$Alex_reads'H2-9K.paired_2.fastq.gz' ref=$Alex_SP'H2-9_spades/H2-9_spades_1k.fasta' out=Spades_map/H2-9K_Spades.mapped.bam t=16 idfilter=90
bbmap.sh in1=$Alex_reads'H2-9N.paired_1.fastq.gz' in2=$Alex_reads'H2-9N.paired_2.fastq.gz' out=Spades_map/H2-9N_spades.mapped.bam t=16 idfilter=90
```

### Coverage

After the mapping of reads to the scaffolds, the next is to calculate the coverage of each of the scaffold! So, that we can map the coverage on the differential coverage plot. Then we use some of the scripts that came with the 'mmgenome' package to make all the necessary files to load it to the R along with the 'mmgenome' package.

``` bash
module load samtools
module load bedtools
#bbmap.sh in1=../../lys1.fq in2=../../lys2.fq ref=./Baja_scaf_1k.fa out=lys_scaf_1k.mapped.bam t=8 idfilter=90
#bbmap.sh in1=../../ut1.fq in2=../../ut2.fq out=ut_scaf_1k.mapped.bam t=8 idfilter=90
ls -1 *map|grep \:|tr -d ':'|while read line; do
    cd $line
    cd MH_map
    ls -1 *bam|while read bam_file; do
        samtools sort -@ 16 -T /tmp/$bam_file'sorted' -o $bam_file'.sorted.bam' $bam_file
          genomeCoverageBed -ibam $bam_file'.sorted.bam' -g ../*MH*genome.txt > $bam_file'.cov'
          bedTool_avg_mh.pl $bam_file'.cov'
    done
    cd ../Spades_map
    ls -1 *bam|while read bam_file; do
        samtools sort -@ 16 -T /tmp/$bam_file'sorted' -o $bam_file'.sorted.bam' $bam_file
        genomeCoverageBed -ibam $bam_file'.sorted.bam' -g ../*spades*genome.txt > $bam_file'.cov'
        bedTool_avg_sp.pl $bam_file'.cov'
    done
    cd ../../
done
```

So, basially what we do here is to take each of the 'bam' mapped assemblies for each of the samples, each of the treatments and each of the assembly and sort them first using the samtools. Then we take the sorted bam file and use the script from bedtools (genomeCoverageBed) with '-g' option that gives the length of each scaffold in a separate file will give us the coverage of each base of the assembly. then I have an in house perl-script (bedTool\_avg\_sp.pl) to calculate the average coverage of each of the scaffold that can be used directly to the 'mmgenome' pipeline in R.

### Essential genes and their taxonomy

Now, following is the code to create other input files for the 'mmgenome' analysis, which are the list of essential genes from each of the scaffold and the corresponding taxonomy of the scaffold! The taxonomy of the scaffolds here are basically done by predicting the single copy genes from the scaffolds. Then 'blasting' them against NR database as the Asgard genomes are not part of the RefSeq database.

Here I use the 'diamond' to blast them as it is faster as well and I do the MEGAn analysis manually instead of the commnd line! Because of the reason that there is a new file for mapping ncbi accession numbers to taxonomy and also it is more intuitive in the graphical user interface.

``` bash
ls -1 *map|grep \:|tr -d ':'|while read line; do
    cd $line
    ls -1 *fasta|while read ass; do
        prodigal -a temp.orfs.faa -i $ass -m -o temp.txt -p meta -q
        cut -f1 -d " " temp.orfs.faa > $ass'.orfs.faa'
        hmmsearch --cpu 16 --tblout $ass'.hmm.orfs.txt' --cut_tc --notextw ~/Diff_Cov_Bin_tools/mmgenome_R_tool/scripts/essential.hmm $ass'.orfs.faa' > hmm.temp.txt
        echo "scaffold orf hmm.id" > $ass'essential.txt'
        tail -n+4  $ass'.hmm.orfs.txt' | sed 's/ * / /g' | cut -f1,4 -d " " | sed 's/_/ /' >> $ass'essential.txt'
          grep -v "#" $ass'.hmm.orfs.txt' | cut -f1 -d " " > list.of.positive.orfs.txt
          perl ~/Diff_Cov_Bin_tools/mmgenome_R_tool/scripts/extract.using.header.list.pl -l list.of.positive.orfs.txt -s $ass'.orfs.faa' -o $ass'.orfs.hmm.faa'

#       blastp -query $ass.orfs.hmm.faa -db refseq_protein -evalue 1e-5 -num_threads 60 -max_target_seqs 5 -outfmt 5 -out $ass.orfs.hmm.blast.xml
        diamond blastp --query $ass'.orfs.hmm.faa' --db /proj/Lokesh/Databases/nr_db/nr -f 5 -e 0.00001 -o $ass'.orfs.hmm.diamond.nr.xml' --threads 16
#       MEGAN +g -x "load taxGIFile= ~/prot_acc2tax-Mar2018X1.abin; import blastfile= $ass'.orfs.hmm.diamond.nr.xml' meganfile=temp.rma;recompute toppercent=5;recompute minsupport=1;update;collapse rank=Species;update;select nodes=all;export what=CSV format=readname_taxonpath separator=tab file=$ass'.orfs.hmm.blast.tax.txt';update;close"
          perl ~/Diff_Cov_Bin_tools/mmgenome_R_tool/scripts/hmm.majority.vote.pl -i $ass'.orfs.hmm.blast.tax.txt' -o $ass'.tax.txt'
#   rm hmm.temp.txt
#   rm list.of.positive.orfs.txt
#   rm $ass'.orfs.hmm.blast.tax.txt'
#   rm temp.orfs.faa
#   rm temp.txt
#   rm temp.rma
# rm $ass.orfs.hmm.blast.xml
#   rm $ass'.orfs.hmm.faa'
#   rm $ass'.hmm.orfs.txt'
    done
    cd ..
done
```

So basically, MEGAN does the LCA analyses of each of the Signal copy gene that went through the diamond blast program and based on the matches it got. then we export the annotations for each of the leaf from MEGAN into a tab separated txt file. The file that comes out of the megan program is not exactly compatible to run the 'hmm.majority.vote.pl' perl script, beacuse of the scaffold names that we use. This perl script likes to have the header format 'scffoldNAME\_ORFid'. The majority vote is definied based on the different 'ORFid's for the same 'scffoldNAME'. But the names of scaffolds from both spades and megahit, they themselves have the '\_'s in their name. So we format the file first na duse the perl script and reformat it for our taste by the following.

``` bash
ls -1 */*tax.txt|while read filename; do
    bin_name=`basename $filename|cut -f 1 -d '.'`
    perl -pi -e 's/\"//g' $filename
    perl -pe 's/(.+)_(.+)\t/\1@\2\t/' $filename|perl -pe 's/_/-/g'|perl -pe 's/@/_/g' > $bin_name'_for_taxa.txt'
    perl ~/Files/Loki_MGNM_Analysis/Diff_Binning_tools/mmgenome_R_tool/scripts/hmm.majority.vote.pl -i $bin_name'_for_taxa.txt' -o $bin_name'.tax.txt'
    perl -pi -e 's/-/_/g' $bin_name'.tax.txt'
done
```

Important Note: We should check if the format of the 'essential genes' text file is also fine with the '\_'s and spaces.

### Quick check on the abundance of our favourite taxa groups

Now we have coverage and taxonomy!! So we can get a quick peek on the relative abundance of our taxa based on the number of reads mapped to the scaffolds that got the taxonomic annotation of the interested groups.

``` bash
# first we only get the headers of the scaffolds that has taxa annotation
tail -n+2 BF3-8_spades_1k.tax.txt|cut -f 1 > H2-2_spades_all_taxa.txt
# Now, the intersting phyla
grep 'TACK' H2-2_spades_1k_TACK.tax.txt|cut -f 1 > H2-2_spades_TACK.txt
# Then further to class (thaumarchaeota)
grep 'Thaum' H2-2_spades_1k_TACK.tax.txt|cut -f 1 > H2-2_spades_Thaum.txt

# get the exact number of reads that mapped to the assembly using samtools
samtools index H2-2K_Spades.mapped.bam.sorted.bam
samtools idxstats H2-2K_Spades.mapped.bam.sorted.bam > H2-2K_Spades.mapped.bam.sorted.bam.idxstats

# Now we use the "grep command perl script" to fish out 'idxstats' info for our interested scaffold 
grep_f_long.pl H2-2_spades_all_taxa.txt ./H2-2K_Spades.mapped.bam.sorted.bam.idxstats H2-2_spades_all_taxa_coverage.txt

#this should give us total number of reads
cut -f 3-4 ./H2-2K_Spades.mapped.bam.sorted.bam.idxstats |perl -pe 's/\s/\+/'|perl -pe 's/\n/\+/'|perl -pe 's/\+$/\n/'|bc

# Then number of reads for the particular cases
cut -f 3-4 H2-2_spades_all_taxa_coverage.txt |perl -pe 's/\s/\+/'|perl -pe 's/\n/\+/'|perl -pe 's/\+$/\n/'|bc
cut -f 3-4 H2-2_spades_TACK_coverage.txt |perl -pe 's/\s/\+/'|perl -pe 's/\n/\+/'|perl -pe 's/\+$/\n/'|bc
cut -f 3-4 H2-2_spades_Thaum_coverage.txt |perl -pe 's/\s/\+/'|perl -pe 's/\n/\+/'|perl -pe 's/\+$/\n/'|bc
```

The ratio of 'TACK' to the 'all\_taxa' will give us the fraction of TACK among the annotations. Keep in mind that this doesn't necessarily reflect the actual relative abundance!! In best cases only about 5% of total reads mapped to scaffolds that had taxonomy annotations based on the Single Copy Genes. This is understable, as not all scaffolds will have SCGs and the primary reason we do this is to find the area in the Differential coverage bi-plot that has our favourite organisms.

Binning
-------

Now it is time for the actual binning! the following code is based on the 'mmgenome' tutorial that is available online.

Here we load the libararies

``` r
library(vegan)
library(plyr)
library(RColorBrewer)
library(alphahull)
library(ggplot2)
library(mmgenome)
```

Here we extract the scaffolds that fall in the region of our interest and we will define the exact location where we would like to get all the scaffolds from.

### Coverage bi-plot

Here we extract the scaffolds using the 'mmgenome' R package and we look at the statistics of the particular scaffolds we extarct

``` r
setwd("~/Files/Alex_MetaGenomes_22_03_2018/Assembly/BF3/")
assembly <- readDNAStringSet("BF3-8_spades_1k.fasta", format = "fasta")
Lysozyme <- read.table("BF3-8L_Spades.mapped.bam.cov_avg.cov", header = T,sep = "\t")  
NS_kit <- read.table("BF3-8N_spades.mapped.bam.cov_avg.cov", header = T,sep = "\t") 
Enz_lys <- read.table("BF3-8L-Enz_Spades.mapped.bam.cov_avg.cov", header = T,sep = "\t") 
Enz_NS <- read.table("BF3-8N-Enz_Spades.mapped.bam.cov_avg.cov", header = T,sep = "\t") 
ess <- read.delim("BF3-8_spades_1k.fastaessential.txt", header = T,sep=" ") 
#colnames(ess) = c("scaffold","orf","hmm.id")
tax <- read.delim("BF3-8_spades_1k.tax.txt", header = T,sep="\t") 
TACK_tax <- read.delim("BF3-8_spades_1k_TACK.tax.txt", header = T,sep="\t") 

m <- mmload(assembly = assembly, 
            pca = T,
            coverage = c("Lysozyme", "NS_kit","Enz_lys","Enz_NS"), 
            essential = ess,
            tax = tax,
            tax.expand = "TACK group")

p <- mmplot(data = m, x = "Lysozyme", y = "NS_kit", log.x = T, log.y = T, color = "essential", minlength = 1000)
p
```

### Visualization of TACK only

``` r
n <- mmload(assembly = assembly, 
            pca = T,
            coverage = c("Lysozyme", "NS_kit"), 
            essential = ess,
            tax = TACK_tax,
            tax.expand = "TACK group")

q <- mmplot(data = n, x = "Lysozyme", y = "NS_kit", log.x = T, log.y = T, color = "essential", minlength = 1000)
q
```

Shearing methods for library prep
---------------------------------

Since we had two different kinds of library preparation, one with usual random shearing using covarice and the other with enzymatic shearing. We would like to look at how the shearing method influenced the taxonomy in these samples

``` r
mmplot(data = m, x = "Lysozyme", y = "Enz_lys", log.x = T, log.y = T, color = "essential", minlength = 1000)
mmplot(data = m, x = "NS_kit", y = "Enz_NS", log.x = T, log.y = T, color = "essential", minlength = 1000)
```

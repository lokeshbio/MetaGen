-   [Sequence processing](#sequence-processing)
    -   [Obtaining reads](#obtaining-reads)
    -   [Quality check](#quality-check)
    -   [Quality trim, filter and check (again)](#quality-trim-filter-and-check-again)
-   [Assembly](#assembly)
    -   [SortMeRNA](#sortmerna)
    -   [Assembly program](#assembly-program)
    -   [Assembly quality](#assembly-quality)
-   [Binning](#binning)
    -   [Mapping](#mapping)

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

Binning
-------

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

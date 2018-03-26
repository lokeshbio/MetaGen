-   [Sequence processing](#sequence-processing)
    -   [Obtaining reads](#obtaining-reads)
    -   [Quality check](#quality-check)
    -   [Quality trim, filter and check (again)](#quality-trim-filter-and-check-again)

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

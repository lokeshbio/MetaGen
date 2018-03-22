-   [Phylogeny](#phylogeny)
    -   [Protein sequences](#protein-sequences)
    -   [Clustering](#clustering)
    -   [Tree](#tree)

Phylogeny
---------

This is an example of how u do the phylogenetic anaysis of certain proteins. So, this page basically contains information starting from obtaining the protein sequence information from the public dtabases and the protein sequences from the Metagenome assembled genomes. Then, clustering them and calculating the tree based on Maximum Liklihood from IQ-tree

### Protein sequences

All the analyses in this section are mainly based on the proteins from the publication [Villaneuva et.al 2016](https://onlinelibrary-wiley-com.uaccess.univie.ac.at/doi/full/10.1111/1462-2920.13361)! We obtained protein sequences based on the accession numbers that are available from this publication using the [NCBI Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) application.

In this particular case, we were interested in the phylogeny of the enzymes coding for Glceraldehyde-3-Phosphate dehydrogenase (G3PDH), so we looked for the protein sequences in our local genomes that were annoated for these enzymes. The archael enzymes were download from the NCBI page separately, as the accession numbers were available for them separately from the other enzymes used in the publication [Villaneuva et.al 2016](https://onlinelibrary-wiley-com.uaccess.univie.ac.at/doi/full/10.1111/1462-2920.13361)! So to see it clearly in the phylogenetic tree, I included 'ARCHAEA' to their headers and then all the sequences are combined into a one single file.

``` bash
perl -pe 's/\>/\>ARCHAEA_/g' Gps.fasta > Gps_e.fasta
cat glp.faa gps.faa ../gene_lists/Baja_G3PDH.faa > G3PDH_all.faa
```

### Clustering

Now, we will cluster all these protein sequences using the program [MAFFT](https://mafft.cbrc.jp/alignment/software/manual/manual.html)! In this we do a local alignment with all the sequences together! My primary aim is to see where does the sequences from our local genomes fall into this phylogenetic tree!? We can also check how the program is running (how many threads is it using, how memory is it consuming) through specialized command called 'htop'

``` bash
mafft-linsi --thread 6 --maxiterate 1000 G3PDH_all.faa > G3PDH_all.mafft.faa
htop
```

Now the file 'G3PDH\_all.mafft.faa' contain the alignment of all the sequences with too many gaps! you can use the tool [AliView](http://ormbunkar.se/aliview/) to see how the alignment looks like! There are two different "alignment trimming" tools, [BMGE](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210) and [TrimAl](http://trimal.cgenomics.org/). These tools help in trimming the alignment to pick only columns that are of evolutionary/phylogenetic significance!! In this way, the tree more meaningful as it will be calculated based on the conserved regions of the enzymes rather than the regions of non-significance!

``` bash
java -jar ~/bin/BMGE-1.12/BMGE.jar -i G3PDH_all.mafft.faa -t AA -of G3PDH_all.mafft.BMGE.faa
trimal -in G3PDH_all.mafft.faa -out G3PDH_all.mafft.gappy.faa -gappyout
trimal -in G3PDH_all.mafft.faa -out G3PDH_all.mafft.automated.faa -automated1
```

### Tree

After trimming the alignment, it is time to build the phylogenetic tree with the trimmed alignemnt, For this, I use [IQTREE](http://www.iqtree.org/)! The good thing about the 'iqtree' is that you can ask the program to test for the different evolutionary models for phylogeny and the program will choose the best model for the alignment! Along with this, we can also add the 'ultrafast-bootstrap' and the 'SH-alrt' values to the branches and their nodes.

``` bash
iqtree-omp -s G3PDH_all.mafft.BMGE.faa -st AA -nt AUTO -quiet -safe -bb 1000 -alrt 1000 -m TEST
iqtree-omp -s G3PDH_all.mafft.gappy.faa -st AA -nt AUTO -quiet -safe -bb 1000 -alrt 1000 -m TEST
iqtree-omp -s G3PDH_all.mafft.automated.faa -st AA -nt AUTO -quiet -safe -bb 1000 -alrt 1000 -m TEST
```

Then, we use [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) for visualization of the tree! Then we can make a tree like the following using the FigTree program!

Inline-style: 
![alt text](https://github.com/lokeshbio/MetaGen/blob/master/GlpK_tree.nexus.jpg "Test-tree")

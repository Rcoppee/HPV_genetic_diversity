# Scripts used for the study "Cervical and anal HPV whole-genome sequencing reveals high variations in APOBEC3-induced mutations among HPV risk categories in a Togolese key-population"
<p>This repository contains all the scripts that were developped to produce the work "Cervical and anal HPV whole-genome sequencing reveals high variations in APOBEC3-induced mutations among HPV risk categories in a Togolese key-population". You can contact the author if you need some helps to use the scripts for your own projects.<br>
 Large metadata can be shared on demand. In the Data/ subdirectory, you will find the GFF annotation file for most HPV types, the list of HPVs (and corresponding GenBank accession numbers) detected in our samples, the list of HPV reference sequences, and the multiple sequence alignments using sequences retrieved on GenBank for HPV types 6, 11, 16 and 18.</p>
 <br>
<h3>1. Prerequisites</h3>
 <p>To use these scripts, you need the following programs:</p>
 <p>
 - R version 4.0.4<br>
 - List of R packages: data.table, dplyr, ggplot2, stringr</p>
 <br>
 <h3>2. Preparing the data</h3>
 <p>Before studying genetic diversity of HPV types, we first applied a common bioinformatics pipeline to produce pileup files (which contain the exhaustive list of bases observed for each position along the genome).</p>
  <p><code> bwa mem hpv_ref_genomes.fasta sample.R1.fastq.gz sample.R2.fastq.gz > sample.sam</code></p>
 <p><code> samtools view -b -S sample.sam > sample.bam</code></p>
 <p><code> samtools sort sample.bam -o sample.sorted.bam </code></p>
 <p><code> samtools index sample.sorted.bam</code></p>
  <p><code> samtools mpileup -a -f reference_genome.fasta file_sorted.bam > file_sorted.pileup</code></p>
  <p>where <code>sample</code> is the name of the sample.
<br>
 For <code>samtools mpileup</code>, <code>-d 1000000</code> indicates that we look a maximum of 1 million of reads for each position of the genome.</p>
 <br>
 <h3>3. Genetic diversity of HPV types</h3>
 <p>xxx</p>
 <br>
  <h3>4. APOBEC3-induced mutations in HPV types at the genome and gene levels using our samples</h3>
 <p>xxx</p>
 <br>
  <h3>5. APOBEC3-induced mutations in HPV types at the genome and gene levels from GenBank sequences</h3>
 <p>xxx</p>
 <br>
 <h3>6. Citation</h3>
 <p>If you use or adapt some scripts for your own work, please cite:</p>
 <p><i>In preparation.</i></p>

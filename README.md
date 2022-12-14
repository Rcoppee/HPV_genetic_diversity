# Scripts used for the study "Cervical and anal HPV whole-genome sequencing reveals high variations in APOBEC3-induced mutations among HPV risk categories in a Togolese key-population"
<p>This repository contains all the scripts that were developed to produce the work "Cervical and anal HPV whole-genome sequencing reveals high variations in APOBEC3-induced mutations among HPV risk categories in a Togolese key-population". You can contact the author if you need some helps to execute or adapt the scripts for your own projects.<br>
 Large metadata can be shared on demand. In the data/ subdirectory, you will find the GFF annotation file for most HPV types, the list of HPVs (and corresponding GenBank accession numbers) detected in our samples, the list of HPV reference sequences, and the multiple sequence alignments from GenBank sequences for HPV types 6, 11, 16 and 18.</p>
 <br>
<h3>1. Prerequisites</h3>
 <p>To use the scripts, you need the following programs/packages:</p>
 <p>
 - R version 4.0.4<br>
 - List of R packages: data.table, dplyr, ggplot2, stringr</p>
 <br>
 <h3>2. Preparing the data</h3>
 <p>Before studying the genetic diversity of HPV types, we first applied a common bioinformatics pipeline to produce pileup files (which contain the exhaustive list of bases observed for each position along the genome).</p>
  <p><code> bwa mem hpv_ref_genomes.fasta sample.R1.fastq.gz sample.R2.fastq.gz > sample.sam</code></p>
 <p><code> samtools view -b -S sample.sam > sample.bam</code></p>
 <p><code> samtools sort sample.bam -o sample.sorted.bam </code></p>
 <p><code> samtools index sample.sorted.bam</code></p>
  <p><code> samtools mpileup -a -f reference_genome.fasta file_sorted.bam > file_sorted.pileup</code></p>
  <p>where <code>sample</code> is the name of the sample.
<br>
 For <code>samtools mpileup</code>, <code>-d 1000000</code> indicates that we consider a maximum of one million reads for each position of the genome.</p>
 <br>
 <h3>3. Genetic diversity of HPV types</h3>
 <p>To explore the genetic diversity of HPV types and per risk category, you developed the script <i><b>genetic_diversity_study.R</b></i>. The script uses as inputs some pileup files and both the <i><b>HPV_ref_type.txt</b></i> and <i><b>reference_observed.txt</b></i> files. The analysis takes from a few minutes to hours, depending on the number of samples investigated). When the analysis is finished, two plots are generated:
<br>- the total number of mutations per HPV type
<br>- the total number of mutations per risk category</p>
 <p> In a second step, the algorithm reads the <i><b> HPV_annotation.gff</b></i> file to perform the analysis at the gene level. After a few minutes, two other plots are produced:
<br>- the ratio(c>t) per gene and per risk category
<br>- a focus of the ratio(c>t) on the E6 gene for low and high risk categories</p>
 <br>
  <h3>4. APOBEC3-induced mutations in HPV types at the genome and gene levels using our samples</h3>
 <p>To explore APOBEC3-induced mutations among HPV types and risk categories using our samples, we developed the scripts <i><b>apobec_genome_study.R</b></i> (genome level) and <i><b>apobec_gene_study.R</b></i> (gene level).</p>
 <p>At the genome level, the pipeline is very similar to the one that evaluates the genetic diversity. The major difference is that we focus on C>T mutations, both in TCW and non-TCW motifs. Systematically, we calculate the ratio(c>t) to identify genome sequences that are enriched in APOBEC3-induced mutations rather than random mutations (we suppose here a ratio(c>t) >= 2 as significant). Two plots are produced after a few minutes:
<br>- the ratio(c>t) according to each HPV type
<br>- the ratio(c>t) according to the risk category</p>
 <p>At the gene level, a part of the script must be executed for each gene investigated (in the deposited script, we focus on the E1 gene as an example). The pipeline is then the same than the one used at the genome level. Once the analysis is done for each gene (E1, E2, E4, E6, E7, L1 and L2), all the results are combined and a plot that shows the distribution of samples with APOBEC3-induced mutations among genes is produced.</p>
 <br>
  <h3>5. APOBEC3-induced mutations in HPV types at the genome and gene levels from GenBank sequences</h3>
 <p>To explore APOBEC3-induced mutations among HPV types 6, 11, 16 and 18, we developed the scripts <i><b>apobec_genome_genbank.R</b></i> (genome level) and <i><b>apobec_gene_genbank.R</b></i> (gene level).</p>
 <p>At the genone level, the script uses one multiple sequence alignment as input (in the deposited script, we focus on HPV 6 as an example, i.e. the <b><i>hpv6_aligned_clean.fasta</i></b> file that can be found in the data/ subdirectory). The pipeline counts the number of C>T mutations in the TCW motif. Once the analysis is done for each type investigated (here, HPV types 6, 11, 16 and 18), 100 replicates comparing the number of APOBEC3-induced mutations among a random draw of 50 HPV sequences per type were performed. The script calculates the number of replicates where a signicant difference between two types (in the deposited script, HPVs 6 vs. 16 as an example) is observed. To illustrate the result, a plot based on the last replicate showing the number of APOBEC3-induced mutations for each HPV type is produced.</p>
 <p>At the gene level, the strategy is similar to the one used at the genome level. A part of the script must be executed for each gene investigated (in the deposited script, we focus on the E1 and E2 genes as an example). The pipeline is then the same than the one used at the genome level. The algorithm counts, for a given HPV type, the proportion of sequences on a gene that has a ratio(c>t) >= 2, suggesting that the gene is enriched in APOBEC3-induced mutations.</p>
 <br>
 <h3>6. Citation</h3>
 <p>If you use or adapt some scripts for your own work, please cite:</p>
 <p><i>In preparation.</i></p>

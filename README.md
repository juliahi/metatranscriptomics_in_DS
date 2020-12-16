# metatranscriptomics_in_DS
## Metatranscriptomics study of the Down’s Syndrome model mice fed high-fat diet

..... 

### Quantification of total-RNA sequences

We run the seed-kraken tool with bash commands:

>for f in /path_to_fasta_files/*_fwd.fasta; do echo $f output/`basename $f`.report; kraken --thread $THREADS --fasta-input --preload -db $KRAKEN_DB_NAME --output /output_folder/`basename $f`.report $f; done

and then: 

>for f in /output_folder/*report; do echo `basename $f`; (kraken-report --db $KRAKEN_DB_NAME $f  >$f.csv); done

The database for seed-kraken was built based on a collection of sequences from RefSeq (Feb. 2016): Bacteria, Viruses, Fungi, H. sapiens and M. musculus. Sequences needed to have version_status=='latest' and assembly_level = "Complete Genome" or assembly_level = "Chromosome".

### Defining differentially expressed genes identified in contigs
R, package DeSeq2 from Bioconductor:

dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)    
dds = DESeq2::estimateSizeFactors( dds ) 
ddsp <- DESeq2::estimateDispersions( dds, fitType="parametric" )
ddsp = DESeq2::nbinomWaldTest( ddsp )    
result <- results(dds)

### Assessing gene expression levels from raw reads 

We run the MetaGeneMark tool with bash commands:

folder_with_MetaGeneMark/mgm/gmhmmp -a -d -m 
folder_with_MetaGeneMark/mgm/MetaGeneMark_v1.mod input_file.fasta -o output_file_mgm

Then for the predicted genes we found orthologs using eggNOGmapper with the default settings using a command-line, also because files were too big to be processed with the online tool. The command used was: 

python folder_with_eggNOGmapper/eggnog-mapper-1.0.3/emapper.py -i 
meta_genemark_predicted_sequences/input_file.faa --output output_file_eggNOG --cpu 10 -m diamond 



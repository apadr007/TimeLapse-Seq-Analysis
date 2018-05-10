# TimeLapse-seq: young vs old RNA

An outline of my code used to analyze Timelapse-seq data

    # download data from ftp server
    mkdir timelapse-seq
    cd timelapse-seq
    wget -m --user=gslftp --password=gsl23ftp!ftp://gslserver.qb3.berkeley.edu/180427_50SR_HS4KA/*
    cd gslserver.qb3.berkeley.edu/L2345678/
    mv Ingolia-RSLiB/ ~/timelapse-seq/
    cd ~/timelapse-seq/
    rm -r gslserver.qb3.berkeley.edu/
    mv Ingolia-RSLiB/ RawFastq
    
    # removing adaptor sequences from library, then splitting by barcode
    cd ~/timelapse-seq
    ./adapter_removal.sh
    
    # 
    cd ~/timelapse-seq/Alignment
    ln -s ~/timelapse-seq/RawFastq/*_clipped.fastq .
    

## Stats can be found here from the Timelapse-seq analysis

[https://docs.google.com/spreadsheets/d/11tny8mEw6agXvnJTfhhDyKMVorB7_Ntgmnz_3YwLaXM/edit#gid=0](https://docs.google.com/spreadsheets/d/11tny8mEw6agXvnJTfhhDyKMVorB7_Ntgmnz_3YwLaXM/edit#gid=0)

## Scripts used to trim adapters and map reads

- **Adapter_removal.sh**

  I made a script to run more easily going forward called `adapter_removal.sh`. I'm placing the Ribosome Profiling libraries and the RNAseq libraries into bash arrays to loop over using different parameters at the `fastx_clipper` (-c) stage. 

      #!/bin/bash
      
      # prep dir
      mkdir -p RawFastq
      mkdir -p Alignment
      mkdir -p Analysis
      
      
      # enter directory name of Ribosome Profilng & RNAseq samples
      cd ./RawFastq
      RPSEQ=(RP*.gz)
      RNASEQ=(RS*.gz)
      DIR='/mnt/ingolialab/apadron/timelapse-seq/RawFastq/'
      ILLUMINA='AGATCGGAAGAGCACACGTCTGAA'
      
      for i in "${RNASEQ[@]}"; do
      			nohup zcat $i | fastx_clipper -Q33 -a ${ILLUMINA} \
      							-v -o ${i}_clipped.fastq &
      done
      
      for i in "${RPSEQ[@]}"; do
      	nohup zcat $i | fastq_illumina_filter --keep N -v | \
      	fastx_clipper -Q33 -a ${ILLUMINA} \
      	-c -v -o ${i}_clipped.fastq &
      done 
      
      # stats from removing barcode
      for i in *_clipped.fastq; do
      		zcat $i | wc -l
      done
      
      wc -l ${i}_clipped.fastq

- **RP_barcode_split.sh**

  This requires having a **pool.csv** file with the names of each file and the corresponding barcode sequences. 

  [pool.csv](https://www.notion.so/file/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2F7f147045-bd68-46fe-8d14-7db4fa2d351e%2Fpool.csv)

      #!/bin/bash
      
      cd Alignment
      RPDIRNAME='RP'
      RPDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/'
      
      RPLIB=(RP*_clipped.fastq)
      rRNADIR='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/rRNARef/'
      
      # splitting reads along with the barcode
      for i in "${RPLIB[@]}"; do
      nohup fastx-split -o ${RPDIRNAME} -p NN -x NNNNNIIIII --min-insert=10 -s pool.csv ${i} &
      done
      
      # bowties index
      #cd /mnt/ingolialab/apadron/Genomes/Homo_sapiens/rRNARef/
      #bowtie2-build -f RefSeqrRNA.fa RefSeqrRNA
      
      # stats on RP barcode splitting
      cd ${RPDIR}${RPDIRNAME}
      cat fates.txt
      
      
      
      # align reads to human rRNA database and capture un-aligned reads
      cd ${RPDIR}${RPDIRNAME}
      for i in *.fastq; do
      bowtie2 -p 2 --very-sensitive --quiet --un ${i}.norrna.fq -x ${rRNADIR}RefSeqrRNA -U ${i} | rrna-stats -o ${i}.stats --tam --maxread 100 --lenrange 5,100 - &
      done

- **~~map_reads.sh~~**

      #!/bin/bash
      
      MAINDIR='/mnt/ingolialab/apadron/timelapse-seq'
      cd ${MAINDIR}
      INDEX='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/GRCh38/hisat/GRCh38'
      # define RNAseq dir and put files into array
      RNASEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment'
      RNASEQ=(`find ${RNASEQDIR} -name RS*_clipped.fastq`)
      
      # define RPseq dir and put files into array
      RPSEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/RP'
      RPSEQ=(`find ${RPSEQDIR} -name *.norrna.fq`)
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RNASEQ[@]}"; do 
      hisat2 -p 32 --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      # Convert SAM to BAM format
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RPSEQ[@]}"; do 
      hisat2 -p 32 --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      # Sort BAM files
      BAMFILES=(`find ./Alignment/* -name *.bam`)
      for i in "${BAMFILES[@]}"; do
      	samtools sort $i -o ${i}.sorted.bam &
      done

- **testing_alignment_param.sh**

  **Note: This is the code I proceeded**

  This script is used to test the `--ignore-quals` option in Hisat2, since the frequency of G>A mutations was low in some of the test samples

      #!/bin/bash
      
      MAINDIR='/mnt/ingolialab/apadron/timelapse-seq'
      cd ${MAINDIR}
      INDEX='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/GRCh38/hisat/GRCh38'
      # define RNAseq dir and put files into array
      RNASEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/testing_other_alignment_params/'
      RNASEQ=(`find ${RNASEQDIR} -name RS*_clipped.fastq`)
      
      # define RPseq dir and put files into array
      RPSEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/testing_other_alignment_params'
      RPSEQ=(`find ${RPSEQDIR} -name *.norrna.fq`)
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RNASEQ[@]}"; do 
      hisat2 -p 32 --ignore-quals --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      # Convert SAM to BAM format
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RPSEQ[@]}"; do 
      hisat2 -p 32 --ignore-quals --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      # Sort BAM files
      BAMFILES=(`find ./Alignment/testing_other_alignment_params/* -name *.bam`)
      for i in "${BAMFILES[@]}"; do
      	samtools sort $i -o ${i}.sorted.bam &
      done

## Finding `A>G` mutations for Timelapse-seq data

After `map_reads.sh` I will run the BAM splitting program, which counts `A>G` mutations. I downloaded his program from his Github in the "synth-data" directory. I'll run it using `python err-stats.py myfile.bam`

- `err-stats.py`

  It will probably take a while to run on your real, large data sets, and so you may want to run it on each of your data sets in parallel. It generates 4 output files:

  1. `myfile-nerr.txt` will indicate how many reads have 0, 1, 2, etc. mismatches against the genomic reference.
  2. `myfile-nalign.txt` will count, for each nucleotide position in the read, how often it encounters each pair of reference nucleotide and read nucleotide. I split these up according to position in the read because the first few bases have more errors.
  3. `myfile-ferr.txt` will indicate how often each kind of change occurs to a certain reference base. That is, the 'AtoC' column is the fraction of reference As that are read as Cs. Like the "nalign.txt" file, these are split up according to read position.
  4. `myfile-ferr-avg.txt` is like the "ferr.txt" file, averaged over bases 6 - 45 (i.e., excluding the ends).

    
    # first I'll make a soft link to all of my bam files into a new dir
    
    cd ~/timelapse-seq/Analysis
    ln -s ~/timelapse-seq/Alignment/*.sorted.bam ./
    ln -s ~/timelapse-seq/Alignment/RP/*.sorted.bam ./
    
    # run err-stats.py program on all of my sorted BAM files 
    for i in *.sorted.bam; do
    `python ~/timelapse-seq/[err-stats](http://err-stats.py/)[.[py](http://err-stats.py/) ${i} &
    done`

## Split BAM file reads based on `A>G` using `split-tl.py`

It will also take a while to run. It generates three output files:

1. `myfile-yestl.bam` contains all reads with TimeLapse mutations (`A>G` mutations on these first-strand cDNA reads).
2. `myfile-notl.bam` contains all reads that don't have TimeLapse mutations
3. `myfile-tl-stats.txt` contains statistics for the different categories of reads.

    #!/bin/bash
    
    # run split-tl.py
    SPLITDIR='/mnt/ingolialab/ingolia/Prog/ribosome-profiling/timelapse/'
    DIR='/mnt/ingolialab/apadron/timelapse-seq/Analysis/'
    
    cd ~/timelapse-seq/Analysis
    
    # for RNAseq samples I'll need to split files by A>G mutations
    for i in RS_*.sorted.bam; do
            python ${SPLITDIR}split-tl.py --mutation A,G ${i} &
    done
    
    # for Ribosome Profiling samples I'll need to split files by T>C mutations
    for i in *.fq.bam.sorted.bam; do
            python ${SPLITDIR}split-tl.py --mutation T,C ${i} &
    done

## Run FeatureCounts on split BAM files and full BAM files

    #!/bin/bash
    
    # run FeatureCounts on ALL BAM files
    cd ~/timelapse-seq/Analysis
    
    FEATURECOUNTSDIR='/mnt/ingolialab/apadron/subread-1.5.3-Linux-x86_64/bin/'
    GTFDIR='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/GRCh38/'
    
    for i in *.bam ; do
    ${FEATURECOUNTSDIR}featureCounts -t exon -M --fraction -g gene_id \
    -a ${GTFDIR}gencode.gtf -o ${i}.featureCounts ${i} &
    done

## Sample replication assessment

    library(data.table)
    require(scales)
    setwd('~/timelapse-seq/Analysis/')
    myFiles <- list.files('.')
    
    featureCounts = list()
    output_list   = list()
    for (i in 1:length(myFiles)){
      featureCounts[[i]] <- fread(myFiles[i], header=T)
      output_list[[i]] <- cbind(featureCounts[[i]]$Geneid, 
                                featureCounts[[i]][,7])
    }
    
    # save count data
    write.csv(samples, '~/timelapse-seq/Ranalysis/gene_counts.csv', quote = F)
    
    length(output_list)
    head(samples)
    
    # calculate pearson correlation
    row.names(samples) <- samples$GeneID
    samples$GeneID <- NULL
    cor(samples)
    
    # make correlation plots
    pdf('~/ZaraZandu/files/correlation_plots.pdf')
    par(pin=c(2,2), tcl=0.25, ps=9, family="Helvetica")
    plot(log10(samples$sample_A), log10(samples$sample_B), cex=0.3, pch=19,
         xlab=c('Sample A [Log Counts]'), ylab=c('Sample B [Log Counts]'),
         col = alpha(1, alpha = 0.1), asp=1)
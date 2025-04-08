
# source package 638df626-d658-40aa-80e5-14a275b7464b # latest samtools package.

# version_root=/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/results/PBMC_MASseq_YL/nf-single_tx_geneBodyCoverage/paper_scripts/
sample=illumina_sample1B
analysis_dir=`pwd`
symlink_dir=./data
source gffread-0.12.7_CBG
source Samtools-1.15.1
source rseqc-3.0.1

# ==== 2. create bam file to used in rseqc =====
if true; then
    cd $analysis_dir/output/
    # list of BC post cellranger filter
    zcat $symlink_dir/filtered_feature_bc_matrix/barcodes.tsv.gz|sort|uniq > $sample.correctedBC.list
    wc -l output/illumina_sample1B.correctedBC.list
    # 4875 output/illumina_sample1B.correctedBC.list

    # extract reads with BC and xf:25: match the numbe in filtered_feature_bc_matrix.  
    # see https://www.10xgenomics.com/analysis-guides/tutorial-navigating-10x-barcoded-bam-files 
    cd $analysis_dir/output && source package 638df626-d658-40aa-80e5-14a275b7464b && \
        samtools view -@4 -d xf:25 $symlink_dir/possorted_genome_bam.bam -b | \
        samtools view -@4 -D CB:$sample.correctedBC.list -b -o $sample.possorted_genome_bam.correctedBC_xf25.bam  && \
        samtools index -@4 $sample.possorted_genome_bam.correctedBC_xf25.bam && \
        samtools flagstat -@4 $sample.possorted_genome_bam.correctedBC_xf25.bam > $sample.possorted_genome_bam.correctedBC_xf25.bam.flagstat
fi 

# ==== 3. make CDS only bed file for each transcript range ====
if true; then
    # 3.1: make CDS only gff3
    # genecode v39 annotation used by pacbio smrtlink
    gff3=$symlink_dir/gencode.v39.annotation.coding.gff3
    # take off exon entries, AND REMOVE COLUMN 8 NUMBER WITH .
    awk '$3!~/exon/' $gff3 |awk -v OFS="\t" '{if ($3=="CDS"){$8="."} print $0 }'| sed -r 's/\sCDS\s/\texon\t/' > $analysis_dir/output/gencode.v39.annotation.coding.CDS.gff3
    gffread $analysis_dir/output/gencode.v39.annotation.coding.CDS.gff3 --bed -o $analysis_dir/output/gencode.v39.annotation.coding.CDS.bed 

    # 3.2: make bed file per transcript range
    for length_range in "0-500" "501-1000" "1001-2000" "2001-3000" "3001-1000000000"; do
        range1=$(echo $length_range | cut -d '-' -f1)
        range2=$(echo $length_range | cut -d '-' -f2)
        
        awk -v range1="$range1" -v range2="$range2" \
            '{split($11,len,",");genelen=0;for (i in len){genelen=genelen+len[i]}; if(genelen<=range2 && genelen>range1){print $0 }}' \
            $analysis_dir/output/gencode.v39.annotation.coding.CDS.bed | sort -k 1,1 -k2,2n > $analysis_dir/output/trans.$length_range.sorted.bed
        echo $length_range done
    done
fi

# ==== 4. run rseqc ====
if true; then
    for length_range in "0-500" "501-1000" "1001-2000" "2001-3000" "3001-1000000000"; do
        if [ ! -d $analysis_dir/output/geneBodyCoverage_output.$length_range ]; then mkdir $analysis_dir/output/geneBodyCoverage_output.$length_range; fi
        cd $analysis_dir/output/geneBodyCoverage_output.$length_range
        cd $analysis_dir/output/geneBodyCoverage_output.$length_range && \
            singularity exec /ei/software/cb/containers/eimap/x86_64/Singularity.img geneBody_coverage.py \
            -i $analysis_dir/output/$sample.possorted_genome_bam.correctedBC_xf25.bam \
            -r $analysis_dir/output/trans.$length_range.sorted.bed -o $sample.CDS_only && \
        sleep 30
    done
fi 


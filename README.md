# Module xpRNApy - X-linked gene eXPRession RNAseq Analysis in PYTHON

Gene expression analysis of RNASeq data sets using genome annotation GTF (genes.gtf) and HTSeq-count read count data *_counts files.

Find and set the two path variables:

`DATAPATH` needs to be set to a path for loading analysis data from and writing output to.
`CHIP_DATA_PATH` needs to be set to a folder containing H3K27me3 ChIPseq data for performing correlation analysis with gene expression.

To load the precomputed analysis data find the following line in the script and set LOAD_ANALYSIS:

`LOAD_ANALYSIS = True
`
The following function will load an existing analysis:

`analysisProject = XpressionAnalysis.from_csv("data.txt")`

`analysisProject = newXAnalysis(DATAPATH="../RNAseq/")`

This version of xpRNApy is based on Python 2.7 and was used with the following packages for performing the analysis:

    pandas.__version__
    u'0.24.2'
    numpy.__version__
    '1.16.6'
    seaborn.__version__
    '0.9.1'
    matplotlib.__version__
    '2.2.5'
    matplotlib_venn.__version__
    '0.11.10'
    scipy.__version__
    '1.2.3'
    HTSeq.__version__
    '0.11.1'


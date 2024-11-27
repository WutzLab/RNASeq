import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import matplotlib.gridspec as gridspec
import os

def import_DSeq2_csv(filename):
    """
    Imports a csv file from DSeq2 differential gene expression analysis.
    The csv file path is provided in filename. Returns a pandas DataFrame.
    """
    return pd.read_csv(filename, index_col=0)

def get_UP_DOWN_genes(df, log2FoldChangeThreshold=1):
    """
    Returns two lists of upregulated and downregulated gene names from a
    pandas DataFrame that countains DSeq2 results. The gene name is from
    the index, the column log2FoldChange is used to test for a fold-change
    threshold, for downregulated genes the negative of threshold is applied.
    """
    UP = [e for e in df[df["log2FoldChange"]>=1].index]    # convert to lists
    DOWN = [e for e in df[df["log2FoldChange"]<=-1].index] # of gene names
    return UP, DOWN

def getDSeq2UpDownGenes(filename, log2FoldChangeThreshold=1):
    """
    Returns two lists of upregulated and downregulated gene names from a
    csv file that countains DSeq2 results. The csv file path is specified
    by filename. Genes must meet or exceed a threshold for log2FoldChange,
    whereby for downregulated genes the negative of threshold is applied.
    """
    return get_UP_DOWN_genes(import_DSeq2_csv(filename), log2FoldChangeThreshold=log2FoldChangeThreshold)

def plotVennDiagrams(difRegGenes, genotypes, fig_filename=None):
    """
    Plots Venn-diagrams of differentially regulated genes specified in a dictionary whose keys are genotypes to be
    considered. For each genotype a dictionary with upregulated and downregulated gene names is accessed via the
    keys up and down. The figure will consiste of two Venn-diagrams side by side showing the number of overlapping
    genes. If fig_filename is specified the figure will be saved to disk else displayed.

    difRegGenes (dictionary) ... contains lists of up- and downregulated genes for genotypes accessed as [gt][reg]
    genotypes (list) ........... a list of strings that must match the keys for genotypes in difRegGenes dictionary
    fig_filename (path) ........ a string that identifies a valid filepath for writing a figure to. If None (default)
                                 then the figure is displayed.
    """
    fig=plt.figure(num="VennDiagrams", facecolor='w', edgecolor='k')   # figsize=(10, 6), dpi=300 # figsize=(10, 6)(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1]*2)
    panel_col=0
    for reg in ["up", "down"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[0, panel_col])
            ax.set_title(reg+"regulated genes")
            sets=[]
            names=[]
            for gt in genotypes:
                sets.append(set(difRegGenes[gt][reg]))
                names.append(gt)
            venn3(sets, names)
        panel_col+=1
    fig.tight_layout()
    if fig_filename is None:
        plt.show()
    else:
        fig.savefig(fig_filename)
    plt.close("VennDiagrams")

def plot_UP_DOWN_Genes(df, log2FoldChangeThreshold=1, padj_cutoff=0.05, filename=None, ax=None, show_figure=False, show_Xist=False):
    fig=None
    if ax is None:
        fig=plt.figure(num="differentialGeneExpression", facecolor='w', edgecolor='k')
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(111)
    if not(filename is None):
        ax.set_title(filename.split(".")[0])
    # ax.set_xscale("log")  # we get log2 change values from DSeq2 therefore no log axis needed
    ax.set_yscale("log")
    filteredDF=df.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    filteredDF=filteredDF.dropna(axis=0, how="any")
    ax.plot(filteredDF["log2FoldChange"], 1/filteredDF["padj"], ".", color="gray", alpha=0.5)
    filteredDF=filteredDF.sort_values(by="log2FoldChange", ascending=False)
    cutoff = padj_cutoff
    upDF=filteredDF[ filteredDF["log2FoldChange"]>=log2FoldChangeThreshold ]    
    upDF=upDF[ upDF["padj"]<=cutoff ]
    downDF=filteredDF[ filteredDF["log2FoldChange"]<=(-float(log2FoldChangeThreshold)) ] # for log neg threshold
    downDF=downDF[ downDF["padj"]<=cutoff ]
    if not (filename is None):
        filteredDF.to_csv(DATAPATH+filename+".csv")
    upGenesList=list(upDF.index)
    ax.plot(upDF["log2FoldChange"], 1/upDF["padj"], ".", color="darkcyan")
    if show_Xist and ("Xist" in filteredDF.index):
        ax.plot(filteredDF.loc["Xist"]["log2FoldChange"], 1/filteredDF.loc["Xist"]["padj"], "o", color="darkcyan")
        ax.text(filteredDF.loc["Xist"]["log2FoldChange"]*0.7, 2/filteredDF.loc["Xist"]["padj"], "Xist", color="black")
    downGenesList=list(downDF.index)
    ax.plot(downDF["log2FoldChange"], 1/downDF["padj"], ".", color="firebrick")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("1/adj. p-value")
    if not (fig is None):
        fig.tight_layout()
    if show_figure:
        plt.show()
    if not (filename is None):
        upDF.to_csv(DATAPATH+filename.split(".")[0]+"UP.csv")
        upDF.to_csv(DATAPATH+filename.split(".")[0]+"DOWN.csv")
        fig.savefig(DATAPATH+filename)    
    if not (fig is None):
        plt.close("differentialGeneExpression")
    print(len(upGenesList),len(upDF),"up and", len(downGenesList), len(downDF), "down regulated")
    return upGenesList, downGenesList


# This function requires import of xpRNApy with a valid analysisProject annotation
def annotate_csv_file(filename):
    """ Reads a csv file with gene names in first column into a pandas DataFrame. Adds gene annotation form an XpressionAnalysis object
    that has been initialized with gene annotation. Columns with gene_annotation information as name are added from the annotation.
    The DataFrame returned has the index set to the gene name.
    """
    df=pd.DataFrame.from_csv(filename)
    for g in df.index:
        for anitem in analysisProject.gene_annotation:
            df.at[g,anitem]=analysisProject.genes[g][anitem]
    return df
 
def annotate_all_Dox_noDox_csv(directory, genotypes=["WT", "Ubn1","Ubn2","Ubn12","Hira"]):
    for gt in genotypes:
        filename = os.path.join(directory, gt + "_Dox_vs_Ctrl_results.csv")
        print("Reading", filename)
        df = annotate_csv_file(filename)
        out_filename = os.path.join(directory,gt + "_Dox_vs_Ctrl.csv")
        print("Writing",out_filename)
        df.to_csv(out_filename)
        

if __name__ == "__main__":
    directory = ''
    genotypes = ["Ubn1","Ubn2","Ubn12","Hira"]
    gs=gridspec.GridSpec(nrows=1+int((len(genotypes)-1.5)/2), ncols=min(len(genotypes)-1,2), width_ratios=[1]*2, height_ratios=[1]*(1+int((len(genotypes)-1.5)/2)))
    panel_row=0
    panel_col=0
    difRegGenes = {}
    fig=plt.figure(num="difGenesExpr", facecolor='w', edgecolor='k', figsize=(10, 10), dpi=300) # figsize=(20, 12)
    for gt in genotypes:
        filename = gt+'_vs_WT_results.csv'
        filename_05 = gt+'_vs_WT_results05.csv'
        fig_filename = "difGeneExprPlot"
        log2FoldChangeThreshold = 1
        UP, DOWN = getDSeq2UpDownGenes(filename_05)
        difRegGenes[gt] = {}
        difRegGenes[gt]["up"] = UP
        difRegGenes[gt]["down"] = DOWN
        print("%d %d-fold upregulated and %d downregulated genes" % (len(UP), 2**log2FoldChangeThreshold, len(DOWN)))
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(gt+" vs. WT")
            ax.set_xlim(left=-4, right=4)
            ax.set_ylim(bottom=1, top=10**7)
            plot_UP_DOWN_Genes(import_DSeq2_csv(filename), log2FoldChangeThreshold=log2FoldChangeThreshold, padj_cutoff=0.05, filename=None, ax=ax, show_figure=False, show_Xist=False)
        panel_col+=1
        if panel_col>=2:
            panel_col=0
            panel_row+=1
        # plot_UP_DOWN_Genes(import_DSeq2_csv(filename), show_figure=True)
    fig.tight_layout()
    for ext in [".svg", ".pdf", "png"]:
        fig.savefig(fig_filename+ext)
    plt.show()
    plt.close("difGenesExpr")
    plotVennDiagrams(difRegGenes, ["Ubn1","Ubn2","Ubn12"], fig_filename = os.path.join(directory, "Venn_Ubn1_Ubn2_Ubn12.svg"))
    plotVennDiagrams(difRegGenes, ["Ubn2","Ubn12","Hira"], fig_filename = os.path.join(directory, "Venn_Ubn2_Ubn12_Hira.svg"))

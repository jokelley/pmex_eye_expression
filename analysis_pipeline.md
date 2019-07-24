# Eye expression analyses *Poecilia mexicana*
Analyses for McGowan *et al.* 2019 "Expression analyses of cave mollies (*Poecilia mexicana*) reveal key genes involved in early-stage eye regression". All code was run on Washington State University's high performance computer.

## Analysis Pipeline

### Raw Reads
Raw paired-end reads sequenced from eye tissues have the NCBI accession PRJNA484191.

### Trimming
Trim raw reads using TrimGalore (v.0.4.2); also requires FastQC (v.0.11.4) and Cutadapt (v.1.9).

#### Round 1
Perform adapter trimming but not quality trimming.

    # trimgalore_fastq_locations.txt is a tab-delimited text file of the locations of the reads and sample output directories
    
    Output_Directory="/path/to/output/directories"

    while read line || [ -n "$line" ];
    do
        read1=$(echo $line | awk '{print $1}')      # path to raw read 1
        read2=$(echo $line | awk '{print $2}')      # path to raw read 2
        dirName=$(echo $line | awk '{print $3}')    # directory names (named by sample)

        mkdir -p $Output_Directory/$dirName

        /path/to/trim_galore --quality 0 --illumina --stringency 6 --fastqc_args "--noextract --nogroup --outdir /path/to/fastqc/output/directory" --output_dir $Output_Directory/$dirName --paired $read1 $read2

    done < trimgalore_fastq_locations.txt
    
#### Round 2
Perform quality trimming and remove poly-A tails.

    # trimgalore_fastq_locations.txt is a tab-delimited text file of the locations of the reads and sample output directories
    
    Output_Directory="/path/to/output/directories"

    while read line || [ -n "$line" ];
    do
	    read1=$(echo $line | awk '{print $1}')      # path to trimmed read 1
	    read2=$(echo $line | awk '{print $2}')      # path to trimmed read 2
	    dirName=$(echo $line | awk '{print $3}')    # directory names (named by sample)

	    mkdir -p $Output_Directory/$dirName

	    /path/to/trim_galore --stringency 6 --adapter "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" --adapter2 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" --fastqc_args "--noextract --nogroup --outdir /path/to/fastqc/output/directory" --output_dir $Output_Directory/$dirName --paired $read1 $read2

    done < trimgalore_fastq_locations.txt

#### Round 3
Remove poly-T tails and discard reads <50 base pairs (bp).

    # trimgalore_fastq_locations.txt is a tab-delimited text file of the locations of the reads and sample output directories
    
    Output_Directory="/path/to/output/directories"

    while read line || [ -n "$line" ];
    do
	    read1=$(echo $line | awk '{print $1}')      # path to trimmed read 1
	    read2=$(echo $line | awk '{print $2}')      # path to trimmed read 2
	    dirName=$(echo $line | awk '{print $3}')    # directory names (named by sample)

	    mkdir -p $Output_Directory/$dirName

	    /path/to/trim_galore --stringency 6 --adapter "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" --adapter2 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" --length 50 --fastqc_args "--noextract --nogroup --outdir /path/to/fastqc/output/directory" --output_dir $Output_Directory/$dirName --paired $read1 $read2

    done < trimgalore_fastq_locations.txt
    
### Mapping
Index the *Poecilia mexicana* reference genome (GenBank: LMXC00000000.1 with an appended mitochondrial genome KC992991) using BWA-mem (v.0.7.12).

    /path/to/bwa index -a bwtsw /path/to/reference/genome/reference.fna

Map trimmed reads to the indexed reference genome using BWA-mem (v.0.7.12).
    
    # sampleID_and_fileNames.txt is a tab-delimited text file of the locations of the reads and sample names
    
    mkdir -p /path/to/output/directory/containing/SAMfiles
    
    while read line || [ -n "$line" ];
    do
        read1=$(echo $line | awk '{print $1}')      # path to trimmed read 1
        read2=$(echo $line | awk '{print $2}')      # path to trimmed read 2
        sampleName=$(echo $line | awk '{print $3}') # sample names
    
        /path/to/bwa mem -t 8 -R '@RG\tID:'${sampleName}'\tSM:'${sampleName}'\tPL:ILLUMINA\tDS:author_name' /path/to/reference/genome/reference.fna $read1 $read2 > /path/to/output/directory/containing/SAMfiles/$sampleName.sam
    
    done < sampleID_and_fileNames.txt
    
### Processing the Mapping Output
Convert the SAM files of the mapped reads to BAM files using SAMtools (v.1.7). Keep the headers.

    while read line || [ -n "$line" ];
    do
        samfile=$(echo $line | awk '{print $1}')      # path to SAM files produced from mapping
        dirName=$(echo $line | awk '{print $2}')      # directory names
    
    /path/to/samtools view -b -h -S $samfile > /path/to/$dirName/$dirName.bam

    done < pathsforSAMtools_view.txt
    
Sort the BAM file by leftmost coordinates using SAMtools (v.1.7).

    while read line || [ -n "$line" ];
    do
        bamfile=$(echo $line | awk '{print $1}')      # path to BAM files
        dirName=$(echo $line | awk '{print $2}')      # directory names

        /path/tosamtools sort $bamfile /data/kelley/projects/kerry/pmex_eye/6_BWA/$dirName/$dirName.sorted.bam

    done < pathsforSAMtools_sort.txt

### Generating a Gene Count Matrix
Produce GTF files for each sample from the BAM files of mapped reads using StringTie (v.1.3.3). GTF files will then be used to create the gene counts matrix.

    # dirnames.txt is a text file with directory names listed in one column
    
    while read line || [ -n "$line" ];
    do
        dirName=$(echo $line | awk '{print $1}')      # directory names

    mkdir -p /path/to/ballgown/$dirName

    /path/to/stringtie /path/to/BAMfiles/$dirName.sorted.bam -o /path/to/ballgown/$dirName/$dirName.gtf -p 4 -G /path/to/GFFfile/reference.gff -B -e

    done < dirnames.txt

Use the Python script prepDE.py (downloaded 2018-08-13) associated with StringTie to extract a gene counts matrix from the GTF files.

    # sample_list_ballgown.txt is a tab-delimited text file listing the sample names and the paths to the GTFfiles in the ballgown directory produced above
    
    module load python/2.7.10

    python prepDE.py -i sample_list_ballgown.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv

### Differential Gene Expression Analyses
All analyses were performed using the R packages EdgeR (v.3.22.5) and limma (v.3.36.5) in RStudio (R v.3.5.1).

Read in gene count matrix.
    
    data <- read.csv("gene_count_matrix.csv", row.names=1)

Install required packages.
    
    source("https://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    biocLite("limma")

Load required packages.
    
    require("limma")
    require("edgeR")

Create grouping factor to be incorporated later into the DGEList.
    
    # Bonita = non-sulfidic surface population
    # Cueva = sulfidic cave population
    # Luna = non-sulfidic cave population
    # PSO = sulfidic surface population
    
    group = c(rep("Bonita",4),rep("Cueva",4),rep("Luna",4),rep("PSO",4))

Remove rows where counts = 0 across all samples.

    total <- data[rowSums(data) > 0,]

Create DGEList object.

    y <- DGEList(counts=total, group=group)

Normalize for RNA composition by creating effective library sizes.

    y <- calcNormFactors(y)

Estimate dispersions using Cox-Reid profile-adjusted likelihood.

    y <- estimateDisp(y, design)

#### Analysis of Cavefish Versus Surface Fish
Perform a quasi-likelihood F-test to determine which genes are significantly differentially expressed between cave populations (Cueva, Luna) compared to surface populations (Bonita, PSO). Write output to a CSV file.

    group <- factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4))

    design <- model.matrix(~0+group)

    fit <- glmQLFit(y, design)

    qlf.1and4vs2and3 <- glmQLFTest(fit, contrast=c(-0.5,0.5,0.5,-0.5))
    
    surface.vs.cave <- topTags(qlf.1and4vs2and3, n = 200, p.value=0.05)
    
    surface.vs.cave <- data.frame(surface.vs.cave)
    
    write.csv(x=surface.vs.cave, file = "surface.vs.cave.csv")
    
#### Analysis of Sulfidic Fish Versus Non-sulfidic Fish
Perform a quasi-likelihood F-test to determine which genes are significantly differentially expressed between sulfidic populations (Cueva, PSO) compared to non-sulfidic populations (Bonita, Luna). Write output to a CSV file.

    qlf.1and3vs2and4 <- glmQLFTest(fit, contrast=c(-0.5,0.5,-0.5,0.5))
    
    sulf.vs.nonsulf <- topTags(qlf.1and3vs2and4, n = 200, p.value=0.05)
    
    write.csv(x=sulf.vs.nonsulf, file = "sulf.vs.nonsulf.csv")
    
### Weighted Gene Co-expression Network Analysis (WGCNA)
Run WGCNA v.1.67 using R v.3.6.0 on Washington State University's high performance computer.
* Code and annotations closely follow https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/, particularly Tutorial I.

Install BiocManager (from https://www.bioconductor.org/install/).
      
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install()

Install WGCNA.
    
    BiocManager::install("WGCNA") 
    require(WGCNA)

Install flashClust.
    
    install.packages("flashClust")
    require(flashClust)

Load edgeR.
    
    BiocManager::install("edgeR")
    require(edgeR)

Load cluster.

    install.packages
    require(cluster)

#### Data Input and Cleaning

The following setting is important, do not omit.

    options(stringsAsFactors = FALSE)

Allow multi-threading.
    
    enableWGCNAThreads()

Read in gene count matrix.

    eye_data <- read.csv("gene_count_matrix.csv", row.names = 1)

Add a grouping factor.
 
    group = c(rep("Bonita",4),rep("Cueva",4),rep("Luna",4),rep("PSO",4))

Filter out gene rows that have 0 reads across all samples and gene rows that have expression in only 1 sample column.

    filtered_data <- eye_data[(rowSums(cpm(eye_data) > 0)) >= 2,]
    dim(filtered_data)

Create a DGElist.

    y <- DGEList(counts=filtered_data, group=group)

Normalize libraries.

    y <- calcNormFactors(y)

Make a data frame of log<sub>2</sub>-counts per million.

    counts_per_mil <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)

Transform the data frame.

    dat_expr_0 = as.data.frame(t(counts_per_mil[,]))

Check data for excessive missing values and identification of outlier samples.
    
    gsg = goodSamplesGenes(dat_expr_0, verbose = 3)
    gsg$allOK

Check for outliers. Clustering by sample not by gene.

    sample_tree = hclust(dist(dat_expr_0), method = "average")
    pdf(file = "1_sample_clustering.pdf", width = 12, height = 9)
    sizeGrWindow(12,9)
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

Plot a line to show the cut.
    
    abline(h = 400, col = "red")
    dev.off()

Determine cluster under the line.

    clust_eye = cutreeStatic(sample_tree, cutHeight = 400, minSize = 1)
    table(clust_eye)

Clust_eye 1 contains the samples we want to keep.

    keep_samples = (clust_eye==1)
    dat_expr = dat_expr_0[keep_samples, ]
    n_genes = ncol(dat_expr)
    n_samples = nrow(dat_expr)

Make a dataframe of traits (order: Bonita, Cueva, Luna, PSO).

    rownames(dat_expr)

    # Non-sulfidic = 0, Sulfidic = 1
    NS_vs_S <- c(rep(0,4), rep(1,4), rep(0,4), rep(1,4))

    # Surface = 0, Cave = 1
    cave_vs_surface <- c(rep(0,4), rep(1,4), rep(1,4), rep(0,4))
    traits <- data.frame(NS_vs_S, cave_vs_surface)

Visualize how the clinical trails relate to the dendrogram.

Re-cluster samples.

    sample_tree_2 = hclust(dist(dat_expr), method = "average")

Convert traits to a color representation: white means low, red means high, grey means missing entry.

    trait_colors = numbers2colors(traits, signed = FALSE)

Plot the sample dendrogram and the colors underneath.
    
    pdf(file = "2_dendrogram_and_trait_heatmap.pdf")
    
    plotDendroAndColors(sample_tree_2, trait_colors,
                    groupLabels = names(traits), 
                    main = "Dendrogram and trait heatmap")
    dev.off()

Save data.

    save(dat_expr, traits, file = "pmex-eye-dataInput.RData")

Load the previously saved data.
    
    l_names = load(file = "pmex-eye-dataInput.RData")

The variable l_names contains the names of the loaded variables.

    l_names

#### Step-by-step (As Opposed to Automatic) Construction of the Gene Network and Identification of Modules

Choose a set of soft-thresholding powers.
     
     powers = c(c(1:10), seq(from = 12, to=20, by=2))

Call the network topology analysis function.

    sft = pickSoftThreshold(dat_expr, powerVector = powers, verbose = 5)

Plot the results.

    pdf(file = 3_scale_indep_and_connectivity.pdf, width = 9, height = 5)
    sizeGrWindow(9, 5)
    par(mfrow = c(1, 2))
    cex1 = 0.9

Scale-free topology fit index as a function of the soft-thresholding power.
    
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
    	xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n",
    	main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    	labels=powers,cex=cex1,col="red")

R<sup>2</sup> cut-off of h.

    abline(h=0.8, col="red")

Mean connectivity as a function of the soft-thresholding power of 4.

    plot(sft$fitIndices[,1], sft$fitIndices[,5],
    	xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
    	main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
    dev.off()

Calculate the adjacencies.

    softPower = 4
    adjacency = adjacency(dat_expr, power = softPower)

Turn adjacency into topological overlap matrix (TOM), then calculate the corresponding dissimilarity.

    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM

Call the hierarchical clustering function.

    geneTree = hclust(as.dist(dissTOM), method = "average")

Plot the resulting clustering tree (dendrogram).

    pdf(file = "4_gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
    sizeGrWindow(12,9)
    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    	labels = FALSE, hang = 0.04)
    dev.off()

Set the minimum module size.

    minModuleSize = 20

Module identification using dynamic tree cut.

    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
    table(dynamicMods)

Now plot the module assignment under the gene dendrogram. Convert numeric labels into colors.

    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)

Plot the dendrogram and colors underneath.

    pdf(file = "5_gene_dendrogram_and_module_colors.pdf", width = 8, height = 6)
    sizeGrWindow(8,6)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
    dev.off()

Save the entire list of environment objects to a file.

    save.image(file='myEnv.RData')
    
Load environment from previous script.

    load('myEnv.RData')

Merge modules whose expression profiles are very similar.

Calculate eigengenes.

    MEList = moduleEigengenes(dat_expr, colors = dynamicColors)
    MEs = MEList$eigengenes

Calculate dissimilarity of module eigengenes.

    MEDiss = 1-cor(MEs)

Cluster module eigengenes.

    METree = hclust(as.dist(MEDiss), method = "average")

Plot the result.

    pdf(file = "6_cluster_module_eigengenes.pdf", width = 7, height = 6)
    sizeGrWindow(7, 6)
    plot(METree, main = "Clustering of module eigengenes",
    	xlab = "", sub = "")

We chose a height cut of 0.25, corresponding to a correlation of 0.75, to merge.

    MEDissThres = 0.25

Plot the cut line into the dendrogram.

    abline(h=MEDissThres, col = "red")
    dev.off()

Call an automatic merging function.

    merge = mergeCloseModules(dat_expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

The merged module colors.

    mergedColors = merge$colors

Eigengenes of the new merged modules.

    mergedMEs = merge$newMEs

To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath.

    pdf(file = "7_geneDendro.pdf", wi = 12, he = 9)
    sizeGrWindow(12, 9)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
    	c("Dynamic Tree Cut", "Merged dynamic"),
    	dendroLabels = FALSE, hang = 0.03,
    	addGuide = TRUE, guideHang = 0.05)
    dev.off()

In the subsequent analysis, we use the merged module colors in mergedColors.

Rename to moduleColors.

    moduleColors = mergedColors

Construct numerical labels corresponding to the colors.

    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(moduleColors, colorOrder)-1
    MEs = mergedMEs

Save module colors and labels for use in subsequent parts.

    save(MEs, moduleLabels, moduleColors, geneTree, file = "pmex-eye-networkConstruction-stepByStep.RData")

#### Relate Modules to External Information and Identify Important Genes

Set up the R session and load results of previous parts.

The following setting is important, do not omit.

    options(stringsAsFactors = FALSE)

Load the expression and trait data saved in the first part.

    lnames = load(file = "pmex-eye-dataInput.RData")

The variable lnames contains the names of loaded variables.

    lnames

Load network data saved in the second part.

    lnames = load(file = "pmex-eye-networkConstruction-stepByStep.RData")
    lnames

Relate modules to external traits.

Quantify module-trait associations, identify modules that are significantly associated with the measured clinical trials.

Define numbers of genes and samples.

    nGenes = ncol(dat_expr)
    nSamples = nrow(dat_expr)

Recalculate MEs with color labels.

    MEs0 = moduleEigengenes(dat_expr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, traits, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

Display correlations and their p-values.

    pdf(file = "8_corelations_and_p_values.pdf", wi = 13, he = 15)
    sizeGrWindow(13,15)
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))

Display the correlation values within a heatmap plot.

    labeledHeatmap(Matrix = moduleTraitCor,
    xLabels = names(traits),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = greenWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1,1),
    main = paste("Module-trait relationships"))
    dev.off()

Gene relationship to trait and important modules: gene significance and module membership.

Define variable enviro containing the environment (cave_vs_surface) column of traits.

    enviro = as.data.frame(traits$cave_vs_surface)
    names(enviro) = "enviro"

Names (colors) of the modules.

    modNames = substring(names(MEs), 3)

    geneModuleMembership = as.data.frame(cor(dat_expr, MEs, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

    names(geneModuleMembership) = paste("MM", modNames, sep="")
    names(MMPvalue) = paste("p.MM", modNames, sep="")

    geneTraitSignificance = as.data.frame(cor(dat_expr, enviro, use = "p"))
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

    names(geneTraitSignificance) = paste("GS.", names(enviro), sep="")
    names(GSPvalue) = paste("p.GS.", names(enviro), sep="")

    module = "coral2"
    column = match(module, modNames)
    moduleGenes = moduleColors==module

    pdf(file = "9_module_membership_vs_gene_significance.pdf", wi = 7, he = 7)
    sizeGrWindow(7, 7)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
    	abs(geneTraitSignificance[moduleGenes, 1]),
    	xlab = paste("Module Membership in", module, "module"),
    	ylab = "Gene significance for environment",
    	main = paste("Module membership vs. gene significance\n"),
    	cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()

    print("These are the gene IDs:")
    names(dat_expr)

    print("These are the gene IDs belonging to the coral2 module:")
    names(dat_expr)[moduleColors=="coral2"]

    annot = read.csv(file = "annotations.csv")
    print("Dimensions of annotation file:")
    dim(annot)
    print("Names in annotation file:")
    names(annot)
    genes = names(dat_expr)
    genes2annot = match(genes, annot$geneID)

The following is the number or probes without annotation. Should return 0.

    print("This should return 0:")
    sum(is.na(genes2annot))

Create the starting data frame.

    geneInfo0 = data.frame(geneName = genes,
    	moduleColor = moduleColors,
    	geneTraitSignificance,
    	GSPvalue)

Order modules by their significance for environment.

    modOrder = order(-abs(cor(MEs, enviro, use = "p")))

Add module membership information in the chosen order.

    for (mod in 1:ncol(geneModuleMembership))
    {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
    	MMPvalue[, modOrder[mod]])
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
    	paste("p.MM.", modNames[modOrder[mod]], sep=""))
    }

Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance.

    geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.enviro))
    geneInfo = geneInfo0[geneOrder, ]

    write.csv(geneInfo, file = "geneInfo.csv")

Save the entire list of environment objects to a file.

    save.image(file='myEnv2.RData')

### Relative Eye Size Statistics

Read in data.

    eye_size <- read.csv("eye_size.csv")

Perform Levene's test for homogeneity of variance.

    install.packages("car")
    library("car")
    leveneTest(RelEye ~ Habitat, data = eye_size)

Data do not have equal variance (Levene's test is significant), therefore a Kruskal-Wallis is required instead of a one-way ANOVA.

Kruskal-Wallis test

    kruskal.test(RelEye ~ Habitat, data = eye_size)

Wilcoxon Rank Sum Test (post-hoc pairwise analyses)

    pairwise.wilcox.test(eye_size$RelEye, eye_size$Habitat,
                     p.adjust.method = "BH")

    

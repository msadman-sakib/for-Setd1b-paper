library(GeneOverlap)

setwd("/Users/sakibs/OneDrive\ -\ stud.uni-goettingen.de/Fischerlab_PhD/0.Lab_Projects201920/Setd1b_final_stuffs_2019/with_neuron_gene_overlap_oddsRatio")

neuron = scan("genelists/neuron_polyA_basemean50_FC5_adj0.01_ruvseq.txt", character(), quote = "")
setd1b = scan("genelists/setd1b_K4me3_down.txt", character(), quote = "")
kmt2a = scan("genelists/Kmt2a_K4me3_down.txt", character(), quote = "")
kmt2b = scan("genelists/Kmt2b_K4me3_down.txt", character(), quote = "")

ChIPseq_gene.list = list("setd1b"=setd1b, "kmt2a"=kmt2a, "kmt2b"=kmt2b)
RNAseq_gene.list = list("neuron"=neuron, "neuron"=neuron)


#gs.RNASeq=11091 ##take the sum of input gene lists...manually o deoa jae..
gs.RNASeq= as.double(sum(length(neuron) + length(setd1b) + length(kmt2a) + length(kmt2b))) ###take the sum of input gene lists...


gom.obj <- newGOM(ChIPseq_gene.list, RNAseq_gene.list,
                  + gs.RNASeq)

print(gom.obj) # to get the values

pdf("odds_ratio.pdf")
drawHeatmap(gom.obj, what="odds.ratio",ncolused=9, grid.col="Oranges", note.col="black", adj.p=T, cutoff=.05, log.scale= F)
dev.off()


save.image("odds_ratio_neuronalgene_setd1bKmt2A2B.Rdata")

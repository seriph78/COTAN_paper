library(Seurat)
library(COTAN)
raw = Read10X("../cotan01/data/2020.02.04/filtered_feature_bc_matrix/")

geni.n = c(20000)
cells.n = c(1000,3000, 6000, 10000)

for(gg in geni.n){
    print(gg)
    genes.kept = sample(c(1:nrow(raw)), gg)

    for (cg in cells.n) {
        print(cg)
        cell.kept = sample(c(1:ncol(raw)), cg)

        data = raw[genes.kept, cell.kept]

        obj= automatic.COTAN.object.creation(df = as.data.frame(data), out_dir = "../Cotan_paper/Data/",
                                             GEO = "10X_neuron_10k_v3", sc.method = "10X",
                                             cond = paste0("T2.genes_",gg,"_cells_",cg))


    }
}


obj = automatic.COTAN.object.creation(df = data,out_dir = "../Cotan_paper/Data/",
                                      GEO = "pbmc33k",sc.method = "10X",cond = "pbmc33k",mt_prefix = "^MT-")

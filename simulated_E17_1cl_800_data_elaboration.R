library(COTAN)
library(data.table)
library(Matrix)
library(ggrepel)
library(latex2exp)
data_dir = "Data/negative_datasets/"
mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))
data = as.data.frame(fread(paste0(data_dir,"simulated.dataset_E17_sqrt.fra.1cl_800cs.csv.gz"),sep = ","))
data = as.data.frame(data)
rownames(data) = data$V1
nu = data[1,]
data = data[2:nrow(data),2:ncol(data)]
out_dir = data_dir

obj = new("scCOTAN",raw = data)
obj = initRaw(obj,GEO="simulated" ,sc.method="dropSeq",cond = "SymE17_1cl_800cs")

genes_to_rem = rownames(obj@raw[grep('^mt', rownames(obj@raw)),]) #genes to remove : mithocondrial
obj@raw = obj@raw[!rownames(obj@raw) %in% genes_to_rem,]
cells_to_rem = colnames(obj@raw[which(colSums(obj@raw) == 0)])
obj@raw = obj@raw[,!colnames(obj@raw) %in% cells_to_rem]

t = "SymE17_1cl_800cs"

print(paste("Condition ",t,sep = ""))
#--------------------------------------
n_cells = length(colnames(obj@raw))
print(paste("n cells", n_cells, sep = " "))

n_it = 1

if(!file.exists(out_dir)){
    dir.create(file.path(out_dir))
}

if(!file.exists(paste(out_dir,"cleaning", sep = ""))){   
    dir.create(file.path(out_dir, "cleaning"))
}
ttm = clean(obj)

obj = ttm$object

ttm$pca.cell.2

pdf(paste(out_dir,"cleaning/",t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = ""))
ttm$pca.cell.2
ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
    geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                    nudge_y      = 0.05,
                    nudge_x      = 0.05,
                    direction    = "x",
                    angle        = 90,
                    vjust        = 0,
                    segment.size = 0.2)+
    ggtitle(label = "B cell group genes mean expression", subtitle = " - B group NOT removed -")+my_theme +
    theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 10,hjust = 0.02 ),
          plot.subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 ))

dev.off()

nu_est = round(obj@nu, digits = 7)

plot.nu <-ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))

plot.nu = plot.nu + geom_point(size = 1,alpha= 0.8)+
    scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                          midpoint = log(mean(nu_est)),name = TeX(" $ln (\\nu) $ "))+
    ggtitle("Cells PCA coloured by cells efficiency") +
    my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                      legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                      legend.text = element_text(color = "#3C5488FF", size = 11),
                      legend.key.width = unit(2, "mm"),
                      legend.position="right")

pdf(paste(out_dir,"cleaning/",t,"_plots_PCA_efficiency_colored.pdf", sep = ""))
plot.nu
dev.off()
plot.nu

nu_df = data.frame("nu"= sort(obj@nu), "n"=c(1:length(obj@nu)))

ggplot(nu_df, aes(x = n, y=nu)) + 
    geom_point(colour = "#8491B4B2", size=1)+
    my_theme #+ ylim(0,1) + xlim(0,70)

pdf(paste(out_dir,"cleaning/",t,"_plots_efficiency.pdf", sep = ""))
ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1) +my_theme + #xlim(0,100)+
    annotate(geom="text", x=50, y=0.25, label="nothing to remove ", color="darkred")
dev.off()


obj = cotan_analysis(obj)

# COEX evaluation and storing 
obj = get.coex(obj)

# saving the structure 
saveRDS(obj,file = paste(out_dir,t,".cotan.RDS", sep = ""))

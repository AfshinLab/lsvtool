
library(ggplot2)
library(reshape2)

save_heatmap <- function(df, save_to) {
    tt_melt <- melt(df)

    g <- ggplot(data=tt_melt,
                aes(x=Var2, y=Var1, fill=value)) + geom_tile() +
                geom_text(aes(label=value), color='white') +
                theme_bw() + scale_fill_gradient(limits=c(0, max(tt_melt$value)),low = "yellow", high = "red",)
    
    ggsave(g, file=save_to , width=6, height=6)
}


# Files
comp_mat_file <- snakemake@input[[1]]
names_file <- snakemake@input[[2]]
output_file <- snakemake@output[[1]]
output_file2 <- snakemake@output[[2]]

samples <- read.table(names_file, header=FALSE)

tt <- as.matrix(read.table(comp_mat_file, header=FALSE))

# If matrix is normalized, round to 2 decimals
if (tt[1,1] == 1 & tt[2,2] == 1) {
    tt = round(tt, 2)
}

colnames(tt) <- as.factor(samples$V1)
rownames(tt) <- as.factor(samples$V1)

save_heatmap(tt, save_to=output_file)

# Jaccard index normalization
tt_norm <- tt
for (i in 1:nrow(tt)) {
    for (j in 1:ncol(tt)) {
        if (i == j) {
            tt_norm[i, j] = 1
        } else {
            #tt_norm[i, j] = round(tt[i, j] / (tt[i, i] + tt[j, j] - tt[i, j]), 2) # Jaccard index
            tt_norm[i, j] = round(tt[i, j] / min(c(tt[i, i], tt[j, j])), 2) # Modified jaccard index
        }
    }
}

save_heatmap(tt_norm, save_to=output_file2)

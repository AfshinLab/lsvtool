library("VennDiagram")
library("tools")
library("ggplot2")

intersection_file <- snakemake@input[[1]]
names_file <- snakemake@input[[2]]
svtype <- snakemake@params$svtype
output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]

sink(file=log_file)

ext = file_ext(output_file)

samples <- read.table(names_file, header=FALSE)
number_of_samples = nrow(samples)

if(number_of_samples > 5) {
    number_of_samples = 5      # plot only the first 5
    print("Warning: More than 5 samples in input, limiting to 5.")
}

#Plotting the intersections
t=read.table(intersection_file, header=F)


#Generate the plotting matrix
list_to_plot=list()
for (len in 1:number_of_samples) {
  list_to_plot <- append(list_to_plot, list(which(t[,len]==1)))
  }

names(list_to_plot) <- samples$V1[1:number_of_samples]

# Form https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
myfill <- c("#E69F00", "#56B4E9" ,"#009E73","#F0E442","#0072B2")
myfill <- myfill[1:number_of_samples]
p <- venn.diagram(
  list_to_plot, 
  filename = NULL,
  disable.logging = TRUE,
  main = gsub(",", " & ", svtype),
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  main.cex = 1.2,
  cat.fontface = "bold",
  height = 3000 , 
  width = 3000 , 
  units = "px",
  resolution = 100,
  compression = "lzw",
  imagetype = ext,
  fill = myfill ,
  cex = 1, #size numbers
  lty = 'blank', #dotted line
  #Set names
  cat.cex = 1.1,
  cat.default.pos = "outer",
  );

ggsave(p, file=output_file, device=ext, units="px")
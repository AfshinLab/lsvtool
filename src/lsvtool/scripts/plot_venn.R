library("VennDiagram")

intersection_file <- snakemake@input[[1]]
names_file <- snakemake@input[[2]]
svtype <- snakemake@params$svtype
output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]

sink(file=log_file)

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


myfill <- c("pink", "orange" ,"green","blue","red")
myfill <- myfill[1:number_of_samples]
venn.diagram(list_to_plot, output=True,
            #image  
            filename =  output_file,
            disable.logging = TRUE,
            main = gsub(",", " & ", svtype),
            #print.mode = "percent",
            main.cex = .7,
            height = 1000 , 
            width = 1000 , 
            resolution = 300,
            compression = "lzw",
            imagetype = "png",
            #cyrcles
            fill = myfill ,
            alpha = rep(0.5,number_of_samples),
            cex = .5, #size numbers
            #lty = 4, #dotted line
            lwd = .5, #thickness
            #Set names
            cat.cex = 0.7,
            cat.default.pos = "outer", #or text
            cat.dist = 0.07
            );

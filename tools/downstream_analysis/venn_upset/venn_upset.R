library(VennDiagram)
library(UpSetR)
library(tidyverse)

# read commandline arguments
args <- commandArgs(trailingOnly=TRUE)
invisible(flog.threshold(ERROR))

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide the MAVISp table as argument.", call.=FALSE)
} 

# read MAVISp table
prt <- read.csv(args[1])
prt <- prt %>% 
  rename(
    Sources = Mutation.sources)
prt <- prt[!is.na(prt$Sources),]

# split the sources in multiple rows duplicating the other columns
prt <- prt%>%
  mutate(Sources = strsplit(as.character(Sources), ",")) %>%
  unnest(Sources) %>%
  filter(Sources != "") %>%
  select(Sources, Mutation)

# upset plot lists
upset_lists <- split(prt$Mutation, prt$Sources)

# longest set of mutations
max_lists <- max(lengths(upset_lists))

# total number of mutlists
n_lists <- length(upset_lists)

# UPSET plot
my_upset <- upset(fromList(upset_lists),
      		  nsets = n_lists,
		  order.by = "freq",
		  
		  # scale text size
		  # c(intersection size title, intersection size tick labels, 
		  # set size title, set size tick labels, 
		  # set names, numbers above bars)
      		  text.scale = c(1.2, 1, 1.2, 1, 1, 0.9),
		  
		  # scale matrix/bar ratio based on 
		  # number of mutation lists
		  mb.ratio = c(0.7 - 0.01*n_lists, 0.3 + 0.01*n_lists),
      		  
		  # customize plot
		  point.size = 2.5,
      		  line.size = 1.2,
      		  sets.x.label = "Mutations by source",
      		  matrix.color = c("royalblue"),
      		  main.bar.color = c("royalblue"),
      		  sets.bar.color = c("royalblue"))

# save upset plot as pdf
pdf(file = 'UpSetPlot.pdf', 
    height = 5, 
    width = 7 + 0.0001*max_lists, 
    onefile = FALSE)
my_upset
invisible(capture.output(dev.off()))

# save upset plot as png
png(file = 'UpSetPlot.png', 
    height = 5, 
    width = 7 + 0.0001*max_lists,
    units = 'in',
    pointsize = 12,
    res = 300)
my_upset
invisible(capture.output(dev.off()))

print('Upset plot done!')

# filter for cBioPortal, COSMIC, and ClinVar as mutation sources
prt_venn <- prt %>%
  filter(str_detect(Sources, "cBioPortal|COSMIC|clinvar"))

# venn diagram lists
venn_lists <- split(prt_venn$Mutation, prt_venn$Sources)

# venn colors
venn_col<- c("cadetblue1", "blue", "chartreuse")

# VENN DIAGRAM
myvenn <- venn.diagram(x=venn_lists, 
                       filename = NULL,
                	     
                       lwd = 3,
                       col = "white",
                       fill = venn_col[1:length(venn_lists)],
                       alpha = 0.7,
                       scaled = FALSE,
		       margin = 0.1,
                       
                       cex = 1.3,
                       fontface = "bold",
                       fontfamily = "sans",
                       
                       cat.cex = 1.3,
                       cat.fontface = "bold",
                       cat.default.pos = "outer")

# make venn diagram and save it as pdf
pdf(file='VennDiagram.pdf', 
    height= 7, 
    width=7, 
    onefile=FALSE)
grid.draw(myvenn)
invisible(capture.output(dev.off()))

# make venn diagram and save it as pdf
png(file = 'VennDiagram.png', 
    height = 7, 
    width = 7,
    unit = 'in',
    pointsize = 12, 
    res = 300)
grid.draw(myvenn)
invisible(capture.output(dev.off()))

print('Venn diagram done!')

invisible(file.exists("Rplots.pdf"))
invisible(file.remove("Rplots.pdf"))

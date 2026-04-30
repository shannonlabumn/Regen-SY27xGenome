library(ape)

tree <- read.tree("SpeciesTree_rooted.txt") # From orthofinder output

# Convert to ultrametric (chronogram)
ultra_tree <- chronos(tree)

write.tree(ultra_tree, file="SpeciesTree_ultrametric.nwk")
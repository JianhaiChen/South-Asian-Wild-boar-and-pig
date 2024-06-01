###Infer the phylogenetic tree with 27 de novo assemblies and 3 outgroup species.

./mash sketch -s 400 -k 19  -o chunk *.fa

./mash dist -t chunk.msh chunk.msh >dist.table

####-------------------------------Note--------------------

##### The following is in R, after changing "#query" to "query" in dist.table. Using Babyrousa as the root  ###

library(ape)

mash_dist <- read.table("dist.table", header = TRUE, row.names = 1, check.names = FALSE)

head(mash_dist)

mash_dist_matrix <- as.matrix(mash_dist)

if (nrow(mash_dist_matrix) != ncol(mash_dist_matrix)) {
  stop("The distance matrix is not square!")
}

dist_matrix <- as.dist(mash_dist_matrix)

tree <- nj(dist_matrix)

rooted_tree <- root(tree, outgroup = "babyrousa.fa", resolve.root = TRUE)

plot(rooted_tree, main = "Phylogenetic Tree of Wild Boar and Domestic Pigs with babyrousa.fa as Outgroup")

write.tree(rooted_tree, file = "rooted_wild_boar_tree.nwk")
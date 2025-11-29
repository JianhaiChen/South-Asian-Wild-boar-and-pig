Phylogenetic inference overview

Phylogenetic relationships of wild boar and domestic pigs were inferred using multiple complementary approaches. Importantly, traditional phylogenetic methods that ignore gene flow can lead to misleading topologies, especially when ancient introgression and inter-lineage admixture are present. Therefore, our analyses combined both classical tree-building approaches and modern recombination- and admixture-aware methods to obtain a robust evolutionary interpretation.

1) Whole-genome assembly–based phylogenetics

We first constructed a phylogeny using k-mer–based distances derived from 27 de novo assemblies and 3 outgroup species.
This approach (script: MASH.dist.phylogenetic_tree.sh) provides a recombination-free, alignment-independent estimate of genome-wide relatedness.
Across assemblies, the South Asian lineage consistently appeared as the basal/sister group to all other Eurasian wild boars and domestic pigs.

2) IBS distance–based phylogenetics using population resequencing data

We further estimated pairwise genetic distances at the population scale using an identity-by-state (IBS) framework
(script: IBS_phylogenetic_tree.sh).
This method, based on large SNP genotype matrices, also recovered South Asia as the earliest diverging lineage.

However, when examining trees inferred from high- vs. low-recombination genomic regions, we found inconsistent topologies.
This discrepancy is expected because low-recombination regions tend to resist introgression, while high-recombination regions often reflect more recent gene flow.
These contrasting signals indicated that a single bifurcating tree cannot fully describe the evolutionary history of Sus.

3) AdmixtureBayes: jointly estimating phylogeny and introgression

To resolve these inconsistencies, we applied AdmixtureBayes, a Bayesian framework capable of jointly inferring species relationships and mapping admixture edges.
Unlike traditional methods, AdmixtureBayes incorporates the possibility of inter-lineage gene flow, allowing the model to disentangle true historical branching from introgressed signals.

The method produced a highly supported topology (posterior ≈0.896) identifying South Asia as the basal/sister lineage, while also highlighting the ancient introgression events responsible for conflicting signals in classical analyses.
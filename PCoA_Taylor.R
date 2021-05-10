library(vegan)
library(phyloseq)
library(ape)
source("https://raw.githubusercontent.com/jesusNPL/Public/master/ordinate_mod.R")

data("varespec")

### Run ENVFIT ###
# Create a distance matrix
dist <- vegdist(varespec,  method = "bray")

# PCoA using APE 
pcoa_ape <- pcoa(dist)
# PCoA using stats
pcoa_stats <- cmdscale(dist, eig = TRUE)

species.envfit <- envfit(pcoa_ape, varespec, choices = c(1, 2), permutations = 999)
#Error in scores.default(ord, display = display, choices = choices, ...) : cannot find scores

species.envfit <- envfit(pcoa_stats, varespec, choices = c(1, 2), permutations = 999)
# Success

### Run PCoA using {phyloseg} ###
data("GlobalPatterns")
GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

pcoa_phyloseq <- ordinate(GP1, method = "PCoA", distance = "bray")

### Run PCoA using modified function ###
pcoa_phyloseq_mod <- ordinate2(GP1, method = "PCoA", distance = "bray")

##### Compare PCoA objects from all four methods #####

# PCoA ape
names(pcoa_ape)
# PCoA phyloseq
names(pcoa_phyloseq)
# PCoA stats
names(pcoa_stats)
# PCoA phyloseq modified
names(pcoa_phyloseq_mod)



################################################################################
# Where_is_my_bug v1.5, 6/8/2022
################################################################################

# This script is designed to carry out a search for one or 
# more taxa in FoodMicrobionet and to perform a number of statistical and
# graphical analyses to support the formulation of hypotheses on their
# ecological distribution in food and food environments.
# The commands are only meant to provide some guidance. For each taxon/food
# classification level combination you probably need specialized filtering steps


# to use this script
# 1. you need to have the script in a folder and make that folder a project
#    folder in RStudio
# 2. the folder must contain a subfolder named FMBN with the FMBN_plus.rds
#    version of FoodMicrobionet (the most recent version can be downloaded from)
#    GitHub (https://github.com/ep142/FoodMicrobionet/tree/master/the_real_thing)
# 

# Install/load packages ---------------------------------------------------

.cran_packages <- c("crayon", "tidyverse", "parallel",  "data.table", "phylotools", 
                    "progress", "beepr", "tictoc")


.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)



# the following command detects the number of cores on UNIX/MacOS
nc <- parallel::detectCores(logical = F) # to detect physical cores in MacOS
# for reproducibility reasons
set.seed(100)
# mainly for reproducibility reasons (this will be saved with the work space)
r_version <- R.Version() # in the future may be check that the version running is compatible
# with the script (it only matters if somebody other than myself is using the script)
session_info <- sessionInfo()

# use this option to play a sound when an operation (for example in the make_tree
# function), is completed
play_sound <- T

# a manual color scale for T/F
TF_colors <- c("#F8766D", "#619CFF")
names(TF_colors) <- c("FALSE", "TRUE")


# load FMBN ---------------------------------------------------------------

FMBN_path <- file.path("data", "FMBN_plus.rds")
FMBN <- readRDS(file = FMBN_path)

# the taxa to search for --------------------------------------------------

# in the future I will transform this in an interactive app
# the taxonomic level for which the search will be performed
# must be %in% c("phylum","class","order","family","genus")
# one can technically search for a species but it is not wise
# the script is designed to look for all distinct taxa at the taxonomic
# level immediately below: so, if I select a family, the script will check 
# for genera in that family (therefore I should probably exclude genus or
# process it in a different way)

taxon_level <- "genus" # perhaps should limit to order, as a maximum?
# choice is between "order", "family", "genus"

# I could set the lower resolution with something like this
go_down_to <- "genus" # "genus" or "species", but many taxa are not identified at species level


match_tax_level <- taxon_level %in% colnames(FMBN$taxa)[3:8]
match_go_down <- go_down_to %in% colnames(FMBN$taxa)[3:9]

# just to have something I will reuse in the app for error trapping
# even if in the app I will use a drop down menu 

if(!(match_tax_level & match_go_down)) cat("\nWARNING: wrong taxonomic level(s)! Must be one of:\n",
                         colnames(FMBN$taxa)[3:9], "\n")


# must be a taxonomic level below taxon_level, so the following must be true
right_tax_order <- which(go_down_to == colnames(FMBN$taxa))[[1]] > which(taxon_level == colnames(FMBN$taxa))[[1]]
# right_tax_order <-T # I am forcing this

# get unique levels for the taxonomic level (to be used in dropdown in the future)
# NOTE FOR SELF I HAVE BENCHMARKED AND THIS IS BETTER THAN USING DISTINCT
# only select bacteria
taxa <- FMBN$taxa
FMBN_bacteria <- taxa |>
  dplyr::filter(domain == "Bacteria")

unique_taxa <- unique(FMBN_bacteria[[taxon_level]])
unique_taxa <- unique_taxa[!is.na(unique_taxa)]

# in a menu in a Shiny app this will autocomplete
# technically can be a vector of length >1

# load the file with taxa
sporeformers <- read_tsv("sporeformers.txt")

taxon_to_search <- sporeformers$Genus


# check if it is in unique_taxa
# %chin% is an oprator provided by data.table and should be faster than %in%
is_valid_taxon <- any(taxon_to_search %chin% unique_taxa)

# let's simplify the list
sporeformers_in_FMBN <- taxon_to_search[which(taxon_to_search %chin% unique_taxa)]
is_valid_taxon <- all(sporeformers_in_FMBN %chin% unique_taxa)
if(!is_valid_taxon) cat("\nWARNING: can't find your taxon in", taxon_level, "\n")

# the procedure for returning the matching taxa is different depending on the length
# of taxon_to_search

if(length(taxon_to_search) == 1){
  my_taxa <- dplyr::filter(FMBN_bacteria, .data[[taxon_level]] == taxon_to_search)
} else {
  my_taxa <- dplyr::filter(FMBN_bacteria, (.data[[taxon_level]] %in% taxon_to_search))
}
# getting the counts of taxa
taxa_count <- my_taxa %>% dplyr::count(.data[[taxon_level]])

# remove NA in the level immediately below taxon_level

tax_level_below <- case_when(
  taxon_level == "order" ~ "family",
  taxon_level == "family" ~ "genus",
  taxon_level == "genus" ~ "species",
)

my_taxa_all <- my_taxa
my_taxa <- dplyr::filter(my_taxa, !is.na(.data[[tax_level_below]]))
# remove Actinobacteriota, only Firmicutes = Bacillota remain
my_taxa <- dplyr::filter(my_taxa, phylum != "Actinobacteriota")



# get the edges -----------------------------------------------------------

edges <- FMBN$edges_B 

# options for removing Chloroplast and/or Mitochondria and and/or poorly identified taxa
remove_chloroplast <- T
remove_mitochondria <- T
remove_no_phylum <- T

if(remove_chloroplast){
  chloroplast_id <- FMBN$taxa %>% dplyr::filter(label == "Chloroplast" | class == "Chloroplast") %>% pull(taxonId)
  edges <- edges <- dplyr::filter(edges, !(taxonId %in% chloroplast_id))
}

if(remove_mitochondria){
  mitochondria_id <- FMBN$taxa %>% dplyr::filter(label == "Mitochondria" | family == "Mitochondria") %>% pull(taxonId)
  edges <- edges %>% dplyr::filter(!(taxonId %in% mitochondria_id))
}
  
if(remove_no_phylum){
  no_phylum_id <- FMBN$taxa %>% dplyr::filter(is.na(phylum)) %>% pull(taxonId)
  edges <- edges %>% dplyr::filter(!(taxonId %in% no_phylum_id))
}
# optionally recalculare weight
if(any(remove_chloroplast, remove_mitochondria, remove_no_phylum)){
  edges <- edges %>% 
    group_by(sampleId) %>%
    mutate(weight_old = weight) %>%
    mutate(weight = weight/sum(weight))
}

# note that now weight is a fraction rather than % and weight_old is the uncorrected weight

n_edges <- nrow(edges)
# this is a filtering join
edges <- semi_join(edges, my_taxa)

prop_edges_my_taxa <- nrow(edges)/n_edges

# the abundance should be summed or not depending on the value of go_down_to

# join taxonomy
edges <- left_join(edges, select(my_taxa, taxonId:species))

# pool, as applicable

pooled_edges <- edges %>%
      group_by(sampleId, .data[[go_down_to]]) %>%
      summarise(weight = sum(weight)) %>%
      ungroup() %>%
      arrange(sampleId, genus)

# join with sample information


# optionally filter samples -----------------------------------------------
# there is only two samples with unbottled water (which are also blanks), I am removing them
filt_samples <- FMBN$samples_B

# I am removing a few samples which are PCR standards or mock community
which(str_detect(filt_samples$L1, "blank") | str_detect(filt_samples$L1, "Mock"))
filt_samples <- filt_samples %>% 
  dplyr::filter(
    !str_detect(filt_samples$L1, "blank") | str_detect(filt_samples$L1, "Mock")
  )
filt_samples <- filt_samples %>% dplyr::filter(s_type == "Sample")
which(str_detect(filt_samples$L1, "^Water"))
filt_samples <- filt_samples %>% dplyr::filter(!str_detect(filt_samples$L1, "^Water"))

# AT THIS STAGE YOU MAY WANT TO PERFORM SOME MORE FILTERING BASED (POSSIBLY)
# ON TAXA, REGION, NUMBER OF SEQUENCES, NUMBER OF ISSUES ETC.
# OR FILTERING BASED ON NATURE (SAMPLE/ENVIRONMENT), FERMENTATION, ETC.
# THIS NEEDS TO BE DONE BEFORE POOLING


# I am adding some diagnostic info from filt_samples just in case
# because this is also a filtering join it removes unwanted samples from edges
pooled_edges_sample <- inner_join(pooled_edges,
                                    select(filt_samples, studyId, sampleId,
                                           s_type, n_reads2, n_issues, foodId,
                                           L1:target2, geo_loc_country)
                                    )



# Optionally filter studies -----------------------------------------------

filt_studies <- FMBN$studies
# add filtering here
# only keep studies with DADA2
filt_studies <- filt_studies %>%
  dplyr::filter(str_detect(bioinf_software, "dada2"))

# join study info (some may be superfluous)
pooled_edges_sample_study <- inner_join(pooled_edges_sample,
                                       select(filt_studies, studyId, read_length_bp,  
                                              Seq_accn, food_group, overlapping,
                                              paired_end)
                                       )

# after filtering
pooled_edges_sample_study %>% distinct(studyId) %>% nrow() # 112 studies from the original 239
pooled_edges_sample_study %>% distinct(sampleId) %>% nrow() # 1990 samples contain sequences belonging to the taxa
                                                            # out of 11253


# calculate prevalence and statistics on abundance  -----------------------

# choose the food grouping variable ---------------------------------------

# food grouping variable

food_grouping_variable <- "L1" 
# should be one of L1, L4, L6; however, if needed, further levels can be 
# joined using foodId from FMBN$foodex2exp

# remove any NA deriving from the left_join of samples
pooled_edges_sample_study <- pooled_edges_sample_study %>%
  dplyr::filter(!is.na(.data[[food_grouping_variable]]))

n_samples <- filt_samples %>%
  dplyr::filter(sampleId %in% unique(pooled_edges_sample_study$sampleId)) %>%
  group_by(.data[[food_grouping_variable]]) %>%
  summarise(n = n_distinct(sampleId))

# join original info and look at proportions (which is the prevalence of Lactobacillaceae)
n_samples_all <- filt_samples %>% dplyr::filter(s_type == "Sample") %>%
  group_by(.data[[food_grouping_variable]]) %>%
  summarise(n = n_distinct(sampleId)) %>%
  dplyr::rename(n_all = n)

n_samples <- left_join(n_samples, n_samples_all) %>% 
  mutate(prev_sporef = n/n_all) %>%
  arrange(desc(prev_sporef))

# calculate prevalence and abundance --------------------------------------

prev_df <- pooled_edges_sample_study %>%
  group_by(genus, .data[[food_grouping_variable]]) %>%
  summarise(abs_prev = n_distinct(sampleId)) %>%
  left_join(., n_samples) %>%
  mutate(prev = abs_prev/n)

# I could use an option for using log(abundance) and/or for choosing the base of the log (2 or 10) 
log_base <- 10
ab_df <- pooled_edges_sample_study %>%
  mutate(log_ab = log(weight, base = log_base)) %>%
  group_by(genus, .data[[food_grouping_variable]]) %>%
  summarise(min_ab = min(weight),
            max_ab = max(weight),
            mean_ab = mean(weight),
            median_ab = median(weight),
            mad_ab = mad(weight),
            median_logab = median(log_ab))

prev_ab_df <- full_join(prev_df, ab_df)
write_tsv(prev_ab_df, "sporeformers_prev_ab_df.txt")

# further filter for samples before plotting ------------------------------

# I am restricting to V1-V3 and V3-V4 because of lack of reliability for V4 at the genus level

pooled_edges_sample_study_all <- pooled_edges_sample_study
pooled_edges_sample_study <- pooled_edges_sample_study %>% 
  dplyr::filter(target2 == "V1-V3" | target2 == "V3-V4")

hist(pooled_edges_sample_study$read_length_bp)
# removing shorter sequences
pooled_edges_sample_study <- pooled_edges_sample_study %>% dplyr::filter(read_length_bp>=350)

# how many studies left?
residual_studies <- unique(pooled_edges_sample_study$studyId) # 68, not too bad
# how many samples left?
residual_samples <- unique(pooled_edges_sample_study$sampleId) # 1007, not too bad

# recalculate prev and ab
prev_df_sel <- pooled_edges_sample_study %>%
  group_by(genus, .data[[food_grouping_variable]]) %>%
  summarise(abs_prev = n_distinct(sampleId)) %>%
  left_join(., n_samples) %>%
  mutate(prev = abs_prev/n)

# I could use an option for using log(abundance) and/or for choosing the base of the log (2 or 10) 
ab_df_sel <- pooled_edges_sample_study %>%
  mutate(log_ab = log(weight, base = log_base)) %>%
  group_by(genus, .data[[food_grouping_variable]]) %>%
  summarise(min_ab = min(weight),
            max_ab = max(weight),
            mean_ab = mean(weight),
            median_ab = median(weight),
            mad_ab = mad(weight),
            median_logab = median(log_ab))

prev_ab_df_sel <- full_join(prev_df_sel, ab_df_sel)

# now save a reordered table
write_tsv(arrange(prev_ab_df_sel, genus, desc(prev), desc(median_ab)), file = "prevab_sporeformers_V1_V4.txt")

# boxplots for prevalence and abundance -----------------------------------

# first: a box plot: gives a better idea of the distribution

pooled_edges_sample_study <- pooled_edges_sample_study %>% 
  mutate(log_ab = log(weight, base = log_base))


# SPECIALISED STEP
# I am transforming genus into an ordered factor, to make plotting more 
# coherent with phylogenetic relatedness. Be aware that not all members of
# family Lactobacillaceae might be present
# Acetilactobacillus and Convivina are missing

genera <- sort(unique(pooled_edges_sample_study$genus))

pooled_edges_sample_study$genus <- factor(pooled_edges_sample_study$genus)

# saveRDS(pooled_edges_sample_study, file = "pooled_edges_sample_study.rds")

weightlabel <- if_else(log_base == 2, 
                       "log2(rel. abundance)",
                       "log10(rel. abundance)")

box_plot <- ggplot(pooled_edges_sample_study, mapping = aes(x = genus, y = log_ab))

if(log_base == 10){
box_plot +
  facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
  geom_boxplot(alpha = I(0.6)) +
  labs(y = weightlabel) +
  annotation_logticks(sides = "l", alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
} else {
  box_plot +
    facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
    geom_boxplot(alpha = I(0.6)) +
    labs(y = weightlabel) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  }

# a color scale for fermentation and spoilage

fermspoil_scale <- c("blue", "red", "cyan", "mediumorchid")
names(fermspoil_scale) <- c("Unspoiled", "Spoiled", "Fermented", "Fermented+Spoiled")
pie(c(1,1,1,1), col = fermspoil_scale, labels = names(fermspoil_scale))

# let's try a variation with colors
if(log_base == 10){
  bplot_jitter_facetfood <- box_plot +
    facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
    geom_boxplot() +
    geom_jitter(aes(colour = spoilage), alpha = I(0.4), width = 0.1, shape = I(16), size = I(1)) +
    labs(y = weightlabel) +
    annotation_logticks(sides = "l", alpha = 0.2) +
    scale_colour_manual(values = fermspoil_scale) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
} else {
  bplot_jitter_facetfood <- box_plot +
    facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
    geom_boxplot() +
    geom_jitter(aes(colour = spoilage), alpha = I(0.4), width = 0.1, shape = I(16), size = I(1)) +
    labs(y = weightlabel) +
    annotation_logticks(sides = "l", alpha = 0.2) +
    scale_colour_manual(values = fermspoil_scale) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
bplot_jitter_facetfood
# not very readable
ggsave(bplot_jitter_facetfood, 
        filename = str_c("bplot_jitter_facetfood", "_Lactobacillaceae", ".tiff", sep = ""), dpi = 300,
        width = 14, height = 10)


# now add prevalence info

if(log_base == 10){
  box_plot +
    facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
    geom_boxplot(alpha = I(0.4)) +
    geom_point(prev_ab_df_sel, mapping = aes(y = median_logab, colour = log(100*prev, base = 10))) +
    geom_hline(yintercept = -3, linetype = "dotted") +
    labs(y = weightlabel, colour = "log(10(prev*100)") +
    annotation_logticks(sides = "l", alpha = 0.2) +
    scale_colour_continuous(type = "viridis", direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"))
} else {
  box_plot +
    facet_wrap(~str_wrap(.data[[food_grouping_variable]],40)) +
    geom_boxplot(alpha = I(0.4)) +
    geom_point(prev_ab_df_sel, mapping = aes(y = median_logab, colour = prev)) +
    geom_hline(yintercept = -3, linetype = "dotted") +
    labs(y = weightlabel, colour = "prevalence") +
    annotation_logticks(sides = "l", alpha = 0.5) +
    scale_colour_continuous(type = "viridis", direction = -1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"))
}
ggsave("Figx_boxplotgenus_sporeformers_facet_food.tiff", width = 14, height = 10, dpi = 300)
ggsave("Figx_boxplotgenus_sporeformers_facet_food.jpg", width = 14, height = 10, dpi = 300)

# same but changing what goes on facets

box_plot_food <- ggplot(pooled_edges_sample_study, 
                        mapping = aes(x = str_trunc(.data[[food_grouping_variable]], 20, "right"), y = log_ab))


box_plot_food +
  facet_wrap(~str_wrap(genus,20)) +
  geom_boxplot(alpha = I(0.4)) +
  labs(y = weightlabel, x = food_grouping_variable) +
  scale_y_continuous(breaks = c(-4, -3, -2,-1, 0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  
box_plot_food +
  facet_wrap(~str_wrap(genus,20)) +
  geom_boxplot(alpha = I(0.6)) +
  geom_point(prev_ab_df_sel, mapping = aes(y = median_logab, colour = prev)) +
  labs(y = weightlabel, x = food_grouping_variable) +
  scale_colour_continuous(type = "viridis", direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# specialized: using str_wrap apparently drops the factor levels
bplot_sporeformers_facet_genera <- box_plot_food +
  facet_wrap(~genus) +
  geom_boxplot(alpha = I(0.6)) +
  geom_point(prev_ab_df_sel, mapping = aes(y = median_logab, colour = prev)) +
  labs(y = weightlabel, x = food_grouping_variable, colour = "prevalence") +
  scale_y_continuous(breaks = seq(-5,0)) +
  scale_colour_continuous(type = "viridis", direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(face = "italic"))
bplot_sporeformers_facet_genera
ggsave(bplot_sporeformers_facet_genera, filename = "boxplotfood_sporeformers.jpg", width = 14, height = 10, dpi = 300)



bplot_coordflip_facet_genus <- ggplot(pooled_edges_sample_study, 
                                      mapping = aes(x = str_trunc(.data[[food_grouping_variable]], 20, "right"), y = log_ab)) +
  geom_boxplot(fatten = 2) +
  geom_point(data = prev_ab_df_sel, mapping = aes(y = median_logab, colour = log(100*prev, base = 10)), shape = I(16)) +
  facet_wrap(~genus) +
  labs(y = weightlabel, x = food_grouping_variable, colour = "log10(prev*100)") +
  coord_flip() +
  geom_hline(yintercept = -3, linetype = "dotted") +
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0)) +
  scale_colour_continuous(type = "viridis", direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, face = "bold"),
        strip.text = element_text(face = "italic"))
bplot_coordflip_facet_genus
# for some reason facet_wrap drops the order of genus and plots in alphabetical order
bplot_coordflip_facet_genus <- ggplot(pooled_edges_sample_study, 
                                      mapping = aes(x = str_trunc(.data[[food_grouping_variable]], 20, "right"), y = log_ab)) +
  geom_boxplot(fatten = 2) +
  geom_point(data = prev_ab_df_sel, mapping = aes(y = median_logab, colour = log(100*prev, base = 10)), shape = I(16)) +
  facet_wrap(~genus) +
  labs(y = weightlabel, x = food_grouping_variable, colour = "log10(prev*100)") +
  coord_flip() +
  geom_hline(yintercept = -3, linetype = "dotted") +
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0)) +
  scale_colour_continuous(type = "viridis", direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, face = "bold"),
        strip.text = element_text(face = "italic"))
bplot_coordflip_facet_genus

ggsave(bplot_coordflip_facet_genus, filename = "FigXboxplotfood_prev_sporeformers_l_LS.tiff", width = 14, height = 10, dpi = 300)

# Credits and copyright ---------------------------------------------------

# The code for sequence alignment and phylogenetic trees is adapted from 
# https://benjjneb.github.io/dada2/tutorial.html
# the code for loading and attaching packages is taken from:
# https://doi.org/10.12688/f1000research.8986.2 


# Assume that this is overall under MIT licence

# Copyright 2022, 2024 Eugenio Parente
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.

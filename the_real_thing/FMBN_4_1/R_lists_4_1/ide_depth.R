
# ide_depth ---------------------------------------------------------------

# a proof of concept script for evaluating the level of taxonomic assignment in FMBN

# the script operates on FMBN_plus, which has to be available in a directory 
# named "FMBN" in the same working directory of the script

# Install/load packages ---------------------------------------------------

.cran_packages <- c("tidyverse", "data.table", "beepr", "tictoc")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(.cran_packages, require, character.only = TRUE)


# for reproducibility reasons
set.seed(100)
# mainly for reproducibility reasons (this will be saved with the work space)
r_version <- R.Version() # in the future may be check that the version running is compatible
# with the script (it only matters if somebody other than myself is using the script)
session_info <- sessionInfo()

# use this option to play a sound when an operation (for example in the make_tree
# function), is completed
play_sound <- T


# open FMBN_plus ----------------------------------------------------------

FMBN_path <- file.path("FMBN","FMBN_plus.rds")

FMBN <- readRDS(FMBN_path)

edges <- FMBN$edges
samples <- FMBN$samples
studies <- FMBN$studies
taxa <- FMBN$taxa


#  get and annotate edges -------------------------------------------------

# get the edges for studies 34 - 180 (those for which dada2 was used)

# get the maximum sample id for samples in study 33
# this will not work in FMBN 4.2 because sampleIds are not in increasing order

max_sample <- samples %>%
  dplyr::filter(studyId == "ST33") %>% 
  summarize(max_sample = max(sampleId)) %>%
  pull()

edges_sel <- edges %>%
  dplyr::filter(sampleId > max_sample)

# annotate with study and region info

edges_sel_ann <- left_join(
  edges_sel, 
  select(samples, sampleId, studyId, n_reads2, n_issues, L1:L6, target1, target2)
)

# annotate with taxonomic info

edges_sel_ann <- left_join(
  edges_sel_ann, 
  select(taxa, taxonId:species, idelevel)
)

# annotate with further info from studies

edges_sel_ann <- left_join(
  edges_sel_ann, 
  select(studies, studyId, overlapping)
)

# create a variable with region/overlap

edges_sel_ann <- edges_sel_ann %>%
  unite(col = "target", target2, overlapping, remove = F)

# tabulations and graphs --------------------------------------------------

# make idelevel an ordered factor

edges_sel_ann <- edges_sel_ann %>%
  mutate(idelevel = factor(idelevel, 
                           levels = c("species", "genus", "family", "order", 
                                      "class", "plylum", "domain"),
                           ordered = T))

# add a columns with number of sequences per edge
edges_sel_ann <- edges_sel_ann %>% 
  mutate(seqs = weight * n_reads2/100)

# get average number of issues and sequence length by study
issues_length <- studies %>%
  dplyr::slice(34:180) %>%
  select(studyId, read_length_bp, target, region, overlapping) 
ave_issues <- samples %>% 
  dplyr::filter(sampleId > 1722) %>%
  group_by(studyId) %>%
  summarize(ave_issues = mean(n_issues))
issues_length <- left_join(issues_length, ave_issues) %>%
  unite(col = "target2", region, overlapping, remove = F)

summary_tab_edges <- edges_sel_ann %>%
  select(studyId, n_issues, L1:seqs) %>%
  group_by(studyId, .drop = F) 

summary_tab_edges_unw <- summary_tab_edges %>%
  count(idelevel) %>%
  mutate(freq = n/sum(n))

summary_tab_edges_w <- summary_tab_edges %>%
  count(idelevel, wt = seqs) %>%
  mutate(freq = n/sum(n))

summary_tab_edges_both <- left_join(summary_tab_edges_unw,
                                    select(summary_tab_edges_w,
                                           studyId, idelevel, nw = n, freqw = freq))

# join region and seq length and issues

summary_tab_edges_both_ann <- left_join(summary_tab_edges_both,
                                        issues_length)



# a box plot, by region, for identifications at the genus level or below
summary_tab_edges_both_ann_sg <- summary_tab_edges_both_ann %>%
  dplyr::filter(idelevel == "genus" | idelevel == "species") %>%
  group_by(studyId, .drop = F) %>%
  summarize(freq_sg = sum(freq),
            freq_sgw = sum(freqw)) %>%
  left_join(., issues_length)

# medians by region, only those including V3 or V4
summaries <- summary_tab_edges_both_ann_sg %>% 
  dplyr::filter(str_detect(target2, "V3") | str_detect(target2, "V4")) %>%
  ungroup() %>%
  group_by(target2) %>% 
  summarize(n = n(),
            medianfreq = median(freq_sg),
            medianfreqw = median(freq_sgw),
            perc90freq = quantile(freq_sg, 0.9),
            perc90freqw = quantile(freq_sgw, 0.9))
# summaries_all
summaries_all <- summary_tab_edges_both_ann_sg %>% 
  ungroup() %>%
  group_by(target2) %>% 
  summarize(n = n(),
            medianfreq = median(freq_sg),
            medianfreqw = median(freq_sgw),
            perc90freq = quantile(freq_sg, 0.9),
            perc90freqw = quantile(freq_sgw, 0.9))
write_tsv(summaries_all, "idefreqgenusspecies.txt")

# now only for species, with overlapping true
summary_species <- summary_tab_edges_both_ann %>%
  dplyr::filter(str_detect(target2, "TRUE")) %>%
  dplyr::filter(idelevel == "species") %>%
  group_by(studyId, .drop = F) %>%
  summarize(freq_sg = sum(freq),
            freq_sgw = sum(freqw)) %>%
  left_join(., issues_length)

summary_species_table <- summary_species %>% 
  ungroup() %>%
  group_by(target2) %>% 
  summarize(n = n(),
            medianfreq = median(freq_sg),
            medianfreqw = median(freq_sgw),
            perc90freq = quantile(freq_sg, 0.9),
            perc90freqw = quantile(freq_sgw, 0.9))
write_tsv(summary_species_table, "idefreqspecies.txt")


# any relationship with length or quality (not a good graph)
ggplot(summary_tab_edges_both_ann_sg, mapping = aes(x = read_length_bp, y = freq_sg, shape = target2)) +
  geom_point(mapping = aes(color = ave_issues)) +
  scale_shape_manual(values = c(1, 2, 0, 5, 6, 16, 17, 15, 18, 3, 4, 10))
  theme_bw()


ggplot(summary_tab_edges_both_ann_sg, mapping = aes(x = target2, y = freq_sg)) +
  geom_boxplot() + 
  geom_jitter(mapping = aes(color = ave_issues)) +
  labs(x = "target region",
    y = "species + genus freq.",
    color = "ave. issues") +
  scale_y_continuous(breaks = seq(0.6, 1, 0.05), minor_breaks = seq(0.6, 1, 0.01)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("idefreqgenspunw.tiff", dpi = 600)

ggplot(summary_tab_edges_both_ann_sg, mapping = aes(x = target2, y = freq_sgw)) +
  geom_boxplot() + 
  geom_jitter(mapping = aes(color = ave_issues)) +
  labs(x = "target region",
       y = "species + genus freq., weighted",
       color = "ave. issues") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("idefreqgenspw.tiff", dpi = 600)

# only those containing V3 or V4, not V5
summary_tab_edges_both_ann_sg %>%
  dplyr::filter(!str_detect(target2, "V5") & !str_detect(target2, "V6")) %>%
  ggplot(mapping = aes(x = target2, y = freq_sgw)) +
  geom_boxplot() + 
  geom_jitter(mapping = aes(color = ave_issues)) +
  labs(x = "target region",
       y = "species + genus freq., weighted",
       color = "ave. issues") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# differences by phylum ---------------------------------------------------

# which are important phyla?
abundant_phyla <- edges_sel_ann %>%
  group_by(phylum) %>%
  summarize(abundance_sum = sum(weight),
            abundance_mean = mean(weight)) %>%
  arrange(desc(abundance_sum))

four_abundant_phyla <- abundant_phyla %>%
  slice(1:4) %>%
  pull(phylum)

summary_tab_edges_4phyla <- edges_sel_ann %>%
  select(studyId, n_issues, L1:seqs) %>%
  dplyr::filter(phylum %in% four_abundant_phyla) %>%
  group_by(studyId, phylum, .drop = F) 

summary_tab_edges_4phyla_unw <- summary_tab_edges_4phyla %>%
  count(idelevel) %>%
  mutate(freq = n/sum(n))

summary_tab_edges_w_4phyla <- summary_tab_edges_4phyla %>%
  count(idelevel, wt = seqs) %>%
  mutate(freq = n/sum(n))

summary_tab_edges_4_phyla_both <- left_join(summary_tab_edges_4phyla_unw,
                                    select(summary_tab_edges_w_4phyla,
                                           studyId, phylum, idelevel, nw = n, freqw = freq))

summary_tab_edges_4_phyla_both_sg <- summary_tab_edges_4_phyla_both %>%
  ungroup() %>%
  dplyr::filter(idelevel == "genus" | idelevel == "species") %>%
  group_by(studyId, phylum, .drop = F) %>%
  summarize(freq_sg = sum(freq),
            freq_sgw = sum(freqw)) %>%
  left_join(., issues_length)

ggplot(summary_tab_edges_4_phyla_both_sg, mapping = aes(x = phylum, y = freq_sg)) +
  facet_wrap(~target2) +
  geom_boxplot() + 
  geom_jitter(alpha = 0.5) +
  labs(x = "phylum",
       y = "species + genus freq.",
       color = "ave. issues") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("idefreqgenspunw_phylum.tiff", dpi = 600)

ggplot(summary_tab_edges_4_phyla_both_sg, mapping = aes(x = phylum, y = freq_sgw)) +
  facet_wrap(~target2) +
  geom_boxplot() + 
  geom_jitter(alpha = 0.5) +
  labs(x = "phylum",
       y = "species + genus freq., weighted",
       color = "ave. issues") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("idefreqgenspw_phylum.tiff", dpi = 600)

# Citations ---------------------------------------------------------------

map(.cran_packages, citation)

# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# with some changes and adaptations

# Assume that this is overall under MIT licence

# Copyright 2021 Eugenio Parente
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

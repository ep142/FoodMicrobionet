# convert_to_microbiomeanalyst

# this script takes as an input the tab delimited tables generated (as a 
# compressed folder) by FoodMicrobionet and generates a otu and metadata table 
# which can be used with MicrobiomeAnalyst (https://www.microbiomeanalyst.ca/)


require(tidyverse)
# load and convert the tax table
tax_table_0 <- read_tsv("tax_table.txt")
# rename and fix NA
tax_table_1 <- tax_table_0 |>
  rename(kingdom = domain)
names(tax_table_1)[2:8] <- str_to_sentence(names(tax_table_1)[2:8])
# handle NA
tax_table_2 <- tax_table_1 |>
  mutate(Phylum = if_else(is.na(Phylum), "", Phylum)) |>
  mutate(Class = if_else(is.na(Class), "", Class)) |>
  mutate(Order = if_else(is.na(Order), "", Order)) |>
  mutate(Family = if_else(is.na(Family), "", Family)) |>
  mutate(Genus = if_else(is.na(Genus), " ", Genus))
tax_table_2$Species <- map2_chr(.x = tax_table_2$Species, 
                                .y = tax_table_2$Genus,
                                \(x,y) str_remove(x, pattern = y)) |> 
  str_squish()
tax_table_2 <- tax_table_2 |>
  mutate(Genus = str_squish(Genus)) |>
  mutate(Species = if_else(is.na(Species),"spp", Species)) |>
  mutate(rown = seq(1:nrow(tax_table_1))) |>
  mutate(NAME = str_c("OTU",rown))

fix_trailing_semicolon <- function(text) {
  # Replace zero or more semicolons at the end ($) with exactly one
  sub(";*$", ";", text)
}
tax_table_3 <- tax_table_2 |>
  unite(col = "lineage", Kingdom:Species, sep = ";") |>
  select(OTU:lineage) |>
  mutate(lineage=fix_trailing_semicolon(lineage))

# handle the OTU table
otu_table_0 <- read_tsv("otutable.txt")
# add the new taxonomy
otu_table_1 <- otu_table_0 |>
  left_join(select(tax_table_3, OTU, `#NAME` = lineage)) |>
  select(-OTU) |>
  relocate(last_col(), .before = 1)

write_tsv(otu_table_1, "otu_table_microbiome_analyst.txt")

# handle the metadata
sample_table_0 <- read_tsv("sampletable.txt")
sample_table_1 <- sample_table_0 |> rename(`#NAME`= NAME)
write_tsv(sample_table_1, "sample_table_microbiome_analyst.txt")


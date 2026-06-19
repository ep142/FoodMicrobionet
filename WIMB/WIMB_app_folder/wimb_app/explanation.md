### Where Is My Bug?

This app queries [FoodMicrobionet](https://github.com/ep142/FoodMicrobionet)
to characterise the distribution of selected microbial taxa across food
categories, using culture-independent (metataxonomic) data from published
studies.

**How to use**

1. Choose a **taxonomic level** to search at (genus or family) and the level
   to **go down to** for plotting.
2. Select one or more **taxa** — type to search, the list updates
   automatically. Up to 12 taxa are supported.
3. Choose a **plot type** and any relevant options.
4. Press **Run** to generate the outputs.

**What the outputs show**

- **Distribution plot**: box or violin plots of log-transformed relative
  abundance across food categories, with optional prevalence colour coding
  and individual sample jitter points.
- **Dot plot**: a matrix of median abundance (point size) vs. relative
  prevalence (fill colour) across food categories and taxa.
- **Bibliography**: the FoodMicrobionet studies contributing data for the
  selected taxa, with DOI links.

**Important caveats**

1. Metataxonomic data are **compositional**: relative abundances within a sample
sum to 1 regardless of the total microbial load, which differs across studies
and food types. FoodMicrobionet aggregates hundreds of studies with variable
sequencing depth (from a few hundred to >100,000 reads per sample) and
different gene targets. A minimum sequencing depth filter is applied to reduce
detection bias. Results should be interpreted as a qualitative ecological
fingerprint rather than precise quantitative estimates.

2. FoodMicrobionet does not represent a random sample of foods: some food
categories are better represented than others, reflecting publication bias.  

3. When you select a large number of taxa boxplots may be crowded; try increasing the size of the saved graph.

**Acknowledgements**  

The app was entirely developed by Claude AI (Sonnet 4.6) with my prompts. Claude was also used to revise the app. 

**Citation**

If you use this app please cite the last version of FoodMicrobionet: 
Parente, E., Ricciardi, A., 2024. A Comprehensive View of Food Microbiota: 
Introducing FoodMicrobionet v5. Foods 13, 1689. [https://doi.org/10.3390/foods13111689]

**Copyright** 

This is version 1 (2026-05-26).

Assume that this is overall under MIT licence.

Copyright 2026 Eugenio Parente.  
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


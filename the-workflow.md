# The workflow

We have written a workflow based on snakemake to cover the end-to-end analysis. Users can provide the raw methylation files and the sample sheet as a starting point, and the pipeline will automatically:

1. Preprocess the raw files to extract the beta values with [`minfi`](https://doi.org/doi:10.18129/B9.bioc.minfi).
2. If reference samples available, call CNVs with [conumee](https://doi.org/doi:10.18129/B9.bioc.conumee) to blacklist potentially disruptive CpGs.
3. Call fluctuating CpG sites (If no fCpG site list is provided).
4. Create the XML file required by BEAST (alternative starting point if beta-values are already given - _trees mode_).
5. Run [`BEAST`](https://beast.community/).
6. Automatic model selection
7. Summarise results and draw the phylogenetic trees.

<figure><img src=".gitbook/assets/imagen.png" alt=""><figcaption><p>PHYFUM's workflow</p></figcaption></figure>




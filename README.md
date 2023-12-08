---
layout:
  title:
    visible: true
  description:
    visible: true
  tableOfContents:
    visible: true
  outline:
    visible: true
  pagination:
    visible: false
---

# Phyfum

Phyfum is a tool for inferring phylogenetic trees on methylation-based studies. We harness fluctuating CpG (fCpG) sites of methylation arrays to study the clonal evolution of samples. You can read more about fCpGs in the [original paper](https://www.nature.com/articles/s41587-021-01109-w).&#x20;

We have implemented the original model described in the paper within BEAST v.1.8.4, which is accesible [here](https://github.com/pbousquets/beast-mcmc-flipflop). In order to make it more accesible, we have designed a snakemake-based pipeline, covering the IDAT preprocessing, fCpG calling, automatic XML generation and BEAST inference. Additionally, if both tumor and normal samples are available, CNVs are called to curate non-fluctuating CpGs.

Before you start using the tool, we recommend you reading [the-workflow.md](the-workflow.md "mention") section to fully understand how our Phyfum works.

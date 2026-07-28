# AGENTS.md

# Repository description
- This repository contains a fork of the phylogenetic software BEAST 1.8.3 with extended models to study somatic evolution using fluctuating methylation. This program is called PHYFUM.
- Most of the codebase is JAVA and none of the new code uses the BEAGLE library.
- The new development starts after tag: v1.8.4.

# Model description
- PHYFUM models the somatic evolution of tissues organized in crypts using fluctuating methylation clocks by using a continuous time markov chain model that models the methylation status of S stem cells as they methylate, demethylate, and replace each other.
- At internal nodes, 4 different likelihood calculation strategies model different biological models of crypt divission. They exend the GeneralCenancestorLikelihoodCore class.

# Coding instructions
- Codebase last modified before tag: v1.8.4 should only be read when it is explicitly necessary
- Write human-readable code and comment it.
- Minimize existing code edits that are not strictly necessary. Ask for my input on improvements that are not necessary but may be beneficial.
- Commit changes after each interaction with me that modifies the code. Write detailed commit messages.
- Save reports in dev/reports/

# Running
- The main executable class is dr.app.beast.BeastMain
- It requires the -beagle_off argument to work
- The -seed argument allows specifying a random generator seed
- Examples/release/flipflop/ contains runnable xml PHYFUM examples

# Other resources
- You can find a sequence simulator under the PHYFUM model in the simulate folder
- You can find a manuscript (under review) presenting the model here: ~/projects/flipFlop/manuscript/PhyfumV3GenomeBiology.docx

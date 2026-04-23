# Statistical Limits of iPSC Study Designs

This repository contains a statistical study that demonstrates the severe limitations of typical small-sample (N=2 to N=10 donors) study designs in induced pluripotent stem cell (iPSC) and organoid research. 

Through a combination of Monte Carlo simulations, exact analytical proofs, and empirical validation on real-world datasets, this codebase evaluates why the field's recommendation of using "3 to 5 donors" is mathematically and biologically insufficient for capturing genetic variance or generalizing treatment effects.

## Repository Structure

The analysis is divided into four main R scripts, each tackling the problem from a distinct methodological angle:

### `01_MCsimulation.R`
**Monte Carlo Simulations of Statistical Power and Error Rates**
* **Statistical Power & ICC:** Evaluates the power to detect donor-level genetic effects across an empirically grounded Intraclass Correlation (ICC) range (5% - 50%).
* **Type M & Type S Errors:** Applies the Gelman & Carlin (2014) framework to demonstrate that significant findings at small N are overwhelmingly likely to be massive exaggerations of true effect sizes (Type M error).
* **Batch-Donor Confounding:** Simulates hierarchical structures (donors -> batches -> organoids) to show that one-way ANOVAs artificially inflate apparent donor variance, and nested ANOVAs fail to separate these components at low N.

### `02_analytical_methods.R`
**Closed-Form Analytical Proofs (No Simulation Assumptions)**
* **Design Effect (Kish, 1965):** Calculates the effective sample size, proving mathematically that adding more organoids per donor yields sharply diminishing returns.
* **Population Genetics Variant Capture:** Calculates the probability of observing minor alleles at various frequencies, demonstrating that small panels miss the vast majority of rare and low-frequency genetic variation.
* **Exact Confidence Intervals:** Uses the Thomas & Hultquist (1978) exact F-distribution method to calculate the theoretical minimum confidence interval widths for ICC, proving that N=3 yields effectively uniform uncertainty bounds.

### `03_validation_real_data.R`
**Empirical Validation using Real RNA-seq Data**
* Utilizes the HipSci iPSC processed RNA-seq dataset (Mirauta et al., 2020; 202 lines, 151 donors).
* Computes true gene-wise ICCs from the full dataset and repeatedly subsamples to typical study sizes (e.g., N=3, N=5, N=10 etc.).
* Evaluates **misclassification rates**: how often a small sample falsely concludes a highly genetic trait is non-genetic, or vice versa.

### `04_donor_by_treatment.R`
**Two-Way Mixed ANOVA for Treatment Interactions**
* Simulates the typical paradigm of testing a treatment across multiple donor lines.
* Demonstrates the inability of small-N studies to distinguish a true, generalizable main treatment effect from underlying donor-by-treatment interaction (heterogeneity).
* Includes a direction-concordance diagnostic showing the high probability that randomly sampled small cohorts will coincidentally respond in the same direction, granting a false sense of generalizability.

## Requirements

The code is written entirely in R and relies on the following standard libraries for data manipulation and visualization:
* `ggplot2`
* `dplyr`
* `tidyr`
* `patchwork`
* `scales`

## Usage

Each script is standalone and can be run interactively or sourced directly in R after adjustment of input/output paths.

*Note: `03_validation_real_data.R` expects local data files (`elife-57390-fig1-data1-v3.xls` and `elife-57390-fig1-data2-v3.xls`) downloaded from Mirauta et al., eLife, 2020. Update the `META_FILE` and `EXPR_FILE` path variables at the top of the script if necessary.*

## Key Takeaways
1. **Power:** Studies with 2-5 donors have less than 37% power to detect a donor-level genetic effect when genetics explains 20% of phenotypic variance. Power drops below 12% at ICC = 5%.
2. **Precision of estimates:** ICC estimates at n = 2-5 are effectively uninformative, with 95% simulation intervals spanning from 0% to 60–80%, regardless of the true value.
3. **'Addressing genetic background' with N=3-5:** Meaningful characterisation of the genetic contribution to phenotypic variance in iPSC-derived organoid models requires donor sample sizes on the order of those used in large-scale studies (n >= 50), consistent with the designs that originally established the ICC estimates we used here. The claim that 3-5 donors "address" or "control for" genetic heterogeneity is not supported by any of these analyses.

## Conclusion
The field's standing recommendation of 3-5 donors is broadly compatible with simulation results that address a single, narrow question: under the assumption that every donor responds identically, how many donors are required to detect a common mean effect? Under this assumption N=3-5 (closer to 5 than 3) is appropriate. The problem arises when this recommendation is extended to questions it was never designed to answer. If the assumption is violated and donors respond differently, N=3-5 is insufficient. Further, the same recommendation is often paired with the practical concern often raised in the field: studies must demonstrate their results generalise across genetic backgrounds. What is overlooked is that this is a fundamentally different statistical question with different answer. With 3-5 donors, the experiment doesn't even sample a substantial fraction of the relevant variants. The question is merely not testable from such data, it requires either much larger samples (n > 50 to capture most of common variants) or external evidence (prior GWAS or eQTL).


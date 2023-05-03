# Uncertainty Propagation in <sup>1</sup>H-MRS
Welcome! This repository replicates the primary analysis from work on **Uncertainty propagation in absolute metabolite quantification for in vivo magnetic resonance spectroscopy of the human brain.** Here, you'll find scripts, test data, and relevant project source code.

**Corresponding Author**: Ronald Instrella, M.S.<br>
**Principal Investigator**: Christoph Juchem, Ph.D.

## Overview
The effect of measurement errors from experimentally-derived parameters in MRS absolute quantification are examined using both analytical derivations and Monte Carlo simulations. The following scripts and utility functions allow users to reproduce published figures, and run the analysis using custom values in the provided JSON configurations.

## Getting Started

Simply invoke the main analysis script to replicate the analysis and produce figures presented in [1].

```
>> CmErrorAnalysis
```

Individual analysis scripts are found in ./analysis, and shared utilities are found in ./utils, including partial derivatives. The employed constants and standard deviations of quantification parameters are defined as JSON files in ./data.

For more detailed information on each script, function or class, please use the "help" command in MATLAB directly, and each file's header comments:

```
>> help CmErrorAnalysis
>> help CmErrorByParameter
```

## Contact
If you have any comments or questions, please email me directly: ronald.instrella@columbia.edu.

## Citation
The primary citation for this work is:

[1] Instrella R, Juchem C. Uncertainty propagation in absolute metabolite quantification
for in vivo magnetic resonance spectroscopy of the human brain. *Magn Reson Med.* 2023 (submitted)

## Selected References

[2] Swanberg KM, Prinsen H, DeStefano K, et al. In vivo evidence of differential frontal cortex metabolic abnormalities in progressive and relapsing-remitting multiple sclerosis. *NMR Biomed.* 2021;34(11). doi:10.1002/nbm.4590

[3] Landheer K, Gajdošík M, Treacy M, Juchem C. Concentration and effective T2 relaxation times of macromolecules at 3T. *Magn Reson Med.* 2020;84(5):2327-2337. doi:10.1002/mrm.28282

# Interference suppression tutorials for OPM data

Tutorial-style scripts for the paper: 
*Seymour et al., (2022). Interference suppression techniques for OPM-based MEG:  Opportunities and challenges. Neuroimage*

## Tutorials:

### 1. Auditory evoked field paradigm during participant movement

- Data: [https://doi.org/10.5281/zenodo.5539414](https://doi.org/10.5281/zenodo.5539414)
- Script: [pipeline_tutorial_1.m](./pipeline_tutorial_1.m)
- *arft.m* : indices of manually-identified bad time segments

### 2. Motor-beta power changes during a finger-tapping paradigm

- Data: [https://doi.org/10.5281/zenodo.5539414](https://doi.org/10.5281/zenodo.5539414)
- Script: [pipeline_tutorial_2.m](./pipeline_tutorial_2.m)

### 3. Auditory evoked field paradigm when participant is sitting still

- Data: Seymour, R. (2022, October 24). WCHN OPM Example Data. Retrieved from [osf.io/vy5n6](osf.io/vy5n6)
- Script: [stationary_AEF.m](./stationary_AEF.m)

## Code Dependencies:

Please download all code dependencies from **Zenodo: https://doi.org/10.5281/zenodo.6599496**. This will allow you to exactly reproduce the analyses presented in the article.

Alternatively, you can download the latest versions of the toolboxes from:
- Fieldtrip:      https://www.fieldtriptoolbox.org/download/
- analyse_OPMEG:  https://github.com/neurofractal/analyse_OPMEG
- optitrack:      https://github.com/FIL-OPMEG/optitrack
- ft_denoise_hfc: email rob.seymour@ucl.ac.uk

## License:

These scripts are shared under the following license:
**Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)**.
Please see **[HERE](https://creativecommons.org/licenses/by-sa/4.0/)** for more information

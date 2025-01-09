# PCAPAM50 1.0.3

## Changes in version 1.03
- the term 'misclassification' has been substituted with 'discordance' for a output diagnostic plot that highlights IHC discordance cases, with PC1-separated cases displayed along the PC1 axis.
- The vignette has been updated to reflect these changes, including the revised figure.
- **Warning Fix**: Resolved an issue where calling `par(no.readonly = TRUE)` in functions created an empty plot if no graphics device was active.
  - Now, the functions check for an active graphics device using `if (!is.null(dev.list()))` before saving and restoring graphical parameters.
  - `on.exit(par(oldpar))` is now conditionally applied only when a device is active, ensuring unintended side effects.
  - Added informative messaging (`message()`) when no active graphics device is detected, making the behavior transparent to users.


## Changes in version 1.0.2
- **Renamed Parameter**: Renamed `inputDir` to `outDir` for clarity and alignment with its intended purpose of storing output files.
- **Heatmap Enhancements**: The height and width of the output heatmap are now dynamically adjusted based on the input sample size. This ensures better visualization, particularly for larger datasets.
- **Citation**: Added citation information.

Repository for analysis scripts and processed data for the paper:
"Regional endemism instead of unlimited gene flow: Phylogeography of nivicolous myxomycetes in the mountains of Eurasia" (Ecological Monographs)

Raw and processed data

- Raw data (sequences, sample metadata, environmental layers):
  - Currently deposited on Zenodo as a review/limited-access record (no public DOI yet).
  - DOI will be added here and in CITATION.cff upon publication.
  - If you need access for review, contact the corresponding author or request the Zenodo preview link (preview links contain tokens and are not posted publicly).
  - Note: large files (e.g. rasters) are hosted on Zenodo and not included in this repository. Checksums and full file list will be available in the Zenodo record.

- Processed input files included in this repository:
  - _data/processed_data/
    - singleton176_mafft.fas
    - unique473_mafft.fas
  - These are the files only used for analyses with python script included.

Reproduce the analyses
1. Place raw data under `data/` and processed inputs under `_data/processed_data/`.
2. Open R in the project root (where `renv.lock` is) and restore environment:
   - `install.packages("renv")` (if needed)
   - `renv::restore()`
3. Run scripts in order (exact filenames):
   - `_script/01_data_cleaning_processing.R`
   - `_script/02_basic_statistics.R`
   - `_script/03_geographic_differentiation.R`
   - `_script/04_haplotype_beta_diversity_analysis.R`
   - `_script/05_species_ecological_analysis.R`
Utility functions: `_script/utils_ecology_myxo.R`
4. Run these scripts at any time to analse the seqeunce data with python:
   - `_script/Double peak detection.py`
   - `_script/Single_substitution_detect_among_singleton.py`

How to cite
- Paper: Inoue et al. 2025, "Regional endemism instead of unlimited gene flow..." Ecological Monographs.
- Repository (Zenodo DOI will be added after release): https://doi.org/10.5281/zenodo.YOUR_ID

License
- Code: MIT (see LICENSE)
- Data: CC-BY-4.0 (see DATA_LICENSE)

Zenodo / release notes
- Create a GitHub release (tag) and Zenodo will mint a DOI. Example:
  ```
  git tag -a v1.0.0 -m "Initial release"
  git push origin v1.0.0
  gh release create v1.0.0 --title "v1.0.0" --notes "Initial release: scripts + renv.lock"
  ```



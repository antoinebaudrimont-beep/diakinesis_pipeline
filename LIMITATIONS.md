# Validation and limitations

This repository is a research analysis tool. Its output requires image-level review before biological interpretation. The historical FIJI/PartSeg measurements remain the quantitative reference; they are not bundled or required at runtime.

## What has been validated

- The prepared LMN mask operation sequence was reconstructed from saved FIJI masks: seven slice-wise dilations, with a slice-wise hole fill after each of the final six dilations. This is tested against historical fixtures when they are available.
- Using saved historical masks, the Python centroid-to-border calculation reproduced FIJI `cen-bor` to floating-point precision. This validates the measurement calculation separately from raw-image segmentation.
- A related 44-stack historical benchmark of the default PartSeg-equivalent DAPI masks found 242 one-to-one matches among 269 eligible historical bodies, with six misses, 21 merged bodies, four over-splits, and 37 extra detections. These figures describe the benchmarked masks, not a guarantee for new experiments.
- A later targeted watershed experiment improved matches to 257, but also raised over-splits from four to seven. It is deliberately **not** part of this package's automatic pipeline. The historical single-body volume rule used in that experiment is not independently validated on new data.

## Limits of the automated workflow

- The complete raw-DV-to-result workflow does not reproduce every historical object or measurement within 5%. DAPI and LMN segmentation, nuclear selection, and gating can all differ from historical manual decisions. Passing the saved-mask measurement test does not establish end-to-end equivalence.
- The program selects one central LMN component per stack. It does not automatically measure every nucleus in a multi-nucleus stack. Check the selected nucleus in each QC image, particularly when LMN selection confidence is below `high`.
- The current `tolerant` DAPI gate includes some objects with centroids just outside the selected LMN when their geometry is near the border. It is an automated candidate-selection rule, not a confirmed reconstruction of every historical manual inclusion decision. `measurements_all.csv` retains excluded and borderline objects for review.
- The same DAPI Moments-bin setting is used for an entire run. There is no automatic per-image optimization. Threshold settings may need expert review for a new acquisition protocol; changing them invalidates direct comparability until revalidated.
- The DV reader supports the audited layout only: two channel blocks, one time point, little-endian float32 image data, embedded emission wavelengths, and physical voxel sizes. Other DeltaVision variants may need a different reader.
- The channel wavelengths in `config.yaml` describe the original dataset. They must be checked against metadata for a new experiment. The default uses metadata voxel sizes; `legacy` calibration is only for comparison with historical FIJI output.
- The optional per-DAPI-body Mann–Whitney test treats bodies as independent observations. Bodies within an image and images within an experiment are nested, so that p-value is exploratory. A biological comparison should use the independent biological replicate as the unit of inference.

## Data and reproducibility

Raw DV files, historical masks and CSVs, and generated example result masks are intentionally excluded from the public repository. The `examples/` tables are synthetic and demonstrate file format and plotting only. Historical regression tests skip unless `DIAKINESIS_REFERENCE_ROOT` and `DIAKINESIS_VALIDATION_ROOT` point to the separate local reference dataset and its mapping file.

The MIT license covers the code and documentation in this repository; it does not grant rights to unpublished microscopy data that are not included here.

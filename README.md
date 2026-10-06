# Diakinesis DAPI–LMN distance pipeline

Deterministic Python analysis of deconvolved, two-channel DeltaVision 3D stacks. For each retained DAPI body, it reports **the distance from the body's voxel-center centroid to the nearest LMN-1 border voxel center in 3D**. The historical FIJI field name is `cen-bor`, and distances are in **µm**. This is a centroid-to-border measurement, not a surface-to-surface distance.

**Status:** research software. Saved-mask distance and LMN morphology steps have strong historical regression checks, but the complete automated raw-image workflow does not reproduce every historical object or measurement. Review segmentation and gating QC before biological interpretation. See [validation and limitations](LIMITATIONS.md).

## Install

Create the Conda environment from the repository root:

```sh
conda env create -f environment.yml
conda activate diakinesis
```

The development environment used Python 3.10, NumPy 1.26, SciPy 1.15, SimpleITK 2.4.1, tifffile 2025.1.10, and Matplotlib 3.10. `environment.yml` declares the direct dependencies. Each run records the installed versions in `analysis_log.txt` and `logs/environment.txt`.

## Analyze a new experiment

1. Put deconvolved `*_R3D_D3D.dv` files in a folder outside this repository. The current program selects **one central LMN nucleus per stack**.
2. Inspect the DV metadata and edit [config.yaml](config.yaml) as needed. The included channel targets (435 nm DAPI, 679 nm LMN-1) describe the original dataset, not a universal channel order. The program checks the embedded wavelengths and refuses missing or ambiguous channels.
3. Choose a DAPI Moments-bin value for the whole run. The default is 20; `--dapi-bins` overrides it. There is no automatic per-image bin selection.
4. Run the analysis into a **new, empty** output folder:

```sh
python run_diakinesis.py \
  --input /path/to/dv_files \
  --output /path/to/new_results \
  --condition auxin \
  --dapi-bins 20
```

If `--input` is omitted, a folder chooser opens. If `--output` is omitted, results go to `<input>/diakinesis_results/`; the program refuses to overwrite a nonempty output folder. If `--condition` is omitted, the input folder name is used. Use `python run_diakinesis.py --help` for all options.

The default calibration reads physical voxel sizes from each DV header. For **historical FIJI comparison only**, `--calibration legacy` uses Z/Y/X = 0.2000/0.0643/0.0643 µm. Audited original DV metadata were approximately 0.200000/0.064346/0.064346 µm. The metadata and used sizes are recorded in the output.

## Processing and selection

- DAPI: Gaussian filtering, whole-stack Moments thresholding, 200-voxel minimum object size, and 26-connected 3D labels. The default does not fill 3D holes.
- LMN-1: Gaussian filtering, Moments thresholding, and 70,000-voxel minimum size. Preparation reproduces the saved FIJI operation sequence: **seven slice-wise 3×3 dilations, with a slice-wise hole fill after each of the final six dilations**.
- Nuclear selection: choose the prepared LMN component nearest the image's XY center, with size as a tie-breaker. Multiple plausible components lower selection confidence.
- DAPI gating: the default `tolerant` rule retains objects whose rounded centroids lie inside the selected LMN and certain near-border objects based on calibrated geometry. It labels detections `inside`, `borderline_retained`, or `outside_excluded`. This is an automated candidate gate; inspect the QC view and review borderline cases.
- Distance: measure each retained object's 3D centroid to the nearest selected LMN border voxel center. The saved-mask distance implementation has been checked separately against historical FIJI values.

All parameters are in [config.yaml](config.yaml). The selected LMN nucleus and accepted DAPI objects should be checked for every new acquisition type. The experimental targeted watershed method is **not** in this automatic pipeline because it created false splits in the historical benchmark.

## Results and QC

| File | Contents |
|:--|:--|
| `measurements_all.csv` | Every detected DAPI object, including excluded and borderline objects, with coordinates, selection reason, calibration, and `cen_bor_um`. |
| `measurements_selected.csv` | Objects retained by the automatic gate. |
| `summary_long.csv` | Compatibility columns `filename,Label1,Label2,cen-bor`. |
| `tables/stack_qc.csv` | Per-stack object counts, LMN selection confidence, and review flags. |
| `tables/input_inventory.csv` | Channel mapping and physical voxel sizes from the DV header. |
| `qc/*_qc.png` | DAPI projection with selected LMN border, individually colored retained bodies, and red excluded bodies. |
| `masks/` | DAPI labels, LMN labels, and selected prepared LMN mask for traceability. |

Inspect every `needs_review` stack, especially uncertain LMN selection, unusual body counts, and borderline gating. `measurements_all.csv` contains blank `manual_selected` and `manual_comment` columns for a separate reviewed copy. Do not erase the original automatic decisions when making manual edits.

## Example plot

[Examples](examples/README.md) includes **synthetic CSVs** and a runnable plotting command. For real results from two separately processed conditions:

```sh
python plot_conditions.py \
  /path/to/auxin/measurements_selected.csv \
  /path/to/ethanol/measurements_selected.csv \
  --output condition_comparison.png
```

The figure shows one point per DAPI body and mean ± SD for each condition. With two conditions it also reports a two-sided Mann–Whitney U test. Bodies within a stack are nested observations, so that per-body test is exploratory; use biological replicates for confirmatory inference.

## Tests

```sh
python -m unittest discover -s tests -v
```

The portable IO tests run without historical data. Historical morphology, distance, and clean-stack checks run only if both `DIAKINESIS_REFERENCE_ROOT` and `DIAKINESIS_VALIDATION_ROOT` point to the separate local reference dataset and validation workspace. These fixtures are not published here and are never runtime inputs to `run_diakinesis.py`.

## License

Code and documentation are released under the [MIT License](LICENSE). Raw microscopy images and historical reference outputs are not included.

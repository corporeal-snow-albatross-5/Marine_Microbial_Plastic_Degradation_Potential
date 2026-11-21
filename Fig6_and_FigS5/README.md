#  Cα backbone SVG plotting

Source: RCSB Pairwise Structure Alignment CIF exports. These CIFs have already been transformed to share an aligned coordinate frame within each **pair**. Do not mix coordinate files between different pair directories.

## Included CSV manifests

- `manifest_all10_full_length.csv` — all 10 full-length reference–marine homolog comparisons from Table S6, for Figure S5.
- `manifest_figure6_main5_full_length.csv` — five representative full-length pairs (archive pairs 2, 6, 7, 8, 10).
- `manifest_CDD_domain5.csv` — five domain-level comparisons (MCL-PHA pairs 5/6, PHB 8/9, P3HB 10), with residue selections exactly as reported in Table S6. This is **not** the same as five distinct families.

The manifest column `reference` is the coordinate file path, not a reference accession. The `homolog` column is the corresponding coordinate file path. Both paths are relative to the CSV directory. Blank residue ranges mean all Cα residues. Ranges such as `53:232` include both endpoints. All coordinate files have chain ID `A`.

## Use

Install Python dependencies:

    python -m pip install biopython numpy matplotlib

From the extracted project folder, run:

    python plot_ca_backbone_overlays.py --manifest manifest_all10_full_length.csv --output Figure_S5_regenerated.svg --columns 5 --elev 0 --azim 0
    python plot_ca_backbone_overlays.py --manifest manifest_figure6_main5_full_length.csv --output Figure_6_full_length_regenerated.svg --columns 5 --elev 0 --azim 0
    python plot_ca_backbone_overlays.py --manifest manifest_CDD_domain5.csv --output Figure_domains_regenerated.svg --columns 5 --elev 0 --azim 0

Both CIF files from one numbered archive are already transformed and should be plotted together, without applying a new TM-align transformation.

## Documentation notes

This script is a documented reconstruction, not the recovered original code for the previously published drafts. Original SVG metadata identifies Matplotlib v3.10.8. The exact earlier viewing orientation is not recoverable from that metadata, so regenerated SVGs may be rotated differently. The current code creates orthographic XY projections at `--elev 0 --azim 0` and permits different common camera angles. It does not compute TM-scores/RMSD and does not establish biological activity.

The _1 alternate ZIP archives were not used: primary numbered exports were selected for consistent replication of Table S6's ten pairs. Sequence identities were assigned by matching chain lengths to accessions in Table S5 and pair order in Table S6; the hashed CIF filenames do not themselves contain accessions.

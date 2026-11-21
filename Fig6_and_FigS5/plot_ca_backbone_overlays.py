#!/usr/bin/env python3
"""Reproducible Matplotlib C-alpha trace visualization of ALREADY ALIGNED proteins.

Input: a CSV manifest containing one row per overlay, with columns
    panel,reference,homolog,reference_chain,homolog_chain,reference_range,homolog_range
All except panel/reference/homolog are optional. PDB or mmCIF supported.
The reference and homolog coordinate files MUST already share a superimposed
coordinate frame, e.g. transformed coordinate exports from a TM-align viewer.
This script DOES NOT run TM-align or compute TM-score/RMSD.

Usage:
    python plot_ca_backbone_overlays.py --manifest pairs.csv --output figure.svg \
      --columns 5 --elev 0 --azim 0

Ranges may be blank, "18:210", or "18:210;240:290" (inclusive residue ids).
Requires: biopython, numpy, matplotlib
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.PDB import PDBParser, MMCIFParser

ORANGE = '#d27845'
TEAL = '#168fa4'


def parse_range(s: str):
    if not s.strip():
        return None
    chunks = []
    for token in s.split(';'):
        lo, hi = token.split(':', 1)
        chunks.append((int(lo), int(hi)))
    return chunks


def get_atoms(path: str, chain_name: str = '', residue_range: str = ''):
    file = Path(path)
    if not file.is_file():
        raise FileNotFoundError(f'Missing coordinate file: {file}')
    if file.suffix.lower() in ('.cif', '.mmcif'):
        parser = MMCIFParser(QUIET=True)
    elif file.suffix.lower() in ('.pdb', '.ent'):
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f'Expected PDB/CIF: {file}')
    model = next(parser.get_structure(file.stem, str(file)).get_models())
    chains = list(model.get_chains())
    if not chains:
        raise ValueError(f'No chains: {file}')
    if chain_name:
        chain = model[chain_name]
    else:
        if len(chains) > 1:
            raise ValueError(f'{file} has {len(chains)} chains; specify a chain in the manifest')
        chain = chains[0]
    intervals = parse_range(residue_range)
    coords, residue_ids = [], []
    for residue in chain:
        if residue.id[0] != ' ' or 'CA' not in residue:
            continue
        pos = residue.id[1]
        if intervals and not any(lo <= pos <= hi for lo, hi in intervals):
            continue
        atom = residue['CA']
        coords.append(atom.get_coord().astype(float))
        residue_ids.append(pos)
    if len(coords) < 2:
        raise ValueError(f'Fewer than two C-alpha coordinates: {file} chain {chain.id}')
    return np.array(coords), np.array(residue_ids)


def view_basis(elev_degrees: float, azim_degrees: float):
    """Construct one common camera frame for all reference/homolog coordinates."""
    elev, azim = np.radians([elev_degrees, azim_degrees])
    # Orthographic camera: elev=azim=0 projects x/y; same camera for both proteins.
    right = np.array([np.cos(azim), np.sin(azim), 0.0])
    up = np.array([-np.cos(elev)*np.sin(azim),
                   np.cos(elev)*np.cos(azim), np.sin(elev)])
    return right, up


def project(coordinates, right, up):
    return np.column_stack([coordinates @ right, coordinates @ up])


def plot_chain(ax, projected, residue_ids, color):
    # Avoid drawing artificial bonds across missing/nonconsecutive residues.
    boundaries = np.flatnonzero(np.diff(residue_ids) != 1) + 1
    chunks = np.split(np.arange(len(residue_ids)), boundaries)
    for idx in chunks:
        if len(idx) < 2:
            continue
        xy = projected[idx]
        # Soft white outline then colored trace, as seen in the supplied SVGs.
        ax.plot(xy[:,0], xy[:,1], color='white', alpha=.38,
                linewidth=2.4, solid_capstyle='round', zorder=1)
        ax.plot(xy[:,0], xy[:,1], color=color, alpha=.80,
                linewidth=1.65, solid_capstyle='round', zorder=2)


def main():
    cli = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    cli.add_argument('--manifest', required=True)
    cli.add_argument('--output', required=True)
    cli.add_argument('--columns', type=int, default=5)
    cli.add_argument('--elev', type=float, default=0,
                     help='View elevation in degrees; default 0 gives an XY projection')
    cli.add_argument('--azim', type=float, default=0)
    cli.add_argument('--panel-width', type=float, default=3.2)
    cli.add_argument('--panel-height', type=float, default=3.2)
    cli.add_argument('--title', default='')
    args = cli.parse_args()
    if args.columns < 1:
        cli.error('--columns must be >= 1')
    manifest_path = Path(args.manifest).resolve()
    with open(manifest_path, newline='') as handle:
        entries = list(csv.DictReader(handle))
    for row in entries:
        for key in ('reference', 'homolog'):
            source = Path(row[key])
            if not source.is_absolute():
                row[key] = str(manifest_path.parent / source)
    if not entries:
        raise ValueError('Manifest is empty')
    for row in entries:
        for key in ('panel','reference','homolog'):
            if not row.get(key):
                raise ValueError(f'Missing {key} in manifest row: {row}')
    right, up = view_basis(args.elev, args.azim)
    ncols = min(args.columns, len(entries))
    nrows = (len(entries) + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows, ncols,
             figsize=(ncols*args.panel_width, nrows*args.panel_height),
             squeeze=False)
    for ax in axs.flat:
        ax.axis('off')
        ax.set_aspect('equal', adjustable='box')
    for i, row in enumerate(entries):
        ax = axs.flat[i]
        reference, rids = get_atoms(row['reference'], row.get('reference_chain',''),
                                      row.get('reference_range',''))
        homolog, hids = get_atoms(row['homolog'], row.get('homolog_chain',''),
                                   row.get('homolog_range',''))
        # IMPORTANT: no fitting occurs here. Both files must already be aligned.
        xy_r = project(reference, right, up)
        xy_h = project(homolog, right, up)
        plot_chain(ax, xy_r, rids, ORANGE)
        plot_chain(ax, xy_h, hids, TEAL)
        combined = np.concatenate((xy_r, xy_h))
        lo = combined.min(axis=0)
        hi = combined.max(axis=0)
        center = (lo + hi)/2
        span = max(float((hi-lo).max()) * 1.12, 1.0)
        ax.set_xlim(center[0] - span/2, center[0] + span/2)
        ax.set_ylim(center[1] - span/2, center[1] + span/2)
        ax.set_title(row['panel'], fontsize=10, pad=7)
    handles = [Line2D([0],[0], color=ORANGE, lw=2, label='Reference'),
               Line2D([0],[0], color=TEAL, lw=2, label='Marine homolog')]
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
               bbox_to_anchor=(.5,.005), fontsize=9)
    if args.title:
        fig.suptitle(args.title, fontsize=13)
    fig.subplots_adjust(left=.02,right=.98,top=.91,bottom=.09,hspace=.16,wspace=.10)
    plt.rcParams['svg.fonttype'] = 'none'
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format='svg', facecolor='white', bbox_inches='tight')
    print(f'Wrote {out} | Matplotlib {matplotlib.__version__}')

if __name__ == '__main__':
    main()

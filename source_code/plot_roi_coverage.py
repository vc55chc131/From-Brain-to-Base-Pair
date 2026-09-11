#!/usr/bin/env python3
"""Plot actual corrected-coordinate coverage distances; no expression inference."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('data/feasibility/results/nearest_samples.csv'))
    parser.add_argument('--output', type=Path, default=Path('outputs/coverage_figure.png'))
    args = parser.parse_args()
    with args.input.open() as f:
        rows = list(csv.DictReader(f))
    values = [min(float(r['distance_to_seed_mm']) for r in rows
                  if r['coordinate_set'] == 'corrected_mni' and int(r['roi_id']) == i)
              for i in (1, 2, 3)]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, ax = plt.subplots(figsize=(7, 3.15), layout='constrained')
    ax.axvspan(0, 2, color='#F0F0F0', zorder=0)
    ax.axvline(2, color='#555555', ls='--', lw=1.1, zorder=1)
    for y, distance, color in zip([2, 1, 0], values, ['#005A9C', '#B2182B', '#B2182B']):
        ax.plot([0, distance], [y, y], color='#D6D6D6', lw=1, zorder=1)
        ax.scatter(distance, y, marker='x', s=115, color=color, linewidths=2, zorder=3)
        ax.annotate(f'{distance:.3f}', (distance, y), xytext=(0, 11),
                    textcoords='offset points', ha='center', color=color, fontsize=10)
    ax.set_yticks([2, 1, 0], ['Left DLPFC target', 'Left caudate target', 'Right caudate target'])
    ax.set_xlim(0, max(10.3, max(values) + 0.9))
    ax.set_ylim(-0.65, 2.6)
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    ax.set_xlabel('Distance from target center to nearest sample center (mm)')
    ax.text(1, -0.53, '2 mm radius', ha='center', fontsize=9, color='#555555')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color('#777777')
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', ls=':', color='#D9D9D9')
    ax.set_axisbelow(True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    main()

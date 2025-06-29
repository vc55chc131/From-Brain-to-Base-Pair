# 🧠 From Brain to Base Pair: Transcriptomic Analysis of Language Control

This repository provides all necessary materials for reproducing the transcriptomic analysis presented in our study:  
**"From Brain to Base Pair: A Dual-Site Transcriptomic Analysis of Language Control in Simultaneous Interpreting."**

The project links functional neuroimaging findings with molecular neuroscience by comparing gene expression profiles between the **dorsolateral prefrontal cortex (DLPFC)** and **caudate nucleus**, two regions critical to simultaneous interpreting cognition.

All analysis scripts used to generate the results and figures in this paper are fully available in the `source_code/` directory of this repository.

---

## 📂 Repository Structure

| Folder / File Name                        | Description                                                                 |
|-------------------------------------------|-----------------------------------------------------------------------------|
| `ROIs/`                                   | MNI coordinates and voxel selection scripts for DLPFC and caudate regions  |
| `Data Sources/`                           | Metadata and links to Allen Human Brain Atlas transcriptome files          |
| `Main packages/`                          | Package list and installation environment (`requirements.txt`)             |
| `Output Highlights/`                      | Visualizations and figure outputs (e.g., volcano plots, gene rankings)     |
| `Reproducibility and FAIR Compliance/`    | FAIR checklist and code reproducibility documentation                      |
| `source_code/`                            | Main pipeline script for Z-score comparison and output visualisation       |
| `Citation` / `Contact`                    | Citation info and project contact email                                    |

---

## 🚀 How to Reproduce

1. **Clone the repository**  
```bash
git clone https://github.com/vc55chc131/From-Brain-to-Base-Pair.git
cd From-Brain-to-Base-Pair


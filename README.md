# 🧠 From Brain to Base Pair: Transcriptomic Analysis of Language Control

This repository provides all necessary materials for reproducing the transcriptomic analysis presented in our study:  
**"From Brain to Base Pair: A Dual-Site Transcriptomic Analysis of Language Control in Simultaneous Interpreting."**

The project links functional neuroimaging findings with molecular neuroscience by comparing gene expression profiles between the **dorsolateral prefrontal cortex (DLPFC)** and **caudate nucleus**, two regions critical to simultaneous interpreting cognition.

---

## 📂 Repository Structure

| Folder / File Name               | Description                                                                 |
|----------------------------------|-----------------------------------------------------------------------------|
| `ROIs/`                          | MNI coordinates and voxel selection scripts for DLPFC and caudate regions  |
| `Data Sources/`                  | Metadata and links to Allen Human Brain Atlas transcriptome files          |
| `Main packages/`                 | Package list and installation environment (`requirements.txt`)             |
| `Output Highlights/`            | Visualizations and figure outputs (e.g., volcano plots, gene rankings)     |
| `Reproducibility and FAIR Compliance/` | FAIR checklist and code reproducibility documentation             |
| `Citation` / `Contact`           | Citation info and project contact email                                    |

---

## 🚀 How to Reproduce

1. **Clone the repository**  
```bash
git clone https://github.com/vc55chc131/From-Brain-to-Base-Pair.git
cd From-Brain-to-Base-Pair
```

2. **Install required packages (Python 3.10+)**  
```bash
pip install -r requirements.txt
```

3. **Run analysis pipeline**  
Follow the step-by-step instructions in the file:  
```
Reproducibility and FAIR Compliance/analysis_steps.md
```

---

## 📊 Data Sources and Licensing

- Gene expression data retrieved from the **Allen Human Brain Atlas** (AHBA):  
  https://human.brain-map.org  
  - License: Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)

- All scripts and results in this repository are shared under the **MIT License** (see `LICENSE` file).

---

## 📜 How to Cite

> Chang, V. C. Y. (2025). *From Brain to Base Pair: A Dual-Site Transcriptomic Analysis of Language Control in Simultaneous Interpreting*. Manuscript submitted for publication.  
> [GitHub Repository](https://github.com/vc55chc131/From-Brain-to-Base-Pair)

If using any portion of the code or dataset parsing procedures, please cite the GitHub repository directly as:

```bibtex
@misc{chang2025brain2base,
  author       = {Vincent Chieh-Ying Chang},
  title        = {From Brain to Base Pair: Transcriptomic Analysis of Language Control},
  year         = 2025,
  url          = {https://github.com/vc55chc131/From-Brain-to-Base-Pair},
  note         = {Version 1.0. Released under MIT License.}
}
```

---

## 🧬 Author and Contact

Vincent Chieh-Ying Chang  
Department of English, Tamkang University  
📧 159738@o365.tku.edu.tw

---

## 🔍 Related Work in Progress

This project is part of a larger research program integrating **neuroimaging**, **transcriptomics**, and **translation studies**, with upcoming publications including:

- Interpreter Genomics I (current manuscript)
- Interpreter Genomics II (in preparation; DLPFC-focused meta-synthesis of gene variants)

- ### Input Example
The script expects a file named `AHBA_expression_matrix.csv` with the following structure:

| GeneID | 1001 | 1002 | 2001 | 2002 |
|--------|------|------|------|------|
| ABC1   | 1.23 | 2.45 | 0.98 | 1.76 |

Sample IDs prefixed with "1xxx" are DLPFC, "2xxx" are caudate samples.


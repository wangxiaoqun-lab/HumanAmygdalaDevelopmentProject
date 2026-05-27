# Single-Cell and Spatial Transcriptomics Joint Analysis Pipeline

This repository contains R and Python scripts supporting the analysis of human fetal amygdala development using single-cell and spatial transcriptomics data.

---

## 📁 Script Overview

| File                                | Description                                                        |
|-------------------------------------|--------------------------------------------------------------------|
| `01.snRNA_data_process.py`          | Processing snRNA-seq data                                          |
| `02.Stereo-seq_process`             | Processing stereo-seq data                                         |
| `03.Progenitor_process.R`           | Progenitor cell analysis and VGAM fitting                          |
| `04.RNA_velocity.py`                | RNA velocity analysis using scVelo                                 |
| `05.Transfer_label.R`               | Label transfer via CCA-based method                                |
| `06.KNN_analysis.R`                 | Integration of InN and GE data with KNN-based label annotation     |
| `07.Pyscenic_upstream.py`           | Upstream part of SCENIC (building gene regulatory networks)        |
| `08.Pyscenic_downstream.R`          | Downstream SCENIC analysis (regulon scoring and visualization)     |

---

## 📂 Data Source

The datasets used in this project from:

- [GSE261153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261153)
- [HRA006461](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA006461)

---

## 🧾 License

This project is released under the MIT License.

Some scripts are adapted from the following excellent open-source tools:

- [Scanpy](https://github.com/scverse/scanpy)
- [Seurat](https://github.com/satijalab/seurat)
- [SCENIC](https://github.com/aertslab/SCENIC)
- [scVelo](https://github.com/theislab/scvelo)
- [Monocle3](https://github.com/cole-trapnell-lab/monocle3)

---

## 📬 Contact

If you encounter issues or have suggestions, feel free to submit an issue or reach out to the maintainer.

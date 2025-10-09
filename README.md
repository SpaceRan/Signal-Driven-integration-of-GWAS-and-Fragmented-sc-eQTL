# Signal-Driven Integration of GWAS and Fragmented sc-eQTL Data

Conventional colocalization methods (e.g., COLOC) assume complete SNP coverage, consistent linkage disequilibrium (LD) structure, and homogeneous variant representation across datasets, which is rarely satisfied in cross-study integration of single-cell eQTL (sc-eQTL) data. To address this, we launched a signal-driven integrative variant prioritization framework designed to operate under the realistic constraints of fragmented single-cell eQTL data.

## 📢 Current Status

- ✅ All core Python scripts have been uploaded and are fully functional as standalone modules.
- ⚠️ The end-to-end workflow is not yet integrated: while each component works independently and the core algorithm is implemented, the input/output interfaces between modules are still being refined.
- ⚠️ Example datasets and toy examples for quick testing are not yet included, but supplementary data is avaliable in zenodo
- ✅ Per-script documentation and usage tutorials will be added shortly.

We are actively preparing standardized test cases and example configurations to facilitate reuse and validation.

## 🛠️ Requirements

See [`requirements.txt`](requirements.txt) for the full list of dependencies.

## 📄 License
MIT © [SpaceRan] 2025

---

- 📬 **Contact**: Feel free to open an issue or reach out via email: yulangxuan@sjtu.edu.cn
- 📚 Related work：(https://submit.medrxiv.org/submission/pdf?msid=MEDRXIV/2025/337327)
- 📚 full supplementary materials：(10.5281/zenodo.17302764)

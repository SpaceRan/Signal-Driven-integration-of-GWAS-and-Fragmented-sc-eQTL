# Signal-Driven Integration of GWAS and Fragmented sc-eQTL Data

Conventional colocalization methods (e.g., COLOC) assume complete SNP coverage across datasets — an assumption rarely met in practice, especially when integrating genome-wide association studies (GWAS) with single-cell eQTL (sc-eQTL) data. Due to technical limitations and study-specific designs, sc-eQTL datasets are often fragmented across cell types, cohorts, and genotyping platforms, leading to substantial missingness in variant overlap.

To address this challenge, we developed a **signal-driven integrative framework** for variant prioritization that does not rely on full SNP concordance. Instead, our method leverages shared regulatory signals and probabilistic imputation strategies to bridge gaps caused by incomplete variant coverage. The approach enables robust prioritization of causal variants and target genes by aggregating sparse sc-eQTL evidence across cell states and datasets.

## 📢 Current Status

- ✅ All core Python scripts are uploaded and fully functional.
- ✅ The overall workflow logic is clearly structured and reproducible.
- ⚠️ Example datasets and toy runs are **not yet included**.
- ✅ Per-script documentation and usage tutorials will be added shortly.

We are actively preparing standardized test cases and example configurations to facilitate reuse and validation.

## 🛠️ Requirements

See [`requirements.txt`](requirements.txt) for the full list of dependencies.

## 📄 License
MIT © [SpaceRan] 2025

---

📬 **Contact**: Feel free to open an issue or reach out via email: yulangxuan@sjtu.edu.cn
📚 Related work: [bioRxiv DOI]

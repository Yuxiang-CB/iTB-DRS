# iTB-DRS: A multi-stage feature selection framework to identify *Mycobacterium tuberculosis* drug-resistance using single-nucleotide polymorphisms

This repository provides the deployment script (`predict.R`) for **iTB-DRS**. It is designed to predict resistance to 13 anti-tuberculosis drugs directly from a single VCF file (or a batch list of VCF files) using pre-trained machine learning models. 

This code is publicly available to support the reproducibility of our research published in *Microbiology Spectrum*.

---

## 🚀 Quick Start

You can use this tool with a single command. Make sure you have the pre-trained models ready.

```bash
# Predict using a single VCF file
Rscript scripts/predict.R data/raw/Example/ERR2184230.vcf.gz

# Predict using a TXT file containing multiple VCF paths
Rscript scripts/predict.R sample_list.txt
```

---

## 📋 Predictable Drugs

The framework predicts resistance (`R`) or susceptibility (`S`) for the following **13 drugs**:

| Abbreviation | Drug Name |
| :--- | :--- |
| **EMB** | ethambutol |
| **INH** | isoniazid |
| **PZA** | pyrazinamide |
| **RIF** | rifampicin |
| **SM** | streptomycin |
| **AMK** | amikacin |
| **CM** | capreomycin |
| **ETO** | ethionamide |
| **KM** | kanamycin |
| **LFX** | levofloxacin |
| **MFX** | moxifloxacin |
| **OFX** | ofloxacin |
| **PAS** | p-aminosalicylic acid |

---

## 🛠️ Requirements & Installation

You need **R** installed. Install the required R packages before running the script:

```R
install.packages(c("randomForest", "vcfR", "ps"))
```
*Note: `randomForest` and `vcfR` are mandatory. `ps` is optional but recommended for memory usage monitoring.*

### Pre-trained Models
The script requires 13 pre-trained model files (`EMB.rds`, `INH.rds`, etc.). By default, the script looks for a `models` folder in the project root:

```text
├── scripts/
│   └── predict.R
├── models/
│   ├── EMB.rds
│   ├── INH.rds
│   └── ... (13 models in total)
```

---

## 📥 Input Data Requirements

To ensure accurate predictions, your VCF files must follow these rules:
1. **Single-sample**: The VCF must contain only one sample (multi-sample VCFs will trigger an error).
2. **Variant-only**: Contains only variant records (e.g., SNPs/Indels).
3. **Single-contig**: The current model uses `rs<POS>` formatting and cannot distinguish between different contigs. 
4. **Consistent Reference Genome**: Coordinates must match the training data reference genome.
5. **Encoding Rule**: Positions with a `FILTER == "PASS"` record are encoded as `1`; absent positions or failing filters are encoded as `0`. (The specific `REF`/`ALT` or `GT` values are not used, exactly as in the training phase).

### Batch Mode (TXT input)
If you provide a `.txt` file, it will be treated as a batch list. Each line should contain one absolute or relative path to a VCF file. Blank lines and lines starting with `#` are ignored.

---

## ⚙️ Usage and Options

```bash
Rscript scripts/predict.R <sample.vcf[.gz] | sample.txt> [options]
```

| Option | Description | Default |
| :--- | :--- | :--- |
| `--input <file>` | Input VCF or TXT file (can also be passed as a positional argument). | **Required** |
| `--vcf <file>` | Alias for `--input`. | |
| `--models <dir>` | Directory containing the 13 `.rds` model files. | `<project_dir>/models` |
| `--outdir <dir>` | Directory to save the prediction CSV. | `./predictions` |
| `--output <file>` | Specify the exact path/name for the output CSV. | `pred_<name>.csv` |
| `--stable <1-5>` | Selects the feature stability model to use. | `3` |
| `--threshold <0-1>` | Threshold for `prob_R`. If `prob_R >= threshold`, it is predicted as `R`; otherwise `S`. | `0.5` |
| `-h, --help` | Show help message. | |

---

## 📤 Output

### Single VCF Output
The output is a CSV file with one row per drug, containing:
- `drug`: Drug abbreviation.
- `prediction`: `R` (Resistant) or `S` (Susceptible).
- `prob_R`: Probability of resistance.
- `prob_S`: Probability of susceptibility.
- `prediction_time_sec`: Time taken to read the VCF and predict all 13 drugs.

### Batch TXT Output
The output is a CSV file with one row per sample. Each sample will have columns for the predicted class and probabilities for all 13 drugs, plus the total `prediction_time_sec` for that sample.

---

## 📈 Performance
- The 13 models are loaded **once** into memory. The loading time is reported separately and is excluded from the per-sample prediction time.
- Memory usage (RSS) and total CPU time are printed at the end of the execution.

---

## 📖 Citation

If you use this tool or the iTB-DRS framework in your research, please cite our paper:

**APA Format:**
> Chen, Y. X., Andrianov, A. M., Gabrielian, A., Rosenthal, A., & Tuzikov, A. V. (Year). iTB-DRS: a multi-stage feature selection framework to identify Mycobacterium tuberculosis drug-resistance using single-nucleotide polymorphisms. *Microbiology Spectrum*, [Volume]([Issue]), [Pages]. https://doi.org/[Insert DOI]

**BibTeX Format:**
```bibtex
@article{chen2023itbdrs,
  title={iTB-DRS: a multi-stage feature selection framework to identify Mycobacterium tuberculosis drug-resistance using single-nucleotide polymorphisms},
  author={Chen, Yu-Xiang and Andrianov, Alexander M and Gabrielian, Andrei and Rosenthal, Alex and Tuzikov, Alexander V},
  journal={Microbiology Spectrum},
  volume={[Volume]},
  number={[Issue]},
  pages={[Pages]},
  year={[Year]},
  publisher={American Society for Microbiology},
  doi={[Insert DOI]}
}
```

---

## 📄 License

This project is licensed under the MIT License. Please see the LICENSE file for details.

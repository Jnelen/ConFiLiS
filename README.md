# ConFiLiS: Consensus Fingerprints for Ligand-based Screening

ConFiLiS provides an efficient and user-friendly approach to consensus ligand-based screening, enabling the comparison of query compounds (SMILES) against fingerprint libraries. It supports parallel processing, Z-score normalization, and the generation of multiple molecular fingerprints to enhance screening accuracy.

## Table of Contents
1. [Features](#features)
2. [Installation & Setup](#installation--setup)
3. [Examples: Getting Started](#examples-getting-started)
    - [Step 0: Obtain a Sample Compound Library](#step-0-obtain-a-sample-compound-library)
    - [Step 1: Generate Fingerprint Libraries](#step-1-generate-fingerprint-libraries)
    - [Step 2: Perform Ligand-Based Screening](#step-2-perform-ligand-based-screening)
    - [Step 3: Visualize Top Results (optional)](#step-3-visualize-top-result-optional)
4. [Script Usage & Options](#script-usage--options)

## Features
- **Generate compound libraries** from CSV files (containing SMILES) or SDF files.
- **Compute molecular fingerprints** for query compounds.
- **Supports multiple fingerprint types** (Avalon, ECFP6, MACCS, PubChem, and more).  
  See [Supported Fingerprints](#supported-fingerprints-detailed-list) for a full list.
- **Batch processing** of multiple libraries.
- **Parallel execution** for increased performance (`--cores <num_workers>`).
- **Outputs similarity scores** with Z-score normalization.
- **Easily visualize key differences** in top hits using [Fingerprint Similarity Maps](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43).


## Installation & Setup
### Clone the repo
Clone the repository and navigate to it:
```bash
git clone https://github.com/Jnelen/ConFiLiS.git
```
```bash
cd ConFiLiS
```
### Install the dependencies using Conda
The required dependencies for ConFiLiS can be installed using a Conda environment.
Run the following command to create the environment:
```bash
conda env create -f environment.yml
```

and activate it using:

```bash
conda activate confilis
```

### Singularity Image (For HPC Systems)
For use on HPC clusters, a Singularity image containing the dependencies is available. You can either download the pre-built image:
```bash
wget --no-check-certificate -r "https://drive.usercontent.google.com/download?id=1W4ggOOjeeC_kMy3gvJVbD7BEPn8x101N&confirm=t" -O singularity/ConFiLiS.sif  
```
Or build the image manually using the provided `.def` file included in the `singularity` directory:

```bash
singularity build singularity/ConFiLiS.sif singularity/ConFiLiS.def 
```

Whenever you run a Python script, use the Singularity container instead of calling `python` directly:
```bash
singularity exec singularity/ConFiLiS.sif python ...
```

## Examples: Getting Started

This section provides a step-by-step tutorial to demonstrate how to:
1. **Generate fingerprint libraries** from a file containing SMILES and compound IDs.
2. **Perform ligand-based screening** using the generated libraries.
3. **Visualize the top hit** using a [Fingerprint Similarity Map](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43).

### Step 0: Obtain a Sample Compound Library

The `convert_library.py` script supports two main formats:
- **SDF files**  
- **Text-based files** (CSV, TSV, ...) containing **both SMILES and compound IDs**  

For text-based files, the script can:
- **Automatically detect the delimiter** (`,`, `;`, tabs, etc.).
- **Adjust for different column orders** (you can specify the SMILES and ID column positions).  

More details are provided in the [Script Usage & Options](#script-usage--options) section.

#### **Download the SCUBIDOO_S Library**
For this example, we will use a **small open-source compound library**, [SCUBIDOO](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00203), as our sample dataset.

Download it using:
```bash
wget -O Scubidoo_S.ism "https://scubidoo.pharmazie.uni-marburg.de/downloads/14_SCU_sample_S.ism"
```
### Step 1: Generate Fingerprint Libraries

Convert the **sample compound library** into a fingerprint database using the `convert_library.py` script. 
You should specify the fingerprint types you want to generate using the `--fingerprint` option. 
The following example will create libraries for five commonly used fingerprints:
```python
python convert_library.py Scubidoo_S.ism libraries/Scubidoo_S --fingerprints avalon circular pubchem rdk-maccs rdkit
```

**Note**: a full list of supported fingerprints is listed in here: [Supported Fingerprints](#supported-fingerprints-detailed-list)

These calculations typically complete within a few minutes and include:
- Reading SMILES data from the input file
- Generating the requested fingerprints
- Saving the results as compressed `.pkl.gz` files in `libraries/Scubidoo_S/`

### Step 2: Perform Ligand-Based Screening

Now, let's **screen a query molecule** against the generated fingerprint library.

For this example, we will **reproduce an experiment** from the original [SCUBIDOO](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00203) paper.  
In the paper, the authors **screened the compound [DB08235](https://go.drugbank.com/drugs/DB08235) against their library.**  

We can replicate this experiment with the following command:
```python
python screen.py -i "CC1=C(CCNC(=O)C2=CC=CS2)C2=CC=CC=C2N1" -n DB08235 -l libraries/Scubidoo_S -f avalon circular pubchem --cores 3 --shortened_rows 0  
```

This will:
- Compare **DB08235** against `libraries/Scubidoo_S/`
- Compute **similarity scores** using **Avalon, circular, and PubChem**
- Run **in parallel** using **3 CPU cores**
- **Skip generating the shortened report** (by setting `--shortened_rows 0`).

#### **Check the Results**
If everything went well, this script should run very quickly (just a few seconds), and you should find the **screening results** in the `outputs/` directory: `outputs/ConFiLiS_DB08235_Scubidoo_S.csv`

Inside this file, the **top hit should be compound `10143065`**, which **matches the results** reported in the [SCUBIDOO](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00203) paper.

### Step 3: Visualize Top Result (optional)
This repo also contains a script to visualise the similarity according to various fingerprints between a query compoud and a library compound.
As an example, we will apply it on result from the last step. You can run the calculation by providing the query smiles, compound hit and the library this hit is from like this:
```python
 python fp_similarity_map.py --query "CC1=C(CCNC(=O)C2=CC=CS2)C2=CC=CC=C2N1" --id 10143065 --library libraries/Scubidoo_S/
```
The output image should be saved in `output/SimilarityMap/`, and look something like this:
![Morgan_1_10143065_Cc1_nH_c2ccccc2c1CCNC__O_c1cccs1](https://github.com/user-attachments/assets/f0f6b038-44af-49ed-8413-d7d633442339)

#### **Understanding the Similarity Map**
The visualization is based on **[Fingerprint Similarity Maps](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43)**, which highlight molecular similarities and differences:  
- **Green outlines** indicate regions of high similarity between the two compounds.  
- **Red outlines** highlight structural differences.  
- **Line density** reflects the magnitude of similarity or difference.  

There are multiple similarity map types available. More details can be found in the [Script Usage & Options](#script-usage--options) section.

## Script Usage & Options

This section provides a detailed overview of the main scripts, their capabilities, and available settings.

### 1. `convert_library.py` - Create Fingerprint Libraries

This script processes **Text or SDF** files and generates fingerprint libraries, storing them as **compressed pickle (`.pkl.gz`) files** for fast lookup during screening.

#### **Arguments**
| Argument | Description |
|----------|-------------|
| `input_file` | Input file (SDF or CSV, TSV, ...) containing compounds and corresponding molecule IDs |
| `output_directory` | Directory where fingerprint files will be saved |
| `-f, --fingerprints, --fps` | List of fingerprints to compute (Required). See the [Supported Fingerprints](#supported-fingerprints-detailed-list) table below. |
| `--delimiter` | (Optional) Delimiter for CSV/TSV (auto-detected if not specified) |
| `--smiles_index` | Column index for SMILES strings (default: `0`) |
| `--name_index` | Column index for compound names (default: inferred from SMILES) |
| `--header` | If set, skips the first row of CSV/TSV files |
| `--num_workers, --cores, -c` | Number of parallel processes (default: `1`) |
| `--silent` | If set, disables progress bars |

### 2. `screen.py` - Ligand-Based Screening

This script **compares a query compound (SMILES format)** against **one or more fingerprint libraries**, computing similarity scores based on Euclidean distance and Z-score normalization.

#### **Arguments**
| Argument | Description |
|----------|-------------|
| `-i, --input, --smiles` | Query compound in SMILES format |
| `-n, --molname, --name` | (Optional) Name of the compound (used in output filenames) |
| `-l, --libraries, --libs` | One or more directories containing fingerprint libraries (`*.pkl.gz` files) |
| `-f, --fingerprints, --fps` | List of fingerprint types (Required). See the [Supported Fingerprints](#supported-fingerprints-detailed-lists) table below. |
| `-o, --output, --output_path` | Output directory or file path prefix (defaults to `"output"`) |
| `-c, --cores, --num_workers` | Number of parallel processes (default: `1`) |
| `--shortened_rows` | Number of rows in the shortened report (default: `1000`, set to `0` to disable). |

### 3. `fp_similarity_map.py` - Fingerprint Similarity Visualization

This script generates **[similarity maps](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43)** to visualize atomic contributions to molecular similarity. It compares a **query molecule (SMILES format)** with a **reference compound** from a molecular library and highlights atoms contributing most to fingerprint-based similarity.

#### **Arguments**
| Argument | Description |
|----------|-------------|
| `--query` | SMILES string of the query molecule (should be between quotes) |
| `--id` | ID of the reference compound from the molecular library |
| `--library` | Path to the fingerprint library |
| `--output` | Output directory for similarity maps (default: `output/SimilarityMap`) |
| `--types` | Fingerprint types for similarity mapping (default: `Morgan_1`). <br> Available options: `Morgan_1`, `Morgan_2`, `Morgan_3`, `ap`, `tt`. Use `'all'` for all. |

---

### Supported Fingerprints (Detailed List)

ConFiLiS supports a variety of fingerprint types for molecular similarity calculations. You can use the identifiers listed here with the `--fingerprints` option in scripts (`convert_library.py`, `screen.py`):

| Fingerprint Name      | Identifier (for `--fingerprints`) |
|-----------------------|--------------------------------|
| Avalon               | `avalon`                          |
| Atom-Pair            | `atom-pair`                       |
| CDK-Substructure     | `cdk-substructure`                |
| Circular             | `circular`                        |
| ECFP4                | `ecfp4`                           |
| ECFP6                | `ecfp6`                           |
| Klekota-Roth         | `klekota-roth`                    |
| Mol2Vec              | `mol2vec`                         |
| PubChem              | `pubchem`                         |
| rdk-MACCS            | `rdk-maccs`                       |
| RDKit                | `rdkit`                           |
| Shortest Path        | `shortest-path`                   |
| Topological-Torsion  | `topological-torsion`             |

**Note:** Fingerprints are computed using [PyFingerprint](https://github.com/hcji/PyFingerprint). Non-default fingerprint options may yield varying results.


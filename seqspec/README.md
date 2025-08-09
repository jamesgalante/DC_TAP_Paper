# seqspec templates

This folder contains Jinja templates for generating seqspec YAML files for DC-TAP-Seq experiments. Templates are modality-specific (e.g., `crispr`, `rna`).

## Environment Setup

It is recommended to use a Python virtual environment.  
To create and activate an environment and install the required libraries:

```bash
python3 -m venv venv
source venv/bin/activate
pip install click jinja2
```

## How to Use

To generate YAML files from a TSV of sample metadata, use the `fill_seqspec_template.py` script from the repository root:

```bash
# To generate YAML files for CRISPR modality, run:
python fill_seqspec_template.py --tsv example/example_crispr.tsv --modality crispr --output-dir output/
# To generate YAML files for RNA modality, run:
python fill_seqspec_template.py --tsv example/example_crispr.tsv --modality rna --output-dir output/
```

- `--tsv`: Path to your input TSV file (see below for structure).
- `--modality`: Choose `crispr` or `rna` to select the template.
- `--output-dir`: Directory for output YAML files.

## TSV Structure

The TSV file should have a header row with columns matching the template variables. For the `crispr` modality, required columns include:

| Column Name                      | Description                                                                 | Example|
|---------------------------------- |-----------------------------------------------------------------------------|-------|
| assay_id                         | Unique identifier for the assay/sample                                      | DC-Tap-Seq |
| library_structure                | Library structure description                                               | https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html |
| library_protocol                 | Library preparation protocol                                                | Custom |
| library_kit                      | Library kit used                                                            | Illumina Truseq Dual Index |
| sequencing_protocol              | Sequencing protocol used                                                    | Illumina NovaSeq 6000 (EFO:0008637) |
| sequencing_kit                   | Sequencing kit used                                                         | Illumina NovaSeq 6000 S4 Reagent Kit v1.5 |
| read1_id                         | Identifier for read 1 file                                                  | IGVFFI4996PJUA.fastq.gz |
| read1_name                       | Name for read 1                                                             | IGVFFI4996PJUA |
| read1_size                       | Size (bytes) of read 1 file                                                 | 123456789 |
| read1_url                        | URL to download read 1 file                                                 | https://api.data.igvf.org/sequence-files/IGVFFI4996PJUA/@@download/IGVFFI4996PJUA.fastq.gz |
| read1_urltype                    | URL type for read 1 (e.g., http, ftp)                                       | https |
| read1_md5                        | MD5 checksum for read 1 file                                                | d41d8cd98f00b204e9800998ecf8427e |
| read2_id                         | Identifier for read 2 file                                                  | IGVFFI4996PJUB.fastq.gz |
| read2_name                       | Name for read 2                                                             | IGVFFI4996PJUB |
| read2_size                       | Size (bytes) of read 2 file                                                 | 123456789 |
| read2_url                        | URL to download read 2 file                                                 | https://api.data.igvf.org/sequence-files/IGVFFI4996PJUB/@@download/IGVFFI4996PJUB.fastq.gz |
| read2_urltype                    | URL type for read 2 (e.g., http, ftp)                                       | https |
| read2_md5                        | MD5 checksum for read 2 file                                                | d41d8cd98f00b204e9800998ecf8427e |
| cell_barcode_file_id      | File ID for cell barcode onlist                                             | IGVFFI4996PJUC.txt |
| cell_barcode_file_name    | Filename for cell barcode onlist                                            | IGVFFI4996PJUC |
| cell_barcode_file_type    | File type for cell barcode onlist (e.g., txt)                               | txt |
| cell_barcode_file_size    | Size (bytes) of cell barcode onlist file                                    | 123456789 |
| cell_barcode_file_url     | URL to download cell barcode onlist file                                    | https://api.data.igvf.org/sequence-files/IGVFFI4996PJUC/@@download/IGVFFI4996PJUC.txt |
| cell_barcode_file_urltype | URL type for cell barcode onlist file (e.g., http, ftp)                     | https |
| cell_barcode_file_md5     | MD5 checksum for cell barcode onlist file                                   | d41d8cd98f00b204e9800998ecf8427e |
| crispr_guides_file_id            | File ID for guides onlist                                                   | IGVFFI4996PJUD.txt |
| crispr_guides_file_name          | Filename for guides onlist                                                  | IGVFFI4996PJUD |
| crispr_guides_file_type          | File type for guides onlist (e.g., txt)                                     | txt |
| crispr_guides_file_size          | Size (bytes) of guides onlist file                                          | 123456789 |
| crispr_guides_file_url           | URL to download guides onlist file                                          | https://api.data.igvf.org/sequence-files/IGVFFI4996PJUD/@@download/IGVFFI4996PJUD.txt |
| crispr_guides_file_urltype       | URL type for guides onlist file (e.g., http, ftp)                           | https |
| crispr_guides_file_md5           | MD5 checksum for guides onlist file                                         | d41d8cd98f00b204e9800998ecf8427e |

**Note:**  
- For the `rna` modality, you can leave the crispr_guide_* columns empty.
- See example/example_crispr.tsv` for a sample TSV.

---

# Evaluation of a silicosis diagostic prediction rule

This repository contains the code and documentation for the evaluation of a diagnostic prediction rule for silicosis, which was developed in the Netherlands [see 2007 research paper](https://oem.bmj.com/lookup/doi/10.1136/oem.2006.027904) to rule out pneumoconiosis and identify workers at high risk for further diagnostic workup. [Lexces](https://www.lexces.nl/) aims to prevent new cases of silicosis in Dutch workers. For this, [a Health Surveillance Program (HSP) for respiratory occupational diseases is being developed](https://www.lexces.nl/en/node/52). The diagnostic prediction rule could be incorporated as part of the HSP. However, the diagnostic rule was developed in the past with chest x-rays (CXR) as the reference standard for the diagnosis of silicosis. Recently, there have been concerns of [suboptimal diagnostic performance of CXR for the diagnosis of silicosis](https://onlinelibrary.wiley.com/doi/10.1111/resp.14755). Thus, the current work aims to assess the potential impact of misclassification error on the diagnostic rule's performance, and to scope alternative diagnostic and prediction models for silicosis/pneumoconiosis. 

## How to use 

This repository is a work in progress and is currently not fully enabled for re-use/collaboration. Further details of how to use it will be provided in the future as the project evolves.

## Project Structure

The project structure distinguishes three kinds of folders:
- read-only (RO): not edited by either code or researcher
- human-writeable (HW): edited by the researcher only.
- project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.

```         
.
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
├── lexces-silicosis-predict.Rproj
├── data                  <- All project data files
│   ├── processed         <- The final, canonical data sets for modeling. (PG)
│   ├── raw               <- The original, immutable data. (RO)
│   └── temp              <- Intermediate data that has been transformed. (PG)
├── docs                  <- Documentation for users (HW)
│   ├── manuscript        <- Manuscript source, docx. (HW)
│   ├── presentations     <- Presentations, pptx, pdf. (HW)
│   └── reports           <- Project reports, pdf. (HW)
├── results
│   ├── output_figures    <- Figures for the manuscript or reports (PG)
│   └── output_tables     <- Output tables for the manuscript (PG)
└── R                     <- Source code for this project (HW)
    ├── scripts           <- Scripts sourced in main R markdown documents (PG)
    └── sessions          <- Text files with information of R sessions (PG)

```

## License

This project is licensed under the terms of the [MIT License](/LICENSE).

This project structure template repository is adapted from the [Good Enough Project](https://github.com/bvreede/good-enough-project) Cookiecutter template by Barbara Vreede (2019).

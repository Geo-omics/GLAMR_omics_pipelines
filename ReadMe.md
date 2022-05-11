# GLAMR: The Great Lakes Atlas of Multi-omics Research
This repository contains the pipleines that power the GLAMR database.

GLAMR is designed to be a centralized resource housing Great Lakes omics datasets analyzed with standardized pipelines and integrated with environmental data. 

## Folder structure
```
GLAMR
├── code
├── config
│   ├── conda_yaml
│   └── slurm_profiles
├── data
│   ├── environment
│   ├── omics
│   │   ├── amplicons
│   │   ├── genomes
│   │   ├── metabolomes
│   │   ├── metagenomes
│   │   └── metatranscriptomes
│   ├── reference
│   │   ├── amplicons
│   │   ├── genomes
│   │   ├── MAGs
│   │   └── UMRAD
│   └── sample_metadata
└── ReadMe.md
```

## Snakemake
Snakemake is used extensively to manage the workflows used in this database. 

Add:
- Instructions for adding samples
- Instructions for running individual pipelines and all pipelines together
- Image summarizing the rule-graph

## Incorporating code from outside repositories

### Examples:

1. Adding an external repository
    ```
    git subtree add --prefix code/Strain-Level_Metagenome_Analysis https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis.git master --squash
    ```

2. Updating an external repository
    ```
    git subtree pull --prefix code/Strain-Level_Metagenome_Analysis https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis.git master --squash
    ```

3. List subtrees
    ```
    git subtree list --resolve
    ```
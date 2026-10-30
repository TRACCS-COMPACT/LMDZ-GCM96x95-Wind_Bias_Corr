# LMDZ-GCM96x95-Wind_Bias_Corr


```DOI TO CREATE```

## Context and Motivation

Purpose of the experiment is to compute the winds bias correction on a realistic atmospheric general circulation configuration. Internal and external computed fluxes are written in an output file with the LMDZ output system (XIOS).

#### Variations
- **K26** : Correction computed with [dummy_ref](https://gitlab.in2p3.fr/ipsl/espri/espri-ia/projects/lmdz_wind_bc) ML model 

## Requirements

- Experiment manager : [IGCM Tools suite](https://ipsl.pages.in2p3.fr/icmc/plateforme/documentation/igcmg_doc) for LMDZOR_v6 config

### Compilation

- LMDZ6 trunk [5891](http://svn.lmd.jussieu.fr/LMDZ/LMDZ6/trunk/) patched with [LMDZ Patch](https://github.com/TRACCS-COMPACT/Patches-LMDZ) and local `CONFIG/my_src` sources.

- ORCHIDEE trunk [7850](https://forge.ipsl.fr/orchidee/browser/branches/ORCHIDEE_3_CMIP6/ORCHIDEE?rev=7850)

- Code Compilation manager : none, use `compile.sh` script provided in `CONFIG`


### Python

- Eophis version : [v1.1.0](https://github.com/alexis-barge/eophis/tree/v1.1.0)
- **K26** dependencies:
  ```bash
  cd GCM96x95-Wind_Bias_Corr.K26/INFERENCES/
  pip install -r requirements.txt
  ```


### Run

- LMDZ Production Manager : [IGCM Tools suite](https://ipsl.pages.in2p3.fr/icmc/plateforme/documentation/igcmg_doc)


### Post-Process

- No Post-Process libraries

- Plotting : Python script `plots_res.py` in `POSTPROCESS`



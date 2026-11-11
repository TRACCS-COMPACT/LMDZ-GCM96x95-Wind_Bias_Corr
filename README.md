# LMDZ-GCM96x95-Wind_Bias_Corr


```DOI TO CREATE```

## Context and Motivation

Purpose of the experiment is to compute the winds bias correction on a realistic atmospheric general circulation configuration. Internal and external computed fluxes are written in an output file with the LMDZ output system (XIOS).

#### Variations
- **K26** : Correction computed with [lmdz_wind_bc](https://gitlab.in2p3.fr/ipsl/espri/espri-ia/projects/lmdz_wind_bc) ML model, [commit dccf81fd](https://gitlab.in2p3.fr/ipsl/espri/espri-ia/projects/lmdz_wind_bc/-/tree/dccf81fdb0a210f93dfa58ea880ebf464352d1ed)

## Requirements

- Experiment manager:
  
   - [IGCM Tools suite](https://ipsl.pages.in2p3.fr/icmc/plateforme/documentation/igcmg_doc) for LMDZOR_v6 config
   - libIGCM trunk [1582](https://forge.ipsl.jussieu.fr/libigcm/browser/trunk/libIGCM?rev=1582)


### Compilation

- Modeles:
 
    - LMDZ6 trunk [5891](http://svn.lmd.jussieu.fr/LMDZ/LMDZ6/trunk/) patched with [LMDZ Patch](https://github.com/TRACCS-COMPACT/Patches-LMDZ) and local `CONFIG/my_src` sources.

    - ORCHIDEE trunk [7850](https://forge.ipsl.fr/orchidee/browser/branches/ORCHIDEE_3_CMIP6/ORCHIDEE?rev=7850)

    - IOIPSL trunk [6339](https://forge.ipsl.jussieu.fr/igcmg/browser/IOIPSL/trunk?rev=6339)

- Code Compilation manager : none, use `compile_MOD.sh` scripts provided in `CONFIG`. Adapt paths in arch files, if necessary.


### Python

- Eophis version : [v1.1.0](https://github.com/alexis-barge/eophis/tree/v1.1.0)
- **K26** dependencies:
  ```bash
  git clone https://gitlab.in2p3.fr/ipsl/espri/espri-ia/projects/lmdz_wind_bc.git GCM96x95-Wind_Bias_Corr.K26/INFERENCES/lmdz_wind_bc
  cd GCM96x95-Wind_Bias_Corr.K26/INFERENCES/lmdz_wind_bc
  git checkout dccf81fdb0a210f93dfa58ea880ebf464352d1ed
  pip install -e .
  ```
- Adapt line 8 in ``DRIVER/eophis.card`` and 7-9 in ``COMP/eophis.card``. Copy ``eophis_wind_bc.py`` in ``modipsl/bin``

### Run

- LMDZ Production Manager : [IGCM Tools suite](https://ipsl.pages.in2p3.fr/icmc/plateforme/documentation/igcmg_doc)
  

### Post-Process

- No Post-Process libraries

- Plotting : Python script `plots_res.py` in `POSTPROCESS`



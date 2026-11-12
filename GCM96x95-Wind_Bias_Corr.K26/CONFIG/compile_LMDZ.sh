#!/bin/bash

./makelmdz_fcm -d 96x95x79 -p lmd -rad rrtm -prod -parallel mpi_omp -cpp key_eophis -c true -io xios -v orchidee2.1 -arch gcc_irene_cpl -j 8 gcm
cd - ; cd modipsl/modeles/LMDZ/bin ; cp gcm_96x95x79_phylmd_rrtm_para_mem_orch_couple.e ../../../bin

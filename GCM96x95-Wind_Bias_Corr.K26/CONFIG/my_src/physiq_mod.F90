
! $Id: physiq_mod.F90 5842 2025-09-25 15:09:22Z rkazeroni $
!
!#define IO_DEBUG
MODULE physiq_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE physiq (nlon,nlev, &
       debut,lafin,pdtphys_, &
       paprs,pplay,pphi,pphis,presnivs, &
       u,v,rot,t,qx, &
       flxmass_w, &
       d_u, d_v, d_t, d_qx, d_ps)

! For clarity, the "USE" section is now arranged in alphabetical order,
! with a separate section for CPP keys
! PLEASE try to follow this rule

    USE ACAMA_GWD_rando_m, only: ACAMA_GWD_rando, ACAMA_GWD_rando_first
    USE add_wake_tend_mod, ONLY: add_wake_tend
    USE aero_mod
    USE add_phys_tend_mod, only : add_pbl_tend, add_phys_tend, diag_phys_tend, prt_enerbil, &
  &      fl_ebil, fl_cor_ebil
    USE ajsec_mod, ONLY: ajsec, ajsec_convv2
    USE alpale_th_mod, ONLY: alpale_th, alpale_th_first
    USE alpale_wk_mod, ONLY: alpale_wk, alpale_wk_first
    USE assert_m, only: assert
    USE change_srf_frac_mod
    USE conf_phys_m, only: conf_phys
    USE carbon_cycle_mod, ONLY : infocfields_init, RCO2_glo, carbon_cycle_rad
    USE CFMIP_point_locations   ! IM stations CFMIP
    USE clouds_bigauss_mod, ONLY: clouds_bigauss
    USE clouds_gno_mod, ONLY: clouds_gno
    USE cmp_seri_mod
    USE diag_slp_mod, ONLY: diag_slp
    USE dimphy
    USE etat0_limit_unstruct_mod
    USE ener_conserv_mod, ONLY: ener_conserv
    USE evappot_mod, ONLY: evappot
    USE FLOTT_GWD_rando_m, only: FLOTT_GWD_rando, FLOTT_GWD_rando_first
    USE fonte_neige_mod, ONLY  : fonte_neige_get_vars
    USE geometry_mod, ONLY: cell_area, latitude_deg, longitude_deg
    USE ioipsl, only: histbeg, histvert, histdef, histend, histsync, &
         histwrite, ju2ymds, ymds2ju, getin
    USE ioipsl_getin_p_mod, ONLY : getin_p
    USE indice_sol_mod
    USE infotrac_phy, ONLY: nqtot, nbtr, nqo, nqtke, tracers, type_trac, addPhase, ivap, iliq, isol, ibs, icf, irvc, itke
    USE strings_mod,  ONLY: strIdx
    USE iophy
    USE limit_read_mod, ONLY : init_limit_read
    USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev, klon_glo, grid1dTo2d_glo, grid_type, unstructured
    USE mod_phys_lmdz_mpi_data, only: is_mpi_root
    USE mod_phys_lmdz_para
    USE netcdf95, only: nf95_close
    USE netcdf, only: nf90_fill_real     ! IM for NMC files
    USE nuage_mod, ONLY: nuage, diagcld1, diagcld2
    USE open_climoz_m, only: open_climoz ! ozone climatology from a file
    USE orbite_mod, ONLY: angle, orbite, zenang
    USE orografi_mod, ONLY: drag_noro, lift_noro, sugwd
    USE orografi_strato_mod, ONLY: drag_noro_strato, lift_noro_strato, sugwd_strato
    USE ozonecm_m, only: ozonecm ! ozone of J.-F. Royer
    USE pbl_surface_mod, ONLY : pbl_surface
    USE phyaqua_mod, only: zenang_an
    USE phyetat0_mod, only: phyetat0
    USE phystokenc_mod, ONLY: offline, phystokenc
    USE phys_cal_mod, only: year_len, mth_len, days_elapsed, jh_1jan, &
         year_cur, mth_cur,jD_cur, jH_cur, jD_ref, day_cur, hour, calend
!!  USE phys_local_var_mod, ONLY : a long list of variables
!!              ==> see below, after "CPP Keys" section
    USE phys_state_var_mod ! Variables sauvegardees de la physique
    USE phys_output_mod
    USE phys_output_ctrlout_mod
    USE print_control_mod, ONLY: mydebug=>debug , lunout, prt_level, &
         alert_first_call, call_alert, prt_alerte
    USE readaerosol_mod, ONLY : init_aero_fromfile
    USE readaerosolstrato_m, ONLY : init_readaerosolstrato
    USE radlwsw_m, only: radlwsw
    USE reevap_mod, ONLY: reevap
    USE regr_horiz_time_climoz_m, ONLY: regr_horiz_time_climoz
    USE regr_pr_time_av_m, only: regr_pr_time_av
    USE stratocu_if_mod, ONLY: stratocu_if
    USE surface_data,     ONLY : type_ocean, ok_veget
    USE tend_to_tke_mod, ONLY: tend_to_tke
    USE time_phylmdz_mod, only: current_time, itau_phy, pdtphys, raz_date, update_time
    USE tracinca_mod, ONLY: config_inca
    USE transp_mod, ONLY: transp
    USE tropopause_m,     ONLY: dyn_tropopause
    USE vampir
    USE water_int_mod, ONLY: water_int
    USE write_field_phy
    use wxios_mod, ONLY: g_ctx, wxios_set_context
    USE lmdz_lscp_main, ONLY : lscp
    USE lmdz_call_cloud_optics_prop, ONLY : call_cloud_optics_prop, call_cloud_optics_prop_post
    USE lmdz_lscp_old, ONLY : fisrtilp, fisrtilp_first
    USE lmdz_call_blowing_snow, ONLY : call_blowing_snow_sublim_sedim
    USE calwake_mod, ONLY : calwake, calwake_first
    USE lmdz_wake_ini, ONLY : wake_ini
    USE lmdz_surf_wind_ini, ONLY : surf_wind_ini
    USE lmdz_surf_wind, ONLY : surf_wind
    USE yamada_ini_mod, ONLY : yamada_ini
    USE lmdz_atke_turbulence_ini, ONLY : atke_ini
    USE lmdz_thermcell_ini, ONLY : thermcell_ini, iflag_thermals_tenv
    USE calltherm_mod, ONLY : calltherm
    USE lmdz_thermcell_dtke, ONLY : thermcell_dtke
    USE lmdz_blowing_snow_ini, ONLY : blowing_snow_ini , qbst_bs
    USE lmdz_lscp_ini, ONLY : lscp_ini
    USE lmdz_ratqs_main, ONLY : ratqs_main, ratqs_main_first
    USE lmdz_ratqs_ini, ONLY : ratqs_ini
    USE lmdz_cloud_optics_prop_ini, ONLY : cloud_optics_prop_ini
    USE phys_output_var_mod, ONLY :      cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv
    USE phys_output_var_mod, ONLY : cloud_cover_sw, cloud_cover_sw_s2, ZFLUX_DIR, ZFLUX_DIR_CLEAR, &
  &      ZFLUX_DIR_s2, ZFLUX_DIR_CLEAR_s2, ZFLUX_DIR_SUN, ZFLUX_DIR_SUN_s2
    USE phys_output_var_mod, ONLY : type_ecri


    !USE cmp_seri_mod
!    USE add_phys_tend_mod, only : add_pbl_tend, add_phys_tend, diag_phys_tend, prt_enerbil, &
!  &      fl_ebil, fl_cor_ebil


    USE phytracr_spl_mod, ONLY: phytracr_spl, phytracr_spl_out_init
    USE phys_output_write_spl_mod, ONLY: phys_output_write_spl
    USE phytrac_mod, ONLY: phytrac_init, phytrac
    USE phys_output_write_mod, ONLY: phys_output_write


    USE geometry_mod,      ONLY: longitude, latitude, boundslon, boundslat, ind_cell_glo
    USE time_phylmdz_mod,  ONLY: ndays
    USE infotrac_phy,      ONLY: nqCO2
    USE lmdz_reprobus_wrappers, ONLY: init_chem_rep_xjour, d_q_rep, d_ql_rep, d_qi_rep, &
         ptrop, ttrop, ztrop, gravit, itroprep, Z1, Z2, fac, B, &
         chemini_rep, rtime, chemtime_rep, coord_hyb_rep
    USE time_phylmdz_mod,    ONLY: annee_ref, day_ini, day_ref, start_time
    USE vertical_layers_mod, ONLY: aps, bps, ap, bp
    USE lmdz_cppkeys_wrapper, ONLY: CPPKEY_INCA, CPPKEY_REPROBUS, CPPKEY_STRATAER, CPPKEY_DUST, CPPKEY_COSP, CPPKEY_COSP2, CPPKEY_COSPV2

#ifdef CPP_RRTM
    USE YOERAD, ONLY : NRADLP
!    USE YOESW, ONLY : RSUN
#endif

    USE phys_local_var_mod, ONLY: d_q_emiss
    USE strataer_local_var_mod, ONLY: strataer_init,flag_emit,flh2o,ok_qemiss,flag_verbose_strataer, &
         year_emit_vol,mth_emit_vol,day_emit_vol,nErupt,nAerErupt,injdur,m_H2O_emiss_vol_daily,m_H2O_emiss_vol, &
         ponde_lonlat_vol,altemiss_vol,sigma_alt_vol,xlat_min_vol,xlat_max_vol,xlon_min_vol,xlon_max_vol
    USE strataer_nuc_mod, ONLY: strataer_nuc_init
    USE strataer_emiss_mod, ONLY: strataer_emiss_init
    
    USE lmdz_xios, ONLY: xios_update_calendar, xios_context_finalize
    USE lmdz_xios, ONLY: xios_get_field_attr, xios_field_is_active, xios_context
    USE lmdz_xios, ONLY: xios_set_current_context
    USE wxios_mod, ONLY: missing_val, using_xios
    USE lmdz_spla_ini, ONLY : spla_ini
    USE pycpl
    USE pyfld
    USE comconst_mod, ONLY: daysec
#ifndef CPP_XIOS
    USE paramLMDZ_phy_mod
#endif
!
!
!!!!!!!!!!!!!!!!!!  END "USE" for CPP keys !!!!!!!!!!!!!!!!!!!!!!

USE physiqex_mod, ONLY : physiqex
USE phys_local_var_mod, ONLY: phys_local_var_init, phys_local_var_end, &
       ! [Variables internes non sauvegardees de la physique]
       ! Variables locales pour effectuer les appels en serie
       p_tropopause, dyntropo, &
       t_seri,q_seri,ql_seri,qs_seri,qbs_seri, &
       u_seri,v_seri,cf_seri,rvc_seri,tr_seri, &
       rhcl, &
       ! Dynamic tendencies (diagnostics)
       d_t_dyn,d_q_dyn,d_ql_dyn,d_qs_dyn,d_qbs_dyn, &
       d_u_dyn,d_v_dyn,d_cf_dyn,d_rvc_dyn,d_tke_dyn,d_tr_dyn, &
       d_q_dyn2d,d_ql_dyn2d,d_qs_dyn2d,d_qbs_dyn2d, &
       ! Physic tendencies
       d_t_con,d_q_con,d_u_con,d_v_con, &
       d_t_con_zmasse,d_q_con_zmasse,d_u_con_zmasse,d_v_con_zmasse, &        
       d_tr, &                              !! to be removed?? (jyg)
       d_t_wake,d_q_wake, &
       d_t_lwr,d_t_lw0,d_t_swr,d_t_sw0, &
       d_t_ajsb,d_q_ajsb, &
       d_t_ajs,d_q_ajs,d_u_ajs,d_v_ajs, &
!       d_t_ajs_w,d_q_ajs_w, &
!       d_t_ajs_x,d_q_ajs_x, &
       !
       d_t_eva,d_q_eva,d_ql_eva,d_qi_eva, &
       d_t_lsc,d_q_lsc,d_ql_lsc,d_qi_lsc, &
       d_t_lscst,d_q_lscst, &
       d_t_lscth,d_q_lscth, &
       plul_st,plul_th, &
       !
       d_t_vdf,d_q_vdf, d_qbs_vdf, d_u_vdf,d_v_vdf,d_t_diss, &
       d_t_vdf_x, d_t_vdf_w, &
       d_q_vdf_x, d_q_vdf_w, &
       d_ts, &
       !
       d_t_bsss,d_q_bsss,d_qbs_bsss, &
       !
!       d_t_oli,d_u_oli,d_v_oli, &
       d_t_oro,d_u_oro,d_v_oro, &
       d_t_oro_gw,d_u_oro_gw,d_v_oro_gw, &
       d_t_lif,d_u_lif,d_v_lif, &
       d_t_ec, &
       !
       du_gwd_hines,dv_gwd_hines,d_t_hin, &
       dv_gwd_rando,dv_gwd_front, &
       east_gwstress,west_gwstress, &
       d_q_ch4, &
       ! proprecip
       qraindiag, qsnowdiag, &
       dqreva, dqssub, &
       dqrauto,dqrcol,dqrmelt,dqrfreez, &
       dqsauto,dqsagg,dqsrim,dqsmelt,dqsfreez, &
       !  Special RRTM
       ZLWFT0_i,ZSWFT0_i,ZFLDN0,  &
       ZFLUP0,ZFSDN0,ZFSUP0,      &
       !
       topswad_aero,solswad_aero,   &
       topswai_aero,solswai_aero,   &
       topswad0_aero,solswad0_aero, &
       !LW additional
       toplwad_aero,sollwad_aero,   &
       toplwai_aero,sollwai_aero,   &
       toplwad0_aero,sollwad0_aero, &
       !pour Ecrad
       topswad_aero_s2, solswad_aero_s2,   &
       topswai_aero_s2, solswai_aero_s2,   &
       topswad0_aero_s2, solswad0_aero_s2, &
       topsw_aero_s2, topsw0_aero_s2,      &
       solsw_aero_s2, solsw0_aero_s2,      &
       topswcf_aero_s2, solswcf_aero_s2,   &
       !LW diagnostics
       toplwad_aero_s2, sollwad_aero_s2,   &
       toplwai_aero_s2, sollwai_aero_s2,   &
       toplwad0_aero_s2, sollwad0_aero_s2, &
       !
       topsw_aero,solsw_aero,       &
       topsw0_aero,solsw0_aero,     &
       topswcf_aero,solswcf_aero,   &
       tausum_aero,tau3d_aero,      &
       drytausum_aero,              &
       !
       !variables CFMIP2/CMIP5
       topswad_aerop, solswad_aerop,   &
       topswai_aerop, solswai_aerop,   &
       topswad0_aerop, solswad0_aerop, &
       topsw_aerop, topsw0_aerop,      &
       solsw_aerop, solsw0_aerop,      &
       topswcf_aerop, solswcf_aerop,   &
       !LW diagnostics
       toplwad_aerop, sollwad_aerop,   &
       toplwai_aerop, sollwai_aerop,   &
       toplwad0_aerop, sollwad0_aerop, &
       !pour Ecrad
       topswad_aero_s2, solswad_aero_s2,   &
       topswai_aero_s2, solswai_aero_s2,   &
       topswad0_aero_s2, solswad0_aero_s2, &
       topsw_aero_s2, topsw0_aero_s2,      &
       solsw_aero_s2, solsw0_aero_s2,      &
       topswcf_aero_s2, solswcf_aero_s2,   &
       !LW diagnostics
       toplwad_aero_s2, sollwad_aero_s2,   &
       toplwai_aero_s2, sollwai_aero_s2,   &
       toplwad0_aero_s2, sollwad0_aero_s2, &
       !
       ptstar, pt0, slp, &
       !
       bils, &
       !
       cldh, cldl,cldm, cldq, cldt,      &
       JrNt,                             &
       dthmin, evap, snowerosion, icesub_lic, fder, plcl, plfc,   &
       prw, prlw, prsw, prbsw, water_budget,         &
       s_lcl, s_pblh, s_pblt, s_therm,   &
       cdragm, cdragh,                   &
       zustar, zu10m, zv10m, rh2m, qsat2m, &
       zq2m, zt2m, zn2mout, weak_inversion, &
       !
       s_pblh_x, s_pblh_w, &
       s_lcl_x, s_lcl_w,   &
       !
       slab_wfbils, tpot, tpote,               &
       ue, uq, ve, vq, zxffonte,               &
       uwat, vwat,                             &
       zxfqcalving, zxfluxlat,                 &
       zxrunofflic,                            &
       zxtsol, snow_lsc, zxfqfonte, zxqsurf,   &
       delta_qsurf,                            &
       rain_lsc, rain_num,                     &
       !
       sens_x, sens_w, &
       zxfluxlat_x, zxfluxlat_w, &
       !
       pbl_tke_input, pbl_eps, l_mix, wprime,&
       t_therm, q_therm, u_therm, v_therm, &
       cdragh_x, cdragh_w, &
       cdragm_x, cdragm_w, &
       kh, kh_x, kh_w, &
       !
       wake_k, &
       alp_wake, &
       wake_h, wake_omg, &
                       ! tendencies of delta T and delta q:
       d_deltat_wk, d_deltaq_wk, &         ! due to wakes
       d_deltat_wk_gw, d_deltaq_wk_gw, &   ! due to wake induced gravity waves
       d_deltat_vdf, d_deltaq_vdf, &       ! due to vertical diffusion
       d_deltat_the, d_deltaq_the, &       ! due to thermals
       d_deltat_ajs_cv, d_deltaq_ajs_cv, & ! due to dry adjustment of (w) before convection
                       ! tendencies of wake fractional area and wake number per unit area:
       d_s_wk, d_s_a_wk, d_dens_wk,  d_dens_a_wk, &  ! due to wakes
!!!       d_s_vdf, d_dens_a_vdf, d_dens_vdf, & ! due to vertical diffusion
!!!       d_s_the, d_dens_a_the, d_dens_the, & ! due to thermals
       !
       ptconv, ratqsc, &
       wbeff, convoccur, zmax_th, &
       sens, flwp, fiwp,  &
       alp_bl_conv,alp_bl_det,  &
       alp_bl_fluct_m,alp_bl_fluct_tke,  &
       alp_bl_stat, n2, s2,  strig, zcong, zlcl_th, &
       proba_notrig, random_notrig,  &
!!       cv_gen,  &  !moved to phys_state_var_mod
       !
       dnwd0,  &
       omega,  &
       epmax_diag,  &
       !    Deep convective variables used in phytrac
       pmflxr, pmflxs,  &
       coef_clos, coef_clos_eff, &
       wdtrainA, wdtrainS, wdtrainM, wdtrainAS,  &
       upwd, dnwd, &
       ep,  &
       da, mp, &
       phi, &
       wght_cvfd, &
       phi2, &
       d1a, dam, &
       ev, &
       elij, &
       qtaa, &
       clw, &
       epmlmMm, eplaMm, &
       sij, &
       !
       rneblsvol, &
       pfraclr, pfracld, & 
       cldfraliq, sigma2_icefracturb, mean_icefracturb,  &
       cldfraliqth, sigma2_icefracturbth, mean_icefracturbth,  &
       distcltop, temp_cltop,  &
       !-- LSCP - condensation and ice supersaturation variables
       qsub, qissr, qcld, subfra, issrfra, gamma_cond, &
       ql_seri_lscp, ratio_ql_qtot, qi_seri_lscp, ratio_qi_qtot, &
       dcf_sub, dcf_con, dcf_mix, dqi_adj, dqi_sub, dqi_con, dqi_mix, &
       dqvc_adj, dqvc_sub, dqvc_con, dqvc_mix, qsatliq, qsatice, &
       !-- LSCP - aviation and contrails variables
       Tcontr, qcontr, qcontr2, fcontrN, fcontrP, &
       dcf_avi, dqi_avi, dqvc_avi, flight_dist, flight_h2o, &
       !
       cldemi,  &
       cldfra, cldtau, fiwc,  &
       fl, re, flwc,  &
       ref_liq, ref_ice, theta,  &
       ref_liq_pi, ref_ice_pi,  &
       zphi, zx_rh, zx_rhl, zx_rhi,  &
       pmfd, pmfu,  &
       !
       t2m, fluxlat,  &
       fsollw, evap_pot,  &
       fsolsw, wfbils, wfevap, &
       prfl, psfl,bsfl, fraca, Vprecip,  &
       zw2,  &
       !
       fluxu, fluxv,  &
       fluxt,  &
       !
       uwriteSTD, vwriteSTD, &                !pour calcul_STDlev.h
       wwriteSTD, phiwriteSTD, &              !pour calcul_STDlev.h
       qwriteSTD, twriteSTD, rhwriteSTD, &    !pour calcul_STDlev.h
       !
       beta_prec,  &
       rneb,  &
       zxsnow,snowhgt,qsnow,to_ice,sissnow,runoff,albsol3_lic, &
!GG       zxfluxt,zxfluxq 
       zxfluxt,zxfluxq, &
       fcds, fcdi, dh_basal_growth, dh_basal_melt, &
       dh_top_melt, dh_snow2sic, &
       dtice_melt, dtice_snow2sic
!GG
       !
       USE phys_local_var_mod, ONLY: zfice, dNovrN, ptconv
       USE phys_output_var_mod, ONLY: scdnc, cldncl, reffclwtop, lcc, reffclws, &
       reffclwc, cldnvi, lcc3d, lcc3dcon, lcc3dstra, icc3dcon, icc3dstra
       USE output_physiqex_mod, ONLY: output_physiqex
       USE yomcst_mod_h
       USE clesphys_mod_h
       USE conema3_mod_h
       USE alpale_mod
       USE yoethf_mod_h
       USE calcul_divers_mod_h, ONLY: calcul_divers
       USE compbl_mod_h
       USE nuage_params_mod_h
       USE dimpft_mod_h, ONLY: nvm_lmdz
       USE radepsi_mod_h
       USE radopt_mod_h
       USE regdim_mod_h
       USE phys_output_write_spl_mod, ONLY: phys_output_write_spl
       USE phytracr_spl_mod, ONLY: phytracr_spl_out_init, phytracr_spl
       USE s2s, ONLY : s2s_initialize, s2s_finalize
    IMPLICIT NONE
    !>======================================================================
    !!
    !! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
    !!
    !! Objet: Moniteur general de la physique du modele
    !!AA      Modifications quant aux traceurs :
    !!AA                  -  uniformisation des parametrisations ds phytrac
    !!AA                  -  stockage des moyennes des champs necessaires
    !!AA                     en mode traceur off-line
    !!======================================================================
    !!   CLEFS CPP POUR LES IO
    !!   =====================
#define histNMC
    !!======================================================================
    !!    modif   ( P. Le Van ,  12/10/98 )
    !!
    !!  Arguments:
    !!
    !! nlon----input-I-nombre de points horizontaux
    !! nlev----input-I-nombre de couches verticales, doit etre egale a klev
    !! debut---input-L-variable logique indiquant le premier passage
    !! lafin---input-L-variable logique indiquant le dernier passage
    !! jD_cur       -R-jour courant a l'appel de la physique (jour julien)
    !! jH_cur       -R-heure courante a l'appel de la physique (jour julien)
    !! pdtphys-input-R-pas d'integration pour la physique (seconde)
    !! paprs---input-R-pression pour chaque inter-couche (en Pa)
    !! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
    !! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
    !! pphis---input-R-geopotentiel du sol
    !! presnivs-input_R_pressions approximat. des milieux couches ( en PA)
    !! u-------input-R-vitesse dans la direction X (de O a E) en m/s
    !! v-------input-R-vitesse Y (de S a N) en m/s
    !! t-------input-R-temperature (K)
    !! qx------input-R-humidite specifique (kg/kg) et d'autres traceurs
    !! d_t_dyn-input-R-tendance dynamique pour "t" (K/s)
    !! d_q_dyn-input-R-tendance dynamique pour "q" (kg/kg/s)
    !! d_ql_dyn-input-R-tendance dynamique pour "ql" (kg/kg/s)
    !! d_qs_dyn-input-R-tendance dynamique pour "qs" (kg/kg/s)
    !! flxmass_w -input-R- flux de masse verticale
    !! d_u-----output-R-tendance physique de "u" (m/s/s)
    !! d_v-----output-R-tendance physique de "v" (m/s/s)
    !! d_t-----output-R-tendance physique de "t" (K/s)
    !! d_qx----output-R-tendance physique de "qx" (kg/kg/s)
    !! d_ps----output-R-tendance physique de la pression au sol
    !!======================================================================
    integer jjmp1
    !  parameter (jjmp1=jjm+1-1/jjm) ! => (jjmp1=nbp_lat-1/(nbp_lat-1))
    !  integer iip1
    !  parameter (iip1=iim+1)

    !======================================================================
    LOGICAL, SAVE :: ok_volcan ! pour activer les diagnostics volcaniques
    !$OMP THREADPRIVATE(ok_volcan)
    INTEGER, SAVE :: flag_volc_surfstrat ! pour imposer le cool/heat rate à la surf/strato
    !$OMP THREADPRIVATE(flag_volc_surfstrat)
    LOGICAL ok_cvl  ! pour activer le nouveau driver pour convection KE
    PARAMETER (ok_cvl=.TRUE.)
    LOGICAL ok_gust ! pour activer l'effet des gust sur flux surface
    PARAMETER (ok_gust=.FALSE.)
    INTEGER, SAVE :: iflag_radia     ! active ou non le rayonnement (MPL)
    !$OMP THREADPRIVATE(iflag_radia)
    !======================================================================
    LOGICAL check ! Verifier la conservation du modele en eau
    PARAMETER (check=.FALSE.)
    LOGICAL ok_stratus ! Ajouter artificiellement les stratus
    PARAMETER (ok_stratus=.FALSE.)
    !======================================================================
    REAL amn, amx
    INTEGER igout
    !======================================================================
    ! Clef iflag_cycle_diurne controlant l'activation du cycle diurne:
    ! en attente du codage des cles par Fred
    ! iflag_cycle_diurne est initialise par conf_phys et se trouve
    ! dans clesphys.h (IM)
    !======================================================================
    ! Modele thermique du sol, a activer pour le cycle diurne:
    !cc      LOGICAL soil_model
    !cc      PARAMETER (soil_model=.FALSE.)
    !======================================================================
    ! Clefs controlant deux parametrisations de l'orographie:
    !c      LOGICAL ok_orodr
    !cc      PARAMETER (ok_orodr=.FALSE.)
    !cc      LOGICAL ok_orolf
    !cc      PARAMETER (ok_orolf=.FALSE.)
    !======================================================================
    LOGICAL ok_journe ! sortir le fichier journalier
    SAVE ok_journe
    !$OMP THREADPRIVATE(ok_journe)
    !
    LOGICAL ok_mensuel ! sortir le fichier mensuel
    SAVE ok_mensuel
    !$OMP THREADPRIVATE(ok_mensuel)
    !
    LOGICAL ok_instan ! sortir le fichier instantane
    SAVE ok_instan
    !$OMP THREADPRIVATE(ok_instan)
    !
    LOGICAL ok_LES ! sortir le fichier LES
    SAVE ok_LES
    !$OMP THREADPRIVATE(ok_LES)
    !
    LOGICAL callstats ! sortir le fichier stats
    SAVE callstats
    !$OMP THREADPRIVATE(callstats)
    !
    LOGICAL ok_region ! sortir le fichier regional
    PARAMETER (ok_region=.FALSE.)
    !======================================================================
    REAL seuil_inversion
    SAVE seuil_inversion
    !$OMP THREADPRIVATE(seuil_inversion)



    real facteur

    REAL wmax_th(klon)
    REAL tau_overturning_th(klon)

    INTEGER lmax_th(klon)
    INTEGER limbas(klon)
    REAL ratqscth(klon,klev)
    REAL ratqsdiff(klon,klev)
    REAL zqsatth(klon,klev)

    !======================================================================
    !
    !
    ! Variables argument:
    !
    INTEGER nlon
    INTEGER nlev
    REAL,INTENT(IN) :: pdtphys_
    ! NB: pdtphys to be used in physics is in time_phylmdz_mod
    LOGICAL debut, lafin
    REAL paprs(klon,klev+1)
    REAL pplay(klon,klev)
    REAL pphi(klon,klev)
    REAL pphis(klon)
    REAL presnivs(klev)
!JLD    REAL znivsig(klev)
!JLD    real pir

    REAL u(klon,klev)
    REAL v(klon,klev)

    REAL, intent(in):: rot(klon, klev)
    ! relative vorticity, in s-1, needed for frontal waves

    REAL t(klon,klev),thetal(klon,klev)
    ! thetal: ligne suivante a decommenter si vous avez les fichiers
    !     MPL 20130625
    ! fth_fonctions.F90 et parkind1.F90
    ! sinon thetal=theta
    !     REAL fth_thetae,fth_thetav,fth_thetal
    REAL qx(klon,klev,nqtot)
    REAL flxmass_w(klon,klev)
    REAL d_u(klon,klev)
    REAL d_v(klon,klev)
    REAL d_t(klon,klev)
    REAL d_qx(klon,klev,nqtot)
    REAL d_ps(klon)
  ! variables pour tend_to_tke
    REAL duadd(klon,klev)
    REAL dvadd(klon,klev)
    REAL dtadd(klon,klev)

!!   Variables moved to phys_local_var_mod
!!    ! Variables pour le transport convectif
!!    real da(klon,klev),phi(klon,klev,klev),mp(klon,klev)
!!    real wght_cvfd(klon,klev)
!!    ! Variables pour le lessivage convectif
!!    ! RomP >>>
!!    real phi2(klon,klev,klev)
!!    real d1a(klon,klev),dam(klon,klev)
!!    real ev(klon,klev)
!!    real clw(klon,klev),elij(klon,klev,klev)
!!    real epmlmMm(klon,klev,klev),eplaMm(klon,klev)
!!    ! RomP <<<
    !IM definition dynamique o_trac dans phys_output_open
    !      type(ctrl_out) :: o_trac(nqtot)

    ! variables a une pression donnee
    !
    include "declare_STDlev.h"

    INTEGER n
    !ym      INTEGER npoints
    !ym      PARAMETER(npoints=klon)
    !
    INTEGER nregISCtot
    PARAMETER(nregISCtot=1)
    !
    ! imin_debut, nbpti, jmin_debut, nbptj : parametres pour sorties
    ! sur 1 region rectangulaire y compris pour 1 point
    ! imin_debut : indice minimum de i; nbpti : nombre de points en
    ! direction i (longitude)
    ! jmin_debut : indice minimum de j; nbptj : nombre de points en
    ! direction j (latitude)
!JLD    INTEGER imin_debut, nbpti
!JLD    INTEGER jmin_debut, nbptj
    !IM: region='3d' <==> sorties en global
    CHARACTER*3 region
    PARAMETER(region='3d')
    LOGICAL ok_hf
    !
    SAVE ok_hf
    !$OMP THREADPRIVATE(ok_hf)

    INTEGER, PARAMETER :: longcles=20
    REAL, SAVE :: clesphy0(longcles)
    !$OMP THREADPRIVATE(clesphy0)
    !
    ! Variables propres a la physique
    INTEGER, SAVE :: itap         ! compteur pour la physique
    !$OMP THREADPRIVATE(itap)

    INTEGER, SAVE :: abortphy=0   ! Reprere si on doit arreter en fin de phys
    !$OMP THREADPRIVATE(abortphy)
    !
    REAL,SAVE ::  solarlong0
    !$OMP THREADPRIVATE(solarlong0)

    !
    !  Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
    !
    !IM 141004     REAL zulow(klon),zvlow(klon),zustr(klon), zvstr(klon)
    REAL zulow(klon),zvlow(klon)
    !
    INTEGER igwd,idx(klon),itest(klon)
    !
    !      REAL,allocatable,save :: run_off_lic_0(:)
    ! !$OMP THREADPRIVATE(run_off_lic_0)
    !ym      SAVE run_off_lic_0
    !KE43
    ! Variables liees a la convection de K. Emanuel (sb):
    !
    REAL, SAVE :: bas, top             ! cloud base and top levels
    !$OMP THREADPRIVATE(bas, top)
    !------------------------------------------------------------------
    ! Upmost level reached by deep convection and related variable (jyg)
    !
!    INTEGER izero
    INTEGER k_upper_cv
    !------------------------------------------------------------------
    ! Compteur de l'occurence de cvpas=1
    INTEGER Ncvpaseq1
    SAVE Ncvpaseq1
    !$OMP THREADPRIVATE(Ncvpaseq1)
    !
    !==========================================================================
    !CR04.12.07: on ajoute les nouvelles variables du nouveau schema
    !de convection avec poches froides
    ! Variables li\'ees \`a la poche froide (jyg)

!!    REAL mipsh(klon,klev)  ! mass flux shed by the adiab ascent at each level
!!      Moved to phys_state_var_mod
    !
    REAL wape_prescr, fip_prescr
    INTEGER it_wape_prescr
    SAVE wape_prescr, fip_prescr, it_wape_prescr
    !$OMP THREADPRIVATE(wape_prescr, fip_prescr, it_wape_prescr)
    !
    ! variables supplementaires de concvl
    REAL Tconv(klon,klev)
!!    variable moved to phys_local_var_mod
!!    REAL sij(klon,klev,klev)
!!    !
!!    ! variables pour tester la conservation de l'energie dans concvl
!!    REAL, DIMENSION(klon,klev)     :: d_t_con_sat
!!    REAL, DIMENSION(klon,klev)     :: d_q_con_sat
!!    REAL, DIMENSION(klon,klev)     :: dql_sat

    REAL, SAVE :: alp_bl_prescr=0.
    REAL, SAVE :: ale_bl_prescr=0.
    REAL, SAVE :: wake_s_min_lsp=0.1
    !$OMP THREADPRIVATE(alp_bl_prescr,ale_bl_prescr)
    !$OMP THREADPRIVATE(wake_s_min_lsp)

    REAL ok_wk_lsp(klon)

    !RC
    ! Variables li\'ees \`a la poche froide (jyg et rr)

    INTEGER,  SAVE               :: iflag_wake_tend  ! wake: if =0, then wake state variables are
                                                     ! updated within calwake
    !$OMP THREADPRIVATE(iflag_wake_tend)
    INTEGER,  SAVE               :: iflag_alp_wk_cond=0 ! wake: if =0, then Alp_wk is the average lifting
                                                        ! power provided by the wakes; else, Alp_wk is the
                                                        ! lifting power conditionned on the presence of a
                                                        ! gust-front in the grid cell.
    !$OMP THREADPRIVATE(iflag_alp_wk_cond)

    REAL t_w(klon,klev),q_w(klon,klev) ! temperature and moisture profiles in the wake region
    REAL t_x(klon,klev),q_x(klon,klev) ! temperature and moisture profiles in the off-wake region

    REAL wake_dth(klon,klev)        ! wake : temp pot difference

    REAL wake_omgbdth(klon,klev)    ! Wake : flux of Delta_Theta
    ! transported by LS omega
    REAL wake_dp_omgb(klon,klev)    ! Wake : vertical gradient of
    ! large scale omega
    REAL wake_dtKE(klon,klev)       ! Wake : differential heating
    ! (wake - unpertubed) CONV
    REAL wake_dqKE(klon,klev)       ! Wake : differential moistening
    ! (wake - unpertubed) CONV
    REAL wake_dp_deltomg(klon,klev) ! Wake : gradient vertical de wake_omg
    REAL wake_spread(klon,klev)     ! spreading term in wake_delt
    !
    !pourquoi y'a pas de save??
    !
!!!    INTEGER, SAVE, DIMENSION(klon)   :: wake_k
!!!    !$OMP THREADPRIVATE(wake_k)
    !
    !jyg<
    !cc      REAL wake_pe(klon)              ! Wake potential energy - WAPE
    !>jyg

    REAL wake_fip_0(klon)           ! Average Front Incoming Power (unconditionned)
    REAL wake_gfl(klon)             ! Gust Front Length
!!!    REAL wake_dens(klon)         ! moved to phys_state_var_mod
    !
    !
    REAL dt_dwn(klon,klev)
    REAL dq_dwn(klon,klev)
    REAL M_dwn(klon,klev)
    REAL M_up(klon,klev)
    REAL dt_a(klon,klev)
    REAL dq_a(klon,klev)
    REAL d_t_adjwk(klon,klev)                !jyg
    REAL d_q_adjwk(klon,klev)                !jyg
    LOGICAL,SAVE :: ok_adjwk=.FALSE.
    !$OMP THREADPRIVATE(ok_adjwk)
    INTEGER,SAVE :: iflag_adjwk=0            !jyg
    !$OMP THREADPRIVATE(iflag_adjwk)         !jyg
    REAL,SAVE :: oliqmax=999.,oicemax=999.
    !$OMP THREADPRIVATE(oliqmax,oicemax)
    REAL, SAVE :: alp_offset
    !$OMP THREADPRIVATE(alp_offset)
    REAL, SAVE :: dtcon_multistep_max=1.e6
    !$OMP THREADPRIVATE(dtcon_multistep_max)
    REAL, SAVE :: dqcon_multistep_max=1.e6
    !$OMP THREADPRIVATE(dqcon_multistep_max)

    !
    !RR:fin declarations poches froides
    !==========================================================================

    REAL ztv(klon,klev),ztva(klon,klev)
    REAL zpspsk(klon,klev)
    REAL ztla(klon,klev),zqla(klon,klev)
    REAL zthl(klon,klev)

    !cc nrlmd le 10/04/2012

    !--------Stochastic Boundary Layer Triggering: ALE_BL--------
    !---Propri\'et\'es du thermiques au LCL
!    real zlcl_th(klon)          ! Altitude du LCL calcul\'e
    ! continument (pcon dans
    ! thermcell_main.F90)
    real fraca0(klon)           ! Fraction des thermiques au LCL
    real w0(klon)               ! Vitesse des thermiques au LCL
    real w_conv(klon)           ! Vitesse verticale de grande \'echelle au LCL
    real tke0(klon,klev+1)      ! TKE au d\'ebut du pas de temps
    real therm_tke_max0(klon)   ! TKE dans les thermiques au LCL
    real env_tke_max0(klon)     ! TKE dans l'environnement au LCL
    INTEGER, SAVE :: iflag_thermcell_tke ! transtport TKE by thermals
    !$OMP THREADPRIVATE(iflag_thermcell_tke)

!JLD    !---D\'eclenchement stochastique
!JLD    integer :: tau_trig(klon)

    REAL,SAVE :: random_notrig_max=1.
    !$OMP THREADPRIVATE(random_notrig_max)

    !--------Statistical Boundary Layer Closure: ALP_BL--------
    !---Profils de TKE dans et hors du thermique
    real therm_tke_max(klon,klev)   ! Profil de TKE dans les thermiques
    real env_tke_max(klon,klev)     ! Profil de TKE dans l'environnement

    !-------Activer les tendances de TKE due a l'orograp??ie---------
     INTEGER, SAVE :: addtkeoro
    !$OMP THREADPRIVATE(addtkeoro)
     REAL, SAVE :: alphatkeoro
    !$OMP THREADPRIVATE(alphatkeoro)
     LOGICAL, SAVE :: smallscales_tkeoro
    !$OMP THREADPRIVATE(smallscales_tkeoro)



    !cc fin nrlmd le 10/04/2012

    ! Variables locales pour la couche limite (al1):
    !
    !Al1      REAL pblh(klon)           ! Hauteur de couche limite
    !Al1      SAVE pblh
    !34EK
    !
    ! Variables locales:
    !
    !AA
    !AA  Pour phytrac
    REAL u1(klon)             ! vents dans la premiere couche U
    REAL v1(klon)             ! vents dans la premiere couche V

    !@$$      LOGICAL offline           ! Controle du stockage ds "physique"
    !@$$      PARAMETER (offline=.false.)
    !@$$      INTEGER physid
    REAL frac_impa(klon,klev) ! fractions d'aerosols lessivees (impaction)
    REAL frac_nucl(klon,klev) ! idem (nucleation)
    ! RomP >>>
    REAL beta_prec_fisrt(klon,klev) ! taux de conv de l'eau cond (fisrt)
    ! RomP <<<

    !IM cf FH pour Tiedtke 080604
    REAL rain_tiedtke(klon),snow_tiedtke(klon)
    !
    !IM 050204 END
    REAL devap(klon) ! evaporation et sa derivee
    REAL dsens(klon) ! chaleur sensible et sa derivee

    !
    ! Conditions aux limites
    !
    !
    REAL :: day_since_equinox
    ! Date de l'equinoxe de printemps
    INTEGER, parameter :: mth_eq=3, day_eq=21
    REAL :: jD_eq

    LOGICAL, parameter :: new_orbit = .TRUE.

    !
    INTEGER lmt_pas
    SAVE lmt_pas                ! frequence de mise a jour
    !$OMP THREADPRIVATE(lmt_pas)
    real zmasse(klon, nbp_lev),exner(klon, nbp_lev)
    !     (column-density of mass of air in a cell, in kg m-2)
    real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

    !IM sorties
    REAL un_jour
    PARAMETER(un_jour=86400.)
    INTEGER itapm1 !pas de temps de la physique du(es) mois precedents
    SAVE itapm1    !mis a jour le dernier pas de temps du mois en cours
    !$OMP THREADPRIVATE(itapm1)
    !======================================================================
    !
    ! Declaration des procedures appelees
    !
    ! EXTERNAL angle     ! calculer angle zenithal du soleil
    EXTERNAL alboc     ! calculer l'albedo sur ocean
    ! EXTERNAL ajsec     ! ajustement sec
    EXTERNAL conlmd    ! convection (schema LMD)
    EXTERNAL conema3  ! convect4.3
    EXTERNAL hgardfou  ! verifier les temperatures
    ! EXTERNAL nuage     ! calculer les proprietes radiatives
    !C      EXTERNAL o3cm      ! initialiser l'ozone
    ! EXTERNAL orbite    ! calculer l'orbite terrestre
    EXTERNAL phyredem  ! ecrire l'etat de redemarrage de la physique
    EXTERNAL suphel    ! initialiser certaines constantes
    ! EXTERNAL transp    ! transport total de l'eau et de l'energie
    !IM
    EXTERNAL haut2bas  !variables de haut en bas
    EXTERNAL ini_undefSTD  !initialise a 0 une variable a 1 niveau de pression
    EXTERNAL undefSTD !somme les valeurs definies d'1 var a 1 niveau de pression
    !     EXTERNAL moy_undefSTD  !moyenne d'1 var a 1 niveau de pression
    ! EXTERNAL moyglo_aire
    ! moyenne globale d'1 var ponderee par l'aire de la maille (moyglo_pondaire)
    ! par la masse/airetot (moyglo_pondaima) et la vraie masse (moyglo_pondmass)
    !
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local variables
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
!    REAL rhcl(klon,klev)    ! humiditi relative ciel clair
    REAL dialiq(klon,klev)  ! eau liquide nuageuse
    REAL diafra(klon,klev)  ! fraction nuageuse
    REAL radocond(klon,klev)  ! eau condensee nuageuse
    !
    !XXX PB
    REAL fluxq(klon,klev, nbsrf)   ! flux turbulent d'humidite
    REAL fluxqbs(klon,klev, nbsrf)   ! flux turbulent de neige soufflee
    !
    !FC    REAL zxfluxt(klon, klev)
    !FC    REAL zxfluxq(klon, klev)
    REAL zxfluxqbs(klon,klev)
    REAL zxfluxu(klon, klev)
    REAL zxfluxv(klon, klev)

    ! Le rayonnement n'est pas calcule tous les pas, il faut donc
    !                      sauvegarder les sorties du rayonnement
    !ym      SAVE  heat,cool,albpla,topsw,toplw,solsw,sollw,sollwdown
    !ym      SAVE  sollwdownclr, toplwdown, toplwdownclr
    !ym      SAVE  topsw0,toplw0,solsw0,sollw0, heat0, cool0
    !
    INTEGER itaprad
    SAVE itaprad
    !$OMP THREADPRIVATE(itaprad)
    !
    REAL conv_q(klon,klev) ! convergence de l'humidite (kg/kg/s)
    REAL conv_t(klon,klev) ! convergence de la temperature(K/s)
    !
    REAL zsav_tsol(klon)
    !
    REAL dist, rmu0(klon), fract(klon)
    REAL zrmu0(klon), zfract(klon)
    REAL zdtime, zdtime1, zdtime2, zlongi
    !
    REAL z_avant(klon), z_apres(klon), z_factor(klon)
    LOGICAL zx_ajustq
    !
    REAL za
    REAL zx_t, zx_qs, zdelta, zcor
    real zqsat(klon,klev)
    !
    INTEGER i, k, iq, nsrf, l, itr
    !
    REAL t_coup
    PARAMETER (t_coup=234.0)

    !ym A voir plus tard !!
    !ym      REAL zx_relief(iim,jjmp1)
    !ym      REAL zx_aire(iim,jjmp1)
    !
    ! Grandeurs de sorties
    REAL s_capCL(klon)
    REAL s_oliqCL(klon), s_cteiCL(klon)
    REAL s_trmb1(klon), s_trmb2(klon)
    REAL s_trmb3(klon)

    ! La convection n'est pas calculee tous les pas, il faut donc
    !                      sauvegarder les sorties de la convection
    !ym      SAVE
    !ym      SAVE
    !ym      SAVE
    !
    INTEGER itapcv, itapwk
    SAVE itapcv, itapwk
    !$OMP THREADPRIVATE(itapcv, itapwk)

    !KE43
    ! Variables locales pour la convection de K. Emanuel (sb):

    REAL tvp(klon,klev)       ! virtual temp of lifted parcel
    CHARACTER*40 capemaxcels  !max(CAPE)

    REAL rflag(klon)          ! flag fonctionnement de convect
    INTEGER iflagctrl(klon)          ! flag fonctionnement de convect

    ! -- convect43:
    INTEGER ntra              ! nb traceurs pour convect4.3
    REAL dtvpdt1(klon,klev), dtvpdq1(klon,klev)
    REAL dplcldt(klon), dplcldr(klon)
    !?     .     condm_con(klon,klev),conda_con(klon,klev),
    !?     .     mr_con(klon,klev),ep_con(klon,klev)
    !?     .    ,sadiab(klon,klev),wadiab(klon,klev)
    ! --
    !34EK
    !
    ! Variables du changement
    !
    ! con: convection
    ! lsc: condensation a grande echelle (Large-Scale-Condensation)
    ! ajs: ajustement sec
    ! eva: evaporation de l'eau liquide nuageuse
    ! vdf: couche limite (Vertical DiFfusion)
    !
    ! tendance nulles
    REAL, dimension(klon,klev):: du0, dv0, dt0, dq0, dql0, dqi0, dqbs0
    REAL, dimension(klon)     :: dsig0, ddens0
    INTEGER, dimension(klon)  :: wkoccur1
    ! tendance buffer pour appel de add_phys_tend
    REAL, DIMENSION(klon,klev)  :: d_q_ch4_dtime
    !
    ! Flag pour pouvoir ne pas ajouter les tendances.
    ! Par defaut, les tendances doivente etre ajoutees et
    ! flag_inhib_tend = 0
    ! flag_inhib_tend > 0 : tendances non ajoutees, avec un nombre
    ! croissant de print quand la valeur du flag augmente
    !!! attention, ce flag doit etre change avec prudence !!!
    INTEGER :: flag_inhib_tend = 0 !  0 is the default value
!!    INTEGER :: flag_inhib_tend = 2
    !
    ! Logical switch to a bug : reseting to 0 convective variables at the
    ! begining of physiq.
    LOGICAL, SAVE :: ok_bug_cv_trac = .TRUE.
    !$OMP THREADPRIVATE(ok_bug_cv_trac)
    !
    ! Logical switch to a bug : changing wake_deltat when thermals are active
    ! even when there are no wakes.
    LOGICAL, SAVE :: ok_bug_split_th = .TRUE.
    !$OMP THREADPRIVATE(ok_bug_split_th)

    ! Logical switch to a bug : modifying directly wake_deltat  by adding
    ! the (w) dry adjustment tendency to wake_deltat
    LOGICAL, SAVE :: ok_bug_ajs_cv = .TRUE.
    !$OMP THREADPRIVATE(ok_bug_ajs_cv)

    !
    !********************************************************
    !     declarations

    !********************************************************
    !IM 081204 END
    !
    REAL pen_u(klon,klev), pen_d(klon,klev)
    REAL pde_u(klon,klev), pde_d(klon,klev)
    INTEGER kcbot(klon), kctop(klon), kdtop(klon)
    !
    REAL ratqsbas,ratqshaut,tau_ratqs
    SAVE ratqsbas,ratqshaut,tau_ratqs
    !$OMP THREADPRIVATE(ratqsbas,ratqshaut,tau_ratqs)
    REAL, SAVE :: ratqsp0=50000., ratqsdp=20000.
    !$OMP THREADPRIVATE(ratqsp0, ratqsdp)

    ! Parametres lies au nouveau schema de nuages (SB, PDF)
    REAL, SAVE :: fact_cldcon
    REAL, SAVE :: facttemps
    !$OMP THREADPRIVATE(fact_cldcon,facttemps)
    LOGICAL, SAVE :: ok_newmicro
    !$OMP THREADPRIVATE(ok_newmicro)

    INTEGER, SAVE :: iflag_cld_th
    !$OMP THREADPRIVATE(iflag_cld_th)
!IM logical ptconv(klon,klev)  !passe dans phys_local_var_mod
    !IM cf. AM 081204 BEG
    LOGICAL ptconvth(klon,klev)

    REAL picefra(klon,klev)
    REAL nm_oro(klon)
    !IM cf. AM 081204 END
    !
    ! Variables liees a l'ecriture de la bande histoire physique
    !
    !======================================================================
    !
    !
!JLD    integer itau_w   ! pas de temps ecriture = itap + itau_phy
    !
    !
    ! Variables locales pour effectuer les appels en serie
    !
    !IM RH a 2m (la surface)
    REAL Lheat

    INTEGER        length
    PARAMETER    ( length = 100 )
    REAL tabcntr0( length       )
    !
!JLD    INTEGER ndex2d(nbp_lon*nbp_lat)
    !IM
    !
    !IM AMIP2 BEG
!JLD    REAL moyglo, mountor
    !IM 141004 BEG
    REAL zustrdr(klon), zvstrdr(klon)
    REAL zustrli(klon), zvstrli(klon)
    REAL zustrph(klon), zvstrph(klon)
    REAL aam, torsfc
    !IM 141004 END
    !IM 190504 BEG
    !  INTEGER imp1jmp1
    !  PARAMETER(imp1jmp1=(iim+1)*jjmp1)
    !ym A voir plus tard
    !  REAL zx_tmp((nbp_lon+1)*nbp_lat)
    !  REAL airedyn(nbp_lon+1,nbp_lat)
    !IM 190504 END
!JLD    LOGICAL ok_msk
!JLD    REAL msk(klon)
    !ym A voir plus tard
    !ym      REAL zm_wo(jjmp1, klev)
    !IM AMIP2 END
    !
    REAL zx_tmp_fi2d(klon)      ! variable temporaire grille physique
    REAL zx_tmp_fi3d(klon,klev) ! variable temporaire pour champs 3D
!JLD    REAL zx_tmp_2d(nbp_lon,nbp_lat)
!JLD    REAL zx_lon(nbp_lon,nbp_lat)
!JLD    REAL zx_lat(nbp_lon,nbp_lat)
    !
    INTEGER nid_ctesGCM
    SAVE nid_ctesGCM
    !$OMP THREADPRIVATE(nid_ctesGCM)
    !
    !IM 280405 BEG
    !  INTEGER nid_bilKPins, nid_bilKPave
    !  SAVE nid_bilKPins, nid_bilKPave
    !  !$OMP THREADPRIVATE(nid_bilKPins, nid_bilKPave)
    !
    REAL ve_lay(klon,klev) ! transport meri. de l'energie a chaque niveau vert.
    REAL vq_lay(klon,klev) ! transport meri. de l'eau a chaque niveau vert.
    REAL ue_lay(klon,klev) ! transport zonal de l'energie a chaque niveau vert.
    REAL uq_lay(klon,klev) ! transport zonal de l'eau a chaque niveau vert.
    !
!JLD    REAL zjulian
!JLD    SAVE zjulian
!JLD!$OMP THREADPRIVATE(zjulian)

!JLD    INTEGER nhori, nvert
!JLD    REAL zsto
!JLD    REAL zstophy, zout

    CHARACTER (LEN=20) :: modname='physiq_mod'
    CHARACTER*80 abort_message
    LOGICAL, SAVE ::  ok_sync, ok_sync_omp
    !$OMP THREADPRIVATE(ok_sync)
    ! ok_sync_omp should not be in a THREADPRIVATE statement
    REAL date0

    ! essai writephys
    INTEGER fid_day, fid_mth, fid_ins
    PARAMETER (fid_ins = 1, fid_day = 2, fid_mth = 3)
    INTEGER prof2d_on, prof3d_on, prof2d_av, prof3d_av
    PARAMETER (prof2d_on = 1, prof3d_on = 2, prof2d_av = 3, prof3d_av = 4)
    REAL ztsol(klon)
    REAL q2m(klon,nbsrf)  ! humidite a 2m
    REAL qbsfra  ! blowing snow fraction
    !IM: t2m, q2m, ustar, u10m, v10m et t2mincels, t2maxcels
    CHARACTER*40 t2mincels, t2maxcels       !t2m min., t2m max
    CHARACTER*40 tinst, tave
    REAL cldtaupi(klon,klev) ! Cloud optical thickness for
    ! pre-industrial (pi) aerosols

    INTEGER :: naero
    ! Aerosol optical properties
    CHARACTER*4, DIMENSION(naero_grp) :: rfname
    REAL, DIMENSION(klon,klev)     :: mass_solu_aero ! total mass
    ! concentration
    ! for all soluble
    ! aerosols[ug/m3]
    REAL, DIMENSION(klon,klev)     :: mass_solu_aero_pi
    ! - " - (pre-industrial value)
    REAL, DIMENSION(klon,klev,naero_tot) :: m_allaer

    ! Parameters
    LOGICAL ok_ade, ok_aie    ! Apply aerosol (in)direct effects or not
    LOGICAL ok_alw            ! Apply aerosol LW effect or not
    LOGICAL ok_cdnc ! ok cloud droplet number concentration (O. Boucher 01-2013)
    REAL bl95_b0, bl95_b1   ! Parameter in Boucher and Lohmann (1995)
    SAVE ok_ade, ok_aie, ok_alw, ok_cdnc, bl95_b0, bl95_b1
    !$OMP THREADPRIVATE(ok_ade, ok_aie, ok_alw, ok_cdnc, bl95_b0, bl95_b1)
    LOGICAL, SAVE :: aerosol_couple ! true  : calcul des aerosols dans INCA
    ! false : lecture des aerosol dans un fichier
    !$OMP THREADPRIVATE(aerosol_couple)
    LOGICAL, SAVE :: chemistry_couple ! true  : use INCA chemistry O3
    ! false : use offline chemistry O3
    !$OMP THREADPRIVATE(chemistry_couple)
    INTEGER, SAVE :: flag_aerosol
    !$OMP THREADPRIVATE(flag_aerosol)
    LOGICAL, SAVE :: flag_bc_internal_mixture
    !$OMP THREADPRIVATE(flag_bc_internal_mixture)
    !
    !--STRAT AEROSOL
    INTEGER, SAVE :: flag_aerosol_strat
    !$OMP THREADPRIVATE(flag_aerosol_strat)
    !
    !--INTERACTIVE AEROSOL FEEDBACK ON RADIATION
    LOGICAL, SAVE :: flag_aer_feedback
    !$OMP THREADPRIVATE(flag_aer_feedback)

    !c-fin STRAT AEROSOL
    !
    ! Declaration des constantes et des fonctions thermodynamiques
    !
    LOGICAL,SAVE :: first=.TRUE.
    !$OMP THREADPRIVATE(first)

    ! VARIABLES RELATED TO OZONE CLIMATOLOGIES ; all are OpenMP shared
    ! Note that pressure vectors are in Pa and in stricly ascending order
    INTEGER,SAVE :: read_climoz                ! Read ozone climatology
    !     (let it keep the default OpenMP shared attribute)
    !     Allowed values are 0, 1 and 2
    !     0: do not read an ozone climatology
    !     1: read a single ozone climatology that will be used day and night
    !     2: read two ozone climatologies, the average day and night
    !     climatology and the daylight climatology
    INTEGER,SAVE :: ncid_climoz                ! NetCDF file identifier
    REAL, ALLOCATABLE, SAVE :: press_cen_climoz(:) ! Pressure levels
    REAL, ALLOCATABLE, SAVE :: press_edg_climoz(:) ! Edges of pressure intervals
    REAL, ALLOCATABLE, SAVE :: time_climoz(:)      ! Time vector

    CHARACTER(LEN=13), PARAMETER :: vars_climoz(2) &
                                  = ["tro3         ","tro3_daylight"]
    ! vars_climoz(1:read_climoz): variables names in climoz file.
    ! vars_climoz(1:read_climoz-2) if read_climoz>2 (temporary)
    REAL :: ro3i ! 0<=ro3i<=360 ; required time index in NetCDF file for
                 ! the ozone fields, old method.

    include "FCTTRE.h"
    !IM 100106 END : pouvoir sortir les ctes de la physique
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Declarations pour Simulateur COSP
    !============================================================
    ! AI 10-22
    include "ini_COSP.h"
    real :: mr_ozone(klon,klev), phicosp(klon,klev)

    !IM stations CFMIP
    INTEGER, SAVE :: nCFMIP
    !$OMP THREADPRIVATE(nCFMIP)
    INTEGER, PARAMETER :: npCFMIP=120
    INTEGER, ALLOCATABLE, SAVE :: tabCFMIP(:)
    REAL, ALLOCATABLE, SAVE :: lonCFMIP(:), latCFMIP(:)
    !$OMP THREADPRIVATE(tabCFMIP, lonCFMIP, latCFMIP)
    INTEGER, ALLOCATABLE, SAVE :: tabijGCM(:)
    REAL, ALLOCATABLE, SAVE :: lonGCM(:), latGCM(:)
    !$OMP THREADPRIVATE(tabijGCM, lonGCM, latGCM)
    INTEGER, ALLOCATABLE, SAVE :: iGCM(:), jGCM(:)
    !$OMP THREADPRIVATE(iGCM, jGCM)
    logical, dimension(nfiles)            :: phys_out_filestations
    logical, parameter :: lNMC=.FALSE.

    !IM betaCRF
    REAL, SAVE :: pfree, beta_pbl, beta_free
    !$OMP THREADPRIVATE(pfree, beta_pbl, beta_free)
    REAL, SAVE :: lon1_beta,  lon2_beta, lat1_beta, lat2_beta
    !$OMP THREADPRIVATE(lon1_beta,  lon2_beta, lat1_beta, lat2_beta)
    LOGICAL, SAVE :: mskocean_beta
    !$OMP THREADPRIVATE(mskocean_beta)
    REAL, dimension(klon, klev) :: beta ! facteur sur cldtaurad et
    ! cldemirad pour evaluer les
    ! retros liees aux CRF
    REAL, dimension(klon, klev) :: cldtaurad   ! epaisseur optique
    ! pour radlwsw pour
    ! tester "CRF off"
    REAL, dimension(klon, klev) :: cldtaupirad ! epaisseur optique
    ! pour radlwsw pour
    ! tester "CRF off"
    REAL, dimension(klon, klev) :: cldemirad   ! emissivite pour
    ! radlwsw pour tester
    ! "CRF off"
    REAL, dimension(klon, klev) :: cldfrarad   ! fraction nuageuse

    REAL :: calday, zxsnow_dummy(klon)
    ! set de variables utilisees pour l'initialisation des valeurs provenant de INCA
    REAL, DIMENSION(klon,klev,naero_grp,nbands) :: init_tauinca
    REAL, DIMENSION(klon,klev,naero_grp,nbands) :: init_pizinca
    REAL, DIMENSION(klon,klev,naero_grp,nbands) :: init_cginca
    REAL, DIMENSION(klon,klev,nbands) :: init_ccminca
    REAL, DIMENSION(klon,nbtr) :: init_source

    !lwoff=y : offset LW CRE for radiation code and other schemes
    REAL, SAVE :: betalwoff
    !$OMP THREADPRIVATE(betalwoff)
!
    INTEGER :: nbtr_tmp ! Number of tracer inside concvl
    REAL, dimension(klon,klev) :: sh_in ! Specific humidity entering in phytrac
    REAL, dimension(klon,klev) :: ch_in ! Condensed humidity entering in phytrac (eau liquide)
    integer iostat

    REAL, dimension(klon,klev+1) :: l_mix_ave, wprime_ave
    REAL zzz
    !albedo SB >>>
    REAL,DIMENSION(6), SAVE :: SFRWL
!$OMP THREADPRIVATE(SFRWL)
    !albedo SB <<<

    !--OB variables for mass fixer (hard coded for now)
    REAL qql1(klon),qql2(klon),corrqql

    REAL, dimension(klon,klev) :: t_env,q_env

    REAL, dimension(klon) :: pr_et
    REAL, dimension(klon) :: w_et, jlr_g_c, jlr_g_s

    REAL pi
    REAL viscom, viscoh
    INTEGER ieru

    ! For Python coupling
    INTEGER :: nstep, nyear, istep

    !AI namelist pour gerer le double appel de Ecrad
    CHARACTER(len=512) :: namelist_ecrad_file


    ! Subgrid scale wind :
    ! Need to be allocatable/save because the number of bin is not known (provided by surf_wind_ini)
    integer, save :: nsurfwind=1
    real, dimension(:,:), allocatable, save :: surf_wind_value, surf_wind_proba ! module and probability of sugrdi wind wind sample
    !$OMP THREADPRIVATE(nsurfwind,surf_wind_value, surf_wind_proba)
    


    !======================================================================!
    ! Bifurcation vers un nouveau moniteur physique pour experimenter      !
    ! des solutions et préparer le couplage avec la physique de MesoNH     !
    ! 14 mai 2023                                                          !
    !======================================================================!
    if (debut) then                                                        !
       iflag_physiq=0
       call getin_p('iflag_physiq', iflag_physiq)                          !
    endif                                                                  !
    if ( iflag_physiq == 2 ) then                                          !
       call physiqex (nlon,nlev, &                                         !
       debut,lafin,pdtphys_, &                                             !
       paprs,pplay,pphi,pphis,presnivs, &                                  !
       u,v,rot,t,qx, &                                                     !
       flxmass_w, &                                                        !
       d_u, d_v, d_t, d_qx, d_ps)                                          !
       return                                                              !
    endif                                                                  !
    !======================================================================!


    pi = 4. * ATAN(1.)

    ! set-up call to alerte function
    call_alert = (alert_first_call .AND. is_master)

    ! Ehouarn: set value of jjmp1 since it is no longer a "fixed parameter"
    jjmp1=nbp_lat

    !======================================================================
    ! Gestion calendrier : mise a jour du module phys_cal_mod
    !
    pdtphys=pdtphys_
    CALL update_time(pdtphys)
    phys_tstep=NINT(pdtphys)
    IF (.NOT. using_xios) missing_val=nf90_fill_real

    IF (using_xios) THEN
      ! switch to XIOS LMDZ physics context
      IF (.NOT. debut .AND. is_omp_master) THEN
        CALL wxios_set_context()
        CALL xios_update_calendar(itap+1)
      ENDIF
    ENDIF

    !======================================================================
    ! Ecriture eventuelle d'un profil verticale en entree de la physique.
    ! Utilise notamment en 1D mais peut etre active egalement en 3D
    ! en imposant la valeur de igout.
    !======================================================================
    IF (prt_level.ge.1) THEN
       igout=klon/2+1/klon
       write(lunout,*) 'DEBUT DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
       write(lunout,*) 'igout, lat, lon ',igout, latitude_deg(igout), &
            longitude_deg(igout)
       write(lunout,*) &
            'nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys'
       write(lunout,*) &
            nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys

       write(lunout,*) 'paprs, play, phi, u, v, t'
       DO k=1,klev
          write(lunout,*) paprs(igout,k),pplay(igout,k),pphi(igout,k), &
               u(igout,k),v(igout,k),t(igout,k)
       ENDDO
       write(lunout,*) 'ovap (g/kg),  oliq (g/kg)'
       DO k=1,klev
          write(lunout,*) qx(igout,k,1)*1000,qx(igout,k,2)*1000.
       ENDDO
    ENDIF

    ! Quick check on pressure levels:
    CALL assert(paprs(:, nbp_lev + 1) < paprs(:, nbp_lev), &
            "physiq_mod paprs bad order")

    IF (first) THEN
        
        CALL s2s_initialize     ! initialization of source to source tools
        
!       CALL init_etat0_limit_unstruct
!       IF (.NOT. create_etat0_limit) CALL init_limit_read(days_elapsed)
       !CR:nvelles variables convection/poches froides

       WRITE(lunout,*) '================================================='
       WRITE(lunout,*) 'Allocation des variables locales et sauvegardees'
       WRITE(lunout,*) '================================================='
       CALL phys_local_var_init
       !
       !     appel a la lecture du run.def physique
       CALL conf_phys(ok_journe, ok_mensuel, &
            ok_instan, ok_hf, &
            ok_LES, &
            callstats, &
            solarlong0,seuil_inversion, &
            fact_cldcon, facttemps,ok_newmicro,iflag_radia, &
            iflag_cld_th,ratqsbas,ratqshaut,tau_ratqs, &
            ok_ade, ok_aie, ok_alw, ok_cdnc, ok_volcan, flag_volc_surfstrat, aerosol_couple, &
            chemistry_couple, flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
            flag_bc_internal_mixture, bl95_b0, bl95_b1, &
                                ! nv flags pour la convection et les
                                ! poches froides
            read_climoz, &
            alp_offset)
       CALL init_etat0_limit_unstruct
       IF (.NOT. create_etat0_limit) CALL init_limit_read(days_elapsed)
       CALL phys_state_var_init(read_climoz)
       CALL phys_output_var_init
       IF (read_climoz>=1 .AND. create_etat0_limit .AND. grid_type==unstructured) &
          CALL regr_horiz_time_climoz(read_climoz,ok_daily_climoz)

       print*, '================================================='
       !
       !CR: check sur le nb de traceurs de l eau
       IF ((iflag_ice_thermo.gt.0).and.(nqo==2)) THEN
          WRITE (lunout, *) ' iflag_ice_thermo==1 requires 3 H2O tracers ', &
               '(H2O_g, H2O_l, H2O_s) but nqo=', nqo, '. Might as well stop here.'
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

       IF (ok_ice_supersat.AND.(iflag_ice_thermo.EQ.0)) THEN
          WRITE (lunout, *) ' ok_ice_supersat=y requires iflag_ice_thermo=1 as well'
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

       IF (ok_ice_supersat.AND.(nqo.LT.5)) THEN
          WRITE (lunout, *) ' ok_ice_supersat=y requires 5 H2O tracers ', &
               '(H2O_g, H2O_l, H2O_s, H2O_f, H2O_c) but nqo=', nqo, '. Might as well stop here.'
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

       IF (ok_plane_h2o.AND..NOT.ok_ice_supersat) THEN
          WRITE (lunout, *) ' ok_plane_h2o=y requires ok_ice_supersat=y '
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

       IF (ok_plane_contrail.AND..NOT.ok_ice_supersat) THEN
          WRITE (lunout, *) ' ok_plane_contrail=y requires ok_ice_supersat=y '
          abort_message='see above'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

        IF (ok_bs) THEN
         IF ((ok_ice_supersat.AND.nqo .LT.6).OR.(.NOT.ok_ice_supersat.AND.nqo.LT.4)) THEN
             WRITE (lunout, *) 'activation of blowing snow needs a specific H2O tracer', &
                               'but nqo=', nqo
             abort_message='see above'
             CALL abort_physic(modname,abort_message, 1)
         ENDIF
        ENDIF

        IF (ok_advtke) THEN
         IF (nqtke .LT. 1) THEN
             WRITE (lunout, *) 'activation of TKE advection need a specific TKE tracer', &
                               'but nqtke=', nqtke
             abort_message='see above'
             CALL abort_physic(modname,abort_message, 1)
         ENDIF
        ENDIF


       Ncvpaseq1 = 0
       dnwd0=0.0
       ftd=0.0
       fqd=0.0
       cin=0.
       !ym Attention pbase pas initialise dans concvl !!!!
       pbase=0
       !IM 180608

       itau_con=0
       first=.FALSE.

    ENDIF  ! first

    !ym => necessaire pour iflag_con != 2
    pmfd(:,:) = 0.
    pen_u(:,:) = 0.
    pen_d(:,:) = 0.
    pde_d(:,:) = 0.
    pde_u(:,:) = 0.
    aam=0.
    d_t_adjwk(:,:)=0
    d_q_adjwk(:,:)=0

    alp_bl_conv(:)=0.

    torsfc=0.
    forall (k=1: nbp_lev) zmasse(:, k) = (paprs(:, k)-paprs(:, k+1)) / rg


    IF (debut) THEN
       CALL suphel ! initialiser constantes et parametres phys.
! tau_gl : constante de rappel de la temperature a la surface de la glace - en
       tau_gl=5.
       CALL getin_p('tau_gl', tau_gl)
! tau_gl : constante de rappel de la temperature a la surface de la glace - en
! secondes
       tau_gl=86400.*tau_gl
       WRITE(lunout,*) 'debut physiq_mod tau_gl=',tau_gl
       iflag_thermcell_tke=0
       call getin_p('iflag_thermcell_tke', iflag_thermcell_tke)                          !

       CALL getin_p('iflag_alp_wk_cond', iflag_alp_wk_cond)
       CALL getin_p('random_notrig_max',random_notrig_max)
       CALL getin_p('ok_adjwk',ok_adjwk)
       IF (ok_adjwk) iflag_adjwk=2  ! for compatibility with older versions
       ! iflag_adjwk: ! 0 = Default: no convective adjustment of w-region
                      ! 1 => convective adjustment but state variables are unchanged
                      ! 2 => convective adjustment and state variables are changed
       CALL getin_p('iflag_adjwk',iflag_adjwk)
       CALL getin_p('dtcon_multistep_max',dtcon_multistep_max)
       CALL getin_p('dqcon_multistep_max',dqcon_multistep_max)
       CALL getin_p('oliqmax',oliqmax)
       CALL getin_p('oicemax',oicemax)
       CALL getin_p('ratqsp0',ratqsp0)
       CALL getin_p('ratqsdp',ratqsdp)
       iflag_wake_tend = 0
       CALL getin_p('iflag_wake_tend',iflag_wake_tend)
       ok_bad_ecmwf_thermo=.TRUE. ! By default thermodynamical constants are set
                                  ! in rrtm/suphec.F90 (and rvtmp2 is set to 0).
       CALL getin_p('ok_bad_ecmwf_thermo',ok_bad_ecmwf_thermo)
       CALL getin_p('ok_bug_cv_trac',ok_bug_cv_trac)
       CALL getin_p('ok_bug_split_th',ok_bug_split_th)
       CALL getin_p('ok_bug_ajs_cv',ok_bug_ajs_cv)
       fl_ebil = 0 ! by default, conservation diagnostics are desactivated
       CALL getin_p('fl_ebil',fl_ebil)
       fl_cor_ebil = 0 ! by default, no correction to ensure energy conservation
       CALL getin_p('fl_cor_ebil',fl_cor_ebil)
       iflag_phytrac = 1 ! by default we do want to call phytrac
       CALL getin_p('iflag_phytrac',iflag_phytrac)

       ok_water_mass_fixer=.FALSE.  ! OB: by default we do not apply the mass fixer
       CALL getin_p('ok_water_mass_fixer',ok_water_mass_fixer)
       iflag_tropopause_height = 1 ! IM: iflag_tropopause_height = 0 => PV tropopause; iflag_tropopause_height=1 => LR tropopause
       CALL getin_p('iflag_tropopause_height', iflag_tropopause_height)
IF (CPPKEY_DUST) THEN
       IF (iflag_phytrac.EQ.0) THEN
         WRITE(lunout,*) 'In order to run with SPLA, iflag_phytrac will be forced to 1'
         iflag_phytrac = 1
       ENDIF
END IF
       nvm_lmdz = 13
       CALL getin_p('NVM',nvm_lmdz)

       WRITE(lunout,*) 'iflag_alp_wk_cond=',  iflag_alp_wk_cond
       WRITE(lunout,*) 'random_ntrig_max=',   random_notrig_max
       WRITE(lunout,*) 'ok_adjwk=',           ok_adjwk
       WRITE(lunout,*) 'iflag_adjwk=',        iflag_adjwk
       WRITE(lunout,*) 'qtcon_multistep_max=',dtcon_multistep_max
       WRITE(lunout,*) 'qdcon_multistep_max=',dqcon_multistep_max
       WRITE(lunout,*) 'ratqsp0=',            ratqsp0
       WRITE(lunout,*) 'ratqsdp=',            ratqsdp
       WRITE(lunout,*) 'iflag_wake_tend=',    iflag_wake_tend
       WRITE(lunout,*) 'ok_bad_ecmwf_thermo=',ok_bad_ecmwf_thermo
       WRITE(lunout,*) 'ok_bug_cv_trac=',     ok_bug_cv_trac
       WRITE(lunout,*) 'ok_bug_split_th=',    ok_bug_split_th
       WRITE(lunout,*) 'fl_ebil=',            fl_ebil
       WRITE(lunout,*) 'fl_cor_ebil=',        fl_cor_ebil
       WRITE(lunout,*) 'iflag_phytrac=',      iflag_phytrac
       WRITE(lunout,*) 'ok_water_mass_fixer=',ok_water_mass_fixer
       WRITE(lunout,*) 'NVM=',                nvm_lmdz

       !--PC: defining fields to be exchanged between LMDz, ORCHIDEE and NEMO
       WRITE(lunout,*) 'Call to infocfields from physiq'
       CALL infocfields_init

       !AI 08 2023
#ifdef CPP_ECRAD
       ok_3Deffect=.false.
       CALL getin_p('ok_3Deffect',ok_3Deffect)
       namelist_ecrad_file='namelist_ecrad'
#endif

    ENDIF

    IF (prt_level.ge.1) print *,'CONVERGENCE PHYSIQUE THERM 1 '

    !======================================================================
    ! Gestion calendrier : mise a jour du module phys_cal_mod
    !
    !     CALL phys_cal_update(jD_cur,jH_cur)

    !
    ! Si c'est le debut, il faut initialiser plusieurs choses
    !          ********
    !
    IF (debut) THEN
       !rv CRinitialisation de wght_th et lalim_conv pour la
       !definition de la couche alimentation de la convection a partir
       !des caracteristiques du thermique
       wght_th(:,:)=1.
       lalim_conv(:)=1
       !RC
       ustar(:,:)=0.
!       u10m(:,:)=0.
!       v10m(:,:)=0.
       rain_con(:)=0.
       snow_con(:)=0.
       topswai(:)=0.
       topswad(:)=0.
       solswai(:)=0.
       solswad(:)=0.

       wmax_th(:)=0.
       tau_overturning_th(:)=0.

       IF (ANY(type_trac == ['inca','inco'])) THEN
          ! jg : initialisation jusqu'au ces variables sont dans restart
          ccm(:,:,:) = 0.
          tau_aero(:,:,:,:) = 0.
          piz_aero(:,:,:,:) = 0.
          cg_aero(:,:,:,:) = 0.
          d_q_ch4(:,:) = 0.

          config_inca='none' ! default
          CALL getin_p('config_inca',config_inca)

       ELSE
          config_inca='none' ! default
       ENDIF

       tau_aero(:,:,:,:) = 1.e-15
       piz_aero(:,:,:,:) = 1.
       cg_aero(:,:,:,:)  = 0.
       d_q_ch4(:,:) = 0.

       IF (aerosol_couple .AND. (config_inca /= "aero" &
            .AND. config_inca /= "aeNP ")) THEN
          abort_message &
               = 'if aerosol_couple is activated, config_inca need to be ' &
               // 'aero or aeNP'
          CALL abort_physic (modname,abort_message,1)
       ENDIF

       rnebcon0(:,:) = 0.0
       clwcon0(:,:) = 0.0
       rnebcon(:,:) = 0.0
       clwcon(:,:) = 0.0

       !
       print*,'iflag_coupl,iflag_clos,iflag_wake', &
            iflag_coupl,iflag_clos,iflag_wake
       print*,'iflag_cycle_diurne', iflag_cycle_diurne
       !
       IF (iflag_con.EQ.2.AND.iflag_cld_th.GT.-1) THEN
          abort_message = 'Tiedtke needs iflag_cld_th=-2 or -1'
          CALL abort_physic (modname,abort_message,1)
       ENDIF
       !
       !
       ! Initialiser les compteurs:
       !
       itap    = 0
       itaprad = 0
       itapcv = 0
       itapwk = 0

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! Un petit travail \`a faire ici.
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       IF (iflag_pbl>1) THEN
          PRINT*, "Using method MELLOR&YAMADA"
       ENDIF

       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! FH 2008/05/02 changement lie a la lecture de nbapp_rad dans
       ! phylmd plutot que dyn3d
       ! Attention : la version precedente n'etait pas tres propre.
       ! Il se peut qu'il faille prendre une valeur differente de nbapp_rad
       ! pour obtenir le meme resultat.
!jyg for fh<
       WRITE(lunout,*) 'Pas de temps phys_tstep pdtphys ',phys_tstep,pdtphys
       IF (abs(phys_tstep-pdtphys)>1.e-10) THEN
          abort_message='pas de temps doit etre entier en seconde pour orchidee et XIOS'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
!>jyg
       IF (MOD(NINT(86400./phys_tstep),nbapp_rad).EQ.0) THEN
          radpas = NINT( 86400./phys_tstep)/nbapp_rad
       ELSE
          WRITE(lunout,*) 'le nombre de pas de temps physique doit etre un ', &
               'multiple de nbapp_rad'
          WRITE(lunout,*) 'changer nbapp_rad ou alors commenter ce test ', &
               'mais 1+1<>2'
          abort_message='nbre de pas de temps physique n est pas multiple ' &
               // 'de nbapp_rad'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
       IF (nbapp_cv .EQ. 0) nbapp_cv=86400./phys_tstep
       IF (nbapp_wk .EQ. 0) nbapp_wk=86400./phys_tstep
       print *,'physiq, nbapp_cv, nbapp_wk ',nbapp_cv,nbapp_wk
       IF (MOD(NINT(86400./phys_tstep),nbapp_cv).EQ.0) THEN
          cvpas_0 = NINT( 86400./phys_tstep)/nbapp_cv
          cvpas = cvpas_0
       print *,'physiq, cvpas ',cvpas
       ELSE
          WRITE(lunout,*) 'le nombre de pas de temps physique doit etre un ', &
               'multiple de nbapp_cv'
          WRITE(lunout,*) 'changer nbapp_cv ou alors commenter ce test ', &
               'mais 1+1<>2'
          abort_message='nbre de pas de temps physique n est pas multiple ' &
               // 'de nbapp_cv'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
       IF (MOD(NINT(86400./phys_tstep),nbapp_wk).EQ.0) THEN
          wkpas = NINT( 86400./phys_tstep)/nbapp_wk
!       print *,'physiq, wkpas ',wkpas
       ELSE
          WRITE(lunout,*) 'le nombre de pas de temps physique doit etre un ', &
               'multiple de nbapp_wk'
          WRITE(lunout,*) 'changer nbapp_wk ou alors commenter ce test ', &
               'mais 1+1<>2'
          abort_message='nbre de pas de temps physique n est pas multiple ' &
               // 'de nbapp_wk'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL init_iophy_new(latitude_deg,longitude_deg)

          !===================================================================
          !IM stations CFMIP
          nCFMIP=npCFMIP
          OPEN(98,file='npCFMIP_param.data',status='old', &
               form='formatted',iostat=iostat)
          IF (iostat == 0) THEN
             READ(98,*,end=998) nCFMIP
998          CONTINUE
             CLOSE(98)
             CONTINUE
             IF(nCFMIP.GT.npCFMIP) THEN
                print*,'nCFMIP > npCFMIP : augmenter npCFMIP et recompiler'
                CALL abort_physic("physiq", "", 1)
             ELSE
                print*,'physiq npCFMIP=',npCFMIP,'nCFMIP=',nCFMIP
             ENDIF

             !
             ALLOCATE(tabCFMIP(nCFMIP))
             ALLOCATE(lonCFMIP(nCFMIP), latCFMIP(nCFMIP))
             ALLOCATE(tabijGCM(nCFMIP))
             ALLOCATE(lonGCM(nCFMIP), latGCM(nCFMIP))
             ALLOCATE(iGCM(nCFMIP), jGCM(nCFMIP))
             !
             ! lecture des nCFMIP stations CFMIP, de leur numero
             ! et des coordonnees geographiques lonCFMIP, latCFMIP
             !
             CALL read_CFMIP_point_locations(nCFMIP, tabCFMIP,  &
                  lonCFMIP, latCFMIP)
             !
             ! identification des
             ! 1) coordonnees lonGCM, latGCM des points CFMIP dans la
             ! grille de LMDZ
             ! 2) indices points tabijGCM de la grille physique 1d sur
             ! klon points
             ! 3) indices iGCM, jGCM de la grille physique 2d
             !
             CALL LMDZ_CFMIP_point_locations(nCFMIP, lonCFMIP, latCFMIP, &
                  tabijGCM, lonGCM, latGCM, iGCM, jGCM)
             !
          ELSE
             ALLOCATE(tabijGCM(0))
             ALLOCATE(lonGCM(0), latGCM(0))
             ALLOCATE(iGCM(0), jGCM(0))
          ENDIF

       !$OMP MASTER
       ! FH : if ok_sync=.true. , the time axis is written at each time step
       ! in the output files. Only at the end in the opposite case
       ok_sync_omp=.FALSE.
       CALL getin('ok_sync',ok_sync_omp)
       CALL phys_output_open(longitude_deg,latitude_deg,nCFMIP,tabijGCM, &
            iGCM,jGCM,lonGCM,latGCM, &
            jjmp1,nlevSTD,clevSTD,rlevSTD, phys_tstep,ok_veget, &
            type_ocean,iflag_pbl,iflag_pbl_split,ok_mensuel,ok_journe, &
            ok_hf,ok_instan,ok_LES,ok_ade,ok_aie, &
            read_climoz, phys_out_filestations, &
            aerosol_couple, &
            flag_aerosol_strat, pdtphys, paprs, pphis,  &
            pplay, lmax_th, ptconv, ptconvth, ivap,  &
            d_u, d_t, qx, d_qx, zmasse, ok_sync_omp)
       !$OMP END MASTER
       !$OMP BARRIER
       ok_sync=ok_sync_omp

!       freq_outNMC(1) = ecrit_files(7)
!       freq_outNMC(2) = ecrit_files(8)
!       freq_outNMC(3) = ecrit_files(9)
!       WRITE(lunout,*)'OK freq_outNMC(1)=',freq_outNMC(1)
!       WRITE(lunout,*)'OK freq_outNMC(2)=',freq_outNMC(2)
!       WRITE(lunout,*)'OK freq_outNMC(3)=',freq_outNMC(3)

#ifndef CPP_XIOS
       CALL ini_paramLMDZ_phy(phys_tstep,nid_ctesGCM)
#endif

       ecrit_reg = ecrit_reg * un_jour
       ecrit_tra = ecrit_tra * un_jour

       !XXXPB Positionner date0 pour initialisation de ORCHIDEE
       date0 = jD_ref
       WRITE(*,*) 'physiq date0 : ',date0
       !

!       CALL create_climoz(read_climoz)
      IF (.NOT. create_etat0_limit) CALL init_aero_fromfile(flag_aerosol, aerosol_couple)  !! initialise aero from file for XIOS interpolation (unstructured_grid)
      IF (.NOT. create_etat0_limit) CALL init_readaerosolstrato(flag_aerosol_strat)  !! initialise aero strato from file for XIOS interpolation (unstructured_grid)

      ! A.I : Initialisations pour le 1er passage a Cosp
      if (ok_cosp) then

IF (CPPKEY_COSP) THEN
        CALL ini_COSP(ref_liq_cosp0,ref_ice_cosp0,pctsrf_cosp0,zu10m_cosp0,zv10m_cosp0, &
               zxtsol_cosp0,zx_rh_cosp0,cldfra_cosp0,rnebcon_cosp0,flwc_cosp0, &
               fiwc_cosp0,prfl_cosp0,psfl_cosp0,pmflxr_cosp0,pmflxs_cosp0, &
               mr_ozone_cosp0,cldtau_cosp0,cldemi_cosp0,JrNt_cosp0)

        CALL phys_cosp(itap,phys_tstep,freq_cosp, &
               ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
               ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
               klon,klev,longitude_deg,latitude_deg,presnivs,overlap, &
               JrNt_cosp0,ref_liq_cosp0,ref_ice_cosp0, &
               pctsrf_cosp0, &
               zu10m_cosp0,zv10m_cosp0,pphis, &
               pphi,paprs(:,1:klev),pplay,zxtsol_cosp0,t, &
               qx(:,:,ivap),zx_rh_cosp0,cldfra_cosp0,rnebcon_cosp0,flwc_cosp0,fiwc_cosp0, &
               prfl_cosp0(:,1:klev),psfl_cosp0(:,1:klev), &
               pmflxr_cosp0(:,1:klev),pmflxs_cosp0(:,1:klev), &
               mr_ozone_cosp0,cldtau_cosp0, cldemi_cosp0)
END IF

IF (CPPKEY_COSPV2) THEN
          CALL ini_COSP(ref_liq_cosp0,ref_ice_cosp0,pctsrf_cosp0,zu10m_cosp0,zv10m_cosp0, &
               zxtsol_cosp0,zx_rh_cosp0,cldfra_cosp0,rnebcon_cosp0,flwc_cosp0, &
               fiwc_cosp0,prfl_cosp0,psfl_cosp0,pmflxr_cosp0,pmflxs_cosp0, &
               mr_ozone_cosp0,cldtau_cosp0,cldemi_cosp0,JrNt_cosp0)

          CALL lmdz_cosp_interface(itap,phys_tstep,freq_cosp, &
               ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
               ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
               klon,klev,longitude_deg,latitude_deg,presnivs,overlap, &
               JrNt_cosp0,ref_liq_cosp0,ref_ice_cosp0, &
               pctsrf_cosp0, &
               zu10m_cosp0,zv10m_cosp0,pphis, &
               pphi,paprs(:,1:klev),pplay,zxtsol_cosp0,t, &
               qx(:,:,ivap),zx_rh_cosp0,cldfra_cosp0,rnebcon_cosp0,flwc_cosp0,fiwc_cosp0, &
               prfl_cosp0(:,1:klev),psfl_cosp0(:,1:klev), &
               pmflxr_cosp0(:,1:klev),pmflxs_cosp0(:,1:klev), &
               mr_ozone_cosp0,cldtau_cosp0, cldemi_cosp0)
END IF
      ENDIF

       !
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Nouvelle initialisation pour le rayonnement RRTM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       CALL iniradia(klon,klev,paprs(1,1:klev+1))

       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL surf_wind_ini(klon,lunout)
       CALL getin_p('nsurfwind',nsurfwind)
       allocate(surf_wind_value(klon,nsurfwind),surf_wind_proba(klon,nsurfwind))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! initialisation des variables tkeoro
       addtkeoro=0
       CALL getin_p('addtkeoro',addtkeoro)

       alphatkeoro=1.
       CALL getin_p('alphatkeoro',alphatkeoro)
       alphatkeoro=min(max(0.,alphatkeoro),1.)

       smallscales_tkeoro=.FALSE.
       CALL getin_p('smallscales_tkeoro',smallscales_tkeoro)
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL wake_ini(iflag_wake,rg,rd,rv,prt_level)
       CALL yamada_ini(klon,lunout,prt_level)
       viscom=1.46E-5
       viscoh=2.06E-5
       CALL atke_ini(RG, RD, RPI, RCPD, RV, viscom, viscoh)
       CALL thermcell_ini(iflag_thermals,prt_level,tau_thermals,lunout, &
   &    RG,RD,RCPD,RKAPPA,RLVTT,RETV)
       CALL ratqs_ini(klon,klev,iflag_thermals,lunout,nbsrf,is_lic,is_ter,RG,RV,RD,RCPD,RLSTT,RLVTT,RTT)
       CALL lscp_ini(pdtphys,lunout,prt_level,ok_ice_supersat,iflag_ratqs,fl_cor_ebil, &
                     RCPD,RLSTT,RLVTT,RLMLT,RVTMP2,RTT,RD,RV,RG,RPI,EPS_W)
       CALL blowing_snow_ini(RCPD, RLSTT, RLVTT, RLMLT, &
                             RVTMP2, RTT,RD,RG, RV, RPI)
       ! Test de coherence sur oc_cdnc utilisé uniquement par cloud_optics_prop
       IF (ok_newmicro) then
          IF (iflag_rrtm.EQ.1) THEN
#ifdef CPP_RRTM
             IF (ok_cdnc.AND.NRADLP.NE.3) THEN
             abort_message='RRTM choix incoherent NRADLP doit etre egal a 3 ' &
                  // 'pour ok_cdnc'
             CALL abort_physic(modname,abort_message,1)
             ENDIF
#else

             abort_message='You should compile with -rrtm if running with '//'iflag_rrtm=1'
             CALL abort_physic(modname,abort_message,1)
#endif
          ENDIF
       ENDIF
       CALL cloud_optics_prop_ini(klon, klev, prt_level, lunout, flag_aerosol, &
                                  & ok_cdnc, bl95_b0, &
                                  & bl95_b1, latitude_deg, rpi, rg, rd, &
                                  & zepsec, novlp, iflag_ice_thermo, ok_new_lscp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Initialisation des champs dans phytrac* qui sont utilises par phys_output_write*
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (CPPKEY_REPROBUS) THEN
       CALL strataer_init
       CALL strataer_emiss_init
END IF

IF (CPPKEY_STRATAER) THEN
       CALL strataer_init
       CALL strataer_nuc_init
       CALL strataer_emiss_init
END IF

IF (CPPKEY_DUST) THEN
       ! Quand on utilise SPLA, on force iflag_phytrac=1
       CALL spla_ini(is_oce,RNAVO,RG,RD,RCPD,RLVTT,RLSTT,RETV,RTT,              &
               R2ES,R3LES,R3IES,R4LES,R4IES,R5LES,R5IES,RVTMP2)
       CALL phytracr_spl_out_init()
       CALL phys_output_write_spl(itap, pdtphys, paprs, pphis,                  &
                                pplay, lmax_th, aerosol_couple,                 &
                                ok_ade, ok_aie, ivap, ok_sync,                  &
                                ptconv, read_climoz, clevSTD,                   &
                                ptconvth, d_t, qx, d_qx, d_tr_dyn, zmasse,      &
                                flag_aerosol, flag_aerosol_strat, ok_cdnc)
ELSE
       ! phys_output_write écrit des variables traceurs seulement si iflag_phytrac == 1
       ! donc seulement dans ce cas on doit appeler phytrac_init()
       IF (iflag_phytrac == 1 ) THEN
          CALL phytrac_init()
       ENDIF
       CALL phys_output_write(itap, pdtphys, paprs, pphis,                    &
                              pplay, lmax_th, aerosol_couple,                 &
                              ok_ade, ok_aie, ok_volcan, ivap, iliq, isol, ibs,  ok_sync,&
                              ptconv, read_climoz, clevSTD,                   &
                              ptconvth, d_u, d_t, qx, d_qx, zmasse,           &
                              flag_aerosol, flag_aerosol_strat, ok_cdnc, t, u1, v1)
END IF


       IF (using_xios) THEN
         IF (is_omp_master) CALL xios_update_calendar(1)
       ENDIF

       IF(read_climoz>=1 .AND. create_etat0_limit) CALL regr_horiz_time_climoz(read_climoz,ok_daily_climoz)
       CALL create_etat0_limit_unstruct
       CALL phyetat0 ("startphy.nc",clesphy0,tabcntr0)

!jyg<
       IF (iflag_pbl<=1) THEN
          ! No TKE for Standard Physics
          pbl_tke(:,:,:)=0.

       ELSE IF (klon_glo==1) THEN
          pbl_tke(:,:,is_ave) = 0.
          pbl_eps(:,:,is_ave) = 0.
          DO nsrf=1,nbsrf
            DO k = 1,klev+1
                 pbl_tke(:,k,is_ave) = pbl_tke(:,k,is_ave) &
                     +pctsrf(:,nsrf)*pbl_tke(:,k,nsrf)
                 pbl_eps(:,k,is_ave) = pbl_eps(:,k,is_ave) &
                     +pctsrf(:,nsrf)*pbl_eps(:,k,nsrf)
            ENDDO
          ENDDO
       ELSE
          pbl_tke(:,:,is_ave) = 0. !ym missing init : maybe must be initialized in the same way that for klon_glo==1 ??
!>jyg
          pbl_eps(:,:,is_ave) = 0.
       ENDIF
       !IM begin
       print*,'physiq: clwcon rnebcon ratqs',clwcon(1,1),rnebcon(1,1) &
            ,ratqs(1,1)
       !IM end


       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       ! on remet le calendrier a zero
       !
       IF (raz_date .eq. 1) THEN
          itau_phy = 0
       ENDIF

!       IF (ABS(phys_tstep-pdtphys).GT.0.001) THEN
!          WRITE(lunout,*) 'Pas physique n est pas correct',phys_tstep, &
!               pdtphys
!          abort_message='Pas physique n est pas correct '
!          !           call abort_physic(modname,abort_message,1)
!          phys_tstep=pdtphys
!       ENDIF
       IF (nlon .NE. klon) THEN
          WRITE(lunout,*)'nlon et klon ne sont pas coherents', nlon,  &
               klon
          abort_message='nlon et klon ne sont pas coherents'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
       IF (nlev .NE. klev) THEN
          WRITE(lunout,*)'nlev et klev ne sont pas coherents', nlev, &
               klev
          abort_message='nlev et klev ne sont pas coherents'
          CALL abort_physic(modname,abort_message,1)
       ENDIF
       !
       IF (phys_tstep*REAL(radpas).GT.21600..AND.iflag_cycle_diurne.GE.1) THEN
          WRITE(lunout,*)'Nbre d appels au rayonnement insuffisant'
          WRITE(lunout,*)"Au minimum 4 appels par jour si cycle diurne"
          abort_message='Nbre d appels au rayonnement insuffisant'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Initialisation pour la convection de K.E. et pour les poches froides
       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       WRITE(lunout,*)"Clef pour la convection, iflag_con=", iflag_con
       WRITE(lunout,*)"Clef pour le driver de la convection, ok_cvl=", ok_cvl
       !
       !KE43
       ! Initialisation pour la convection de K.E. (sb):
       IF (iflag_con.GE.3) THEN

          WRITE(lunout,*)"*** Convection de Kerry Emanuel 4.3  "
          WRITE(lunout,*) &
               "On va utiliser le melange convectif des traceurs qui"
          WRITE(lunout,*)"est calcule dans convect4.3"
          WRITE(lunout,*)" !!! penser aux logical flags de phytrac"

          DO i = 1, klon
             ema_cbmf(i) = 0.
             ema_pcb(i)  = 0.
             ema_pct(i)  = 0.
             !          ema_workcbmf(i) = 0.
          ENDDO
          !IM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>BEG
          DO i = 1, klon
             ibas_con(i) = 1
             itop_con(i) = 1
          ENDDO
          !IM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>END
          !================================================================
          !CR:04.12.07: initialisations poches froides
          ! Controle de ALE et ALP pour la fermeture convective (jyg)
          IF (iflag_wake>=1) THEN
             CALL ini_wake(0.,0.,it_wape_prescr,wape_prescr,fip_prescr &
                  ,alp_bl_prescr, ale_bl_prescr)
             ! 11/09/06 rajout initialisation ALE et ALP du wake et PBL(YU)
             !        print*,'apres ini_wake iflag_cld_th=', iflag_cld_th
             !
             ! Initialize tendencies of wake state variables (for some flag values
             ! they are not computed).
             d_deltat_wk(:,:) = 0.
             d_deltaq_wk(:,:) = 0.
             d_deltat_wk_gw(:,:) = 0.
             d_deltaq_wk_gw(:,:) = 0.
             d_deltat_vdf(:,:) = 0.
             d_deltaq_vdf(:,:) = 0.
             d_deltat_the(:,:) = 0.
             d_deltaq_the(:,:) = 0.
             d_deltat_ajs_cv(:,:) = 0.
             d_deltaq_ajs_cv(:,:) = 0.
             d_s_wk(:) = 0.
             d_s_a_wk(:) = 0.
             d_dens_wk(:) = 0.
             d_dens_a_wk(:) = 0.
          ENDIF  !  (iflag_wake>=1)

          !        do i = 1,klon
          !           Ale_bl(i)=0.
          !           Alp_bl(i)=0.
          !        enddo

       !ELSE
       !   ALLOCATE(tabijGCM(0))
       !   ALLOCATE(lonGCM(0), latGCM(0))
       !   ALLOCATE(iGCM(0), jGCM(0))
       ENDIF  !  (iflag_con.GE.3)
       !
       DO i=1,klon
          rugoro(i) = f_rugoro * MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
       ENDDO

       !34EK
       IF (ok_orodr) THEN

          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! FH sans doute a enlever de finitivement ou, si on le
          ! garde, l'activer justement quand ok_orodr = false.
          ! ce rugoro est utilise par la couche limite et fait double emploi
          ! avec les param\'etrisations sp\'ecifiques de Francois Lott.
          !           DO i=1,klon
          !             rugoro(i) = MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
          !           ENDDO
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (ok_strato) THEN
             CALL SUGWD_strato(klon,klev,paprs,pplay)
          ELSE
             CALL SUGWD(klon,klev,paprs,pplay)
          ENDIF

          DO i=1,klon
             zuthe(i)=0.
             zvthe(i)=0.
             IF (zstd(i).gt.10.) THEN
                zuthe(i)=(1.-zgam(i))*cos(zthe(i))
                zvthe(i)=(1.-zgam(i))*sin(zthe(i))
             ENDIF
          ENDDO
       ENDIF
       !
       !
       lmt_pas = NINT(86400./phys_tstep * 1.0)   ! tous les jours
       WRITE(lunout,*)'La frequence de lecture surface est de ',  &
            lmt_pas
       !
       capemaxcels = 't_max(X)'
       t2mincels = 't_min(X)'
       t2maxcels = 't_max(X)'
       tinst = 'inst(X)'
       tave = 'ave(X)'
       !IM cf. AM 081204 BEG
       write(lunout,*)'AVANT HIST IFLAG_CON=',iflag_con
       !IM cf. AM 081204 END
       !
       !=============================================================
       !   Initialisation des sorties
       !=============================================================

       IF (using_xios) THEN
         ! Get "missing_val" value from XML files (from temperature variable)
         IF (is_omp_master) CALL xios_get_field_attr("temp",default_value=missing_val)
         CALL bcast_omp(missing_val)
       ENDIF

       IF (using_xios) THEN
         ! Need to put this initialisation after phyetat0 as in the coupled model the XIOS context is only
         ! initialised at that moment
         ! Get "missing_val" value from XML files (from temperature variable)
         IF (is_omp_master) CALL xios_get_field_attr("temp",default_value=missing_val)
         CALL bcast_omp(missing_val)
       !
       ! Now we activate some double radiation call flags only if some
       ! diagnostics are requested, otherwise there is no point in doing this
         IF (is_master) THEN
           !--setting up swaero_diag to TRUE in XIOS case
           IF (xios_field_is_active("topswad").OR.xios_field_is_active("topswad0").OR. &
              xios_field_is_active("solswad").OR.xios_field_is_active("solswad0").OR. &
              xios_field_is_active("topswai").OR.xios_field_is_active("solswai").OR.  &
                (iflag_rrtm==1.AND.(xios_field_is_active("toplwad").OR.xios_field_is_active("toplwad0").OR. &
                                    xios_field_is_active("sollwad").OR.xios_field_is_active("sollwad0"))))  &
              !!!--for now these fields are not in the XML files so they are omitted
              !!!  xios_field_is_active("toplwai").OR.xios_field_is_active("sollwai") !))) &
              swaero_diag=.TRUE.

           !--setting up swaerofree_diag to TRUE in XIOS case
           IF (xios_field_is_active("SWdnSFCcleanclr").OR.xios_field_is_active("SWupSFCcleanclr").OR. &
              xios_field_is_active("SWupTOAcleanclr").OR.xios_field_is_active("rsucsaf").OR.   &
              xios_field_is_active("rsdcsaf") .OR. xios_field_is_active("LWdnSFCcleanclr").OR. &
              xios_field_is_active("LWupTOAcleanclr")) &
              swaerofree_diag=.TRUE.

           !--setting up dryaod_diag to TRUE in XIOS case
           DO naero = 1, naero_tot-1
             IF (xios_field_is_active("dryod550_"//name_aero_tau(naero))) dryaod_diag=.TRUE.
           ENDDO
           !
          !--setting up ok_4xCO2atm to TRUE in XIOS case
           IF (xios_field_is_active("rsut4co2").OR.xios_field_is_active("rlut4co2").OR. &
              xios_field_is_active("rsutcs4co2").OR.xios_field_is_active("rlutcs4co2").OR. &
              xios_field_is_active("rsu4co2").OR.xios_field_is_active("rsucs4co2").OR. &
              xios_field_is_active("rsd4co2").OR.xios_field_is_active("rsdcs4co2").OR. &
              xios_field_is_active("rlu4co2").OR.xios_field_is_active("rlucs4co2").OR. &
              xios_field_is_active("rld4co2").OR.xios_field_is_active("rldcs4co2")) &
              ok_4xCO2atm=.TRUE.
           ENDIF
           !$OMP BARRIER
           CALL bcast(swaero_diag)
           CALL bcast(swaerofree_diag)
           CALL bcast(dryaod_diag)
           CALL bcast(ok_4xCO2atm)
         ENDIF !using_xios
       !
       CALL printflag( tabcntr0,radpas,ok_journe, &
            ok_instan, ok_region )
       !
       !
       ! Prescrire l'ozone dans l'atmosphere
       !
       !c         DO i = 1, klon
       !c         DO k = 1, klev
       !c            CALL o3cm (paprs(i,k)/100.,paprs(i,k+1)/100., wo(i,k),20)
       !c         ENDDO
       !c         ENDDO
       !
       IF (ANY(type_trac == ['inca','inco'])) THEN ! ModThL
IF (CPPKEY_INCA) THEN
          CALL VTe(VTphysiq)
          CALL VTb(VTinca)
          calday = REAL(days_elapsed) + jH_cur
          WRITE(lunout,*) 'initial time chemini', days_elapsed, calday

          call init_const_lmdz( &
          ndays, nbsrf, is_oce,is_sic, is_ter,is_lic, calend, &
          config_inca)

          CALL init_inca_geometry( &
               longitude, latitude, &
               boundslon, boundslat, &
               cell_area, ind_cell_glo)

          if (grid_type==unstructured) THEN
             CALL chemini(  pplay, &
                  nbp_lon, nbp_lat, &
                  latitude_deg, &
                  longitude_deg, &
                  presnivs, &
                  calday, &
                  klon, &
                  nqtot, &
                  nqo+nqCO2, &
                  pdtphys, &
                  annee_ref, &
                  year_cur, &
                  day_ref,  &
                  day_ini, &
                  start_time, &
                  itau_phy, &
                  date0, &
                  chemistry_couple, &
                  init_source, &
                  init_tauinca, &
                  init_pizinca, &
                  init_cginca, &
                  init_ccminca)
          ELSE
             CALL chemini(  pplay, &
                  nbp_lon, nbp_lat, &
                  latitude_deg, &
                  longitude_deg, &
                  presnivs, &
                  calday, &
                  klon, &
                  nqtot, &
                  nqo+nqCO2, &
                  pdtphys, &
                  annee_ref, &
                  year_cur, &
                  day_ref,  &
                  day_ini, &
                  start_time, &
                  itau_phy, &
                  date0, &
                  chemistry_couple, &
                  init_source, &
                  init_tauinca, &
                  init_pizinca, &
                  init_cginca, &
                  init_ccminca, &
                  io_lon, &
                  io_lat)
          ENDIF


          ! initialisation des variables depuis le restart de inca
          ccm(:,:,:) = init_ccminca
          tau_aero(:,:,:,:) = init_tauinca
          piz_aero(:,:,:,:) = init_pizinca
          cg_aero(:,:,:,:) = init_cginca
!


          CALL VTe(VTinca)
          CALL VTb(VTphysiq)
END IF
       ENDIF
       !
       IF (type_trac == 'repr') THEN
IF (CPPKEY_REPROBUS) THEN
          CALL chemini_rep(  &
               presnivs, &
               pdtphys, &
               annee_ref, &
               day_ref,  &
               day_ini, &
               start_time, &
               itau_phy, &
               io_lon, &
               io_lat)
END IF
       ENDIF

       !$omp single
       IF (read_climoz >= 1) CALL open_climoz(ncid_climoz, press_cen_climoz,   &
           press_edg_climoz, time_climoz, ok_daily_climoz, adjust_tropopause)
       !$omp end single
       !
       !IM betaCRF
       pfree=70000. !Pa
       beta_pbl=1.
       beta_free=1.
       lon1_beta=-180.
       lon2_beta=+180.
       lat1_beta=90.
       lat2_beta=-90.
       mskocean_beta=.FALSE.

       !albedo SB >>>
       SELECT CASE(nsw)
       CASE(2)
          SFRWL(1)=0.45538747
          SFRWL(2)=0.54461211
       CASE(4)
          SFRWL(1)=0.45538747
          SFRWL(2)=0.32870591
          SFRWL(3)=0.18568763
          SFRWL(4)=3.02191470E-02
       CASE(6)
          SFRWL(1)=1.28432794E-03
          SFRWL(2)=0.12304168
          SFRWL(3)=0.33106142
          SFRWL(4)=0.32870591
          SFRWL(5)=0.18568763
          SFRWL(6)=3.02191470E-02
       END SELECT
       !albedo SB <<<

       OPEN(99,file='beta_crf.data',status='old', &
            form='formatted',err=9999)
       READ(99,*,end=9998) pfree
       READ(99,*,end=9998) beta_pbl
       READ(99,*,end=9998) beta_free
       READ(99,*,end=9998) lon1_beta
       READ(99,*,end=9998) lon2_beta
       READ(99,*,end=9998) lat1_beta
       READ(99,*,end=9998) lat2_beta
       READ(99,*,end=9998) mskocean_beta
9998   Continue
       CLOSE(99)
9999   Continue
       WRITE(*,*)'pfree=',pfree
       WRITE(*,*)'beta_pbl=',beta_pbl
       WRITE(*,*)'beta_free=',beta_free
       WRITE(*,*)'lon1_beta=',lon1_beta
       WRITE(*,*)'lon2_beta=',lon2_beta
       WRITE(*,*)'lat1_beta=',lat1_beta
       WRITE(*,*)'lat2_beta=',lat2_beta
       WRITE(*,*)'mskocean_beta=',mskocean_beta

      !lwoff=y : offset LW CRE for radiation code and other schemes
      !lwoff=y : betalwoff=1.
      betalwoff=0.
      IF (ok_lwoff) THEN
         betalwoff=1.
      ENDIF
      WRITE(*,*)'ok_lwoff=',ok_lwoff
      !
      !lwoff=y to begin only sollw and sollwdown are set up to CS values
      sollw = sollw + betalwoff * (sollw0 - sollw)
      sollwdown(:)= sollwdown(:) + betalwoff *(-1.*ZFLDN0(:,1) - &
                    sollwdown(:))

      !--Init for LSCP - condensation
      ratio_ql_qtot(:,:) = 0.
      ratio_qi_qtot(:,:) = 0.


    ENDIF
    !
    !   ****************     Fin  de   IF ( debut  )   ***************
    !
    !
    ! Incrementer le compteur de la physique
    !
    itap   = itap + 1
    IF (is_master .OR. prt_level > 9) THEN
      IF (prt_level > 5 .or. MOD(itap,5) == 0) THEN
         WRITE(LUNOUT,*)'Entering physics elapsed seconds since start ', current_time
         WRITE(LUNOUT,100)year_cur,mth_cur,day_cur,hour/3600.
 100     FORMAT('Date = ',i4.4,' / ',i2.2, ' / ',i2.2,' : ',f20.17)
      ENDIF
    ENDIF
    !
    !
    ! Update fraction of the sub-surfaces (pctsrf) and
    ! initialize, where a new fraction has appeared, all variables depending
    ! on the surface fraction.
    !
    CALL change_srf_frac(itap, phys_tstep, days_elapsed+1,  &
         pctsrf, fevap, z0m, z0h, agesno,              &
         falb_dir, falb_dif, ftsol, ustar, u10m, v10m, pbl_tke)

    ! Update time and other variables in Reprobus
    IF (type_trac == 'repr') THEN
IF (CPPKEY_REPROBUS) THEN
       CALL Init_chem_rep_xjour(jD_cur-jD_ref+day_ref)
       print*,'xjour equivalent rjourvrai',jD_cur-jD_ref+day_ref
       CALL Rtime(debut)
END IF
    ENDIF

    ! Tendances bidons pour les processus qui n'affectent pas certaines
    ! variables.
    du0(:,:)=0.
    dv0(:,:)=0.
    dt0 = 0.
    dq0(:,:)=0.
    dql0(:,:)=0.
    dqi0(:,:)=0.
    dqbs0(:,:)=0.
    dsig0(:) = 0.
    ddens0(:) = 0.
    wkoccur1(:)=1
    !
    ! Mettre a zero des variables de sortie (pour securite)
    !
    DO i = 1, klon
       d_ps(i) = 0.0
    ENDDO
    DO k = 1, klev
       DO i = 1, klon
          d_t(i,k) = 0.0
          d_u(i,k) = 0.0
          d_v(i,k) = 0.0
       ENDDO
    ENDDO
    DO iq = 1, nqtot
       DO k = 1, klev
          DO i = 1, klon
             d_qx(i,k,iq) = 0.0
          ENDDO
       ENDDO
    ENDDO
    beta_prec_fisrt(:,:)=0.
    beta_prec(:,:)=0.
    !
    !   Output variables from the convective scheme should not be set to 0
    !   since convection is not always called at every time step.
    IF (ok_bug_cv_trac) THEN
      da(:,:)=0.
      mp(:,:)=0.
      phi(:,:,:)=0.
      ! RomP >>>
      phi2(:,:,:)=0.
      epmlmMm(:,:,:)=0.
      eplaMm(:,:)=0.
      d1a(:,:)=0.
      dam(:,:)=0.
      pmflxr(:,:)=0.
      pmflxs(:,:)=0.
      ! RomP <<<
    ENDIF
    !
    ! Ne pas affecter les valeurs entrees de u, v, h, et q
    !
    DO k = 1, klev
       DO i = 1, klon
          t_seri(i,k)  = t(i,k)
          u_seri(i,k)  = u(i,k)
          v_seri(i,k)  = v(i,k)
          q_seri(i,k)  = qx(i,k,ivap)
          ql_seri(i,k) = qx(i,k,iliq)
          qbs_seri(i,k)= 0.
          cf_seri(i,k) = 0.
          rvc_seri(i,k)= 0.
          !CR: ATTENTION, on rajoute la variable glace
          IF (nqo.EQ.2) THEN             !--vapour and liquid only
             qs_seri(i,k) = 0.
          ELSE IF (nqo.EQ.3) THEN        !--vapour, liquid and ice
             qs_seri(i,k) = qx(i,k,isol)
          ELSE IF (nqo.GE.4) THEN        !--vapour, liquid, ice, blowing snow, cloud fraction and cloudy water vapor to total water vapor ratio
             qs_seri(i,k) = qx(i,k,isol)
             IF (ok_ice_supersat) THEN
               cf_seri(i,k) = qx(i,k,icf)
               rvc_seri(i,k) = qx(i,k,irvc)
             ENDIF
             IF (ok_bs) THEN
               qbs_seri(i,k)= qx(i,k,ibs)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    ! +++ Python coupling: Send inputs to NN  +++
    IF (lk_pycpl) THEN
       nstep = NINT( daysec / phys_tstep )
       nyear = year_len * nstep
       istep = MOD( days_elapsed*nstep + INT(jh_cur*nstep), nyear ) + 1
       nn_cosday = COS( 2.*pi*istep/nyear ) * SQRT(2.)
       nn_sinday = SIN( 2.*pi*istep/nyear ) * SQRT(2.)
       !
       ! Send
       CALL send_to_python_phys('u', u_seri, itap)
       CALL send_to_python_phys('v', v_seri, itap)
       CALL send_to_python_phys('tpot', theta, itap)
       CALL send_to_python_phys('pres', pplay, itap)
       CALL send_to_python_phys('psol', paprs, itap)
       CALL send_to_python_phys('topo', pphis / 9.81, itap)
       CALL send_to_python_cpl('cos_day', nn_cosday, itap)
       CALL send_to_python_cpl('sin_day', nn_sinday, itap)
    END IF
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++

    ! in case of TKE advection, we interpolate vertically
    ! since TKE is defined at the bottom interface of layers but
    ! it is interpolated onto middle layers for advection

    IF (ok_advtke) THEN
      DO k=2,klev
        DO i=1,klon
           pbl_tke(i,k,:)=(qx(i,k-1,itke)*zmasse(i,k-1)+qx(i,k,itke)*zmasse(i,k))/(zmasse(i,k-1)+zmasse(i,k))
        ENDDO
      ENDDO
    ENDIF

    tke0(:,:)=pbl_tke(:,:,is_ave)


    !--OB water mass fixer
    IF (ok_water_mass_fixer) THEN
    !--store initial water burden
    qql1(:)=0.0
    DO k = 1, klev
      qql1(:)=qql1(:)+(q_seri(:,k)+ql_seri(:,k))*zmasse(:,k)
      IF (nqo >= 3) THEN
        qql1(:)=qql1(:)+qs_seri(:,k)*zmasse(:,k)
      ENDIF
      IF (ok_bs) THEN
        qql1(:)=qql1(:)+qbs_seri(:,k)*zmasse(:,k)
      ENDIF
    ENDDO
    ENDIF
    !--fin mass fixer

    IF (nqtot > (nqo+nqtke)) THEN
       ! water isotopes are not included in tr_seri
       itr = 0
       DO iq = 1, nqtot
         IF(.NOT.tracers(iq)%isInPhysics) CYCLE
         itr = itr+1
          DO  k = 1, klev
             DO  i = 1, klon
                tr_seri(i,k,itr) = qx(i,k,iq)
             ENDDO
          ENDDO
       ENDDO
    ELSE
! DC: make sure the final "1" index was meant for 1st H2O phase (vapor) !!!
       tr_seri(:,:,strIdx(tracers(:)%name,addPhase('H2O','g'))) = 0.0
    ENDIF
!
! Temporary solutions adressing ticket #104 and the non initialisation of tr_ancien
! LF
    IF (debut) THEN
      WRITE(lunout,*)' WARNING: tr_ancien initialised to tr_seri'
       itr = 0
       do iq = 1, nqtot
         IF(.NOT.tracers(iq)%isInPhysics) CYCLE
         itr = itr+1
         tr_ancien(:,:,itr)=tr_seri(:,:,itr)
       enddo
    ENDIF
    !
    DO i = 1, klon
       ztsol(i) = 0.
    ENDDO
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          ztsol(i) = ztsol(i) + ftsol(i,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO
    ! Initialize variables used for diagnostic purpose
    IF (flag_inhib_tend .ne. 0) CALL init_cmp_seri

    ! Diagnostiquer la tendance dynamique
    !
    IF (ancien_ok) THEN
    !
       d_u_dyn(:,:)  = (u_seri(:,:)-u_ancien(:,:))/phys_tstep
       d_v_dyn(:,:)  = (v_seri(:,:)-v_ancien(:,:))/phys_tstep
       d_t_dyn(:,:)  = (t_seri(:,:)-t_ancien(:,:))/phys_tstep
       d_q_dyn(:,:)  = (q_seri(:,:)-q_ancien(:,:))/phys_tstep
       d_ql_dyn(:,:) = (ql_seri(:,:)-ql_ancien(:,:))/phys_tstep
       d_qs_dyn(:,:) = (qs_seri(:,:)-qs_ancien(:,:))/phys_tstep
       d_qbs_dyn(:,:)= (qbs_seri(:,:)-qbs_ancien(:,:))/phys_tstep
       d_cf_dyn(:,:) = (cf_seri(:,:)-cf_ancien(:,:))/phys_tstep
       d_rvc_dyn(:,:)= (rvc_seri(:,:)-rvc_ancien(:,:))/phys_tstep
       d_tke_dyn(:,:)= (pbl_tke(:,:,is_ave)-tke_ancien(:,:))/phys_tstep
       CALL water_int(klon,klev,q_seri,zmasse,zx_tmp_fi2d)
       d_q_dyn2d(:)=(zx_tmp_fi2d(:)-prw_ancien(:))/phys_tstep
       CALL water_int(klon,klev,ql_seri,zmasse,zx_tmp_fi2d)
       d_ql_dyn2d(:)=(zx_tmp_fi2d(:)-prlw_ancien(:))/phys_tstep
       CALL water_int(klon,klev,qs_seri,zmasse,zx_tmp_fi2d)
       d_qs_dyn2d(:)=(zx_tmp_fi2d(:)-prsw_ancien(:))/phys_tstep
       CALL water_int(klon,klev,qbs_seri,zmasse,zx_tmp_fi2d)
       d_qbs_dyn2d(:)=(zx_tmp_fi2d(:)-prbsw_ancien(:))/phys_tstep
       ! !! RomP >>>   td dyn traceur
       IF (nqtot > (nqo+nqtke)) d_tr_dyn(:,:,:)=(tr_seri(:,:,:)-tr_ancien(:,:,:))/phys_tstep
       ! !! RomP <<<
    ELSE
       d_u_dyn(:,:)  = 0.0
       d_v_dyn(:,:)  = 0.0
       d_t_dyn(:,:)  = 0.0
       d_q_dyn(:,:)  = 0.0
       d_ql_dyn(:,:) = 0.0
       d_qs_dyn(:,:) = 0.0
       d_qbs_dyn(:,:)= 0.0
       d_cf_dyn(:,:) = 0.0
       d_rvc_dyn(:,:)= 0.0
       d_tke_dyn(:,:)= 0.0
       d_q_dyn2d(:)  = 0.0
       d_ql_dyn2d(:) = 0.0
       d_qs_dyn2d(:) = 0.0
       d_qbs_dyn2d(:)= 0.0
       ! !! RomP >>>   td dyn traceur
       IF (nqtot > (nqo+nqtke)) d_tr_dyn(:,:,:)= 0.0
       ! !! RomP <<<
       ancien_ok = .TRUE.
    ENDIF
    !
    ! Ajouter le geopotentiel du sol:
    !
    DO k = 1, klev
       DO i = 1, klon
          zphi(i,k) = pphi(i,k) + pphis(i)
       ENDDO
    ENDDO
    !
    ! Verifier les temperatures
    !
    !IM BEG
    IF (check) THEN
       amn=MIN(ftsol(1,is_ter),1000.)
       amx=MAX(ftsol(1,is_ter),-1000.)
       DO i=2, klon
          amn=MIN(ftsol(i,is_ter),amn)
          amx=MAX(ftsol(i,is_ter),amx)
       ENDDO
       !
       PRINT*,' debut avant hgardfou min max ftsol',itap,amn,amx
    ENDIF !(check) THEN
    !IM END
    !
    CALL hgardfou(t_seri,ftsol,'debutphy',abortphy)
    IF (abortphy==1) Print*,'ERROR ABORT hgardfou debutphy'

    !
    !IM BEG
    IF (check) THEN
       amn=MIN(ftsol(1,is_ter),1000.)
       amx=MAX(ftsol(1,is_ter),-1000.)
       DO i=2, klon
          amn=MIN(ftsol(i,is_ter),amn)
          amx=MAX(ftsol(i,is_ter),amx)
       ENDDO
       !
       PRINT*,' debut apres hgardfou min max ftsol',itap,amn,amx
    ENDIF !(check) THEN
    !IM END
    !
    ! Mettre en action les conditions aux limites (albedo, sst, etc.).
    ! Prescrire l'ozone et calculer l'albedo sur l'ocean.
    !
    ! Update ozone if day change
    IF (MOD(itap-1,lmt_pas) == 0) THEN
       IF (read_climoz <= 0) THEN
          ! Once per day, update ozone from Royer:
          IF (solarlong0<-999.) then
             ! Generic case with evolvoing season
             zzz=real(days_elapsed+1)
          ELSE IF (abs(solarlong0-1000.)<1.e-4) then
             ! Particular case with annual mean insolation
             zzz=real(90) ! could be revisited
             IF (read_climoz/=-1) THEN
                abort_message ='read_climoz=-1 is recommended when ' &
                     // 'solarlong0=1000.'
                CALL abort_physic (modname,abort_message,1)
             ENDIF
          ELSE
             ! Case where the season is imposed with solarlong0
             zzz=real(90) ! could be revisited
          ENDIF

          wo(:,:,1)=ozonecm(latitude_deg, paprs,read_climoz,rjour=zzz)
IF (CPPKEY_REPROBUS) THEN
          ptrop=dyn_tropopause(t_seri, ztsol, paprs, pplay, rot)/100.
          DO i = 1, klon
             Z1=t_seri(i,itroprep(i)+1)
             Z2=t_seri(i,itroprep(i))
             fac=(Z1-Z2)/alog(pplay(i,itroprep(i)+1)/pplay(i,itroprep(i)))
             B=Z2-fac*alog(pplay(i,itroprep(i)))
             ttrop(i)= fac*alog(ptrop(i))+B
!
             Z1= 1.e-3 * ( pphi(i,itroprep(i)+1)+pphis(i) ) / gravit
             Z2= 1.e-3 * ( pphi(i,itroprep(i))  +pphis(i) ) / gravit
             fac=(Z1-Z2)/alog(pplay(i,itroprep(i)+1)/pplay(i,itroprep(i)))
             B=Z2-fac*alog(pplay(i,itroprep(i)))
             ztrop(i)=fac*alog(ptrop(i))+B
          ENDDO
END IF
       ELSE
          !--- ro3i = elapsed days number since current year 1st january, 0h
          ro3i=days_elapsed+jh_cur-jh_1jan
          !--- scaling for old style files (360 records)
          IF(SIZE(time_climoz)==360.AND..NOT.ok_daily_climoz) ro3i=ro3i*360./year_len
          IF(adjust_tropopause) THEN
             CALL stratosphere_mask(missing_val, pphis, t_seri, pplay, latitude_deg)
             IF(iflag_tropopause_height.EQ.0) THEN
               dyntropo=dyn_tropopause(t_seri, ztsol, paprs, pplay, rot)
             ELSEIF(iflag_tropopause_height.EQ.1) THEN
               dyntropo=p_tropopause
             ENDIF
             CALL regr_pr_time_av(ncid_climoz, vars_climoz(1:read_climoz),   &
                      ro3i, 'C', press_cen_climoz, pplay, wo, paprs(:,1),    &
                      time_climoz ,  longitude_deg,   latitude_deg,          &
                      dyntropo)
          ELSE
             CALL regr_pr_time_av(ncid_climoz,  vars_climoz(1:read_climoz),  &
                      ro3i, 'C', press_cen_climoz, pplay, wo, paprs(:,1),    &
                      time_climoz )
          ENDIF
          ! Convert from mole fraction of ozone to column density of ozone in a
          ! cell, in kDU:
          FORALL (l = 1: read_climoz) wo(:, :, l) = wo(:, :, l) * rmo3 / rmd &
               * zmasse / dobson_u / 1e3
          ! (By regridding ozone values for LMDZ only once a day, we
          ! have already neglected the variation of pressure in one
          ! day. So do not recompute "wo" at each time step even if
          ! "zmasse" changes a little.)
       ENDIF
    ENDIF

    !-- Needed for LSCP - condensation and ice supersaturation
    IF (ok_ice_supersat) THEN
      DO k = 1, klev
        DO i = 1, klon
          IF ( ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) ) .GT. 0. ) THEN
            ratio_ql_qtot(i,k) = ql_seri(i,k) / ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) )
            ratio_qi_qtot(i,k) = qs_seri(i,k) / ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) )
            rvc_seri(i,k) = rvc_seri(i,k) * q_seri(i,k) / ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) )
          ELSE
            ratio_ql_qtot(i,k) = 0.
            ratio_qi_qtot(i,k) = 0.
            rvc_seri(i,k) = 0.
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO k = 1, klev
        DO i = 1, klon
          IF ( ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) ) .GT. 0. ) THEN
            ratio_ql_qtot(i,k) = ql_seri(i,k) / ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) )
            ratio_qi_qtot(i,k) = qs_seri(i,k) / ( q_seri(i,k) + ql_seri(i,k) + qs_seri(i,k) )
          ELSE
            ratio_ql_qtot(i,k) = 0.
            ratio_qi_qtot(i,k) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    !
    ! Re-evaporer l'eau liquide nuageuse
    !
     CALL reevap (klon,klev,iflag_ice_thermo,t_seri,q_seri,ql_seri,qs_seri, &
   &         d_t_eva,d_q_eva,d_ql_eva,d_qi_eva)

     CALL add_phys_tend &
            (du0,dv0,d_t_eva,d_q_eva,d_ql_eva,d_qi_eva,dqbs0,paprs,&
               'eva',abortphy,flag_inhib_tend,itap,0)
    CALL prt_enerbil('eva',itap)

    !=========================================================================
    ! Calculs de l'orbite.
    ! Necessaires pour le rayonnement et la surface (calcul de l'albedo).
    ! doit donc etre plac\'e avant radlwsw et pbl_surface

    ! !!   jyg 17 Sep 2010 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL ymds2ju(year_cur, mth_eq, day_eq,0., jD_eq)
    day_since_equinox = (jD_cur + jH_cur) - jD_eq
    !
    !   choix entre calcul de la longitude solaire vraie ou valeur fixee a
    !   solarlong0
    IF (solarlong0<-999.) THEN
       IF (new_orbit) THEN
          ! calcul selon la routine utilisee pour les planetes
          CALL solarlong(day_since_equinox, zlongi, dist)
       ELSE
          ! calcul selon la routine utilisee pour l'AR4
          CALL orbite(REAL(days_elapsed+1),zlongi,dist)
       ENDIF
    ELSE
       zlongi=solarlong0  ! longitude solaire vraie
       dist=1.            ! distance au soleil / moyenne
    ENDIF

    IF (prt_level.ge.1) write(lunout,*)'Longitude solaire ',zlongi,solarlong0,dist


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calcul de l'ensoleillement :
    ! ============================
    ! Pour une solarlong0=1000., on calcule un ensoleillement moyen sur
    ! l'annee a partir d'une formule analytique.
    ! Cet ensoleillement est sym\'etrique autour de l'\'equateur et
    ! non nul aux poles.
    IF (abs(solarlong0-1000.)<1.e-4) THEN
       CALL zenang_an(iflag_cycle_diurne.GE.1,jH_cur, &
            latitude_deg,longitude_deg,rmu0,fract)
       swradcorr(:) = 1.0
       JrNt(:) = 1.0
       zrmu0(:) = rmu0(:)
    ELSE
       ! recode par Olivier Boucher en sept 2015
       SELECT CASE (iflag_cycle_diurne)
       CASE(0)
          !  Sans cycle diurne
          CALL angle(zlongi, latitude_deg, fract, rmu0)
          swradcorr = 1.0
          JrNt = 1.0
          zrmu0 = rmu0
       CASE(1)
          !  Avec cycle diurne sans application des poids
          !  bit comparable a l ancienne formulation cycle_diurne=true
          !  on integre entre gmtime et gmtime+radpas
          zdtime=phys_tstep*REAL(radpas) ! pas de temps du rayonnement (s)
          CALL zenang(zlongi,jH_cur,0.0,zdtime, &
               latitude_deg,longitude_deg,rmu0,fract)
          zrmu0 = rmu0
          swradcorr = 1.0
          ! Calcul du flag jour-nuit
          JrNt = 0.0
          WHERE (fract.GT.0.0) JrNt = 1.0
       CASE(2)
          !  Avec cycle diurne sans application des poids
          !  On integre entre gmtime-pdtphys et gmtime+pdtphys*(radpas-1)
          !  Comme cette routine est appele a tous les pas de temps de
          !  la physique meme si le rayonnement n'est pas appele je
          !  remonte en arriere les radpas-1 pas de temps
          !  suivant. Petite ruse avec MOD pour prendre en compte le
          !  premier pas de temps de la physique pendant lequel
          !  itaprad=0
          zdtime1=phys_tstep*REAL(-MOD(itaprad,radpas)-1)
          zdtime2=phys_tstep*REAL(radpas-MOD(itaprad,radpas)-1)
          CALL zenang(zlongi,jH_cur,zdtime1,zdtime2, &
               latitude_deg,longitude_deg,rmu0,fract)
          !
          ! Calcul des poids
          !
          zdtime1=-phys_tstep !--on corrige le rayonnement pour representer le
          zdtime2=0.0    !--pas de temps de la physique qui se termine
          CALL zenang(zlongi,jH_cur,zdtime1,zdtime2, &
               latitude_deg,longitude_deg,zrmu0,zfract)
          swradcorr = 0.0
          WHERE (rmu0.GE.1.e-10 .OR. fract.GE.1.e-10) &
               swradcorr=zfract/fract*zrmu0/rmu0
          ! Calcul du flag jour-nuit
          JrNt = 0.0
          WHERE (zfract.GT.0.0) JrNt = 1.0
       END SELECT
    ENDIF
    sza_o = ACOS (rmu0) *180./pi

    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Appel au pbl_surface : Planetary Boudary Layer et Surface
    ! Cela implique tous les interactions des sous-surfaces et la
    ! partie diffusion turbulent du couche limit.
    !
    ! Certains varibales de sorties de pbl_surface sont utiliser que pour
    ! ecriture des fihiers hist_XXXX.nc, ces sont :
    !   qsol,      zq2m,      s_pblh,  s_lcl,
    !   s_capCL,   s_oliqCL,  s_cteiCL,s_pblT,
    !   s_therm,   s_trmb1,   s_trmb2, s_trmb3,
    !   zu10m,     zv10m,   fder,
    !   zxqsurf,   delta_qsurf,
    !   rh2m,      zxfluxu, zxfluxv,
    !   frugs,     agesno,    fsollw,  fsolsw,
    !   d_ts,      fevap,     fluxlat, t2m,
    !   wfbils,    fluxt,   fluxu, fluxv,
    !
    ! Certains ne sont pas utiliser du tout :
    !   dsens, devap, zxsnow, zxfluxt, zxfluxq, q2m, fluxq
    !

    ! Calcul de l'humidite de saturation au niveau du sol

! Tests Fredho, instensibilite au pas de temps -------------------------------
! A detruire en 2024 une fois les tests documentes et les choix faits        !
! Conservation des variables avant l'appel à l a diffusion pour les tehrmic  !
    if (iflag_thermals_tenv / 10 == 1 ) then                                 !
        do k=1,klev                                                          !
           do i=1,klon                                                       !
              t_env(i,k)=t_seri(i,k)                                         !
              q_env(i,k)=q_seri(i,k)                                         !
           enddo                                                             !
        enddo                                                                !
    else if (iflag_thermals_tenv / 10 == 2 ) then                            !
        do k=1,klev                                                          !
           do i=1,klon                                                       !
              t_env(i,k)=t_seri(i,k)                                         !
           enddo                                                             !
        enddo                                                                !
    endif                                                                    !
! Tests Fredho, instensibilite au pas de temps -------------------------------


    IF (iflag_pbl/=0) THEN

       !jyg+nrlmd<
!!jyg       IF (prt_level .ge. 2 .and. mod(iflag_pbl_split,2) .eq. 1) THEN
       IF (prt_level .ge. 2 .and. mod(iflag_pbl_split,10) .ge. 1) THEN
          print *,'debut du splitting de la PBL, wake_s = ', wake_s(:)
          print *,'debut du splitting de la PBL, wake_deltat = ', wake_deltat(:,1)
          print *,'debut du splitting de la PBL, wake_deltaq = ', wake_deltaq(:,1)
       ENDIF
       ! !!
       !>jyg+nrlmd
       !
       !-------gustiness calculation-------!
       !ym : Warning gustiness non inialized for iflag_gusts=2 & iflag_gusts=3
       gustiness=0  !ym missing init

       IF (iflag_gusts==0) THEN
          gustiness(1:klon)=0
       ELSE IF (iflag_gusts==1) THEN
          gustiness(1:klon)=f_gust_bl*ale_bl(1:klon)+f_gust_wk*ale_wake(1:klon)
       ELSE IF (iflag_gusts==2) THEN
          gustiness(1:klon)=f_gust_bl*ale_bl_stat(1:klon)+f_gust_wk*ale_wake(1:klon)
       !!!! modif olivier torres
       ELSE IF (iflag_gusts==3) THEN
          w_et=wstar(1,3)
          jlr_g_s=(0.65*w_et)**2
          pr_et=rain_con*8640
          jlr_g_c = (((19.8*(pr_et(1:klon)**2))/(1.5+pr_et(1:klon)+pr_et(1:klon)**2))**(0.4))**2
          gustiness(1:klon)=jlr_g_c+jlr_g_s
!!       write(*,*) "rain ",pr_et
!!       write(*,*) "jlr_g_c",jlr_g_c
!!       write(*,*) "wstar",wstar(1,3)
!!       write(*,*) "jlr_g_s",jlr_g_s
          ! ELSE IF (iflag_gusts==2) THEN
          !    do i = 1, klon
          !       gustiness(i)=f_gust_bl*ale_bl(i)+sigma_wk(i)*f_gust_wk&
          !           *ale_wake(i) !! need to make sigma_wk accessible here
          !    enddo
          ! ELSE IF (iflag_gusts==3) THEN
          !    do i = 1, klon
          !       gustiness(i)=f_gust_bl*alp_bl(i)+f_gust_wk*alp_wake(i)
           !    enddo
        ENDIF


        CALL pbl_surface(  &
            phys_tstep,     date0,     itap,    days_elapsed+1, &
            debut,     lafin, &
            longitude_deg, latitude_deg, rugoro,  zrmu0,      &
        !GG    sollwdown,    cldt,      &
            sollwdown, pphi,   cldt,      &
        !GG
            rain_fall, snow_fall, bs_fall, solsw,   solswfdiff, sollw,     &
            gustiness,                                &
            t_seri,    q_seri,   qbs_seri,  u_seri,  v_seri,    &
                                !nrlmd+jyg<
            wake_deltat, wake_deltaq, wake_cstar, wake_s, &
                                !>nrlmd+jyg
            pplay,     paprs,     pctsrf,             &
            ftsol,SFRWL,falb_dir,falb_dif,ustar,u10m,v10m,wstar, &
                                !albedo SB <<<
            cdragh,    cdragm,  u1,    v1,            &
            beta_aridity, &
            albsol_dir,   albsol_dif,   sens,    evap, snowerosion, icesub_lic, &
            albsol3_lic,runoff,   snowhgt,   qsnow, to_ice, sissnow, &
            zxtsol,    zxfluxlat, zt2m,    qsat2m,  zn2mout, &
            d_t_vdf,   d_q_vdf, d_qbs_vdf,  d_u_vdf, d_v_vdf, d_t_diss, &
                                !nrlmd<
                                !jyg<
            d_t_vdf_w, d_q_vdf_w, &
            d_t_vdf_x, d_q_vdf_x, &
            sens_x, zxfluxlat_x, sens_w, zxfluxlat_w, &
                                !>jyg
            delta_tsurf,wake_dens, &
            cdragh_x,cdragh_w,cdragm_x,cdragm_w, &
            kh,kh_x,kh_w, &
                                !>nrlmd
            coefh(1:klon,1:klev,1:nbsrf+1), coefm(1:klon,1:klev,1:nbsrf+1), &
            slab_wfbils,                 &
            qsol,      zq2m,      s_pblh,  s_lcl, &
                                !jyg<
            s_pblh_x, s_lcl_x, s_pblh_w, s_lcl_w, &
                                !>jyg
            s_capCL,   s_oliqCL,  s_cteiCL,s_pblT, &
            s_therm,   s_trmb1,   s_trmb2, s_trmb3, &
            zustar, zu10m,     zv10m,   fder, &
            zxqsurf, delta_qsurf,   rh2m,      zxfluxu, zxfluxv, &
            z0m, z0h,     agesno,    fsollw,  fsolsw, &
            d_ts,      fevap,     fluxlat, t2m, &
            wfbils, wfevap, &
            fluxt,   fluxu,  fluxv, &
            dsens,     devap,     zxsnow, &
            zxfluxt,   zxfluxq,  zxfluxqbs,  q2m, fluxq, fluxqbs, pbl_tke, pbl_eps,  &
                                !nrlmd+jyg<
            wake_delta_pbl_TKE, &
                                !>nrlmd+jyg
!GG             treedrg )
            treedrg,hice, tice, bilg_cumul, &
            fcds, fcdi, dh_basal_growth, dh_basal_melt, &
            dh_top_melt, dh_snow2sic, &
            dtice_melt, dtice_snow2sic , &
!GG
!FC
!AM 
            tsurf_tersrf, tsoil_tersrf, qsurf_tersrf, tsurf_new_tersrf, &
            cdragm_tersrf, cdragh_tersrf, &
            swnet_tersrf, lwnet_tersrf, fluxsens_tersrf, fluxlat_tersrf)
       !
       !  Add turbulent diffusion tendency to the wake difference variables
!!jyg       IF (mod(iflag_pbl_split,2) .NE. 0) THEN
       IF (mod(iflag_pbl_split,10) .NE. 0) THEN
!jyg<
          d_deltat_vdf(:,:) = d_t_vdf_w(:,:)-d_t_vdf_x(:,:)
          d_deltaq_vdf(:,:) = d_q_vdf_w(:,:)-d_q_vdf_x(:,:)
          CALL add_wake_tend &
             (d_deltat_vdf, d_deltaq_vdf, dsig0, dsig0, ddens0, ddens0, wkoccur1, 'vdf', abortphy)
       ELSE
          d_deltat_vdf(:,:) = 0.
          d_deltaq_vdf(:,:) = 0.
!>jyg
       ENDIF

       !---------------------------------------------------------------------
       ! ajout des tendances de la diffusion turbulente
       IF (klon_glo==1) THEN
          CALL add_pbl_tend &
               (d_u_vdf,d_v_vdf,d_t_vdf+d_t_diss,d_q_vdf,dql0,dqi0,d_qbs_vdf,paprs,&
               'vdf',abortphy,flag_inhib_tend,itap)
       ELSE
          CALL add_phys_tend &
               (d_u_vdf,d_v_vdf,d_t_vdf+d_t_diss,d_q_vdf,dql0,dqi0,d_qbs_vdf,paprs,&
               'vdf',abortphy,flag_inhib_tend,itap,0)
       ENDIF
       CALL prt_enerbil('vdf',itap)

       !--------------------------------------------------------------------

       IF (mydebug) THEN
          CALL writefield_phy('u_seri',u_seri,nbp_lev)
          CALL writefield_phy('v_seri',v_seri,nbp_lev)
          CALL writefield_phy('t_seri',t_seri,nbp_lev)
          CALL writefield_phy('q_seri',q_seri,nbp_lev)
       ENDIF

       !albedo SB >>>
       albsol1=0.
       albsol2=0.
       falb1=0.
       falb2=0.
       SELECT CASE(nsw)
       CASE(2)
          albsol1=albsol_dir(:,1)
          albsol2=albsol_dir(:,2)
          falb1=falb_dir(:,1,:)
          falb2=falb_dir(:,2,:)
       CASE(4)
          albsol1=albsol_dir(:,1)
          albsol2=albsol_dir(:,2)*SFRWL(2)+albsol_dir(:,3)*SFRWL(3) &
               +albsol_dir(:,4)*SFRWL(4)
          albsol2=albsol2/(SFRWL(2)+SFRWL(3)+SFRWL(4))
          falb1=falb_dir(:,1,:)
          falb2=falb_dir(:,2,:)*SFRWL(2)+falb_dir(:,3,:)*SFRWL(3) &
               +falb_dir(:,4,:)*SFRWL(4)
          falb2=falb2/(SFRWL(2)+SFRWL(3)+SFRWL(4))
       CASE(6)
          albsol1=albsol_dir(:,1)*SFRWL(1)+albsol_dir(:,2)*SFRWL(2) &
               +albsol_dir(:,3)*SFRWL(3)
          albsol1=albsol1/(SFRWL(1)+SFRWL(2)+SFRWL(3))
          albsol2=albsol_dir(:,4)*SFRWL(4)+albsol_dir(:,5)*SFRWL(5) &
               +albsol_dir(:,6)*SFRWL(6)
          albsol2=albsol2/(SFRWL(4)+SFRWL(5)+SFRWL(6))
          falb1=falb_dir(:,1,:)*SFRWL(1)+falb_dir(:,2,:)*SFRWL(2) &
               +falb_dir(:,3,:)*SFRWL(3)
          falb1=falb1/(SFRWL(1)+SFRWL(2)+SFRWL(3))
          falb2=falb_dir(:,4,:)*SFRWL(4)+falb_dir(:,5,:)*SFRWL(5) &
               +falb_dir(:,6,:)*SFRWL(6)
          falb2=falb2/(SFRWL(4)+SFRWL(5)+SFRWL(6))
       END SELECt
       !albedo SB <<<


       CALL evappot(klon,nbsrf,ftsol,pplay(:,1),cdragh, &
            t_seri(:,1),q_seri(:,1),u_seri(:,1),v_seri(:,1),evap_pot)

    ENDIF

    ! ==================================================================
    ! Blowing snow sublimation and sedimentation

    d_t_bsss(:,:)=0.
    d_q_bsss(:,:)=0.
    d_qbs_bsss(:,:)=0.
    bsfl(:,:)=0.
    bs_fall(:)=0.
    IF (ok_bs) THEN

     CALL call_blowing_snow_sublim_sedim(klon,klev,phys_tstep,t_seri,q_seri,qbs_seri,pplay,paprs, &
                                        d_t_bsss,d_q_bsss,d_qbs_bsss,bsfl,bs_fall)

     CALL add_phys_tend &
               (du0,dv0,d_t_bsss,d_q_bsss,dql0,dqi0,d_qbs_bsss,paprs,&
               'bsss',abortphy,flag_inhib_tend,itap,0)

    ENDIF

    ! =================================================================== c
    !   Calcul de Qsat

    DO k = 1, klev
       DO i = 1, klon
          zx_t = t_seri(i,k)
          IF (thermcep) THEN
             zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
             zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
             zx_qs  = MIN(0.5,zx_qs)
             zcor   = 1./(1.-retv*zx_qs)
             zx_qs  = zx_qs*zcor
          ELSE
             !!           IF (zx_t.LT.t_coup) THEN             !jyg
             IF (zx_t.LT.rtt) THEN                  !jyg
                zx_qs = qsats(zx_t)/pplay(i,k)
             ELSE
                zx_qs = qsatl(zx_t)/pplay(i,k)
             ENDIF
          ENDIF
          zqsat(i,k)=zx_qs
       ENDDO
    ENDDO

    IF (prt_level.ge.1) THEN
       write(lunout,*) 'L   qsat (g/kg) avant clouds_gno'
       write(lunout,'(i4,f15.4)') (k,1000.*zqsat(igout,k),k=1,klev)
    ENDIF
    !
    ! Appeler la convection (au choix)
    !
    DO k = 1, klev
       DO i = 1, klon
          conv_q(i,k) = d_q_dyn(i,k)  &
               + d_q_vdf(i,k)/phys_tstep
          conv_t(i,k) = d_t_dyn(i,k)  &
               + d_t_vdf(i,k)/phys_tstep
       ENDDO
    ENDDO

    ! Calcule de vitesse verticale a partir de flux de masse verticale
    DO k = 1, klev
       DO i = 1, klon
          omega(i,k) = RG*flxmass_w(i,k) / cell_area(i)
       ENDDO
    ENDDO

    IF (prt_level.ge.1) write(lunout,*) 'omega(igout, :) = ', &
         omega(igout, :)
    !
    ! Appel de la convection tous les "cvpas"
    !
!!jyg    IF (MOD(itapcv,cvpas).EQ.0) THEN
!!    print *,' physiq : itapcv, cvpas, itap-1, cvpas_0 ', &
!!                       itapcv, cvpas, itap-1, cvpas_0
    IF (MOD(itapcv,cvpas).EQ.0 .OR. MOD(itap-1,cvpas_0).EQ.0) THEN

    !
    ! Mettre a zero des variables de sortie (pour securite)
    !
    pmflxr(:,:) = 0.
    pmflxs(:,:) = 0.
    wdtrainA(:,:) = 0.
    wdtrainS(:,:) = 0.
    wdtrainM(:,:) = 0.
    wdtrainAS(:,:) = 0.
    upwd(:,:) = 0.
    dnwd(:,:) = 0.
    ep(:,:) = 0.
    da(:,:)=0.
    mp(:,:)=0.
    wght_cvfd(:,:)=0.
    phi(:,:,:)=0.
    phi2(:,:,:)=0.
    epmlmMm(:,:,:)=0.
    eplaMm(:,:)=0.
    d1a(:,:)=0.
    dam(:,:)=0.
    elij(:,:,:)=0.
    ev(:,:)=0.
    qtaa(:,:)=0.
    clw(:,:)=0.
    sij(:,:,:)=0.
    !
    IF (iflag_con.EQ.1) THEN
       abort_message ='reactiver le call conlmd dans physiq.F'
       CALL abort_physic (modname,abort_message,1)
       !     CALL conlmd (phys_tstep, paprs, pplay, t_seri, q_seri, conv_q,
       !    .             d_t_con, d_q_con,
       !    .             rain_con, snow_con, ibas_con, itop_con)
    ELSE IF (iflag_con.EQ.2) THEN
       CALL conflx(phys_tstep, paprs, pplay, t_seri, q_seri, &
            conv_t, conv_q, -evap, omega, &
            d_t_con, d_q_con, rain_con, snow_con, &
            pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
            kcbot, kctop, kdtop, pmflxr, pmflxs)
       d_u_con = 0.
       d_v_con = 0.

       WHERE (rain_con < 0.) rain_con = 0.
       WHERE (snow_con < 0.) snow_con = 0.
       DO i = 1, klon
          ibas_con(i) = klev+1 - kcbot(i)
          itop_con(i) = klev+1 - kctop(i)
       ENDDO
    ELSE IF (iflag_con.GE.3) THEN
       ! nb of tracers for the KE convection:
       ! MAF la partie traceurs est faite dans phytrac
       ! on met ntra=1 pour limiter les appels mais on peut
       ! supprimer les calculs / ftra.
       ntra = 1

       !=======================================================================
       !ajout pour la parametrisation des poches froides: calcul de
       !t_w et t_x: si pas de poches froides, t_w=t_x=t_seri
       IF (iflag_wake>=1) THEN
         DO k=1,klev
            DO i=1,klon
                t_w(i,k) = t_seri(i,k) + (1-wake_s(i))*wake_deltat(i,k)
                q_w(i,k) = q_seri(i,k) + (1-wake_s(i))*wake_deltaq(i,k)
                t_x(i,k) = t_seri(i,k) - wake_s(i)*wake_deltat(i,k)
                q_x(i,k) = q_seri(i,k) - wake_s(i)*wake_deltaq(i,k)
            ENDDO
         ENDDO
       ELSE
                t_w(:,:) = t_seri(:,:)
                q_w(:,:) = q_seri(:,:)
                t_x(:,:) = t_seri(:,:)
                q_x(:,:) = q_seri(:,:)
       ENDIF
       !
       !jyg<
       ! Perform dry adiabatic adjustment on wake profile
       ! The corresponding tendencies are added to the convective tendencies
       ! after the call to the convective scheme.
       IF (iflag_wake>=1) then
          IF (iflag_adjwk >= 1) THEN
             limbas(:) = 1
             CALL ajsec(paprs, pplay, t_w, q_w, limbas, &
                  d_t_adjwk, d_q_adjwk)
             !
             DO k=1,klev
                DO i=1,klon
                   IF (wake_s(i) .GT. 1.e-3) THEN
                      t_w(i,k) = t_w(i,k) + d_t_adjwk(i,k)
                      q_w(i,k) = q_w(i,k) + d_q_adjwk(i,k)
                      d_deltat_ajs_cv(i,k) = d_t_adjwk(i,k)
                      d_deltaq_ajs_cv(i,k) = d_q_adjwk(i,k)
                   ELSE
                      d_deltat_ajs_cv(i,k) = 0.
                      d_deltaq_ajs_cv(i,k) = 0.
                   ENDIF
                ENDDO
             ENDDO
             IF (iflag_adjwk == 2 .AND. OK_bug_ajs_cv) THEN
               CALL add_wake_tend &
                 (d_deltat_ajs_cv, d_deltaq_ajs_cv, dsig0, dsig0, ddens0, ddens0, wkoccur1, 'ajs_cv', abortphy)
             ENDIF  ! (iflag_adjwk == 2 .AND. OK_bug_ajs_cv)
          ENDIF  ! (iflag_adjwk >= 1)
       ENDIF ! (iflag_wake>=1)
       !>jyg
       !

!!      print *,'physiq. q_w(1,k), q_x(1,k) ', &
!!             (k, q_w(1,k), q_x(1,k),k=1,25)

!jyg<
       CALL alpale( debut, itap, phys_tstep, paprs, omega, t_seri,   &
                    alp_offset, it_wape_prescr,  wape_prescr, fip_prescr, &
                    ale_bl_prescr, alp_bl_prescr, &
                    wake_pe, wake_fip,  &
                    Ale_bl, Ale_bl_trig, Alp_bl, &
                    Ale, Alp , Ale_wake, Alp_wake)
!>jyg
!
       ! sb, oct02:
       ! Schema de convection modularise et vectorise:
       ! (driver commun aux versions 3 et 4)
       !
       IF (ok_cvl) THEN ! new driver for convectL
          !
          !jyg<
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculate the upmost level of deep convection loops: k_upper_cv
          !  (near 22 km)
          k_upper_cv = klev
          !izero = klon/2+1/klon
          !DO k = klev,1,-1
          !   IF (pphi(izero,k) > 22.e4) k_upper_cv = k
          !ENDDO
          ! FH : nouveau calcul base sur un profil global sans quoi
          ! le modele etait sensible au decoupage de domaines
          DO k = klev,1,-1
             IF (-7*log(presnivs(k)/presnivs(1)) > 25.) k_upper_cv = k
          ENDDO
          IF (prt_level .ge. 5) THEN
             Print *, 'upmost level of deep convection loops: k_upper_cv = ', &
                  k_upper_cv
          ENDIF
          !
          !>jyg
          IF (type_trac == 'repr') THEN
             nbtr_tmp=ntra
          ELSE
             nbtr_tmp=nbtr
          ENDIF
          !jyg   iflag_con est dans clesphys
          !c          CALL concvl (iflag_con,iflag_clos,
          CALL concvl (iflag_clos, &
               phys_tstep, paprs, pplay, k_upper_cv, t_x,q_x, &
               t_w,q_w,wake_s, &
               u_seri,v_seri,tr_seri,nbtr_tmp, &
               ALE,ALP, &
               sig1,w01, &
               d_t_con,d_q_con,fqcomp,d_u_con,d_v_con,d_tr, &
               rain_con, snow_con, ibas_con, itop_con, sigd, &
               ema_cbmf,plcl,plfc,wbeff,convoccur,upwd,dnwd,dnwd0, &
               Ma,mipsh,Vprecip,cape,cin,tvp,Tconv,iflagctrl, &
               pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr,qcondc,wd, &
                                ! RomP >>>
                                !!     .        pmflxr,pmflxs,da,phi,mp,
                                !!     .        ftd,fqd,lalim_conv,wght_th)
               pmflxr,pmflxs, &
               coef_clos, coef_clos_eff, &
               da,phi,mp,phi2,d1a,dam,sij,qtaa,clw,elij, &
               ftd,fqd,lalim_conv,wght_th, &
               ev, ep,epmlmMm,eplaMm, &
               wdtrainA, wdtrainS, wdtrainM,wght_cvfd,qtc_cv,sigt_cv,detrain_cv, &
               tau_cld_cv,coefw_cld_cv,epmax_diag)

          ! RomP <<<

          !IM begin
          !       print*,'physiq: cin pbase dnwd0 ftd fqd ',cin(1),pbase(1),
          !    .dnwd0(1,1),ftd(1,1),fqd(1,1)
          !IM end
          !IM cf. FH
          clwcon0=qcondc
          pmfu(:,:)=upwd(:,:)+dnwd(:,:)
          fm_cv(:,:)=upwd(:,:)+dnwd(:,:)+dnwd0(:,:)
          !
          !jyg<
          ! If convective tendencies are too large, then call convection
          !  every time step
          cvpas = cvpas_0
          DO k=1,k_upper_cv
             DO i=1,klon
               IF (d_t_con(i,k) > 6.721 .AND. d_t_con(i,k) < 6.722 .AND.&
                   d_q_con(i,k) > -.0002171 .AND. d_q_con(i,k) < -.0002170) THEN
                     dtcon_multistep_max = 3.
                     dqcon_multistep_max = 0.02
               ENDIF
             ENDDO
          ENDDO
!
          DO k=1,k_upper_cv
             DO i=1,klon
!!               IF (abs(d_t_con(i,k)) > 0.24 .OR. &
!!                   abs(d_q_con(i,k)) > 2.e-2) THEN
               IF (abs(d_t_con(i,k)) > dtcon_multistep_max .OR. &
                   abs(d_q_con(i,k)) > dqcon_multistep_max) THEN
                 cvpas = 1
!!                 print *,'physiq1, i,k,d_t_con(i,k),d_q_con(i,k) ', &
!!                                   i,k,d_t_con(i,k),d_q_con(i,k)
               ENDIF
             ENDDO
          ENDDO
!!!   Ligne a ne surtout pas remettre sans avoir murement reflechi (jyg)
!!!          call bcast(cvpas)
!!!   ------------------------------------------------------------
          !>jyg
          !
          DO i = 1, klon
             IF (iflagctrl(i).le.1) itau_con(i)=itau_con(i)+cvpas
          ENDDO
          !
          !jyg<
          !    Add the tendency due to the dry adjustment of the wake profile
          IF (iflag_wake>=1) THEN
            IF (iflag_adjwk == 2) THEN
              DO k=1,klev
                 DO i=1,klon
                    ftd(i,k) = ftd(i,k) + wake_s(i)*d_t_adjwk(i,k)/phys_tstep
                    fqd(i,k) = fqd(i,k) + wake_s(i)*d_q_adjwk(i,k)/phys_tstep
                    d_t_con(i,k) = d_t_con(i,k) + wake_s(i)*d_t_adjwk(i,k)
                    d_q_con(i,k) = d_q_con(i,k) + wake_s(i)*d_q_adjwk(i,k)
                 ENDDO
              ENDDO
            ENDIF  ! (iflag_adjwk = 2)
          ENDIF   ! (iflag_wake>=1)
          !>jyg
          !
       ELSE ! ok_cvl

          ! MAF conema3 ne contient pas les traceurs
          CALL conema3 (phys_tstep, &
               paprs,pplay,t_seri,q_seri, &
               u_seri,v_seri,tr_seri,ntra, &
               sig1,w01, &
               d_t_con,d_q_con,d_u_con,d_v_con,d_tr, &
               rain_con, snow_con, ibas_con, itop_con, &
               upwd,dnwd,dnwd0,bas,top, &
               Ma,cape,tvp,rflag, &
               pbase &
               ,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr &
               ,clwcon0)

       ENDIF ! ok_cvl

       !
       ! Correction precip
       rain_con = rain_con * cvl_corr
       snow_con = snow_con * cvl_corr
       !

       IF (.NOT. ok_gust) THEN
          do i = 1, klon
             wd(i)=0.0
          enddo
       ENDIF

       ! =================================================================== c
       ! Calcul des proprietes des nuages convectifs
       !

       !   calcul des proprietes des nuages convectifs
       clwcon0(:,:)=fact_cldcon*clwcon0(:,:)
       IF (iflag_cld_cv == 0) THEN
          CALL clouds_gno &
               (klon,klev,q_seri,zqsat,clwcon0,ptconv,ratqsc,rnebcon0)
       ELSE
          CALL clouds_bigauss &
               (klon,klev,q_seri,zqsat,qtc_cv,sigt_cv,ptconv,ratqsc,rnebcon0)
       ENDIF


       ! =================================================================== c

       DO i = 1, klon
          itop_con(i) = min(max(itop_con(i),1),klev)
          ibas_con(i) = min(max(ibas_con(i),1),itop_con(i))
       ENDDO

       DO i = 1, klon
          ! C Risi modif: pour éviter pb de dépassement d'indice dans les cas
          ! où i n'est pas un point convectif et donc ibas_con(i)=0
          ! c'est un pb indépendant des isotopes
          if (ibas_con(i) > 0) then
             ema_pcb(i)  = paprs(i,ibas_con(i))
          else
             ema_pcb(i)  = 0.0
          endif
       ENDDO
       DO i = 1, klon
          ! L'idicage de itop_con peut cacher un pb potentiel
          ! FH sous la dictee de JYG, CR
          ema_pct(i)  = paprs(i,itop_con(i)+1)

          IF (itop_con(i).gt.klev-3) THEN
             IF (prt_level >= 9) THEN
                write(lunout,*)'La convection monte trop haut '
                write(lunout,*)'itop_con(,',i,',)=',itop_con(i)
             ENDIF
          ENDIF
       ENDDO
    ELSE IF (iflag_con.eq.0) THEN
       write(lunout,*) 'On n appelle pas la convection'
       clwcon0=0.
       rnebcon0=0.
       d_t_con=0.
       d_q_con=0.
       d_u_con=0.
       d_v_con=0.
       rain_con=0.
       snow_con=0.
       bas=1
       top=1
    ELSE
       WRITE(lunout,*) "iflag_con non-prevu", iflag_con
       CALL abort_physic("physiq", "", 1)
    ENDIF

    !--saving d_X_con * zmass for next timestep if convection is not called every timestep
    IF (ok_mass_dqcon) THEN
      d_q_con_zmasse(:,:) = d_q_con(:,:) * zmasse(:,:)
    ENDIF

    IF (ok_mass_dtcon) THEN
      d_t_con_zmasse(:,:) = d_t_con(:,:) * zmasse(:,:)
    ENDIF

    IF (ok_mass_duvcon) THEN
      d_u_con_zmasse(:,:) = d_u_con(:,:) * zmasse(:,:)
      d_v_con_zmasse(:,:) = d_v_con(:,:) * zmasse(:,:)
    ENDIF


    !     CALL homogene(paprs, q_seri, d_q_con, u_seri,v_seri,
    !    .              d_u_con, d_v_con)

!jyg    Reinitialize proba_notrig and itapcv when convection has been called
    proba_notrig(:) = 1.
    itapcv = 0
    ENDIF !  (MOD(itapcv,cvpas).EQ.0 .OR. MOD(itapcv,cvpas_0).EQ.0)
!
    itapcv = itapcv+1
    !
    ! Compter les steps ou cvpas=1
    IF (cvpas == 1) THEN
      Ncvpaseq1 = Ncvpaseq1+1
    ENDIF
    IF (mod(itap,1000) == 0) THEN
      print *,' physiq, nombre de steps ou cvpas = 1 : ', Ncvpaseq1
    ENDIF

!!!jyg  Appel diagnostique a add_phys_tend pour tester la conservation de
!!!     l'energie dans les courants satures.
!!    d_t_con_sat(:,:) = d_t_con(:,:) - ftd(:,:)*dtime
!!    d_q_con_sat(:,:) = d_q_con(:,:) - fqd(:,:)*dtime
!!    dql_sat(:,:) = (wdtrainA(:,:)+wdtrainM(:,:))*dtime/zmasse(:,:)
!!    CALL add_phys_tend(d_u_con, d_v_con, d_t_con_sat, d_q_con_sat, dql_sat,   &
!!                     dqi0, paprs, 'convection_sat', abortphy, flag_inhib_tend,&
!!                     itap, 1)
!!    call prt_enerbil('convection_sat',itap)
!!
!!

    !--recompute d_X_con with zmasse from new timestep
    IF (ok_mass_dqcon) THEN
      d_q_con(:,:)=d_q_con_zmasse(:,:)/zmasse(:,:)
    ENDIF

    IF (ok_mass_dtcon) THEN
      d_t_con(:,:)=d_t_con_zmasse(:,:)/zmasse(:,:)
    ENDIF

    IF (ok_mass_duvcon) THEN
      d_u_con(:,:)=d_u_con_zmasse(:,:)/zmasse(:,:)
      d_v_con(:,:)=d_v_con_zmasse(:,:)/zmasse(:,:)
    ENDIF



    CALL add_phys_tend(d_u_con, d_v_con, d_t_con, d_q_con, dql0, dqi0, dqbs0, paprs, &
         'convection',abortphy,flag_inhib_tend,itap,0)
    CALL prt_enerbil('convection',itap)

    !-------------------------------------------------------------------------

    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    !
    !==========================================================================
    !RR:Evolution de la poche froide: on ne fait pas de separation wake/env
    !pour la couche limite diffuse pour l instant
    !
    !
    ! nrlmd le 22/03/2011---Si on met les poches hors des thermiques
    ! il faut rajouter cette tendance calcul\'ee hors des poches
    ! froides
    !
    IF (iflag_wake>=1) THEN
       !
       !
       ! Call wakes every "wkpas" step
       !
       IF (MOD(itapwk,wkpas).EQ.0) THEN
          !
          DO k=1,klev
             DO i=1,klon
                dt_dwn(i,k)  = ftd(i,k)
                dq_dwn(i,k)  = fqd(i,k)
                M_dwn(i,k)   = dnwd0(i,k)
                M_up(i,k)    = upwd(i,k)
                dt_a(i,k)    = d_t_con(i,k)/phys_tstep - ftd(i,k)
                dq_a(i,k)    = d_q_con(i,k)/phys_tstep - fqd(i,k)
             ENDDO
          ENDDO

          IF (mod(iflag_wake,10)==2) THEN
             ok_wk_lsp(:)=max(sign(1.,wake_s(:)-wake_s_min_lsp),0.)
             DO k = 1,klev
                dt_dwn(:,k)= dt_dwn(:,k)+ &
                     ok_wk_lsp(:)*(d_t_eva(:,k)+d_t_lsc(:,k))/phys_tstep
                dq_dwn(:,k)= dq_dwn(:,k)+ &
                     ok_wk_lsp(:)*(d_q_eva(:,k)+d_q_lsc(:,k))/phys_tstep
             ENDDO
          ELSEIF (mod(iflag_wake,10)==3) THEN
             ok_wk_lsp(:)=max(sign(1.,wake_s(:)-wake_s_min_lsp),0.)
             DO k = 1,klev
                DO i=1,klon
                   IF (rneb(i,k)==0.) THEN
                      ! On ne tient compte des tendances qu'en dehors des
                      ! nuages (c'est-\`a-dire a priri dans une region ou
                      ! l'eau se reevapore).
                      dt_dwn(i,k)= dt_dwn(i,k)+ &
                           ok_wk_lsp(i)*d_t_lsc(i,k)/phys_tstep
                      dq_dwn(i,k)= dq_dwn(i,k)+ &
                           ok_wk_lsp(i)*d_q_lsc(i,k)/phys_tstep
                   ENDIF
                ENDDO
             ENDDO
          ENDIF

          !
          !calcul caracteristiques de la poche froide
          CALL calWAKE_first(phys_tstep)
          CALL calWAKE (iflag_wake_tend, paprs, pplay, phys_tstep, &
               t_seri, q_seri, omega,  &
               dt_dwn, dq_dwn, M_dwn, M_up,  &
               dt_a, dq_a, cv_gen,  &
               sigd, cin,  &
               wake_deltat, wake_deltaq, wake_s, awake_s, wake_dens, awake_dens,  &
               wake_dth, wake_h,  &
!!               wake_pe, wake_fip, wake_gfl,  &
               wake_pe, wake_fip_0, wake_gfl,  &   !! jyg
               d_t_wake, d_q_wake,  &
               wake_k, t_x, q_x,  &
               wake_omgbdth, wake_dp_omgb,  &
               wake_dtKE, wake_dqKE,  &
               wake_omg, wake_dp_deltomg,  &
               wake_spread, wake_Cstar, d_deltat_wk_gw,  &
               d_deltat_wk, d_deltaq_wk, d_s_wk, d_s_a_wk, d_dens_wk, d_dens_a_wk)
          !
          !jyg    Reinitialize itapwk when wakes have been called
          itapwk = 0
       ENDIF !  (MOD(itapwk,wkpas).EQ.0)
       !
       itapwk = itapwk+1
       !
       !-----------------------------------------------------------------------
       ! ajout des tendances des poches froides
       CALL add_phys_tend(du0,dv0,d_t_wake,d_q_wake,dql0,dqi0,dqbs0,paprs,'wake', &
            abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('wake',itap)
       !------------------------------------------------------------------------

       ! Increment Wake state variables
       IF (iflag_wake_tend .GT. 0.) THEN

         CALL add_wake_tend &
            (d_deltat_wk, d_deltaq_wk, d_s_wk, d_s_a_wk, d_dens_wk, d_dens_a_wk, wake_k, &
             'wake', abortphy)
          CALL prt_enerbil('wake',itap)
       ENDIF   ! (iflag_wake_tend .GT. 0.)
       !
       IF (prt_level .GE. 10) THEN
         print *,' physiq, after calwake, wake_s: ',wake_s(:)
         print *,' physiq, after calwake, wake_deltat: ',wake_deltat(:,1)
         print *,' physiq, after calwake, wake_deltaq: ',wake_deltaq(:,1)
       ENDIF

       IF (iflag_alp_wk_cond .GT. 0.) THEN

         CALL alpale_wk_first(cell_area)
         CALL alpale_wk(phys_tstep, cell_area, wake_k, wake_s, wake_dens, wake_fip_0, &
                        wake_fip)
       ELSE
         wake_fip(:) = wake_fip_0(:)
       ENDIF   ! (iflag_alp_wk_cond .GT. 0.)

    ENDIF  ! (iflag_wake>=1)
    !
    !===================================================================
    ! Convection seche (thermiques ou ajustement)
    !===================================================================
    !
    CALL stratocu_if(klon,klev,pctsrf,paprs, pplay,t_seri &
         ,seuil_inversion,weak_inversion,dthmin)



    d_t_ajsb(:,:)=0.
    d_q_ajsb(:,:)=0.
    d_t_ajs(:,:)=0.
    d_u_ajs(:,:)=0.
    d_v_ajs(:,:)=0.
    d_q_ajs(:,:)=0.
    clwcon0th(:,:)=0.
    !
    !      fm_therm(:,:)=0.
    !      entr_therm(:,:)=0.
    !      detr_therm(:,:)=0.
    !
    IF (prt_level>9) WRITE(lunout,*) &
         'AVANT LA CONVECTION SECHE , iflag_thermals=' &
         ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals
    IF (iflag_thermals<0) THEN
       !  Rien
       !  ====
       IF (prt_level>9) WRITE(lunout,*)'pas de convection seche'
       WRITE(lunout,*) 'WARNING : running without dry convection. Somme intermediate variables are not properly defined in physiq_mod.F90'
       ! Reprendre proprement les initialisation ci dessouds si on veut vraiment utiliser l'option (FH)
          fraca(:,:)=0.
          fm_therm(:,:)=0.
          ztv(:,:)=t_seri(:,:)
          zqasc(:,:)=q_seri(:,:)
          ztla(:,:)=0.
          zthl(:,:)=0.
          zpspsk(:,:)=(pplay(:,:)/100000.)**RKAPPA



    ELSE

       !  Thermiques
       !  ==========
       IF (prt_level>9) WRITE(lunout,*)'JUSTE AVANT , iflag_thermals=' &
            ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals


       !cc nrlmd le 10/04/2012
       DO k=1,klev+1
          DO i=1,klon
             pbl_tke_input(i,k,is_oce)=pbl_tke(i,k,is_oce)
             pbl_tke_input(i,k,is_ter)=pbl_tke(i,k,is_ter)
             pbl_tke_input(i,k,is_lic)=pbl_tke(i,k,is_lic)
             pbl_tke_input(i,k,is_sic)=pbl_tke(i,k,is_sic)
          ENDDO
       ENDDO
       !cc fin nrlmd le 10/04/2012

       IF (iflag_thermals>=1) THEN

! Tests Fredho, instensibilite au pas de temps -------------------------------
! A detruire en 2024 une fois les tests documentes et les choix faits        !
          if (iflag_thermals_tenv /10 == 0 ) then                            !
            do k=1,klev                                                      !
               do i=1,klon                                                   !
                  t_env(i,k)=t_seri(i,k)                                     !
                  q_env(i,k)=q_seri(i,k)                                     !
               enddo                                                         !
            enddo                                                            !
          else if (iflag_thermals_tenv / 10 == 2 ) then                      !
            do k=1,klev                                                      !
               do i=1,klon                                                   !
                  q_env(i,k)=q_seri(i,k)                                     !
               enddo                                                         !
            enddo                                                            !
          else if (iflag_thermals_tenv / 10 == 3 ) then                      !
            do k=1,klev                                                      !
               do i=1,klon                                                   !
                  t_env(i,k)=t(i,k)                                          !
                  q_env(i,k)=qx(i,k,1)                                       !
               enddo                                                         !
            enddo                                                            !
          endif                                                              !
! Tests Fredho, instensibilite au pas de temps ------------------------------

          !jyg<
!!       IF (mod(iflag_pbl_split/2,2) .EQ. 1) THEN
          IF (mod(iflag_pbl_split/10,10) .GE. 1) THEN
             !  Appel des thermiques avec les profils exterieurs aux poches
             DO k=1,klev
                DO i=1,klon
                   t_therm(i,k) = t_seri(i,k) - wake_s(i)*wake_deltat(i,k)
                   q_therm(i,k) = q_seri(i,k) - wake_s(i)*wake_deltaq(i,k)
                   t_env(i,k)   = t_env(i,k) - wake_s(i)*wake_deltat(i,k)
                   q_env(i,k)   = q_env(i,k) - wake_s(i)*wake_deltaq(i,k)
                   u_therm(i,k) = u_seri(i,k)
                   v_therm(i,k) = v_seri(i,k)
                ENDDO
             ENDDO
          ELSE
             !  Appel des thermiques avec les profils moyens
             DO k=1,klev
                DO i=1,klon
                   t_therm(i,k) = t_seri(i,k)
                   q_therm(i,k) = q_seri(i,k)
                   u_therm(i,k) = u_seri(i,k)
                   v_therm(i,k) = v_seri(i,k)
                ENDDO
             ENDDO
          ENDIF
          !>jyg
          CALL calltherm(itap, pdtphys &
               ,pplay,paprs,pphi,weak_inversion &
                        ! ,u_seri,v_seri,t_seri,q_seri,zqsat,debut & !jyg
               ,u_therm,v_therm,t_therm,q_therm,t_env,q_env,zqsat,debut &  !jyg
               ,d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs &
               ,fm_therm,entr_therm,detr_therm &
               ,zqasc,clwcon0th,lmax_th,ratqscth &
               ,ratqsdiff,zqsatth &
                                !on rajoute ale et alp, et les
                                !caracteristiques de la couche alim
               ,Ale_bl,Alp_bl,lalim_conv,wght_th, zmax0, f0, zw2,fraca &
               ,ztv,zpspsk,ztla,zthl &
                                !cc nrlmd le 10/04/2012
               ,pbl_tke_input,pctsrf,omega,cell_area &
               ,zlcl_th,fraca0,w0,w_conv,therm_tke_max0,env_tke_max0 &
               ,n2,s2,strig,zcong,ale_bl_stat &
               ,therm_tke_max,env_tke_max &
               ,alp_bl_det,alp_bl_fluct_m,alp_bl_fluct_tke &
               ,alp_bl_conv,alp_bl_stat &
                                !cc fin nrlmd le 10/04/2012
               ,zqla,ztva )
          !
          !jyg<
!!jyg          IF (mod(iflag_pbl_split/2,2) .EQ. 1) THEN
          IF (mod(iflag_pbl_split/10,10) .GE. 1) THEN
             !  Si les thermiques ne sont presents que hors des
             !  poches, la tendance moyenne associ\'ee doit etre
             !  multipliee par la fraction surfacique qu'ils couvrent.
             IF (mod(iflag_pbl_split/10,10) == 1) THEN
                ! On tient compte du splitting pour modifier les profils deltatq/T des poches
                DO k=1,klev
                   DO i=1,klon
                      d_deltat_the(i,k) = - d_t_ajs(i,k)
                      d_deltaq_the(i,k) = - d_q_ajs(i,k)
                   ENDDO
                ENDDO
             ELSE
                d_deltat_the(:,:) = 0.
                d_deltaq_the(:,:) = 0.
             ENDIF

             DO k=1,klev
                DO i=1,klon
                   d_u_ajs(i,k) = d_u_ajs(i,k)*(1.-wake_s(i))
                   d_v_ajs(i,k) = d_v_ajs(i,k)*(1.-wake_s(i))
                   d_t_ajs(i,k) = d_t_ajs(i,k)*(1.-wake_s(i))
                   d_q_ajs(i,k) = d_q_ajs(i,k)*(1.-wake_s(i))
                ENDDO
             ENDDO
          !
             IF (ok_bug_split_th) THEN
               CALL add_wake_tend &
                   (d_deltat_the, d_deltaq_the, dsig0, dsig0, ddens0, ddens0, wkoccur1, 'the', abortphy)
             ELSE
               CALL add_wake_tend &
                   (d_deltat_the, d_deltaq_the, dsig0, dsig0, ddens0, ddens0, wake_k, 'the', abortphy)
             ENDIF
             CALL prt_enerbil('the',itap)
          !
          ENDIF  ! (mod(iflag_pbl_split/10,10) .GE. 1)
          !
          CALL add_phys_tend(d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs,  &
                             dql0,dqi0,dqbs0,paprs,'thermals', abortphy,flag_inhib_tend,itap,0)
          CALL prt_enerbil('thermals',itap)
          !
!
          CALL alpale_th_first()
          CALL alpale_th( phys_tstep, lmax_th, t_seri, cell_area,  &
                          cin, s2, n2, strig, &
                          ale_bl_trig, ale_bl_stat, ale_bl,  &
                          alp_bl, alp_bl_stat, &
                          proba_notrig, random_notrig, cv_gen)
          !>jyg

          ! ------------------------------------------------------------------
          ! Transport de la TKE par les panaches thermiques.
          ! FH : 2010/02/01
               if (iflag_thermcell_tke==1) then
               call thermcell_dtke(klon,klev,nbsrf,pdtphys,fm_therm,entr_therm,rg,paprs,pbl_tke)
               endif
          ! -------------------------------------------------------------------

          DO i=1,klon
             !           zmax_th(i)=pphi(i,lmax_th(i))/rg
             !CR:04/05/12:correction calcul zmax
             zmax_th(i)=zmax0(i)
          ENDDO

       ENDIF

       !  Ajustement sec
       !  ==============

       ! Dans le cas o\`u on active les thermiques, on fait partir l'ajustement
       ! a partir du sommet des thermiques.
       ! Dans le cas contraire, on demarre au niveau 1.

       IF (iflag_thermals>=13.or.iflag_thermals<=0) THEN

          IF (iflag_thermals.eq.0) THEN
             IF (prt_level>9) WRITE(lunout,*)'ajsec'
             limbas(:)=1
          ELSE
             limbas(:)=lmax_th(:)
          ENDIF

          ! Attention : le call ajsec_convV2 n'est maintenu que momentanneement
          ! pour des test de convergence numerique.
          ! Le nouveau ajsec est a priori mieux, meme pour le cas
          ! iflag_thermals = 0 (l'ancienne version peut faire des tendances
          ! non nulles numeriquement pour des mailles non concernees.

          IF (iflag_thermals==0) THEN
             ! Calling adjustment alone (but not the thermal plume model)
             CALL ajsec_convV2(paprs, pplay, t_seri,q_seri &
                  , d_t_ajsb, d_q_ajsb)
          ELSE IF (iflag_thermals>0) THEN
             ! Calling adjustment above the top of thermal plumes
             CALL ajsec(paprs, pplay, t_seri,q_seri,limbas &
                  , d_t_ajsb, d_q_ajsb)
          ENDIF

          !--------------------------------------------------------------------
          ! ajout des tendances de l'ajustement sec ou des thermiques
          CALL add_phys_tend(du0,dv0,d_t_ajsb,d_q_ajsb,dql0,dqi0,dqbs0,paprs, &
               'ajsb',abortphy,flag_inhib_tend,itap,0)
          CALL prt_enerbil('ajsb',itap)
          d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_ajsb(:,:)
          d_q_ajs(:,:)=d_q_ajs(:,:)+d_q_ajsb(:,:)

          !---------------------------------------------------------------------

       ENDIF

    ENDIF
    !
    !===================================================================
    ! Computation of subrgid scale near-surface wind distribution
    ! Developed for dust lifting. Could be extended to coupling with ocean and others
    ! by default : 1 bin equal to the mean wind

     call surf_wind(klon,nsurfwind,zu10m,zv10m,wake_s,wake_Cstar,zustar,ale_bl,surf_wind_value,surf_wind_proba)

    !===================================================================
    ! Computation of ratqs, the width (normalized) of the subrid scale
    ! water distribution

    l_mix_ave(:,:)=0.
    wprime_ave(:,:)=0.

    DO nsrf = 1, nbsrf
       DO i = 1, klon
          l_mix_ave(i,:) = l_mix_ave(i,:) + l_mix(i,:,nsrf)*pctsrf(i,nsrf)
          wprime_ave(i,:) = wprime_ave(i,:) + wprime(i,:,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO

    CALL ratqs_main_first(klon, cell_area)
    CALL ratqs_main(klon,klev,nbsrf,prt_level,lunout,        &
         iflag_ratqs,iflag_con,iflag_cld_th,pdtphys,  &
         ratqsbas,ratqshaut,ratqsp0, ratqsdp, &
         pctsrf,s_pblh,zstd, &
         tau_ratqs,fact_cldcon,wake_s, wake_deltaq,   &
         ptconv,ptconvth,clwcon0th, rnebcon0th,     &
         paprs,pplay,t_seri,q_seri, &
         qtc_cv, sigt_cv,detrain_cv,fm_cv,fqd,fqcomp,sigd,zqsat, &
         omega,pbl_tke(:,:,is_ave),pbl_eps(:,:,is_ave),l_mix_ave,wprime_ave, &
         t2m,q2m,fm_therm,entr_therm,detr_therm,cell_area, &
         ratqs,ratqsc,ratqs_inter_,sigma_qtherm)

    !
    ! Appeler le processus de condensation a grande echelle
    ! et le processus de precipitation
    !-------------------------------------------------------------------------
    IF (prt_level .GE.10) THEN
       print *,'itap, ->fisrtilp ',itap
    ENDIF
    !

    picefra(:,:)=0.

    IF (ok_new_lscp) THEN

 
    DO k = 1, klev
      DO i = 1, klon
        ql_seri_lscp(i,k) = ratio_ql_qtot(i,k) * q_seri(i,k)
        qi_seri_lscp(i,k) = ratio_qi_qtot(i,k) * q_seri(i,k)
      ENDDO
    ENDDO


    !--mise à jour de flight_m et flight_h2o dans leur module
    !IF (ok_plane_h2o .OR. ok_plane_contrail) THEN
    !  CALL airplane(debut,pphis,pplay,paprs,t_seri)
    !ENDIF

    CALL lscp(klon,klev,phys_tstep,missing_val,paprs,pplay,omega, &
         t_seri, q_seri, ql_seri_lscp, qi_seri_lscp, ptconv, ratqs, sigma_qtherm, &
         d_t_lsc, d_q_lsc, d_ql_lsc, d_qi_lsc, rneb, rneblsvol, &
         pfraclr, pfracld, cldfraliq, cldfraliqth,              &
         sigma2_icefracturb, sigma2_icefracturbth,              &
         mean_icefracturb,  mean_icefracturbth,                 &
         radocond, picefra, rain_lsc, snow_lsc, &
         frac_impa, frac_nucl, beta_prec_fisrt, &
         prfl, psfl, rhcl,  &
         zqasc, fraca,ztv,zpspsk,ztla,zthl,zw2,iflag_cld_th, &
         iflag_ice_thermo, distcltop, temp_cltop,   &
         pbl_tke(:,:,is_ave), pbl_eps(:,:,is_ave), &
         entr_therm, detr_therm, &
         cell_area, &
         cf_seri, rvc_seri, u_seri, v_seri, &
         qsub, qissr, qcld, subfra, issrfra, gamma_cond,  &
         dcf_sub, dcf_con, dcf_mix, dqi_adj, dqi_sub, dqi_con, dqi_mix, &
         dqvc_adj, dqvc_sub, dqvc_con, dqvc_mix, qsatliq, qsatice, &
         Tcontr, qcontr, qcontr2, fcontrN, fcontrP, &
         dcf_avi, dqi_avi, dqvc_avi, flight_dist, flight_h2o, &
         cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv, &
         qraindiag, qsnowdiag, dqreva, dqssub, dqrauto, dqrcol, dqrmelt, &
         dqrfreez, dqsauto, dqsagg, dqsrim, dqsmelt, dqsfreez)


    ELSE
    
    CALL fisrtilp_first(klon, klev, phys_tstep, pfrac_impa, pfrac_nucl, pfrac_1nucl)
    CALL fisrtilp(klon,klev,phys_tstep,paprs,pplay, &
         t_seri, q_seri,ptconv,ratqs,sigma_qtherm, &
         d_t_lsc, d_q_lsc, d_ql_lsc, d_qi_lsc, rneb, rneblsvol, radocond, &
         rain_lsc, snow_lsc, &
         pfrac_impa, pfrac_nucl, pfrac_1nucl, &
         frac_impa, frac_nucl, beta_prec_fisrt, &
         prfl, psfl, rhcl,  &
         zqasc, fraca,ztv,zpspsk,ztla,zthl,iflag_cld_th, &
         iflag_ice_thermo, &
         cloudth_sth,cloudth_senv,cloudth_sigmath,cloudth_sigmaenv)

    ENDIF
    !
    WHERE (rain_lsc < 0) rain_lsc = 0.
    WHERE (snow_lsc < 0) snow_lsc = 0.

!+JLD
!    write(*,9000) 'phys lsc',"enerbil: bil_q, bil_e,",rain_lsc+snow_lsc &
!        & ,((rcw-rcpd)*rain_lsc + (rcs-rcpd)*snow_lsc)*t_seri(1,1)-rlvtt*rain_lsc+rlstt*snow_lsc &
!        & ,rain_lsc,snow_lsc
!    write(*,9000) "rcpv","rcw",rcpv,rcw,rcs,t_seri(1,1)
!-JLD
    CALL add_phys_tend(du0,dv0,d_t_lsc,d_q_lsc,d_ql_lsc,d_qi_lsc,dqbs0,paprs, &
         'lsc',abortphy,flag_inhib_tend,itap,0)
    CALL prt_enerbil('lsc',itap)
    rain_num(:)=0.
    DO k = 1, klev
       DO i = 1, klon
          IF (ql_seri(i,k)>oliqmax) THEN
             rain_num(i)=rain_num(i)+(ql_seri(i,k)-oliqmax)*zmasse(i,k)/pdtphys
             ql_seri(i,k)=oliqmax
          ENDIF
       ENDDO
    ENDDO
    IF (nqo >= 3) THEN
    DO k = 1, klev
       DO i = 1, klon
          IF (qs_seri(i,k)>oicemax) THEN
             rain_num(i)=rain_num(i)+(qs_seri(i,k)-oicemax)*zmasse(i,k)/pdtphys
             qs_seri(i,k)=oicemax
          ENDIF
       ENDDO
    ENDDO
    ENDIF


!---------------------------------------------------------------------------
    DO k = 1, klev
       DO i = 1, klon
          cldfra(i,k) = rneb(i,k)
          ! keep only liquid droplets in radocond if not liqice_in_radocond
          IF (.NOT.liqice_in_radocond) radocond(i,k) = ql_seri(i,k)
       ENDDO
    ENDDO


    ! Option to activate the radiative effect of blowing snow (ok_rad_bs)
    ! makes sense only if the new large scale condensation scheme is active
    ! with the ok_icefra_lscp flag active as well

    IF (ok_bs .AND. ok_rad_bs) THEN
       IF (ok_new_lscp .AND. ok_icefra_lscp) THEN
           DO k=1,klev
             DO i=1,klon
                radocond(i,k)=radocond(i,k)+qbs_seri(i,k)
                picefra(i,k)=(radocond(i,k)*picefra(i,k)+qbs_seri(i,k))/(radocond(i,k))
                qbsfra=min(qbs_seri(i,k)/qbst_bs,1.0)
                cldfra(i,k)=max(cldfra(i,k),qbsfra)
             ENDDO
           ENDDO
       ELSE
          WRITE(lunout,*)"PAY ATTENTION, you try to activate the radiative effect of blowing snow"
          WRITE(lunout,*)"with ok_new_lscp=false and/or ok_icefra_lscp=false"
          abort_message='inconsistency in cloud phase for blowing snow'
          CALL abort_physic(modname,abort_message,1)
       ENDIF

    ENDIF

    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    !
    !-------------------------------------------------------------------
    !  PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT
    !-------------------------------------------------------------------

    ! 1. NUAGES CONVECTIFS
    !
    !IM cf FH
    !     IF (iflag_cld_th.eq.-1) THEN ! seulement pour Tiedtke
    IF (iflag_cld_th.le.-1) THEN ! seulement pour Tiedtke
       snow_tiedtke=0.
       !     print*,'avant calcul de la pseudo precip '
       !     print*,'iflag_cld_th',iflag_cld_th
       IF (iflag_cld_th.eq.-1) THEN
          rain_tiedtke=rain_con
       ELSE
          !       print*,'calcul de la pseudo precip '
          rain_tiedtke=0.
          !         print*,'calcul de la pseudo precip 0'
          DO k=1,klev
             DO i=1,klon
                IF (d_q_con(i,k).lt.0.) THEN
                   rain_tiedtke(i)=rain_tiedtke(i)-d_q_con(i,k)/pdtphys &
                        *(paprs(i,k)-paprs(i,k+1))/rg
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
       !     call dump2d(iim,jjm,rain_tiedtke(2:klon-1),'PSEUDO PRECIP ')
       !

       ! Nuages diagnostiques pour Tiedtke
       CALL diagcld1(paprs,pplay, &
                                !IM cf FH. rain_con,snow_con,ibas_con,itop_con,
            rain_tiedtke,snow_tiedtke,ibas_con,itop_con, &
            diafra,dialiq)
       DO k = 1, klev
          DO i = 1, klon
             IF (diafra(i,k).GT.cldfra(i,k)) THEN
                radocond(i,k) = dialiq(i,k)
                cldfra(i,k) = diafra(i,k)
             ENDIF
          ENDDO
       ENDDO

    ELSE IF (iflag_cld_th.ge.3) THEN
       !  On prend pour les nuages convectifs le max du calcul de la
       !  convection et du calcul du pas de temps precedent diminue d'un facteur
       !  facttemps
       facteur = pdtphys *facttemps
       DO k=1,klev
          DO i=1,klon
             rnebcon(i,k)=rnebcon(i,k)*facteur
             IF (rnebcon0(i,k)*clwcon0(i,k).GT.rnebcon(i,k)*clwcon(i,k)) THEN
                rnebcon(i,k)=rnebcon0(i,k)
                clwcon(i,k)=clwcon0(i,k)
             ENDIF
          ENDDO
       ENDDO

       !   On prend la somme des fractions nuageuses et des contenus en eau

       IF (iflag_cld_th>=5) THEN

          DO k=1,klev
             ptconvth(:,k)=fm_therm(:,k+1)>0.
          ENDDO

          IF (iflag_coupl==4) THEN

             ! Dans le cas iflag_coupl==4, on prend la somme des convertures
             ! convectives et lsc dans la partie des thermiques
             ! Le controle par iflag_coupl est peut etre provisoire.
             DO k=1,klev
                DO i=1,klon
                   IF (ptconv(i,k).AND.ptconvth(i,k)) THEN
                      radocond(i,k)=radocond(i,k)+rnebcon(i,k)*clwcon(i,k)
                      cldfra(i,k)=min(cldfra(i,k)+rnebcon(i,k),1.)
                   ELSE IF (ptconv(i,k)) THEN
                      cldfra(i,k)=rnebcon(i,k)
                      radocond(i,k)=rnebcon(i,k)*clwcon(i,k)
                   ENDIF
                ENDDO
             ENDDO

          ELSE IF (iflag_coupl==5) THEN
             DO k=1,klev
                DO i=1,klon
                   cldfra(i,k)=min(cldfra(i,k)+rnebcon(i,k),1.)
                   radocond(i,k)=radocond(i,k)+rnebcon(i,k)*clwcon(i,k)
                ENDDO
             ENDDO

          ELSE

             ! Si on est sur un point touche par la convection
             ! profonde et pas par les thermiques, on prend la
             ! couverture nuageuse et l'eau nuageuse de la convection
             ! profonde.

             !IM/FH: 2011/02/23
             ! definition des points sur lesquels ls thermiques sont actifs

             DO k=1,klev
                DO i=1,klon
                   IF (ptconv(i,k).AND. .NOT.ptconvth(i,k)) THEN
                      cldfra(i,k)=rnebcon(i,k)
                      radocond(i,k)=rnebcon(i,k)*clwcon(i,k)
                   ENDIF
                ENDDO
             ENDDO

          ENDIF

       ELSE

          ! Ancienne version
          cldfra(:,:)=min(max(cldfra(:,:),rnebcon(:,:)),1.)
          radocond(:,:)=radocond(:,:)+rnebcon(:,:)*clwcon(:,:)
       ENDIF

    ENDIF

    !     plulsc(:)=0.
    !     do k=1,klev,-1
    !        do i=1,klon
    !              zzz=prfl(:,k)+psfl(:,k)
    !           if (.not.ptconvth.zzz.gt.0.)
    !        enddo prfl, psfl,
    !     enddo
    !
    ! 2. NUAGES STARTIFORMES
    !
    IF (ok_stratus) THEN
       CALL diagcld2(paprs,pplay,t_seri,q_seri, diafra,dialiq)
       DO k = 1, klev
          DO i = 1, klon
             IF (diafra(i,k).GT.cldfra(i,k)) THEN
                radocond(i,k) = dialiq(i,k)
                cldfra(i,k) = diafra(i,k)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    ! Precipitation totale
    !
    DO i = 1, klon
       rain_fall(i) = rain_con(i) + rain_lsc(i)
       snow_fall(i) = snow_con(i) + snow_lsc(i)
    ENDDO
    !
    ! Calculer l'humidite relative pour diagnostique
    !
    DO k = 1, klev
       DO i = 1, klon
          zx_t = t_seri(i,k)
          IF (thermcep) THEN
             !!           if (iflag_ice_thermo.eq.0) then                 !jyg
             zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
             !!           else                                            !jyg
             !!           zdelta = MAX(0.,SIGN(1.,t_glace_min-zx_t))      !jyg
             !!           endif                                           !jyg
             zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
             zx_qs  = MIN(0.5,zx_qs)
             zcor   = 1./(1.-retv*zx_qs)
             zx_qs  = zx_qs*zcor
          ELSE
             !!           IF (zx_t.LT.t_coup) THEN             !jyg
             IF (zx_t.LT.rtt) THEN                  !jyg
                zx_qs = qsats(zx_t)/pplay(i,k)
             ELSE
                zx_qs = qsatl(zx_t)/pplay(i,k)
             ENDIF
          ENDIF
          zx_rh(i,k) = q_seri(i,k)/zx_qs
            IF (iflag_ice_thermo .GT. 0) THEN
          zx_rhl(i,k) = q_seri(i,k)/(qsatl(zx_t)/pplay(i,k))
          zx_rhi(i,k) = q_seri(i,k)/(qsats(zx_t)/pplay(i,k))
            ENDIF
          zqsat(i,k)=zx_qs
       ENDDO
    ENDDO

    !IM Calcul temp.potentielle a 2m (tpot) et temp. potentielle
    !   equivalente a 2m (tpote) pour diagnostique
    !
    DO i = 1, klon
       tpot(i)=zt2m(i)*(100000./paprs(i,1))**RKAPPA
       IF (thermcep) THEN
          IF(zt2m(i).LT.RTT) then
             Lheat=RLSTT
          ELSE
             Lheat=RLVTT
          ENDIF
       ELSE
          IF (zt2m(i).LT.RTT) THEN
             Lheat=RLSTT
          ELSE
             Lheat=RLVTT
          ENDIF
       ENDIF
       tpote(i) = tpot(i)*      &
            EXP((Lheat *qsat2m(i))/(RCPD*zt2m(i)))
    ENDDO

    IF (ANY(type_trac == ['inca','inco'])) THEN ! ModThL
IF (CPPKEY_INCA) THEN
       CALL VTe(VTphysiq)
       CALL VTb(VTinca)
       calday = REAL(days_elapsed + 1) + jH_cur

       CALL chemtime(itap+itau_phy-1, date0, phys_tstep, itap)
       CALL AEROSOL_METEO_CALC( &
            calday,pdtphys,pplay,paprs,t,pmflxr,pmflxs, &
            prfl,psfl,pctsrf,cell_area, &
            latitude_deg,longitude_deg,u10m,v10m)

       zxsnow_dummy(:) = 0.0

       CALL chemhook_begin (calday, &
            days_elapsed+1, &
            jH_cur, &
            pctsrf(1,1), &
            latitude_deg, &
            longitude_deg, &
            cell_area, &
            paprs, &
            pplay, &
            coefh(1:klon,1:klev,is_ave), &
            pphi, &
            t_seri, &
            u, &
            v, &
            rot, &
            wo(:, :, 1), &
            q_seri, &
            zxtsol, &
            zt2m, &
            zxsnow_dummy, &
            solsw, &
            albsol1, &
            rain_fall, &
            snow_fall, &
            itop_con, &
            ibas_con, &
            cldfra, &
            nbp_lon, &
            nbp_lat-1, &
            tr_seri(:,:,1+nqCO2:nbtr), &
            ftsol, &
            paprs, &
            cdragh, &
            cdragm, &
            pctsrf, &
            pdtphys, &
            itap)

       CALL VTe(VTinca)
       CALL VTb(VTphysiq)
END IF
    ENDIF !type_trac = inca or inco
    IF (type_trac == 'repr') THEN
IF (CPPKEY_REPROBUS) THEN
    !CALL chemtime_rep(itap+itau_phy-1, date0, dtime, itap)
    CALL chemtime_rep(itap+itau_phy-1, date0, phys_tstep, itap)
END IF
    ENDIF

    !
    ! Appeler le rayonnement mais calculer tout d'abord l'albedo du sol.
    !
    IF (MOD(itaprad,radpas).EQ.0) THEN

       !
       !jq - introduce the aerosol direct and first indirect radiative forcings
       !jq - Johannes Quaas, 27/11/2003 (quaas@lmd.jussieu.fr)
       IF (flag_aerosol .GT. 0) THEN
          IF (iflag_rrtm .EQ. 0) THEN !--old radiation
             IF (.NOT. aerosol_couple) THEN
                !
                CALL readaerosol_optic( &
                     debut, flag_aerosol, itap, jD_cur-jD_ref, &
                     pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
                     mass_solu_aero, mass_solu_aero_pi,  &
                     tau_aero, piz_aero, cg_aero,  &
                     tausum_aero, tau3d_aero)
             ENDIF
          ELSE IF (iflag_rrtm .EQ.1) THEN  ! RRTM radiation
             IF (aerosol_couple .AND. config_inca == 'aero' ) THEN
                abort_message='config_inca=aero et rrtm=1 impossible'
                CALL abort_physic(modname,abort_message,1)
             ELSE
                !
#ifdef CPP_RRTM
                IF (NSW.EQ.6) THEN
                   !--new aerosol properties SW and LW
                   !
IF (CPPKEY_DUST) THEN
                   !--SPL aerosol model
                   CALL splaerosol_optic_rrtm( ok_alw, pplay, paprs, t_seri, rhcl, &
                        tr_seri, mass_solu_aero, mass_solu_aero_pi,  &
                        tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm,  &
                        tausum_aero, tau3d_aero)
ELSE
                   !--climatologies or INCA aerosols
                   CALL readaerosol_optic_rrtm( debut, aerosol_couple, ok_alw, ok_volcan, &
                        flag_aerosol, flag_bc_internal_mixture, itap, jD_cur-jD_ref, &
                        pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
                        tr_seri, mass_solu_aero, mass_solu_aero_pi,  &
                        tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm,  &
                        tausum_aero, drytausum_aero, tau3d_aero)
END IF

                   IF (flag_aerosol .EQ. 7) THEN
                      CALL MACv2SP(pphis,pplay,paprs,longitude_deg,latitude_deg,  &
                                   tau_aero_sw_rrtm,piz_aero_sw_rrtm,cg_aero_sw_rrtm)
                   ENDIF

                   !
                ELSE IF (NSW.EQ.2) THEN
                   !--for now we use the old aerosol properties
                   !
                   CALL readaerosol_optic( &
                        debut, flag_aerosol, itap, jD_cur-jD_ref, &
                        pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
                        mass_solu_aero, mass_solu_aero_pi,  &
                        tau_aero, piz_aero, cg_aero,  &
                        tausum_aero, tau3d_aero)
                   !
                   !--natural aerosols
                   tau_aero_sw_rrtm(:,:,1,:)=tau_aero(:,:,3,:)
                   piz_aero_sw_rrtm(:,:,1,:)=piz_aero(:,:,3,:)
                   cg_aero_sw_rrtm (:,:,1,:)=cg_aero (:,:,3,:)
                   !--all aerosols
                   tau_aero_sw_rrtm(:,:,2,:)=tau_aero(:,:,2,:)
                   piz_aero_sw_rrtm(:,:,2,:)=piz_aero(:,:,2,:)
                   cg_aero_sw_rrtm (:,:,2,:)=cg_aero (:,:,2,:)
                   !
                   !--no LW optics
                   tau_aero_lw_rrtm(:,:,:,:) = 1.e-15
                   !
                ELSE
                   abort_message='Only NSW=2 or 6 are possible with ' &
                        // 'aerosols and iflag_rrtm=1'
                   CALL abort_physic(modname,abort_message,1)
                ENDIF
#else
                abort_message='You should compile with -rrtm if running ' &
                     // 'with iflag_rrtm=1'
                CALL abort_physic(modname,abort_message,1)
#endif
                !
             ENDIF
          ELSE IF (iflag_rrtm .EQ.2) THEN    ! ecrad RADIATION
#ifdef CPP_ECRAD
             !--climatologies or INCA aerosols
             CALL readaerosol_optic_ecrad( debut, aerosol_couple, ok_alw, ok_volcan, &
                  flag_aerosol, flag_bc_internal_mixture, itap, jD_cur-jD_ref, &
                  pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
                  tr_seri, mass_solu_aero, mass_solu_aero_pi, m_allaer)
#else
                abort_message='You should compile with -rad ecrad if running with iflag_rrtm=2'
                CALL abort_physic(modname,abort_message,1)
#endif
          ENDIF

       ELSE   !--flag_aerosol = 0
          tausum_aero(:,:,:) = 0.
          drytausum_aero(:,:) = 0.
          mass_solu_aero(:,:) = 0.
          mass_solu_aero_pi(:,:) = 0.
          IF (iflag_rrtm .EQ. 0) THEN !--old radiation
             tau_aero(:,:,:,:) = 1.e-15
             piz_aero(:,:,:,:) = 1.
             cg_aero(:,:,:,:)  = 0.
          ELSE
             tau_aero_sw_rrtm(:,:,:,:) = 1.e-15
             tau_aero_lw_rrtm(:,:,:,:) = 1.e-15
             piz_aero_sw_rrtm(:,:,:,:) = 1.0
             cg_aero_sw_rrtm(:,:,:,:)  = 0.0
          ENDIF
       ENDIF
       !
       !--WMO criterion to determine tropopause
       CALL stratosphere_mask(missing_val, pphis, t_seri, pplay, latitude_deg)
       !
       !--STRAT AEROSOL
       !--updates tausum_aero,tau_aero,piz_aero,cg_aero
       IF (flag_aerosol_strat.GT.0) THEN
          IF (prt_level .GE.10) THEN
             PRINT *,'appel a readaerosolstrat', mth_cur
          ENDIF
          IF (iflag_rrtm.EQ.0) THEN
           IF (flag_aerosol_strat.EQ.1) THEN
             CALL readaerosolstrato(debut)
           ELSE
             abort_message='flag_aerosol_strat must equal 1 for rrtm=0'
             CALL abort_physic(modname,abort_message,1)
           ENDIF
          ELSE
#ifdef CPP_RRTM
IF (.NOT. CPPKEY_STRATAER) THEN
          !--prescribed strat aerosols
          !--only in the case of non-interactive strat aerosols
            IF (flag_aerosol_strat.EQ.1) THEN
             CALL readaerosolstrato1_rrtm(debut)
            ELSEIF (flag_aerosol_strat.EQ.2) THEN
             CALL readaerosolstrato2_rrtm(debut, ok_volcan)
            ELSE
             abort_message='flag_aerosol_strat must equal 1 or 2 for rrtm=1'
             CALL abort_physic(modname,abort_message,1)
            ENDIF
END IF
#else
             abort_message='You should compile with -rrtm if running ' &
                  // 'with iflag_rrtm=1'
             CALL abort_physic(modname,abort_message,1)
#endif
          ENDIF
       ELSE
          tausum_aero(:,:,id_STRAT_phy) = 0.
       ENDIF
!
#ifdef CPP_RRTM
IF (CPPKEY_STRATAER) THEN
       !--compute stratospheric mask
       CALL stratosphere_mask(missing_val, pphis, t_seri, pplay, latitude_deg)
       !--interactive strat aerosols
       CALL calcaerosolstrato_rrtm(pplay,t_seri,paprs,debut)
END IF
#endif
       !--fin STRAT AEROSOL
       !

       ! Calculer les parametres optiques des nuages et quelques
       ! parametres pour diagnostiques:
       !
       IF (aerosol_couple.AND.config_inca=='aero') THEN
          mass_solu_aero(:,:)    = ccm(:,:,1)
          mass_solu_aero_pi(:,:) = ccm(:,:,2)
       ENDIF

       !Rajout appel a interface calcul proprietes optiques des nuages
       CALL call_cloud_optics_prop(klon, klev, ok_newmicro, &
               paprs, pplay, t_seri, radocond, picefra, cldfra, &
               cldtau, cldemi, cldh, cldl, cldm, cldt, cldq, &
               flwp, fiwp, flwc, fiwc, ok_aie, &
               mass_solu_aero, mass_solu_aero_pi, &
               cldtaupi, distcltop, temp_cltop, re, fl, ref_liq, ref_ice, &
               ref_liq_pi, ref_ice_pi, scdnc, cldncl, reffclwtop, lcc, reffclws, &
               reffclwc, cldnvi, lcc3d, lcc3dcon, lcc3dstra, icc3dcon, icc3dstra,  &
               zfice, dNovrN, ptconv, rnebcon, clwcon)
       CALL call_cloud_optics_prop_post(ok_newmicro)

       !
       !IM betaCRF
       !
       cldtaurad   = cldtau
       cldtaupirad = cldtaupi
       cldemirad   = cldemi
       cldfrarad   = cldfra

       !
       IF (lon1_beta.EQ.-180..AND.lon2_beta.EQ.180..AND. &
           lat1_beta.EQ.90..AND.lat2_beta.EQ.-90.) THEN
          !
          ! global
          !
!IM 251017 begin
!               print*,'physiq betaCRF global zdtime=',zdtime
!IM 251017 end
          DO k=1, klev
             DO i=1, klon
                IF (pplay(i,k).GE.pfree) THEN
                   beta(i,k) = beta_pbl
                ELSE
                   beta(i,k) = beta_free
                ENDIF
                IF (mskocean_beta) THEN
                   beta(i,k) = beta(i,k) * pctsrf(i,is_oce)
                ENDIF
                cldtaurad(i,k)   = cldtau(i,k) * beta(i,k)
                cldtaupirad(i,k) = cldtaupi(i,k) * beta(i,k)
                cldemirad(i,k)   = cldemi(i,k) * beta(i,k)
                cldfrarad(i,k)   = cldfra(i,k) * beta(i,k)
             ENDDO
          ENDDO
          !
       ELSE
          !
          ! regional
          !
          DO k=1, klev
             DO i=1,klon
                !
                IF (longitude_deg(i).ge.lon1_beta.AND. &
                    longitude_deg(i).le.lon2_beta.AND. &
                    latitude_deg(i).le.lat1_beta.AND.  &
                    latitude_deg(i).ge.lat2_beta) THEN
                   IF (pplay(i,k).GE.pfree) THEN
                      beta(i,k) = beta_pbl
                   ELSE
                      beta(i,k) = beta_free
                   ENDIF
                   IF (mskocean_beta) THEN
                      beta(i,k) = beta(i,k) * pctsrf(i,is_oce)
                   ENDIF
                   cldtaurad(i,k)   = cldtau(i,k) * beta(i,k)
                   cldtaupirad(i,k) = cldtaupi(i,k) * beta(i,k)
                   cldemirad(i,k)   = cldemi(i,k) * beta(i,k)
                   cldfrarad(i,k)   = cldfra(i,k) * beta(i,k)
                ENDIF
             !
             ENDDO
          ENDDO
       !
       ENDIF

       !lecture de la chlorophylle pour le nouvel albedo de Sunghye Baek
       IF (ok_chlorophyll) THEN
          print*,"-- reading chlorophyll"
          CALL readchlorophyll(debut)
       ENDIF

!--if ok_suntime_rrtm we use ancillay data for RSUN
!--previous values are therefore overwritten
!--this is needed for CMIP6 runs
!--and only possible for new radiation scheme
       IF (iflag_rrtm.EQ.1.AND.ok_suntime_rrtm) THEN
#ifdef CPP_RRTM
         CALL read_rsun_rrtm(debut)
#endif
       ENDIF

       IF (mydebug) THEN
          CALL writefield_phy('u_seri',u_seri,nbp_lev)
          CALL writefield_phy('v_seri',v_seri,nbp_lev)
          CALL writefield_phy('t_seri',t_seri,nbp_lev)
          CALL writefield_phy('q_seri',q_seri,nbp_lev)
       ENDIF

       !
       !sonia : If Iflag_radia >=2, pertubation of some variables
       !input to radiation (DICE)
       !
       IF (iflag_radia .ge. 2) THEN
          zsav_tsol (:) = zxtsol(:)
          CALL perturb_radlwsw(zxtsol,iflag_radia)
       ENDIF

       IF (aerosol_couple.AND.config_inca=='aero') THEN
IF (CPPKEY_INCA) THEN
          CALL radlwsw_inca  &
               (chemistry_couple, kdlon,kflev,dist, rmu0, fract, solaire, &
               paprs, pplay,zxtsol,albsol1, albsol2, t_seri,q_seri, &
               size(wo,3), wo, &
               cldfrarad, cldemirad, cldtaurad, &
               heat,heat0,cool,cool0,albpla, &
               topsw,toplw,solsw,sollw, &
               sollwdown, &
               topsw0,toplw0,solsw0,sollw0, &
               lwdn0, lwdn, lwup0, lwup,  &
               swdn0, swdn, swup0, swup, &
               ok_ade, ok_aie, &
               tau_aero, piz_aero, cg_aero, &
               topswad_aero, solswad_aero, &
               topswad0_aero, solswad0_aero, &
               topsw_aero, topsw0_aero, &
               solsw_aero, solsw0_aero, &
               cldtaupirad, &
               topswai_aero, solswai_aero)
END IF
       ELSE
          !
          !IM calcul radiatif pour le cas actuel
          !
          RCO2 = RCO2_act
          RCH4 = RCH4_act
          RN2O = RN2O_act
          RCFC11 = RCFC11_act
          RCFC12 = RCFC12_act
          !
          !--interactive CO2 in ppm from carbon cycle
          IF (carbon_cycle_rad) RCO2=RCO2_glo
          !
          IF (prt_level .GE.10) THEN
             print *,' ->radlwsw, number 1 '
          ENDIF
          !
          ! AI namelist utilise pour l appel principal de radlwsw (ecrad)
          namelist_ecrad_file='namelist_ecrad'
          !
          CALL radlwsw &
               (debut, dist, rmu0, fract,  &
                                !albedo SB >>>
                                !      paprs, pplay,zxtsol,albsol1, albsol2,  &
               paprs, pplay,zxtsol,SFRWL,albsol_dir, albsol_dif,  &
                                !albedo SB <<<
               t_seri,q_seri,wo, &
               cldfrarad, cldemirad, cldtaurad, &
               ok_ade.OR.flag_aerosol_strat.GT.0, ok_aie,  ok_volcan, flag_volc_surfstrat, &
               flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
               tau_aero, piz_aero, cg_aero, &
               tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm, &
               ! Rajoute par OB pour RRTM
               tau_aero_lw_rrtm, &
               cldtaupirad, m_allaer, &
!              zqsat, flwcrad, fiwcrad, &
               zqsat, flwc, fiwc, &
               ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
               namelist_ecrad_file, &
               heat,heat0,cool,cool0,albpla, &
               heat_volc,cool_volc, &
               topsw,toplw,solsw,solswfdiff,sollw, &
               sollwdown, &
               topsw0,toplw0,solsw0,sollw0, &
               lwdnc0, lwdn0, lwdn, lwupc0, lwup0, lwup,  &
               lwtoa0b, lwtoab , &  !FC
               swdnc0, swdn0, swdn, swupc0, swup0, swup, &
               topswad_aero, solswad_aero, &
               topswai_aero, solswai_aero, &
               topswad0_aero, solswad0_aero, &
               topsw_aero, topsw0_aero, &
               solsw_aero, solsw0_aero, &
               topswcf_aero, solswcf_aero, &
                                !-C. Kleinschmitt for LW diagnostics
               toplwad_aero, sollwad_aero,&
               toplwai_aero, sollwai_aero, &
               toplwad0_aero, sollwad0_aero,&
                                !-end
               ZLWFT0_i, ZFLDN0, ZFLUP0, &
               ZSWFT0_i, ZFSDN0, ZFSUP0, &
               ZFLUX_DIR, ZFLUX_DIR_CLEAR, ZFLUX_DIR_SUN, &
               cloud_cover_sw)

          !lwoff=y, betalwoff=1. : offset LW CRE for radiation code and other
          !schemes
          toplw = toplw + betalwoff * (toplw0 - toplw)
          sollw = sollw + betalwoff * (sollw0 - sollw)
          lwdn = lwdn + betalwoff * (lwdn0 - lwdn)
          lwup = lwup + betalwoff * (lwup0 - lwup)
          sollwdown(:)= sollwdown(:) + betalwoff *(-1.*ZFLDN0(:,1) - &
                        sollwdown(:))
          cool = cool + betalwoff * (cool0 - cool)

          IF (.NOT. using_xios) THEN
            !
            !IM 2eme calcul radiatif pour le cas perturbe ou au moins un
            !IM des taux doit etre different du taux actuel
            !IM Par defaut on a les taux perturbes egaux aux taux actuels
            !
            IF (RCO2_per.NE.RCO2_act.OR. &
                RCH4_per.NE.RCH4_act.OR. &
                RN2O_per.NE.RN2O_act.OR. &
                RCFC11_per.NE.RCFC11_act.OR. &
                RCFC12_per.NE.RCFC12_act) ok_4xCO2atm =.TRUE.
          ENDIF
   !
          IF (ok_4xCO2atm) THEN
                !
                RCO2 = RCO2_per
                RCH4 = RCH4_per
                RN2O = RN2O_per
                RCFC11 = RCFC11_per
                RCFC12 = RCFC12_per
                !
                IF (prt_level .GE.10) THEN
                   print *,' ->radlwsw, number 2 '
                ENDIF
                !
                namelist_ecrad_file='namelist_ecrad'
                !
                CALL radlwsw &
                     (debut, dist, rmu0, fract,  &
                                !albedo SB >>>
                                !      paprs, pplay,zxtsol,albsol1, albsol2,  &
                     paprs, pplay,zxtsol,SFRWL,albsol_dir, albsol_dif, &
                                !albedo SB <<<
                     t_seri,q_seri,wo, &
                     cldfrarad, cldemirad, cldtaurad, &
                     ok_ade.OR.flag_aerosol_strat.GT.0, ok_aie,  ok_volcan, flag_volc_surfstrat, &
                     flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
                     tau_aero, piz_aero, cg_aero, &
                     tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm, &
                                ! Rajoute par OB pour RRTM
                     tau_aero_lw_rrtm, &
                     cldtaupi, m_allaer, &
!                    zqsat, flwcrad, fiwcrad, &
                     zqsat, flwc, fiwc, &
                     ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
                     namelist_ecrad_file, &
                     heatp,heat0p,coolp,cool0p,albplap, &
                     heat_volc,cool_volc, &
                     topswp,toplwp,solswp,solswfdiffp,sollwp, &
                     sollwdownp, &
                     topsw0p,toplw0p,solsw0p,sollw0p, &
                     lwdnc0p, lwdn0p, lwdnp, lwupc0p, lwup0p, lwupp,  &
                     lwtoa0bp, lwtoabp , &  !FC
                     swdnc0p, swdn0p, swdnp, swupc0p, swup0p, swupp, &
                     topswad_aerop, solswad_aerop, &
                     topswai_aerop, solswai_aerop, &
                     topswad0_aerop, solswad0_aerop, &
                     topsw_aerop, topsw0_aerop, &
                     solsw_aerop, solsw0_aerop, &
                     topswcf_aerop, solswcf_aerop, &
                                !-C. Kleinschmitt for LW diagnostics
                     toplwad_aerop, sollwad_aerop,&
                     toplwai_aerop, sollwai_aerop, &
                     toplwad0_aerop, sollwad0_aerop,&
                                !-end
                     ZLWFT0_i, ZFLDN0, ZFLUP0, &
                     ZSWFT0_i, ZFSDN0, ZFSUP0, &
                     ZFLUX_DIR, ZFLUX_DIR_CLEAR, ZFLUX_DIR_SUN, &
                     cloud_cover_sw)
          ENDIF !ok_4xCO2atm

! A.I aout 2023
! Effet 3D des nuages Ecrad
! a passer : nom du ficher namelist et cles ok_3Deffect
! a declarer comme iflag_rrtm et a lire dans physiq.def
#ifdef CPP_ECRAD
          IF (ok_3Deffect) then
!                print*,'ok_3Deffect = ',ok_3Deffect
                namelist_ecrad_file='namelist_ecrad_s2'
                CALL radlwsw &
                     (debut, dist, rmu0, fract,  &
                     paprs, pplay,zxtsol,SFRWL,albsol_dir, albsol_dif, &
                     t_seri,q_seri,wo, &
                     cldfrarad, cldemirad, cldtaurad, &
                     ok_ade.OR.flag_aerosol_strat.GT.0, ok_aie,  ok_volcan, flag_volc_surfstrat, &
                     flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
                     tau_aero, piz_aero, cg_aero, &
                     tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm, &
                     tau_aero_lw_rrtm, &
                     cldtaupi, m_allaer, &
                     zqsat, flwc, fiwc, &
                     ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
                     namelist_ecrad_file, &
! A modifier
                     heat_s2,heat0_s2,cool_s2,cool0_s2,albpla_s2, &
                     heat_volc,cool_volc, &
                     topsw_s2,toplw_s2,solsw_s2,solswfdiff_s2,sollw_s2, &
                     sollwdown_s2, &
                     topsw0_s2,toplw0_s2,solsw0_s2,sollw0_s2, &
                     lwdnc0_s2, lwdn0_s2, lwdn_s2, lwupc0_s2, lwup0_s2, lwup_s2,  &
                     lwtoa0b_s2, lwtoab_s2 , &  !FC
                     swdnc0_s2, swdn0_s2, swdn_s2, swupc0_s2, swup0_s2, swup_s2, &
                     topswad_aero_s2, solswad_aero_s2, &
                     topswai_aero_s2, solswai_aero_s2, &
                     topswad0_aero_s2, solswad0_aero_s2, &
                     topsw_aero_s2, topsw0_aero_s2, &
                     solsw_aero_s2, solsw0_aero_s2, &
                     topswcf_aero_s2, solswcf_aero_s2, &
                                !-C. Kleinschmitt for LW diagnostics
                     toplwad_aero_s2, sollwad_aero_s2,&
                     toplwai_aero_s2, sollwai_aero_s2, &
                     toplwad0_aero_s2, sollwad0_aero_s2,&
                                !-end
                     ZLWFT0_i, ZFLDN0, ZFLUP0, &
                     ZSWFT0_i, ZFSDN0, ZFSUP0, &
                     ZFLUX_DIR, ZFLUX_DIR_CLEAR, ZFLUX_DIR_SUN, &
                     cloud_cover_sw_s2)
          ENDIF ! ok_3Deffect
#endif

       ENDIF ! aerosol_couple
       itaprad = 0
       !
       !  If Iflag_radia >=2, reset pertubed variables
       !
       IF (iflag_radia .ge. 2) THEN
          zxtsol(:) = zsav_tsol (:)
       ENDIF
    ENDIF ! MOD(itaprad,radpas)
    itaprad = itaprad + 1

    IF (iflag_radia.eq.0) THEN
       IF (prt_level.ge.9) THEN
          PRINT *,'--------------------------------------------------'
          PRINT *,'>>>> ATTENTION rayonnement desactive pour ce cas'
          PRINT *,'>>>>           heat et cool mis a zero '
          PRINT *,'--------------------------------------------------'
       ENDIF
       heat=0.
       cool=0.
       sollw=0.   ! MPL 01032011
       solsw=0.
       radsol=0.
       swup=0.    ! MPL 27102011 pour les fichiers AMMA_profiles et AMMA_scalars
       swup0=0.
       lwup=0.
       lwup0=0.
       lwdn=0.
       lwdn0=0.
    ENDIF

    !
    ! Calculer radsol a l'exterieur de radlwsw
    ! pour prendre en compte le cycle diurne
    ! recode par Olivier Boucher en sept 2015
    !
    radsol=solsw*swradcorr+sollw

    IF (ok_4xCO2atm) THEN
       radsolp=solswp*swradcorr+sollwp
    ENDIF

    !
    ! Ajouter la tendance des rayonnements (tous les pas)
    ! avec une correction pour le cycle diurne dans le SW
    !

    DO k=1, klev
       d_t_swr(:,k)=swradcorr(:)*heat(:,k)*phys_tstep/RDAY
       d_t_sw0(:,k)=swradcorr(:)*heat0(:,k)*phys_tstep/RDAY
       d_t_lwr(:,k)=-cool(:,k)*phys_tstep/RDAY
       d_t_lw0(:,k)=-cool0(:,k)*phys_tstep/RDAY
    ENDDO

    CALL add_phys_tend(du0,dv0,d_t_swr,dq0,dql0,dqi0,dqbs0,paprs,'SW',abortphy,flag_inhib_tend,itap,0)
    CALL prt_enerbil('SW',itap)
    CALL add_phys_tend(du0,dv0,d_t_lwr,dq0,dql0,dqi0,dqbs0,paprs,'LW',abortphy,flag_inhib_tend,itap,0)
    CALL prt_enerbil('LW',itap)

    !
    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    ! Calculer l'hydrologie de la surface
    !
    !      CALL hydrol(dtime,pctsrf,rain_fall, snow_fall, zxevap,
    !     .            agesno, ftsol,fqsurf,fsnow, ruis)
    !

    !
    ! Calculer le bilan du sol et la derive de temperature (couplage)
    !
    DO i = 1, klon
       !         bils(i) = radsol(i) - sens(i) - evap(i)*RLVTT
       ! a la demande de JLD
       bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
    ENDDO
    !
    !moddeblott(jan95)
    ! Appeler le programme de parametrisation de l'orographie
    ! a l'echelle sous-maille:
    !

    ! calculation of nm_oro
    DO i=1,klon
          ! nm_oro is a proxy for the number of subgrid scale mountains
          ! -> condition on nm_oro can deactivate the lifting on tilted planar terrains
          !    such as ice sheets (work by V. Wiener)
          ! in such a case, the SSO scheme should activate only where nm_oro>0 i.e. by setting
          ! nm_oro_t=0.
          nm_oro(i)=zsig(i)*sqrt(cell_area(i)*(pctsrf(i,is_ter)+pctsrf(i,is_lic)))/(4.*MAX(zstd(i),1.e-8))-1.
    ENDDO

    IF (prt_level .GE.10) THEN
       print *,' call orography ? ', ok_orodr
    ENDIF
    !
    IF (ok_orodr) THEN
       !
       !  selection des points pour lesquels le shema est actif:
       igwd=0
       DO i=1,klon
          itest(i)=0
          ! zpmm_orodr_t and zstd_orodr_t are activation thresholds set by F. Lott to
          ! earn computation time but they are not physical.
          IF (((zpic(i)-zmea(i)).GT.zpmm_orodr_t).AND.(zstd(i).GT.zstd_orodr_t).AND.(nm_oro(i).GT.nm_oro_t)) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
          ENDIF
       ENDDO
       !        igwdim=MAX(1,igwd)
       !
       IF (ok_strato) THEN

          CALL drag_noro_strato(0,klon,klev,phys_tstep,paprs,pplay, &
               zmea,zstd, zsig, zgam, zthe,zpic,zval, &
               igwd,idx,itest, &
               t_seri, u_seri, v_seri, &
               zulow, zvlow, zustrdr, zvstrdr, &
               d_t_oro, d_u_oro, d_v_oro)

       ELSE
          CALL drag_noro(klon,klev,phys_tstep,paprs,pplay, &
               zmea,zstd, zsig, zgam, zthe,zpic,zval, &
               igwd,idx,itest, &
               t_seri, u_seri, v_seri, &
               zulow, zvlow, zustrdr, zvstrdr, &
               d_t_oro, d_u_oro, d_v_oro)
       ENDIF
       !
       !  ajout des tendances
       !-----------------------------------------------------------------------
       ! ajout des tendances de la trainee de l'orographie
       CALL add_phys_tend(d_u_oro,d_v_oro,d_t_oro,dq0,dql0,dqi0,dqbs0,paprs,'oro', &
            abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('oro',itap)
       !----------------------------------------------------------------------
       !
    ENDIF ! fin de test sur ok_orodr
    !
    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    IF (ok_orolf) THEN
       !
       !  selection des points pour lesquels le shema est actif:
       igwd=0
       DO i=1,klon
          itest(i)=0
          IF (((zpic(i)-zmea(i)).GT.zpmm_orolf_t).AND.(nm_oro(i).GT.nm_oro_t)) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
          ENDIF
       ENDDO
       !        igwdim=MAX(1,igwd)
       !
       IF (ok_strato) THEN

          CALL lift_noro_strato(klon,klev,phys_tstep,paprs,pplay, &
               latitude_deg,zmea,zstd,zpic,zgam,zthe,zpic,zval, &
               igwd,idx,itest, &
               t_seri, u_seri, v_seri, &
               zulow, zvlow, zustrli, zvstrli, &
               d_t_lif, d_u_lif, d_v_lif               )

       ELSE
          CALL lift_noro(klon,klev,phys_tstep,paprs,pplay, &
               latitude_deg,zmea,zstd,zpic, &
               itest, &
               t_seri, u_seri, v_seri, &
               zulow, zvlow, zustrli, zvstrli, &
               d_t_lif, d_u_lif, d_v_lif)
       ENDIF

       ! ajout des tendances de la portance de l'orographie
       CALL add_phys_tend(d_u_lif, d_v_lif, d_t_lif, dq0, dql0, dqi0, dqbs0,paprs, &
            'lif', abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('lif',itap)
    ENDIF ! fin de test sur ok_orolf

    IF (ok_hines) then
       !  HINES GWD PARAMETRIZATION
       east_gwstress=0.
       west_gwstress=0.
       du_gwd_hines=0.
       dv_gwd_hines=0.
       CALL hines_gwd(klon, klev, phys_tstep, paprs, pplay, latitude_deg, t_seri, &
            u_seri, v_seri, zustr_gwd_hines, zvstr_gwd_hines, d_t_hin, &
            du_gwd_hines, dv_gwd_hines)
       zustr_gwd_hines=0.
       zvstr_gwd_hines=0.
       DO k = 1, klev
          zustr_gwd_hines(:)=zustr_gwd_hines(:)+ du_gwd_hines(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
          zvstr_gwd_hines(:)=zvstr_gwd_hines(:)+ dv_gwd_hines(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
       ENDDO

       d_t_hin(:, :)=0.
       CALL add_phys_tend(du_gwd_hines, dv_gwd_hines, d_t_hin, dq0, dql0, &
            dqi0, dqbs0, paprs, 'hin', abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('hin',itap)
    ENDIF

    IF (.not. ok_hines .and. ok_gwd_rando) then
       ! ym missing init for east_gwstress & west_gwstress -> added in phys_local_var_mod
       CALL acama_GWD_rando_first()
       CALL acama_GWD_rando(PHYS_TSTEP, pplay, latitude_deg, t_seri, u_seri, &
            v_seri, rot, zustr_gwd_front, zvstr_gwd_front, du_gwd_front, &
            dv_gwd_front, east_gwstress, west_gwstress)
       zustr_gwd_front=0.
       zvstr_gwd_front=0.
       DO k = 1, klev
          zustr_gwd_front(:)=zustr_gwd_front(:)+ du_gwd_front(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
          zvstr_gwd_front(:)=zvstr_gwd_front(:)+ dv_gwd_front(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
       ENDDO

       CALL add_phys_tend(du_gwd_front, dv_gwd_front, dt0, dq0, dql0, dqi0, dqbs0, &
            paprs, 'front_gwd_rando', abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('front_gwd_rando',itap)
    ENDIF

    IF (ok_gwd_rando) THEN
       CALL FLOTT_GWD_rando_first()
       CALL FLOTT_GWD_rando(PHYS_TSTEP, pplay, t_seri, u_seri, v_seri, &
            rain_fall + snow_fall, zustr_gwd_rando, zvstr_gwd_rando, &
            du_gwd_rando, dv_gwd_rando, east_gwstress, west_gwstress)
       CALL add_phys_tend(du_gwd_rando, dv_gwd_rando, dt0, dq0, dql0, dqi0, dqbs0, &
            paprs, 'flott_gwd_rando', abortphy,flag_inhib_tend,itap,0)
       CALL prt_enerbil('flott_gwd_rando',itap)
       zustr_gwd_rando=0.
       zvstr_gwd_rando=0.
       DO k = 1, klev
          zustr_gwd_rando(:)=zustr_gwd_rando(:)+ du_gwd_rando(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
          zvstr_gwd_rando(:)=zvstr_gwd_rando(:)+ dv_gwd_rando(:, k)/phys_tstep &
               * (paprs(:, k)-paprs(:, k+1))/rg
       ENDDO
    ENDIF

    ! STRESS NECESSAIRES: TOUTE LA PHYSIQUE

    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    DO i = 1, klon
       zustrph(i)=0.
       zvstrph(i)=0.
    ENDDO
    DO k = 1, klev
       DO i = 1, klon
          zustrph(i)=zustrph(i)+(u_seri(i,k)-u(i,k))/phys_tstep* &
               (paprs(i,k)-paprs(i,k+1))/rg
          zvstrph(i)=zvstrph(i)+(v_seri(i,k)-v(i,k))/phys_tstep* &
               (paprs(i,k)-paprs(i,k+1))/rg
       ENDDO
    ENDDO
    !
    !IM calcul composantes axiales du moment angulaire et couple des montagnes
    !
    IF (is_sequential .and. ok_orodr) THEN
       CALL aaam_bud (27,klon,klev,jD_cur-jD_ref,jH_cur, &
            ra,rg,romega, &
            latitude_deg,longitude_deg,pphis, &
            zustrdr,zustrli,zustrph, &
            zvstrdr,zvstrli,zvstrph, &
            paprs,u,v, &
            aam, torsfc)
    ENDIF
    !IM cf. FLott END
    !DC Calcul de la tendance due au methane
    IF (ok_qch4) THEN
!      d_q_ch4: H2O source from CH4 in MMR/s (mass mixing ratio/s or kg H2O/kg air/s)
    IF (CPPKEY_STRATAER) THEN

       CALL stratH2O_methox(debut,paprs,d_q_ch4)
    ELSE
!      ECMWF routine METHOX
       CALL METHOX(1,klon,klon,klev,q_seri,d_q_ch4,pplay)
    END IF
       ! add humidity tendency due to methane
       d_q_ch4_dtime(:,:) = d_q_ch4(:,:)*phys_tstep
       CALL add_phys_tend(du0, dv0, dt0, d_q_ch4_dtime, dql0, dqi0, dqbs0, paprs, &
            'q_ch4', abortphy,flag_inhib_tend,itap,0)
       d_q_ch4(:,:) = d_q_ch4_dtime(:,:)/phys_tstep ! update with H2O conserv done in add_phys_tend
    ENDIF
    !
    !
IF (CPPKEY_STRATAER) THEN
    IF (ok_qemiss) THEN
       flh2o=1
       IF(flag_verbose_strataer) THEN
          print *,'IN physiq_mod: ok_qemiss =yes (',ok_qemiss,'), flh2o=',flh2o
          print *,'IN physiq_mod: flag_emit=',flag_emit,', nErupt=',nErupt
          print *,'IN physiq_mod: nAerErupt=',nAerErupt
       ENDIF

       SELECT CASE(flag_emit)
       CASE(1) ! emission volc H2O in LMDZ
          DO ieru=1, nErupt
             IF (year_cur==year_emit_vol(ieru).AND.&
                  mth_cur==mth_emit_vol(ieru).AND.&
                  day_cur>=day_emit_vol(ieru).AND.&
                  day_cur<(day_emit_vol(ieru)+injdur)) THEN

                IF(flag_verbose_strataer) print *,'IN physiq_mod: date=',year_cur,mth_cur,day_cur
                ! initialisation of q tendency emission
                d_q_emiss(:,:)=0.
                ! daily injection mass emission - NL
                m_H2O_emiss_vol_daily = m_H2O_emiss_vol(ieru)/(REAL(injdur)&
                     *REAL(ponde_lonlat_vol(ieru)))
                !
                CALL STRATEMIT(pdtphys,pdtphys,latitude_deg,longitude_deg,t_seri,&
                    pplay,paprs,tr_seri,&
                    m_H2O_emiss_vol_daily,&
                    xlat_min_vol(ieru),xlat_max_vol(ieru),&
                    xlon_min_vol(ieru),xlon_max_vol(ieru),&
                    altemiss_vol(ieru),sigma_alt_vol(ieru),1,1.,&
                    nAerErupt+1,0)

                IF(flag_verbose_strataer) print *,'IN physiq_mod: min max d_q_emiss=',&
                     minval(d_q_emiss),maxval(d_q_emiss)

                CALL add_phys_tend(du0, dv0, dt0, d_q_emiss, dql0, dqi0, dqbs0, paprs, &
                     'q_emiss',abortphy,flag_inhib_tend,itap,0)
                IF (abortphy==1) Print*,'ERROR ABORT TEND EMISS'
             ENDIF
          ENDDO
          flh2o=0
       END SELECT ! emission scenario (flag_emit)
    ENDIF
END IF

!===============================================================
!            Additional tendency of TKE due to orography
!===============================================================

       dtadd(:,:)=0.
       duadd(:,:)=0.
       dvadd(:,:)=0.

! Choices for addtkeoro:
!      ** 0 no TKE tendency from orography
!      ** 1 we include a fraction alphatkeoro of the whole tendency duoro
!      ** 2 we include a fraction alphatkeoro of the gravity wave part of duoro
!

       IF (addtkeoro .GT. 0 .AND. ok_orodr ) THEN
!      -------------------------------------------


       !  selection des points pour lesquels le schema est actif:


  IF (addtkeoro .EQ. 1 ) THEN

            duadd(:,:)=alphatkeoro*d_u_oro(:,:)
            dvadd(:,:)=alphatkeoro*d_v_oro(:,:)

  ELSE IF (addtkeoro .EQ. 2) THEN

     IF (smallscales_tkeoro) THEN
       igwd=0
       DO i=1,klon
          itest(i)=0
! Etienne: ici je prends en compte plus de relief que la routine drag_noro_strato
! car on peut s'attendre a ce que les petites echelles produisent aussi de la TKE
! Mais attention, cela ne va pas dans le sens de la conservation de l'energie!
          IF ((zstd(i).GT.1.0) .AND.(nm_oro(i).GT.nm_oro_t)) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
          ENDIF
       ENDDO

     ELSE

       igwd=0
       DO i=1,klon
          itest(i)=0
        IF (((zpic(i)-zmea(i)).GT.zpmm_orodr_t).AND.(zstd(i).GT.zstd_orodr_t).AND.(nm_oro(i).GT.nm_oro_t)) THEN
             itest(i)=1
             igwd=igwd+1
             idx(igwd)=i
        ENDIF
       ENDDO

     ENDIF

     CALL drag_noro_strato(addtkeoro,klon,klev,phys_tstep,paprs,pplay, &
               zmea,zstd, zsig, zgam, zthe,zpic,zval, &
               igwd,idx,itest, &
               t_seri, u_seri, v_seri, &
               zulow, zvlow, zustrdr, zvstrdr, &
               d_t_oro_gw, d_u_oro_gw, d_v_oro_gw)

     zustrdr(:)=0.
     zvstrdr(:)=0.
     zulow(:)=0.
     zvlow(:)=0.

     duadd(:,:)=alphatkeoro*d_u_oro_gw(:,:)
     dvadd(:,:)=alphatkeoro*d_v_oro_gw(:,:)
  ENDIF


   ! TKE update from subgrid temperature and wind tendencies
   !----------------------------------------------------------
    forall (k=1: nbp_lev) exner(:, k) = (pplay(:, k)/paprs(:,1))**RKAPPA


    CALL tend_to_tke(pdtphys,paprs,exner,t_seri,u_seri,v_seri,dtadd,duadd,dvadd,pctsrf,pbl_tke)
   !
   ! Prevent pbl_tke_w from becoming negative
    wake_delta_pbl_tke(:,:,:) = max(wake_delta_pbl_tke(:,:,:), -pbl_tke(:,:,:))
   !

       ENDIF
!      -----
!===============================================================


    !====================================================================
    ! Interface Simulateur COSP (Calipso, ISCCP, MISR, ..)
    !====================================================================
    ! Abderrahmane 24.08.09

    IF (ok_cosp) THEN
       ! adeclarer
IF (CPPKEY_COSP) THEN
       IF (itap.eq.1.or.MOD(itap,NINT(freq_cosp/phys_tstep)).EQ.0) THEN

          IF (prt_level .GE.10) THEN
             print*,'freq_cosp',freq_cosp
          ENDIF
          mr_ozone=wo(:, :, 1) * dobson_u * 1e3 / zmasse
          !       print*,'Dans physiq.F avant appel cosp ref_liq,ref_ice=',
          !     s        ref_liq,ref_ice
          CALL phys_cosp(itap,phys_tstep,freq_cosp, &
               ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
               ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
               klon,klev,longitude_deg,latitude_deg,presnivs,overlap, &
               JrNt,ref_liq,ref_ice, &
               pctsrf(:,is_ter)+pctsrf(:,is_lic), &
               zu10m,zv10m,pphis, &
               zphi,paprs(:,1:klev),pplay,zxtsol,t_seri, &
               qx(:,:,ivap),zx_rh,cldfra,rnebcon,flwc,fiwc, &
               prfl(:,1:klev),psfl(:,1:klev), &
               pmflxr(:,1:klev),pmflxs(:,1:klev), &
               mr_ozone,cldtau, cldemi)

          !     L         calipso2D,calipso3D,cfadlidar,parasolrefl,atb,betamol,
          !     L          cfaddbze,clcalipso2,dbze,cltlidarradar,
          !     M          clMISR,
          !     R          clisccp2,boxtauisccp,boxptopisccp,tclisccp,ctpisccp,
          !     I          tauisccp,albisccp,meantbisccp,meantbclrisccp)

       ENDIF
END IF

IF (CPPKEY_COSP2) THEN
       IF (itap.eq.1.or.MOD(itap,NINT(freq_cosp/phys_tstep)).EQ.0) THEN

          IF (prt_level .GE.10) THEN
             print*,'freq_cosp',freq_cosp
          ENDIF
          mr_ozone=wo(:, :, 1) * dobson_u * 1e3 / zmasse
                 print*,'Dans physiq.F avant appel '
          !     s        ref_liq,ref_ice
          CALL phys_cosp2(itap,phys_tstep,freq_cosp, &
               ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
               ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
               klon,klev,longitude_deg,latitude_deg,presnivs,overlap, &
               JrNt,ref_liq,ref_ice, &
               pctsrf(:,is_ter)+pctsrf(:,is_lic), &
               zu10m,zv10m,pphis, &
               zphi,paprs(:,1:klev),pplay,zxtsol,t_seri, &
               qx(:,:,ivap),zx_rh,cldfra,rnebcon,flwc,fiwc, &
               prfl(:,1:klev),psfl(:,1:klev), &
               pmflxr(:,1:klev),pmflxs(:,1:klev), &
               mr_ozone,cldtau, cldemi)
       ENDIF
END IF

IF (CPPKEY_COSPV2) THEN
       IF (itap.eq.1.or.MOD(itap,NINT(freq_cosp/phys_tstep)).EQ.0) THEN
!        IF (MOD(itap,NINT(freq_cosp/phys_tstep)).EQ.0) THEN

          IF (prt_level .GE.10) THEN
             print*,'freq_cosp',freq_cosp
          ENDIF
           DO k = 1, klev
             DO i = 1, klon
               phicosp(i,k) = pphi(i,k) + pphis(i)
             ENDDO
           ENDDO
          mr_ozone=wo(:, :, 1) * dobson_u * 1e3 / zmasse
                 print*,'Dans physiq.F avant appel '
          !     s        ref_liq,ref_ice
          CALL lmdz_cosp_interface(itap,phys_tstep,freq_cosp, &
               ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
               ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, missing_val, &
               klon,klev,longitude_deg,latitude_deg,presnivs,overlap, &
               JrNt,ref_liq,ref_ice, &
               pctsrf(:,is_ter)+pctsrf(:,is_lic), &
               zu10m,zv10m,pphis, &
               zphi,paprs(:,1:klev),pplay,zxtsol,t_seri, &
               qx(:,:,ivap),zx_rh,cldfra,rnebcon,flwc,fiwc, &
               prfl(:,1:klev),psfl(:,1:klev), &
               pmflxr(:,1:klev),pmflxs(:,1:klev), &
               mr_ozone,cldtau, cldemi)
       ENDIF
END IF

    ENDIF  !ok_cosp


! Marine

  IF (ok_airs) then

  IF (itap.eq.1.or.MOD(itap,NINT(freq_airs/phys_tstep)).EQ.0) THEN
     write(*,*) 'je vais appeler simu_airs, ok_airs, freq_airs=', ok_airs, freq_airs
     CALL simu_airs(itap,rneb, t_seri, cldemi, fiwc, ref_ice, pphi, pplay, paprs,&
        & map_prop_hc,map_prop_hist,&
        & map_emis_hc,map_iwp_hc,map_deltaz_hc,map_pcld_hc,map_tcld_hc,&
        & map_emis_Cb,map_pcld_Cb,map_tcld_Cb,&
        & map_emis_ThCi,map_pcld_ThCi,map_tcld_ThCi,&
        & map_emis_Anv,map_pcld_Anv,map_tcld_Anv,&
        & map_emis_hist,map_iwp_hist,map_deltaz_hist,map_rad_hist,&
        & map_ntot,map_hc,map_hist,&
        & map_Cb,map_ThCi,map_Anv,&
        & alt_tropo )
  ENDIF

  ENDIF  ! ok_airs


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !AA
    !AA Installation de l'interface online-offline pour traceurs
    !AA
    !====================================================================
    !   Calcul  des tendances traceurs
    !====================================================================
    !

    IF (type_trac == 'repr') THEN
!MM pas d'impact, car on recupere q_seri,tr_seri,t_seri via phys_local_var_mod
!MM                               dans Reprobus
       sh_in(:,:) = q_seri(:,:)
IF (CPPKEY_REPROBUS) THEN
       d_q_rep(:,:) = 0.
       d_ql_rep(:,:) = 0.
       d_qi_rep(:,:) = 0.
END IF
    ELSE
       sh_in(:,:) = qx(:,:,ivap)
       IF (nqo >= 3) THEN
          ch_in(:,:) = qx(:,:,iliq) + qx(:,:,isol)
       ELSE
          ch_in(:,:) = qx(:,:,iliq)
       ENDIF
    ENDIF

    ! Merge wdtrainA and wdtrainS in the total source of precipitation due to
    ! adiabatic updraughts.
    !
    wdtrainAS(:,:) = wdtrainA(:,:) + wdtrainS(:,:)

IF (CPPKEY_DUST) THEN
    ! Avec SPLA, iflag_phytrac est forcé =1

    CALL       phytracr_spl ( debut,lafin , jD_cur,jH_cur,iflag_con,       &
                      pdtphys,ftsol,                                       &
                      t,q_seri,paprs,pplay,RHcl,                           &
                      pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,              &
                      coefh(1:klon,1:klev,is_ave), cdragh, cdragm, u1, v1, &  
                      u_seri, v_seri, latitude_deg, longitude_deg,         &
                      pphis,pctsrf,pmflxr,pmflxs,prfl,psfl,                &
                      da,phi,phi2,d1a,dam,mp,ep,sigd,sij,clw,elij,         &
                      epmlmMm,eplaMm,upwd,dnwd,itop_con,ibas_con,          &
                      ev,wdtrainAS,  wdtrainM,wght_cvfd,                    &
                      fm_therm, entr_therm, rneb,                          &
                      beta_prec_fisrt,beta_prec,                           &
                      zu10m,zv10m,wstar,ale_bl,ale_wake,                   &
                      nsurfwind,surf_wind_value, surf_wind_proba,          &
                      d_tr_dyn,tr_seri)

ELSE
    IF (iflag_phytrac == 1 ) THEN
      CALL phytrac ( &
         itap,     days_elapsed+1,    jH_cur,   debut, &
         lafin,    phys_tstep,     u, v,     t, &
         paprs,    pplay,     pmfu,     pmfd, &
         pen_u,    pde_u,     pen_d,    pde_d, &
         cdragh,   coefh(1:klon,1:klev,is_ave),   fm_therm, entr_therm, &
         u1,       v1,        ftsol,    pctsrf, &
         zustar,   zu10m,     zv10m, &
         wstar(:,is_ave),    ale_bl,         ale_wake, &
         latitude_deg, longitude_deg, &
         frac_impa,frac_nucl, beta_prec_fisrt,beta_prec, &
         presnivs, pphis,     pphi,     albsol1, &
         sh_in,   ch_in,    rhcl,      cldfra,   rneb, &
         diafra,   radocond,    itop_con, ibas_con, &
         pmflxr,   pmflxs,    prfl,     psfl, &
         da,       phi,       mp,       upwd, &
         phi2,     d1a,       dam,      sij, wght_cvfd, &        !<<RomP+RL
         wdtrainAS, wdtrainM,  sigd,     clw,elij, &   !<<RomP
         ev,       ep,        epmlmMm,  eplaMm, &     !<<RomP
         dnwd,     aerosol_couple,      flxmass_w, &
         tau_aero, piz_aero,  cg_aero,  ccm, &
         rfname, &
         d_tr_dyn, &                                 !<<RomP
         tr_seri, init_source)
IF (CPPKEY_REPROBUS) THEN


          print*,'avt add phys rep',abortphy

     CALL add_phys_tend &
            (du0,dv0,dt0,d_q_rep,d_ql_rep,d_qi_rep,dqbs0,paprs,&
             'rep',abortphy,flag_inhib_tend,itap,0)
        IF (abortphy==1) Print*,'ERROR ABORT REP'

          print*,'apr add phys rep',abortphy

END IF
    ENDIF    ! (iflag_phytrac=1)

END IF
    !ENDIF    ! (iflag_phytrac=1)

    IF (offline) THEN

       IF (prt_level.ge.9) &
            print*,'Attention on met a 0 les thermiques pour phystoke'
       CALL phystokenc ( &
            nlon,klev,pdtphys,longitude_deg,latitude_deg, &
            t,pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
            fm_therm,entr_therm, &
            cdragh,coefh(1:klon,1:klev,is_ave),u1,v1,ftsol,pctsrf, &
            frac_impa, frac_nucl, &
            pphis,cell_area,phys_tstep,itap, &
            qx(:,:,ivap),da,phi,mp,upwd,dnwd)


    ENDIF

    !
    ! Calculer le transport de l'eau et de l'energie (diagnostique)
    !
    CALL transp (paprs,zxtsol, t_seri, q_seri, ql_seri, qs_seri, u_seri, v_seri, zphi, &
                 ue, ve, uq, vq, uwat, vwat)
    !
    !IM global posePB BEG
    IF(1.EQ.0) THEN
       !
       CALL transp_lay (paprs,zxtsol, t_seri, q_seri, u_seri, v_seri, zphi, &
            ve_lay, vq_lay, ue_lay, uq_lay)
       !
    ENDIF !(1.EQ.0) THEN
    !IM global posePB END
    !
    ! Accumuler les variables a stocker dans les fichiers histoire:
    !

    !================================================================
    ! Conversion of kinetic and potential energy into heat, for
    ! parameterisation of subgrid-scale motions
    !================================================================

    d_t_ec(:,:)=0.
    forall (k=1: nbp_lev) exner(:, k) = (pplay(:, k)/paprs(:,1))**RKAPPA
    CALL ener_conserv(klon,klev,pdtphys,u,v,t,qx,ivap,iliq,isol, &
         u_seri,v_seri,t_seri,q_seri,ql_seri,qs_seri,pbl_tke(:,:,is_ave)-tke0(:,:), &
         zmasse,exner,d_t_ec)
    t_seri(:,:)=t_seri(:,:)+d_t_ec(:,:)

    !==================================================================
    !--OB water mass fixer for the physics
    !--water profiles are corrected to force mass conservation of water
    !--currently flag is turned off
    !==================================================================
    IF (ok_water_mass_fixer) THEN
    qql2(:)=0.0
    DO k = 1, klev
      qql2(:)=qql2(:)+(q_seri(:,k)+ql_seri(:,k))*zmasse(:,k)
      IF (nqo >= 3) THEN
        qql2(:)=qql2(:)+qs_seri(:,k)*zmasse(:,k)
      ENDIF
      IF (ok_bs) THEN
        qql2(:)=qql2(:)+qbs_seri(:,k)*zmasse(:,k)
      ENDIF
    ENDDO

IF (CPPKEY_STRATAER) THEN
    IF (ok_qemiss) THEN
       DO k = 1, klev
          qql1(:) = qql1(:)+d_q_emiss(:,k)*zmasse(:,k)
       ENDDO
    ENDIF
END IF
    IF (ok_qch4) THEN
       DO k = 1, klev
          qql1(:) = qql1(:)+d_q_ch4_dtime(:,k)*zmasse(:,k)
       ENDDO
    ENDIF

    DO i = 1, klon
      !--compute ratio of what q+ql should be with conservation to what it is
      IF (ok_bs) THEN
        corrqql=(qql1(i)+(evap(i)-snowerosion(i)-rain_fall(i)-snow_fall(i)-bs_fall(i))*pdtphys)/qql2(i)
      ELSE
        corrqql=(qql1(i)+(evap(i)-rain_fall(i)-snow_fall(i))*pdtphys)/qql2(i)
      ENDIF
      DO k = 1, klev
        q_seri(i,k) =q_seri(i,k)*corrqql
        ql_seri(i,k)=ql_seri(i,k)*corrqql
        IF (nqo >= 3) THEN
          qs_seri(i,k)=qs_seri(i,k)*corrqql
        ENDIF
        IF (ok_bs) THEN
          qbs_seri(i,k)=qbs_seri(i,k)*corrqql
        ENDIF
      ENDDO
    ENDDO
    ENDIF
    !--fin mass fixer

    !cc prw  = eau precipitable
    !   prlw = colonne eau liquide
    !   prlw = colonne eau solide
    !   prbsw = colonne neige soufflee
    !   water_budget = non-conservation residual from the LMDZ physics
    !                  (should be equal to machine precision if mass fixer is activated)
    prw(:) = 0.
    prlw(:) = 0.
    prsw(:) = 0.
    prbsw(:) = 0.
    water_budget(:) = 0.0
    DO k = 1, klev
       prw(:)  = prw(:)  + q_seri(:,k)*zmasse(:,k)
       prlw(:) = prlw(:) + ql_seri(:,k)*zmasse(:,k)
       water_budget(:) = water_budget(:) + (q_seri(:,k)-qx(:,k,ivap)+ql_seri(:,k)-qx(:,k,iliq))*zmasse(:,k)
       IF (nqo >= 3) THEN
         prsw(:) = prsw(:) + qs_seri(:,k)*zmasse(:,k)
         water_budget(:) = water_budget(:) + (qs_seri(:,k)-qx(:,k,isol))*zmasse(:,k)
       ENDIF
       IF (nqo >= 4 .AND. ok_bs) THEN
         prbsw(:)= prbsw(:) + qbs_seri(:,k)*zmasse(:,k)
         water_budget(:) = water_budget(:) + (qbs_seri(:,k)-qx(:,k,ibs))*zmasse(:,k)
       ENDIF
    ENDDO
    water_budget(:)=water_budget(:)+(rain_fall(:)+snow_fall(:)-evap(:))*pdtphys
    IF (ok_bs) THEN
      water_budget(:)=water_budget(:)+bs_fall(:)*pdtphys
    ENDIF

    !=======================================================================
    !   SORTIES
    !=======================================================================
    !
    !IM initialisation + calculs divers diag AMIP2
    CALL calcul_divers(itap, itapm1, un_jour)
    !
    !IM Interpolation sur les niveaux de pression du NMC
    !   -------------------------------------------------
    !
    include "calcul_STDlev.h"
    !
    ! slp sea level pressure derived from Arpege-IFS : CALL ctstar + CALL pppmer
    CALL diag_slp(klon,t_seri,paprs,pplay,pphis,ptstar,pt0,slp)
    !
    !
    IF (ANY(type_trac == ['inca','inco'])) THEN
IF (CPPKEY_INCA) THEN
       CALL VTe(VTphysiq)
       CALL VTb(VTinca)

       CALL chemhook_end ( &
            phys_tstep, &
            pplay, &
            t_seri, &
            tr_seri(:,:,1+nqCO2:nbtr), &
            nbtr, &
            paprs, &
            q_seri, &
            cell_area, &
            pphi, &
            pphis, &
            zx_rh, &
            aps, bps, ap, bp, lafin)

       CALL VTe(VTinca)
       CALL VTb(VTphysiq)
END IF
    ENDIF

    IF (type_trac == 'repr') THEN
IF (CPPKEY_REPROBUS) THEN
        CALL coord_hyb_rep(paprs, pplay, aps, bps, ap, bp, cell_area)
END IF
    ENDIF

    !
    ! Convertir les incrementations en tendances
    !
    IF (prt_level .GE.10) THEN
       print *,'Convertir les incrementations en tendances '
    ENDIF
    !
    IF (mydebug) THEN
       CALL writefield_phy('u_seri',u_seri,nbp_lev)
       CALL writefield_phy('v_seri',v_seri,nbp_lev)
       CALL writefield_phy('t_seri',t_seri,nbp_lev)
       CALL writefield_phy('q_seri',q_seri,nbp_lev)
    ENDIF

    DO k = 1, klev
       DO i = 1, klon
          d_u(i,k) = ( u_seri(i,k) - u(i,k) ) / phys_tstep
          d_v(i,k) = ( v_seri(i,k) - v(i,k) ) / phys_tstep
          d_t(i,k) = ( t_seri(i,k)-t(i,k) ) / phys_tstep
          d_qx(i,k,ivap) = ( q_seri(i,k) - qx(i,k,ivap) ) / phys_tstep
          d_qx(i,k,iliq) = ( ql_seri(i,k) - qx(i,k,iliq) ) / phys_tstep
          !CR: on ajoute le contenu en glace
          IF (nqo >= 3) THEN
             d_qx(i,k,isol) = ( qs_seri(i,k) - qx(i,k,isol) ) / phys_tstep
          ENDIF
          !--ice_supersat: nqo=5, we add cloud fraction and cloudy water vapor to total water vapor ratio
          IF (nqo.ge.5 .and. ok_ice_supersat) THEN
             d_qx(i,k,icf) = ( cf_seri(i,k) - qx(i,k,icf) ) / phys_tstep
             d_qx(i,k,irvc) = ( rvc_seri(i,k) - qx(i,k,irvc) ) / phys_tstep
          ENDIF

           IF (nqo.ge.4 .and. ok_bs) THEN
             d_qx(i,k,ibs) = ( qbs_seri(i,k) - qx(i,k,ibs) ) / phys_tstep
          ENDIF

       ENDDO
    ENDDO

    ! in case of advection of TKE
    IF (ok_advtke) THEN
      DO k=1,klev
         DO i=1,klon
           d_qx(i,k,itke)=((pbl_tke(i,k,is_ave)+pbl_tke(i,k+1,is_ave))/2. - qx(i,k,itke)) / phys_tstep
         ENDDO
      ENDDO
    ENDIF
    !
    ! DC: All iterations are cycled if nqtot==nqo, so no nqtot>nqo condition required
    itr = 0
    DO iq = 1, nqtot
       IF(.NOT.tracers(iq)%isInPhysics) CYCLE
       itr = itr+1
       DO  k = 1, klev
          DO  i = 1, klon
             d_qx(i,k,iq) = ( tr_seri(i,k,itr) - qx(i,k,iq) ) / phys_tstep
          ENDDO
       ENDDO
    ENDDO
    !
    !IM rajout diagnostiques bilan KP pour analyse MJO par Jun-Ichi Yano
    !IM global posePB      include "write_bilKP_ins.h"
    !IM global posePB      include "write_bilKP_ave.h"
    !
    ! Sauvegarder les valeurs de t et q a la fin de la physique:
    !
    u_ancien(:,:)  = u_seri(:,:)
    v_ancien(:,:)  = v_seri(:,:)
    t_ancien(:,:)  = t_seri(:,:)
    q_ancien(:,:)  = q_seri(:,:)
    ql_ancien(:,:) = ql_seri(:,:)
    qs_ancien(:,:) = qs_seri(:,:)
    qbs_ancien(:,:)= qbs_seri(:,:)
    cf_ancien(:,:) = cf_seri(:,:)
    rvc_ancien(:,:)= rvc_seri(:,:)
    tke_ancien(:,:)= pbl_tke(:,:,is_ave)

    CALL water_int(klon,klev,q_ancien,zmasse,prw_ancien)
    CALL water_int(klon,klev,ql_ancien,zmasse,prlw_ancien)
    CALL water_int(klon,klev,qs_ancien,zmasse,prsw_ancien)
    CALL water_int(klon,klev,qbs_ancien,zmasse,prbsw_ancien)
    ! !! RomP >>>
    IF (nqtot > nqo) tr_ancien(:,:,:) = tr_seri(:,:,:)
    ! !! RomP <<<
    !==========================================================================
    ! Sorties des tendances pour un point particulier
    ! a utiliser en 1D, avec igout=1 ou en 3D sur un point particulier
    ! pour le debug
    ! La valeur de igout est attribuee plus haut dans le programme
    !==========================================================================

    IF (prt_level.ge.1) THEN
       write(lunout,*) 'FIN DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
       write(lunout,*) &
            'nlon,klev,nqtot,debut,lafin,jD_cur, jH_cur, pdtphys pct tlos'
       write(lunout,*) &
            nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur ,pdtphys, &
            pctsrf(igout,is_ter), pctsrf(igout,is_lic),pctsrf(igout,is_oce), &
            pctsrf(igout,is_sic)
       write(lunout,*) 'd_t_dyn,d_t_con,d_t_lsc,d_t_ajsb,d_t_ajs,d_t_eva'
       DO k=1,klev
          write(lunout,*) d_t_dyn(igout,k),d_t_con(igout,k), &
               d_t_lsc(igout,k),d_t_ajsb(igout,k),d_t_ajs(igout,k), &
               d_t_eva(igout,k)
       ENDDO
       write(lunout,*) 'cool,heat'
       DO k=1,klev
          write(lunout,*) cool(igout,k),heat(igout,k)
       ENDDO

       !jyg<     (En attendant de statuer sur le sort de d_t_oli)
       !jyg!     write(lunout,*) 'd_t_oli,d_t_vdf,d_t_oro,d_t_lif,d_t_ec'
       !jyg!     do k=1,klev
       !jyg!        write(lunout,*) d_t_oli(igout,k),d_t_vdf(igout,k), &
       !jyg!             d_t_oro(igout,k),d_t_lif(igout,k),d_t_ec(igout,k)
       !jyg!     enddo
       write(lunout,*) 'd_t_vdf,d_t_oro,d_t_lif,d_t_ec'
       DO k=1,klev
          write(lunout,*) d_t_vdf(igout,k), &
               d_t_oro(igout,k),d_t_lif(igout,k),d_t_ec(igout,k)
       ENDDO
       !>jyg

       write(lunout,*) 'd_ps ',d_ps(igout)
       write(lunout,*) 'd_u, d_v, d_t, d_qx1, d_qx2 '
       DO k=1,klev
          write(lunout,*) d_u(igout,k),d_v(igout,k),d_t(igout,k), &
               d_qx(igout,k,1),d_qx(igout,k,2)
       ENDDO
    ENDIF

    !============================================================
    !   Calcul de la temperature potentielle
    !============================================================
    DO k = 1, klev
       DO i = 1, klon
          !JYG/IM theta en debut du pas de temps
          !JYG/IM       theta(i,k)=t(i,k)*(100000./pplay(i,k))**(RD/RCPD)
          !JYG/IM theta en fin de pas de temps de physique
          theta(i,k)=t_seri(i,k)*(100000./pplay(i,k))**(RD/RCPD)
          ! thetal: 2 lignes suivantes a decommenter si vous avez les fichiers
          !     MPL 20130625
          ! fth_fonctions.F90 et parkind1.F90
          ! sinon thetal=theta
          !       thetal(i,k)=fth_thetal(pplay(i,k),t_seri(i,k),q_seri(i,k),
          !    :         ql_seri(i,k))
          thetal(i,k)=theta(i,k)
       ENDDO
    ENDDO
    !

    ! 22.03.04 BEG
    !=============================================================
    !   Ecriture des sorties
    !=============================================================

    ! Recupere des varibles calcule dans differents modules
    ! pour ecriture dans histxxx.nc

    ! Get some variables from module fonte_neige_mod
    CALL fonte_neige_get_vars(pctsrf,  &
         zxfqcalving, zxfqfonte, zxffonte, zxrunofflic)


    !=============================================================
    ! Separation entre thermiques et non thermiques dans les sorties
    ! de fisrtilp
    !=============================================================

    IF (iflag_thermals>=1) THEN
       d_t_lscth=0.
       d_t_lscst=0.
       d_q_lscth=0.
       d_q_lscst=0.
       DO k=1,klev
          DO i=1,klon
             IF (ptconvth(i,k)) THEN
                d_t_lscth(i,k)=d_t_eva(i,k)+d_t_lsc(i,k)
                d_q_lscth(i,k)=d_q_eva(i,k)+d_q_lsc(i,k)
             ELSE
                d_t_lscst(i,k)=d_t_eva(i,k)+d_t_lsc(i,k)
                d_q_lscst(i,k)=d_q_eva(i,k)+d_q_lsc(i,k)
             ENDIF
          ENDDO
       ENDDO

       DO i=1,klon
          plul_st(i)=prfl(i,lmax_th(i)+1)+psfl(i,lmax_th(i)+1)
          plul_th(i)=prfl(i,1)+psfl(i,1)
       ENDDO
    ENDIF

    !On effectue les sorties:

IF (CPPKEY_DUST) THEN
  CALL phys_output_write_spl(itap, pdtphys, paprs, pphis,  &
       pplay, lmax_th, aerosol_couple,                 &
       ok_ade, ok_aie, ivap, ok_sync,                  &
       ptconv, read_climoz, clevSTD,                   &
       ptconvth, d_t, qx, d_qx, d_tr_dyn, zmasse,      &
       flag_aerosol, flag_aerosol_strat, ok_cdnc)
ELSE
    CALL phys_output_write(itap, pdtphys, paprs, pphis,  &
         pplay, lmax_th, aerosol_couple,                 &
         ok_ade, ok_aie, ok_volcan, ivap, iliq, isol, ibs,   &
         ok_sync, ptconv, read_climoz, clevSTD,          &
         ptconvth, d_u, d_t, qx, d_qx, zmasse,           &
         flag_aerosol, flag_aerosol_strat, ok_cdnc,t, u1, v1)
END IF

#ifndef CPP_XIOS
      CALL write_paramLMDZ_phy(itap,nid_ctesGCM,ok_sync)
#endif

    ! Petit appelle de sorties pour accompagner le travail sur phyex
    if ( iflag_physiq == 1 ) then
        call output_physiqex(debut,jD_eq,pdtphys,presnivs,paprs,u,v,t,qx,cldfra,0.*t,0.*t,0.*t,pbl_tke,theta)
    endif

    ! +++ Python coupling: Receive results from NN +++
    IF (lk_pycpl) THEN
       CALL receive_from_python_phys('du', py_du_phys, itap)
       CALL receive_from_python_phys('dv', py_dv_phys, itap)
       !
       ! Swap physics band -> caldyn grid
       IF (MOD(itap-1,24) == 0) THEN
          CALL phys_to_dyn(py_du_phys, py_dv_phys, py_du_dyn, py_dv_dyn)
          CALL pyfld_swap_ph2cd(py_du_dyn, py_dv_dyn, py_du_caldyn, py_dv_caldyn)
       END IF
    END IF
    ! +++++++++++++

    !====================================================================
    ! Arret du modele apres hgardfou en cas de detection d'un
    ! plantage par hgardfou
    !====================================================================

    IF (abortphy==1) THEN
       abort_message ='Plantage hgardfou'
       CALL abort_physic (modname,abort_message,1)
    ENDIF

    ! 22.03.04 END
    !
    !====================================================================
    ! Si c'est la fin, il faut conserver l'etat de redemarrage
    !====================================================================
    !

    ! Disabling calls to the prt_alerte function
    alert_first_call = .FALSE.


    IF (lafin) THEN
       CALL s2s_finalize     ! finalization of source to source tools
       itau_phy = itau_phy + itap
       CALL phyredem ("restartphy.nc")
       !         open(97,form="unformatted",file="finbin")
       !         write(97) u_seri,v_seri,t_seri,q_seri
       !         close(97)

       IF (is_omp_master) THEN

         IF (read_climoz >= 1) THEN
           IF (is_mpi_root) CALL nf95_close(ncid_climoz)
            DEALLOCATE(press_edg_climoz)
            DEALLOCATE(press_cen_climoz)
         ENDIF

       ENDIF

       IF (using_xios) THEN

IF (CPPKEY_INCA) THEN
          IF (ANY(type_trac == ['inca','inco'])) THEN
             IF (is_omp_master) THEN
                CALL finalize_inca
             ENDIF
          ENDIF
END IF

! close xios physiq context (call LMDZ)
          IF (is_omp_master) CALL xios_context_finalize
       ENDIF

       WRITE(lunout,*) ' physiq fin, nombre de steps ou cvpas = 1 : ', Ncvpaseq1

    ENDIF

    !      first=.false.

  END SUBROUTINE physiq

END MODULE physiq_mod

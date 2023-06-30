authorData_SpecificMixed <- function(x,
                                     consent_1L_chr = "",
                                     depnt_var_min_val_1L_dbl = numeric(0),
                                     title_1L_chr = "An R model object",
                                     what_1L_chr = "Shareable",
                                     ...){
  if(what_1L_chr == "Shareable"){
    results_ls <- purrr::map(manufacture(x@c_SpecificResults,
                                         what_1L_chr = "indexed_shareable"),
                             ~{
                               outp_smry_ls <- append(procureSlot(x,
                                                                  "c_SpecificResults@b_SpecificPrivate@private_outp_ls"),
                                                      .x)
                               outp_smry_ls <- outp_smry_ls %>%
                                 write_shareable_mdls(new_dir_nm_1L_chr = "G_Shareable",
                                                      consent_1L_chr = consent_1L_chr,
                                                      depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                      shareable_title_detail_1L_chr = title_1L_chr)
                               
                               outp_smry_ls[-1]
                             })
    x <- renewSlot(x,
                   "c_SpecificResults@a_SpecificShareable@shareable_outp_ls",
                   append(results_ls[[1]],
                          results_ls[-1]))
  }
  return(x)
}
depict_SpecificSynopsis <- function(x,
                                    axis_text_sclg_1L_dbl = 1.5,
                                    axis_title_sclg_1L_dbl = 2,
                                    base_height_1L_dbl = 13,
                                    base_size_1L_dbl = 30,
                                    consent_1L_chr = "",
                                    depnt_var_desc_1L_chr = NA_character_,
                                    depnt_var_min_val_1L_dbl = numeric(0),
                                    dim_plot_heights_int = c(10L, 1L),
                                    dim_plot_log_log_tfmn_1L_lgl = F,
                                    dim_plot_rows_cols_pair_int = c(3L,2L),
                                    labels_chr = c("A","B","C","D"),
                                    label_x_1L_dbl = 0.2,
                                    label_y_1L_dbl = 0.9,
                                    label_size_1L_dbl = 30,
                                    legend_sclg_1L_dbl = 2,
                                    mdl_indcs_int = 1:2,
                                    rel_heights_dbl = c(4,10,1),
                                    scale_dbl = c(0.9,0.9,0.9),
                                    timepoint_old_nms_chr = NA_character_,
                                    timepoint_new_nms_chr = NA_character_,
                                    use_png_fls_1L_lgl = F,
                                    utl_plot_label_1L_chr = " ",
                                    utl_by_rnd_plots_params_ls = list(width_1L_dbl = 6,
                                                                      height_1L_dbl = 4),
                                    what_1L_chr = "composite_mdl",
                                    write_1L_lgl = F,
                                    x_labels_chr = character(0),
                                    y_label_1L_chr = " ",
                                    ...){
  plt <- NULL
  outp_smry_ls <- append(x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                         x@b_SpecificResults@b_SpecificPrivate@private_outp_ls)
  if(!is.na(timepoint_new_nms_chr[1])){
    correspondences_lup <- ready4show::ready4show_correspondences() %>%
      renew(old_nms_chr = timepoint_old_nms_chr,
            new_nms_chr = timepoint_new_nms_chr)
  }else{
    correspondences_lup <- NULL
  }
  if(what_1L_chr == "composite_mdl"){
    plt <- make_cmpst_sctr_and_dnst_plt(outp_smry_ls,
                                        base_size_1L_dbl =  base_size_1L_dbl,
                                        correspondences_lup = correspondences_lup,
                                        depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                        depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                        labels_chr = labels_chr,
                                        label_x_1L_dbl = label_x_1L_dbl,
                                        label_y_1L_dbl = label_y_1L_dbl,
                                        label_size_1L_dbl = label_size_1L_dbl,
                                        mdl_indcs_int = mdl_indcs_int,
                                        use_png_fls_1L_lgl = use_png_fls_1L_lgl)
    write_path_1L_chr <- paste0(outp_smry_ls$path_to_write_to_1L_chr, "/dens_and_sctr.png")
  }
  if(what_1L_chr == "composite_utl"){
    ds_descvs_ls <- manufacture(x,
                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                what_1L_chr = "ds_descvs_ls")
    outp_smry_ls <- append(x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                           x@b_SpecificResults@b_SpecificPrivate@private_outp_ls)
    maui_domains_col_nms_chr <- x@c_SpecificParameters@domain_labels_chr
    first_plt <- rlang::exec(youthvars::make_var_by_round_plt,
                             !!!list(data_tb = outp_smry_ls$scored_data_tb,
                                     legend_sclg_1L_dbl = legend_sclg_1L_dbl,
                                     var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                     round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                     x_label_1L_chr = ds_descvs_ls$dictionary_tb %>%
                                       ready4::get_from_lup_obj(match_value_xx = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                                match_var_nm_1L_chr = "var_nm_chr",
                                                                target_var_nm_1L_chr = "var_desc_chr",
                                                                evaluate_1L_lgl = F) %>% as.vector(),
                                     label_fill_1L_chr = utl_plot_label_1L_chr,
                                     axis_text_sclg_1L_dbl = axis_text_sclg_1L_dbl,
                                     axis_title_sclg_1L_dbl = axis_title_sclg_1L_dbl,
                                     y_label_1L_chr = y_label_1L_chr))
    second_plt <- rlang::exec(youthvars::make_sub_tot_plts,
                              !!!list(data_tb = outp_smry_ls$scored_data_tb,
                                      add_legend_1L_lgl = F,
                                      axis_text_sclg_1L_dbl = axis_text_sclg_1L_dbl,
                                      axis_title_sclg_1L_dbl = axis_title_sclg_1L_dbl,
                                      col_nms_chr = maui_domains_col_nms_chr,
                                      legend_sclg_1L_dbl = legend_sclg_1L_dbl,
                                      plot_rows_cols_pair_int = dim_plot_rows_cols_pair_int,
                                      round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                      heights_int = dim_plot_heights_int,
                                      make_log_log_tfmn_1L_lgl = dim_plot_log_log_tfmn_1L_lgl,
                                      x_labels_chr = x_labels_chr,
                                      y_label_1L_chr = y_label_1L_chr))
    legend_ls <- cowplot::get_legend(first_plt)
    plt <- cowplot::plot_grid(first_plt +
                                ggplot2::theme(legend.position="none"),
                              second_plt,
                              legend_ls,
                              nrow = 3L,
                              rel_heights = rel_heights_dbl,
                              scale = scale_dbl)
    write_path_1L_chr <- paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr, "/Output/_Descriptives/combined_utl.png")
  }
  if(write_1L_lgl){
    ready4::write_with_consent(consented_fn = cowplot::save_plot,
                               prompt_1L_chr = paste0("Do you confirm that you want to write the file ",
                                                      write_path_1L_chr,
                                                      "?"),
                               consent_1L_chr = consent_1L_chr,
                               consented_args_ls = list(filename = write_path_1L_chr,
                                                        plot = plt,
                                                        base_height = base_height_1L_dbl),
                               consented_msg_1L_chr = paste0("File ",
                                                             write_path_1L_chr,
                                                             " has been written."),
                               declined_msg_1L_chr = "Write request cancelled - no new files have been written.")
  }
  return(plt)
}

make_fake_ts_data <- function (outp_smry_ls, depnt_var_min_val_1L_dbl = numeric(0), 
          depnt_vars_are_NA_1L_lgl = T) {
  data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                                                     depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                                                                     predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% 
                                                                       purrr::flatten_chr() %>% unique(), 
                                                                     id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                                                                     round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                                                                     round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr,
                                                                     scaling_fctr_dbl = outp_smry_ls$predr_vars_nms_ls %>% 
                                                                       purrr::flatten_chr() %>% unique() %>%
                                                                       purrr::map_dbl(~ifelse(.x %in% outp_smry_ls$predictors_lup$short_name_chr,
                                                                                              ready4::get_from_lup_obj(outp_smry_ls$predictors_lup,
                                                                                                                       match_var_nm_1L_chr = "short_name_chr",
                                                                                                                       match_value_xx = .x,
                                                                                                                       target_var_nm_1L_chr = "mdl_scaling_dbl"),
                                                                                              1)),
                                                                     tidy_1L_lgl = T)
  # if (identical(outp_smry_ls$round_var_nm_1L_chr, character(0)) | 
  #     ifelse(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)), 
  #            T, is.na(outp_smry_ls$round_var_nm_1L_chr))) {
  #   data_tb <- data_tb %>% dplyr::select(-(outp_smry_ls$predr_vars_nms_ls %>% 
  #                                            purrr::flatten_chr() %>% unique() %>% paste0("_change")))
  # }
  fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) != 
                                                                         outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
  fk_data_tb <- fk_data_ls$syn
  if (identical(outp_smry_ls$round_var_nm_1L_chr, character(0)) | 
      ifelse(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)), 
             T, is.na(outp_smry_ls$round_var_nm_1L_chr))) {
    fk_data_tb <- fk_data_tb %>% dplyr::ungroup()
  }
  else {
    fk_data_tb <- fk_data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr), 
                                                    as.character(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr)))) %>% 
      dplyr::group_by(!!rlang::sym(outp_smry_ls$id_var_nm_1L_chr)) %>% 
      dplyr::mutate(`:=`(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr), 
                         !!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) %>% 
                           transform_timepoint_vals(timepoint_levels_chr = outp_smry_ls$scored_data_tb %>% 
                                                      dplyr::pull(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr)) %>% 
                                                      unique(), bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr))) %>% 
      dplyr::ungroup()
  }
  if (depnt_vars_are_NA_1L_lgl) {
    depnt_vars_chr <- names(fk_data_tb)[names(fk_data_tb) %>% 
                                          purrr::map_lgl(~startsWith(.x, outp_smry_ls$depnt_var_nm_1L_chr))]
    fk_data_tb <- fk_data_tb %>% dplyr::mutate(dplyr::across(dplyr::all_of(depnt_vars_chr), 
                                                             ~NA_real_))
  }
  return(fk_data_tb)
}
make_shareable_mdl <- function (fake_ds_tb, mdl_smry_tb, control_1L_chr = NA_character_, 
          depnt_var_nm_1L_chr = "utl_total_w", id_var_nm_1L_chr = "fkClientID", 
          mdl_type_1L_chr = "OLS_CLL", mdl_types_lup = NULL, seed_1L_int = 12345L, 
          start_1L_chr = NA_character_, tfmn_1L_chr = "CLL") 
{
  if (is.null(mdl_types_lup)) 
    utils::data(mdl_types_lup, envir = environment())
  if (is.na(tfmn_1L_chr)) 
    tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                                            match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                                            target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
  predr_var_nms_chr <- mdl_smry_tb$Parameter[!mdl_smry_tb$Parameter %in% 
                                               c("SD (Intercept)", "Intercept", "R2", "RMSE", "Sigma")] %>% 
    stringi::stri_replace_last_fixed(" baseline", "_baseline") %>% 
    stringi::stri_replace_last_fixed(" change", "_change") %>% 
    stringi::stri_replace_last_fixed(" scaled", "_scaled") %>% 
    stringi::stri_replace_last_fixed(" unscaled", "_unscaled")
  X <- ready4use::Ready4useDyad(ds_tb = outp_smry_ls$scored_data_tb, dictionary_r3 = outp_smry_ls$dictionary_tb)
  dummys_chr <- manufacture(X, flatten_1L_lgl = T) 
  predr_var_nms_chr <- predr_var_nms_chr %>% purrr::map_chr(~ifelse(.x %in% dummys_chr,
                                                                    manufacture(X, flatten_1L_lgl = T, what_1L_chr = "factors-d", match_1L_chr = .x),
                                                                    .x)) %>% unique()
  tfd_depnt_var_nm_1L_chr <- transform_depnt_var_nm(depnt_var_nm_1L_chr, 
                                                    tfmn_1L_chr = tfmn_1L_chr)
  if (length(predr_var_nms_chr) > 1) {
    covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
  } else {
    covar_var_nms_chr <- NA_character_
  }
  model_mdl <- make_mdl(fake_ds_tb %>% dplyr::select(tidyselect::all_of(c(id_var_nm_1L_chr, 
                                                                          tfd_depnt_var_nm_1L_chr, predr_var_nms_chr))), depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                        predr_var_nm_1L_chr = predr_var_nms_chr[1], covar_var_nms_chr = covar_var_nms_chr, 
                        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
                        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr, 
                        start_1L_chr = start_1L_chr)
  if (ready4::get_from_lup_obj(mdl_types_lup, match_value_xx = mdl_type_1L_chr, 
                               match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "fn_chr", 
                               evaluate_1L_lgl = F) == "betareg::betareg") {
    model_coeffs_dbl <- model_mdl$coefficients$mean
  }
  else {
    model_coeffs_dbl <- model_mdl$coefficients
  }
  param_nms_chr <- model_coeffs_dbl %>% names()
  mdl_smry_tb <- mdl_smry_tb %>% dplyr::mutate(Parameter = dplyr::case_when(Parameter == 
                                                                              "Intercept" ~ "(Intercept)", TRUE ~ purrr::map_chr(Parameter, 
                                                                                                                                 ~stringr::str_replace_all(.x, " ", "_")))) %>% dplyr::filter(Parameter %in% 
                                                                                                                                                                                                param_nms_chr) %>% dplyr::slice(match(param_nms_chr, 
                                                                                                                                                                                                                                      Parameter))
  assertthat::assert_that(all(param_nms_chr == mdl_smry_tb$Parameter), 
                          msg = "Parameter names mismatch between data and model summary table")
  model_coeffs_dbl <- mdl_smry_tb$Estimate
  names(model_coeffs_dbl) <- param_nms_chr
  if (ready4::get_from_lup_obj(mdl_types_lup, match_value_xx = mdl_type_1L_chr, 
                               match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = "fn_chr", 
                               evaluate_1L_lgl = F) == "betareg::betareg") {
    model_mdl$coefficients$mean <- model_coeffs_dbl
  }
  else {
    model_mdl$coefficients <- model_coeffs_dbl
  }
  return(model_mdl)
}
make_smry_of_brm_mdl <- function (mdl_ls, data_tb, depnt_var_nm_1L_chr = "utl_total_w", 
                                  predr_vars_nms_chr, mdl_nm_1L_chr = NA_character_, seed_1L_dbl = 23456, 
                                  tfmn_1L_chr) 
{
  if (is.na(mdl_nm_1L_chr)) 
    mdl_nm_1L_chr <- predr_vars_nms_chr[1]
  set.seed(seed_1L_dbl)
  predictions <- stats::predict(mdl_ls, summary = F) %>% calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, 
                                                                                  tfmn_is_outp_1L_lgl = T)
  sd_intcpt_df <- summary(mdl_ls, digits = 4)$random[[1]]
  sd_intcpt_df <- sd_intcpt_df[1:nrow(sd_intcpt_df), 1:4] %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
  coef <- summary(mdl_ls, digits = 4)$fixed
  coef <- coef[1:nrow(coef), 1:4] %>% dplyr::mutate(dplyr::across(dplyr::everything(), 
                                                                  as.numeric))
  R2 <- brms::bayes_R2(mdl_ls) %>% as.vector()
  RMSE <- psych::describe(apply(predictions, 1, calculate_rmse, 
                                y_dbl = data_tb %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), 
                          quant = c(0.25, 0.75), skew = F, ranges = F)
  RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75) %>% 
    as.vector()
  Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
  smry_of_brm_mdl_tb <- data.frame(round(rbind(sd_intcpt_df, 
                                               coef, R2, RMSE, Sigma), 3)) %>% 
    dplyr::mutate(Parameter = c("SD (Intercept)", 
                                "Intercept",
                                purrr::map(predr_vars_nms_chr, 
                                           ~{
                                             possibilities_chr <- paste0(.x, c("", " baseline", " change", " scaled", " unscaled"))
                                             if(possibilities_chr[1] %in% names(mdl_ls$data)){
                                               values_xx <- mdl_ls$data %>% dplyr::pull(.x)
                                               if(is.factor(values_xx)){
                                                 possibilities_chr <- c(possibilities_chr[1], paste0(.x,levels(values_xx)[2:length(levels(values_xx))]))
                                               }
                                             }
                                             possibilities_chr
                                           }) %>% 
                                  purrr::flatten_chr() %>% intersect(purrr::map(names(mdl_ls$data), 
                                                                                ~{
                                                                                  values_xx <- mdl_ls$data %>% dplyr::pull(.x)
                                                                                  if(is.factor(values_xx)){
                                                                                    paste0(.x,levels(values_xx)[2:length(levels(values_xx))])
                                                                                  }else{
                                                                                    stringi::stri_replace_last_fixed(.x, "_baseline", " baseline") %>% 
                                                                                      stringi::stri_replace_last_fixed("_change",  " change") %>% 
                                                                                      stringi::stri_replace_last_fixed("_scaled"," scaled") %>% 
                                                                                      stringi::stri_replace_last_fixed("_unscaled"," unscaled")
                                                                                  }
                                                                                }) %>%
                                                                       purrr::flatten_chr()
                                  ),
                                "R2", "RMSE", "Sigma"), Model = mdl_nm_1L_chr) %>% 
    dplyr::mutate(`95% CI` = paste(l.95..CI, ",", u.95..CI)) %>% 
    dplyr::rename(SE = Est.Error) %>% dplyr::select(Model, 
                                                    Parameter, Estimate, SE, `95% CI`)
  rownames(smry_of_brm_mdl_tb) <- NULL
  return(smry_of_brm_mdl_tb)
}
transform_tb_to_mdl_inp <- function (data_tb, depnt_var_min_val_1L_dbl = numeric(0), depnt_var_max_val_1L_dbl = 0.99999, 
          depnt_var_nm_1L_chr = "utl_total_w", predr_vars_nms_chr, 
          id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
          round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, 
          scaling_fctr_dbl = 1, tfmn_1L_chr = "NTF", tidy_1L_lgl = F, ungroup_1L_lgl = F) 
{
  if (length(scaling_fctr_dbl) != length(predr_vars_nms_chr)) {
    scaling_fctr_dbl <- rep(scaling_fctr_dbl[1], length(predr_vars_nms_chr))
  }
  data_tb <- data.frame(data_tb) %>% ready4use::remove_labels_from_ds()
  tfd_for_mdl_inp_tb <- data_tb %>% dplyr::select(dplyr::all_of(id_var_nm_1L_chr), 
                                                  dplyr::all_of(round_var_nm_1L_chr), dplyr::all_of(predr_vars_nms_chr), 
                                                  dplyr::all_of(depnt_var_nm_1L_chr)) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr))
  tfd_for_mdl_inp_tb <- if (!identical(round_var_nm_1L_chr, 
                                       character(0)) && ifelse(identical(round_var_nm_1L_chr, 
                                                                         character(0)), T, !is.na(round_var_nm_1L_chr))) {
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), 
                                                                !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr), 
                                        .init = tfd_for_mdl_inp_tb, ~{
                                          idx_1L_int <- as.integer(.y)
                                          .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]), 
                                                                             .fns = list(baseline = ~if (!is.numeric(.)) {
                                                                               .
                                                                             } else {
                                                                               dplyr::first(.) * scaling_fctr_dbl[idx_1L_int]
                                                                             }, change = ~ifelse(!!rlang::sym(round_var_nm_1L_chr) == 
                                                                                                   round_bl_val_1L_chr, 0, if (!is.numeric(.)) {
                                                                                                     .
                                                                                                   } else {
                                                                                                     (. - dplyr::lag(.)) * scaling_fctr_dbl[idx_1L_int]
                                                                                                   }))))
                                        })
  }  else {
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr), 
                                        .init = tfd_for_mdl_inp_tb, ~{
                                          idx_1L_int <- as.integer(.y)
                                          table_tb <- .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]), 
                                                                                         .fns = list(baseline = ~if (!is.numeric(.)) {
                                                                                           .
                                                                                         } else {
                                                                                           dplyr::first(.) * scaling_fctr_dbl[idx_1L_int]
                                                                                         }, change = ~0)))
                                          old_name_1L_chr <- paste0(predr_vars_nms_chr[idx_1L_int], 
                                                                    "_baseline")
                                          new_name_1L_chr <- paste0(predr_vars_nms_chr[idx_1L_int], 
                                                                    ifelse(scaling_fctr_dbl[idx_1L_int] == 1, "_unscaled", 
                                                                           "_scaled"))
                                          table_tb <- table_tb %>% dplyr::rename(`:=`(!!rlang::sym(new_name_1L_chr), 
                                                                                      !!rlang::sym(old_name_1L_chr)))
                                        })
  }
  if (!identical(depnt_var_min_val_1L_dbl, numeric(0))) {
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr), 
                                                                    !!rlang::sym(depnt_var_nm_1L_chr) %>% purrr::map_dbl(~max(.x, 
                                                                                                                              depnt_var_min_val_1L_dbl))))
  }
  tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                                                                 tfmn_1L_chr = tfmn_1L_chr, depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl)
  if (drop_all_msng_1L_lgl) {
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% stats::na.omit()
  }
  if (ungroup_1L_lgl) {
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::ungroup()
  }
  tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr)
  if(tidy_1L_lgl){
    if (identical(round_var_nm_1L_chr, character(0)) | 
        ifelse(identical(round_var_nm_1L_chr, character(0)), 
               T, is.na(round_var_nm_1L_chr))) {
      tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% 
        dplyr::select(-(predr_vars_nms_chr %>% paste0("_change"))) %>%
        dplyr::select(-(intersect(predr_vars_nms_chr %>% paste0("_unscaled"), names(tfd_for_mdl_inp_tb))))
    }
  }
  return(tfd_for_mdl_inp_tb)
}
write_shareable_mdls <- function (outp_smry_ls, consent_1L_chr = "", consent_indcs_int = 1L, 
                                  depnt_var_min_val_1L_dbl = numeric(0), new_dir_nm_1L_chr = "G_Shareable", 
                                  options_chr = c("Y", "N"), shareable_title_detail_1L_chr = "", 
                                  write_mdls_to_dv_1L_lgl = F) {
  output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls, 
                                        consent_1L_chr = consent_1L_chr, consent_indcs_int = consent_indcs_int, 
                                        new_dir_nm_1L_chr = new_dir_nm_1L_chr, options_chr = options_chr)
  incld_mdl_paths_chr <- make_incld_mdl_paths(outp_smry_ls)
  fake_ds_tb <- make_fake_ts_data(outp_smry_ls, depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                  depnt_vars_are_NA_1L_lgl = F)
  mdl_types_lup <- outp_smry_ls$mdl_types_lup
  shareable_mdls_ls <- outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr() %>% 
    purrr::map2(incld_mdl_paths_chr, ~{
      model_mdl <- readRDS(paste0(outp_smry_ls$path_to_write_to_1L_chr, "/", .y))
      mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model == .x)
      mdl_nm_1L_chr <- .x
      mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr, 
                                              mdl_types_lup = mdl_types_lup)
      tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                                              match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                                              target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
      predn_type_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                                                    match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                                                    target_var_nm_1L_chr = "predn_type_chr", evaluate_1L_lgl = F)
      if (is.na(predn_type_1L_chr)) 
        predn_type_1L_chr <- NULL
      control_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup, 
                                                 match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                                                 target_var_nm_1L_chr = "control_chr", evaluate_1L_lgl = F)
      sd_dbl <- mdl_smry_tb %>% dplyr::filter(Parameter == 
                                                "SD (Intercept)") %>% dplyr::select(Estimate, 
                                                                                    SE) %>% t() %>% as.vector()
      mdl_fake_ds_tb <- fake_ds_tb %>% add_tfd_var_to_ds(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                                                         tfmn_1L_chr = tfmn_1L_chr, depnt_var_max_val_1L_dbl = 0.999) %>% 
        dplyr::select(names(model_mdl$data))
      model_mdl$data <- mdl_fake_ds_tb
      table_predn_mdl <- make_shareable_mdl(fake_ds_tb = mdl_fake_ds_tb, 
                                            mdl_smry_tb = mdl_smry_tb, 
                                            depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                                            id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                                            tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr, 
                                            mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr, 
                                            start_1L_chr = NA_character_, seed_1L_int = outp_smry_ls$seed_1L_int)
      c(4, 3) %>% purrr::walk2(list(table_predn_mdl, model_mdl), 
                               ~{
                                 ready4::write_with_consent(consented_fn = saveRDS, 
                                                            prompt_1L_chr = paste0("Do you confirm that you want to write the file ", 
                                                                                   paste0(mdl_nm_1L_chr, ".RDS"), " to ", 
                                                                                   output_dir_chr[.x], "?"), consent_1L_chr = consent_1L_chr, 
                                                            consent_indcs_int = consent_indcs_int, consented_args_ls = list(object = .y, 
                                                                                                                            file = paste0(output_dir_chr[.x], "/", 
                                                                                                                                          mdl_nm_1L_chr, ".RDS")), consented_msg_1L_chr = paste0("File ", 
                                                                                                                                                                                                 paste0(mdl_nm_1L_chr, ".RDS"), " has been written to ", 
                                                                                                                                                                                                 output_dir_chr[.x], "."), declined_msg_1L_chr = "Write request cancelled - no new files have been written.", 
                                                            options_chr = options_chr)
                               })
      scaling_fctr_dbl <- make_scaling_fctr_dbl(outp_smry_ls)
      write_ts_mdl_plts(brms_mdl = model_mdl, consent_1L_chr = consent_1L_chr, 
                        consent_indcs_int = consent_indcs_int, depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                        mdl_nm_1L_chr = mdl_nm_1L_chr, options_chr = options_chr, 
                        path_to_write_to_1L_chr = output_dir_chr[3], 
                        predn_type_1L_chr = predn_type_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                        sd_dbl = sd_dbl, sfx_1L_chr = " from table", 
                        table_predn_mdl = table_predn_mdl, tfd_data_tb = outp_smry_ls$scored_data_tb %>% 
                          transform_tb_to_mdl_inp(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                                  depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                                                  predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% 
                                                    purrr::flatten_chr() %>% unique(), id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                                                  round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                                                  round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr, 
                                                  scaling_fctr_dbl = scaling_fctr_dbl), tfmn_1L_chr = tfmn_1L_chr, 
                        utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl), 
                                                    outp_smry_ls$utl_min_val_1L_dbl, -1))
      table_predn_mdl
    }) %>% stats::setNames(outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr())
  outp_smry_ls$shareable_mdls_ls <- shareable_mdls_ls
  outp_smry_ls$shareable_mdls_tb <- NULL
  ingredients_ls <- list(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, 
                         dictionary_tb = outp_smry_ls$dictionary_tb %>% dplyr::filter(var_nm_chr %in% 
                                                                                        names(fake_ds_tb)), id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, 
                         fake_ds_tb = fake_ds_tb, mdls_lup = outp_smry_ls$shareable_mdls_ls %>% 
                           purrr::map2_dfr(names(outp_smry_ls$shareable_mdls_ls), 
                                           ~{
                                             if (inherits(.x, "betareg")) {
                                               coeffs_dbl <- .x$coefficients$mean
                                             } else {
                                               coeffs_dbl <- .x$coefficients
                                             }
                                             mdl_type_1L_chr = get_mdl_type_from_nm(.y, 
                                                                                    mdl_types_lup = outp_smry_ls$mdl_types_lup)
                                             tibble::tibble(mdl_nms_chr = .y) %>% dplyr::mutate(predrs_ls = list(coeffs_dbl %>% 
                                                                                                                   names() %>% stringr::str_remove_all("_change") %>% 
                                                                                                                   stringr::str_remove_all("_baseline") %>% 
                                                                                                                   stringr::str_remove_all("_scaled") %>% stringr::str_remove_all("_unscaled") %>% 
                                                                                                                   unique() %>% purrr::discard(~.x == "(Intercept)")), 
                                                                                                mdl_type_chr = mdl_type_1L_chr, tfmn_chr = ready4::get_from_lup_obj(outp_smry_ls$mdl_types_lup, 
                                                                                                                                                                    match_value_xx = mdl_type_1L_chr, match_var_nm_1L_chr = "short_name_chr", 
                                                                                                                                                                    target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F))
                                           }), mdls_smry_tb = outp_smry_ls$mdls_smry_tb, 
                         mdl_types_lup = mdl_types_lup, predictors_lup = outp_smry_ls$predictors_lup, 
                         round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr, 
                         seed_1L_int = outp_smry_ls$seed_1L_int, utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl), 
                                                                                             outp_smry_ls$utl_min_val_1L_dbl, -1))
  ready4::write_with_consent(consented_fn = saveRDS, prompt_1L_chr = paste0("Do you confirm that you want to write the file ", 
                                                                            paste0("mdl_ingredients", ".RDS"), " to ", output_dir_chr[2], 
                                                                            "?"), consent_1L_chr = consent_1L_chr, consent_indcs_int = consent_indcs_int, 
                             consented_args_ls = list(object = ingredients_ls, file = paste0(output_dir_chr[2], 
                                                                                             "/", "mdl_ingredients", ".RDS")), consented_msg_1L_chr = paste0("File ", 
                                                                                                                                                             paste0("mdl_ingredients", ".RDS"), " has been written to ", 
                                                                                                                                                             output_dir_chr[2], "."), declined_msg_1L_chr = "Write request cancelled - no new files have been written.", 
                             options_chr = options_chr)
  outp_smry_ls <- write_mdls_to_dv(outp_smry_ls, consent_1L_chr = consent_1L_chr, 
                                   consent_indcs_int = consent_indcs_int, new_dir_nm_1L_chr = new_dir_nm_1L_chr, 
                                   options_chr = options_chr, output_dir_chr = output_dir_chr, 
                                   shareable_title_detail_1L_chr = shareable_title_detail_1L_chr)
  return(outp_smry_ls)
}
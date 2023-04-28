## For stand alone program
make_csnl_example_dict <- function(ds_tb){
  dictionary_r3 <- Ready4useRepos(dv_nm_1L_chr = "TTU",
                                  dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                  dv_server_1L_chr = "dataverse.harvard.edu") %>% 
    ingest(fls_to_ingest_chr = c("dictionary_r3"), 
           metadata_1L_lgl = F) 
  dictionary_r3 <- dictionary_r3 %>%
    dplyr::filter(var_nm_chr %in% names(ds_tb)) %>%
    renew.ready4use_dictionary(var_nm_chr = setdiff(names(ds_tb),dictionary_r3$var_nm_chr) %>% sort(),
                               var_ctg_chr = c(rep("clinical symptom",2),
                                               rep("multi-attribute utility instrument question",9),
                                               "health_utility",
                                               rep("demographic",3),
                                               "psychological distress",
                                               "quality of life",
                                               rep("spatial",2)),
                               var_desc_chr = c("Days Unable To Perform Usual Activities",
                                                "Days Cut Back On Usual Activities",
                                                paste0("Child Health Utility (9 Dimension) question ",1:9),
                                                "Child Health Utility (9 Dimension) total score",
                                                "Employed",
                                                "Employment Type",
                                                "Studying",
                                                "Kessler Psychological Distress Scale (10 Item)",
                                                "My Life Tracker", 
                                                "Area Index of Relative Social Disadvantage",
                                                "Area Remoteness"
                               ),
                               var_type_chr = setdiff(names(ds_tb),dictionary_r3$var_nm_chr) %>% 
                                 sort() %>% purrr::map_chr(~{
                                   classes_chr <- ds_tb[,.x][[1]] %>% class()
                                   ifelse("numeric" %in% classes_chr,
                                          ifelse(is.integer(ds_tb[,.x][[1]]),
                                                 "integer",
                                                 "double"),
                                          classes_chr[1])
                                 })) %>% dplyr::arrange(var_ctg_chr, var_nm_chr)
  dictionary_r3
}
make_csnl_example_predrs <- function(){
  predictors_r3 <- Ready4useRepos(dv_nm_1L_chr = "TTU", 
                                  dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                  dv_server_1L_chr = "dataverse.harvard.edu") %>%
    ingest(fls_to_ingest_chr = c("predictors_r3"),
           metadata_1L_lgl = F)
  predictors_r3 <- renew.specific_predictors(predictors_r3,
                                             filter_cdn_1L_chr = "short_name_chr == 'SOFAS'") %>%
    renew.specific_predictors(short_name_chr = c("K10", "MLT"),
                              long_name_chr = c("K10 total score", "MLT total score"),
                              min_val_dbl = c(10,0),
                              max_val_dbl = c(50,100),
                              class_chr = c("integer","numeric"),
                              increment_dbl = 1,
                              class_fn_chr = c("as.integer","as.double"), # update when new youthvars classes are created
                              mdl_scaling_dbl = 0.01,
                              covariate_lgl = F) %>%
    dplyr::mutate(covariate_lgl = dplyr::case_when(short_name_chr == "SOFAS" ~ F,
                                                   T ~ covariate_lgl))
  return(predictors_r3)
}
## From ready4use
renew_Ready4useDyad <- function(x,
                                remove_old_lbls_1L_lgl = T,
                                tfmn_1L_chr = "capitalise",
                                type_1L_chr = "label"){
  if(type_1L_chr %in% c("label","case")){
    dictionary_tb <- x@dictionary_r3
    if(tfmn_1L_chr == "capitalise")
      dictionary_tb$var_desc_chr <- dictionary_tb$var_desc_chr %>%
        Hmisc::capitalize()
    if(tfmn_1L_chr == "title")
      dictionary_tb$var_desc_chr <- dictionary_tb$var_desc_chr %>%
        stringr::str_to_title()
  }
  if(type_1L_chr == "case"){
    x@dictionary_r3 <- dictionary_tb
  }
  if(type_1L_chr == "label"){
    tfd_ds_tb <- add_labels_from_dictionary(x@ds_tb,
                                            dictionary_tb = dictionary_tb %>% ready4::remove_lbls_from_df(),
                                            remove_old_lbls_1L_lgl = remove_old_lbls_1L_lgl)
    x@ds_tb <- tfd_ds_tb
  }
  if(type_1L_chr == "unlabel"){
    x@ds_tb <- remove_labels_from_ds(x@ds_tb)
  }
  return(x)
}
## From youthvars
# make_descv_stats_tbl <- function(data_tb,
#                                  key_var_nm_1L_chr = "round",
#                                  key_var_vals_chr = NULL,
#                                  variable_nms_chr,
#                                  dictionary_tb = NULL,
#                                  test_1L_lgl = F,
#                                  sections_as_row_1L_lgl = F,
#                                  nbr_of_digits_1L_int = NA_integer_){
#   if(!identical(key_var_nm_1L_chr, character(0))){
#     if(is.null(key_var_vals_chr)){
#       key_var_vals_chr <- data_tb %>% dplyr::pull(key_var_nm_1L_chr) %>% unique() %>% as.character()
#     }
#   }
#   if(identical(key_var_vals_chr, character(0))){
#     key_var_vals_chr <- "Overall"
#   }
#   if(length(key_var_vals_chr)<2 & test_1L_lgl){ ### EDIT
#     descv_stats_tbl_tb <- NULL
#   }else{
#     descv_stats_tbl_tb <- make_tableby_ls(data_tb,
#                                           key_var_nm_1L_chr = key_var_nm_1L_chr,
#                                           variable_nms_chr = variable_nms_chr,
#                                           test_1L_lgl = test_1L_lgl) %>%
#       as.data.frame() 
#     descv_stats_tbl_tb <- descv_stats_tbl_tb %>%
#       dplyr::select(c("variable","label",
#                       key_var_vals_chr,
#                       ifelse(test_1L_lgl,
#                              "p.value",
#                              character(0))) %>% purrr::discard(is.na))
#     if(!is.null(dictionary_tb)){
#       descv_stats_tbl_tb <- descv_stats_tbl_tb %>%
#         dplyr::mutate(variable = variable %>% purrr::map_chr(~ready4::get_from_lup_obj(dictionary_tb,
#                                                                                        target_var_nm_1L_chr = "var_desc_chr",
#                                                                                        match_var_nm_1L_chr = "var_nm_chr",
#                                                                                        match_value_xx = .x,
#                                                                                        evaluate_1L_lgl = F) %>% as.vector()))
#       
#     }
#     vars_with_mdns_chr <- descv_stats_tbl_tb %>% dplyr::filter(label == "Median (Q1, Q3)") %>% dplyr::pull(variable)
#     descv_stats_tbl_tb <- descv_stats_tbl_tb %>%
#       dplyr::mutate(dplyr::across(tidyselect::all_of(key_var_vals_chr),
#                                   ~ list(.x) %>% purrr::pmap_dbl(~{
#                                     ifelse(..1[[1]][[1]] =="",
#                                            NA_real_,
#                                            ..1[[1]][[1]])
#                                   }),
#                                   .names = "{col}_val_1_dbl"),
#                     dplyr::across(tidyselect::all_of(key_var_vals_chr),
#                                   ~ list(.x,variable,label) %>%
#                                     purrr::pmap(~ {
#                                       if(..2 %in% vars_with_mdns_chr){
#                                         if(..3 == "Median (Q1, Q3)"){
#                                           return_dbl <- c(..1[[2]],..1[[3]])
#                                         }else{
#                                           return_dbl <- ifelse(length(..1) == 1,
#                                                                NA_real_,
#                                                                ..1[[2]])
#                                         }
#                                       }else{
#                                         return_dbl <- ifelse(length(..1) == 1,
#                                                              NA_real_,
#                                                              ifelse(..1[[2]]=="",
#                                                                     NA_real_,
#                                                                     ..1[[2]]))
#                                       }
#                                     }
#                                     ),
#                                   .names = "{col}_val_2_ls")) %>%
#       dplyr::select(variable,
#                     label,
#                     key_var_vals_chr %>% purrr::map(~c(paste0(.x, c("_val_1_dbl","_val_2_ls")))) %>%
#                       purrr::flatten_chr(),
#                     ifelse(test_1L_lgl,"p.value",character(0)) %>% purrr::discard(is.na)
#       )
#     if(sections_as_row_1L_lgl){
#       descv_stats_tbl_tb <- descv_stats_tbl_tb %>%
#         dplyr::select(-variable)
#     }else{
#       descv_stats_tbl_tb <- descv_stats_tbl_tb %>%
#         dplyr::filter(label != variable)
#     }
#     if(!is.na(nbr_of_digits_1L_int)){
#       descv_stats_tbl_tb <- c(key_var_vals_chr %>%
#                                 purrr::map(~c(paste0(.x, c("_val_1_dbl","_val_2_ls"
#                                 )))) %>%
#                                 purrr::flatten_chr(),
#                               ifelse(test_1L_lgl,"p.value",character(0)) %>% purrr::discard(is.na)) %>%
#         purrr::reduce(.init = descv_stats_tbl_tb,
#                       ~ .x %>% dplyr::mutate(!!rlang::sym(.y) := !!rlang::sym(.y) %>%
#                                                purrr::map_chr(~ {
#                                                  ifelse(length(.x) == 1,
#                                                         ifelse(is.na(.x),
#                                                                "",
#                                                                paste0("",
#                                                                       format(round(.x, nbr_of_digits_1L_int),
#                                                                              nsmall = nbr_of_digits_1L_int),
#                                                                       ""
#                                                                )
#                                                         ),
#                                                         paste0("",
#                                                                .x %>%
#                                                                  purrr::map_chr(~format(round(.x,
#                                                                                               nbr_of_digits_1L_int),
#                                                                                         nsmall = nbr_of_digits_1L_int)) %>%
#                                                                  paste0(collapse = ", "),
#                                                                ""
#                                                         )
#                                                  )
#                                                })))
#       descv_stats_tbl_tb <- paste0(key_var_vals_chr, "_val_2_ls") %>%
#         purrr::reduce(.init = descv_stats_tbl_tb,
#                       ~ .x %>%
#                         dplyr::mutate(!!rlang::sym(.y) := !!rlang::sym(.y) %>%
#                                         purrr::map2_chr(label,
#                                                         ~ ifelse(.x=="" | .y== "Min - Max",
#                                                                  .x,
#                                                                  paste0("(",
#                                                                         .x,
#                                                                         ifelse(.y %in% c("Mean (SD)","Median (Q1, Q3)","Missing"),
#                                                                                "",
#                                                                                "%"),
#                                                                         ")")))
#                                       
#                         ))
#     }
#   }
#   return(descv_stats_tbl_tb)
# }
# make_tableby_ls <- function(data_tb,
#                             key_var_nm_1L_chr,
#                             variable_nms_chr,
#                             test_1L_lgl = F){
#   forumla_fml <- make_formula(key_var_nm_1L_chr,
#                               predictors_chr = variable_nms_chr)
#   tableby_ls <- arsenal::tableby(forumla_fml,
#                                  data = data_tb,
#                                  control = make_tableby_cntrls(test_1L_lgl))
#   return(tableby_ls)
# }
# make_cors_with_utl_tbl <- function(data_tb,
#                                    ds_descvs_ls,
#                                    dictionary_tb = NULL,
#                                    cor_type_1L_chr = "pearson"){
#   cors_with_utl_tb <- purrr::map(ds_descvs_ls$round_vals_chr, ~data_tb %>% dplyr::filter(!!rlang::sym(ds_descvs_ls$round_var_nm_1L_chr) ==
#                                                                                            .x) %>% dplyr::select(!!!rlang::syms(c(ds_descvs_ls$utl_wtd_var_nm_1L_chr,ds_descvs_ls$candidate_predrs_chr))) %>% as.matrix() %>% Hmisc::rcorr(type = cor_type_1L_chr)) %>%
#     purrr::map2_dfc(ds_descvs_ls$round_vals_chr,
#                     ~tibble::tibble(!!rlang::sym(paste0(.y,"_cor_dbl")) := .x[[1]][2:(length(ds_descvs_ls$candidate_predrs_chr)+1)],
#                                     !!rlang::sym(paste0(.y,"_sig_dbl")) := .x[[3]][2:(length(ds_descvs_ls$candidate_predrs_chr)+1)])) %>%
#     dplyr::mutate(variable_chr = ds_descvs_ls$candidate_predrs_chr) %>%
#     dplyr::select(variable_chr, dplyr::everything())
#   if(!is.null(dictionary_tb)){
#     cors_with_utl_tb <- cors_with_utl_tb %>%
#       dplyr::mutate(variable_chr = variable_chr %>% purrr::map_chr(~ready4::get_from_lup_obj(dictionary_tb,
#                                                                                              target_var_nm_1L_chr = "var_desc_chr",
#                                                                                              match_var_nm_1L_chr = "var_nm_chr",
#                                                                                              match_value_xx = .x,
#                                                                                              evaluate_1L_lgl = F)))
#     
#   }
#   return(cors_with_utl_tb)
# }
transform_csnl_example_ds <- function(ds_df){
  ds_tb <- ds_df %>% 
    tibble::as_tibble() %>%
    dplyr::rename(c_days_cut_back = K12_DaysCutDown,
                  c_days_unable = K11_DaysTotallyUnable,
                  c_p_diag_s = DiagnosisPrimary,
                  d_age = Age,
                  d_ATSI = ATSI,
                  #d_CALD = CALD, # Uncomment when dictionary is updated
                  d_employed = Working,
                  d_employment_type = EmploymentType,
                  d_gender = Gender,
                  d_studying = Studying,
                  K10 = K10_total,
                  MLT = MLT_mean, 
                  s_IRSD = IRSD,
                  s_remoteness = Remoteness) %>%
    dplyr::select(-c("aqol6d_total_c",
                     "aqol6d_total_w"))
  ds_tb <- ds_tb %>%
    dplyr::mutate(dplyr::across(c(dplyr::starts_with("aqol6d_q"),
                                  dplyr::starts_with("chu9_q"),
                                  dplyr::starts_with("c_days_"),
                                  d_age,
                                  K10,
                                  s_IRSD,
                                  SOFAS,
    ), ~as.integer(.x)))
  ds_tb <- youthvars::add_uids_to_tbs_ls(list(ds_tb),"Participant_") %>% purrr::pluck(1)
  return(ds_tb)
}
write_descv_tbls <- function(data_tb,
                             ds_descvs_ls,
                             predictors_lup,
                             descv_outp_dir_1L_chr,
                             nbr_of_digits_1L_int = 2,
                             participation_var_1L_chr = "participation"){
  
  descv_tbl_ls <- list(cohort_desc_tb = make_descv_stats_tbl(data_tb = data_tb, #
                                                             key_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                             key_var_vals_chr = ds_descvs_ls$round_vals_chr,
                                                             dictionary_tb = ds_descvs_ls$dictionary_tb,
                                                             variable_nms_chr = ds_descvs_ls$cohort_descv_var_nms_chr,
                                                             nbr_of_digits_1L_int = nbr_of_digits_1L_int),
                       main_outc_tbl_tb = make_descv_stats_tbl(data_tb = data_tb, #
                                                               key_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                               key_var_vals_chr = ds_descvs_ls$round_vals_chr,#if(identical(ds_descvs_ls$round_vals_chr, character(0))){"Overall"}else{ds_descvs_ls$round_vals_chr},
                                                               dictionary_tb = ds_descvs_ls$dictionary_tb,
                                                               variable_nms_chr = c(ds_descvs_ls$candidate_predrs_chr,
                                                                                    ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                                                    ds_descvs_ls$utl_unwtd_var_nm_1L_chr),
                                                               test_1L_lgl = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){F}else{T},
                                                               nbr_of_digits_1L_int = nbr_of_digits_1L_int),
                       outc_by_partcn_tbl_tb = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NULL}else{make_descv_stats_tbl(data_tb = data_tb %>% dplyr::filter(!!rlang::sym(ds_descvs_ls$round_var_nm_1L_chr) == ds_descvs_ls$round_vals_chr[1]), #
                                                                    key_var_nm_1L_chr = participation_var_1L_chr,
                                                                    key_var_vals_chr = data_tb %>% dplyr::pull(participation_var_1L_chr) %>% unique(),
                                                                    dictionary_tb = ds_descvs_ls$dictionary_tb,
                                                                    variable_nms_chr = c(ds_descvs_ls$candidate_predrs_chr,
                                                                                         ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                                                         ds_descvs_ls$utl_unwtd_var_nm_1L_chr),
                                                                    test_1L_lgl = T,
                                                                    nbr_of_digits_1L_int = nbr_of_digits_1L_int)},
                       bl_cors_tb = transform_ds_for_tstng(data_tb,
                                                           depnt_var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                           depnt_var_max_val_1L_dbl = Inf,
                                                           candidate_predrs_chr = ds_descvs_ls$candidate_predrs_chr,
                                                           round_var_nm_1L_chr = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_descvs_ls$round_var_nm_1L_chr},
                                                           round_val_1L_chr = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_descvs_ls$round_vals_chr[1]}) %>%
                         make_corstars_tbl_xx(result_chr = "none"),##
                       fup_cors_tb = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NULL}else{transform_ds_for_tstng(data_tb,
                                                            depnt_var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                            depnt_var_max_val_1L_dbl = Inf,
                                                            candidate_predrs_chr = ds_descvs_ls$candidate_predrs_chr,
                                                            round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                            round_val_1L_chr = ds_descvs_ls$round_vals_chr[2]) %>%
                         make_corstars_tbl_xx(result_chr = "none")},
                       cors_with_utl_tb = make_cors_with_utl_tbl(data_tb %>%
                                                                   dplyr::mutate(catch_all_round_chr = ds_descvs_ls$round_vals_chr[1]),
                                                                 ds_descvs_ls = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){ds_descvs_ls %>% purrr::list_modify(round_var_nm_1L_chr = "catch_all_round_chr")}else{ds_descvs_ls}),
                       ds_descvs_ls = ds_descvs_ls)
  descv_tbl_ls$predr_pars_and_cors_tb <- make_predr_pars_and_cors_tbl(data_tb, ## ERROR HERE
                                                                      ds_descvs_ls = ds_descvs_ls,
                                                                      descv_tbl_ls = descv_tbl_ls,
                                                                      dictionary_tb = ds_descvs_ls$dictionary_tb,
                                                                      nbr_of_digits_1L_int = nbr_of_digits_1L_int,
                                                                      predictors_lup = predictors_lup)
  saveRDS(descv_tbl_ls,paste0(descv_outp_dir_1L_chr,"/descv_tbls_ls.RDS"))
  return(descv_tbl_ls)
}
# From specific
author_SpecificModels <- function(x,
                                  prefd_mdl_types_chr = NULL,
                                  what_1L_chr = "all",
                                  digits_1L_int = 3L,
                                  reference_1L_int = NULL){
  if(what_1L_chr %in% c("all","descriptives","models","workspace")){
    session_data_ls <- sessionInfo()
    if(what_1L_chr %in% c("workspace","all")){
      if(!is.null(reference_1L_int)){
        transform_paths_ls <- list(fn = transform_paths_ls_for_scndry,
                                   args_ls = list(reference_1L_int = reference_1L_int)) # FOR SECONDARY
      }else{
        transform_paths_ls <- NULL
      }
      path_params_ls <- make_path_params_ls()
      path_params_ls$path_from_top_level_1L_chr <- x@paths_chr
      path_params_ls$use_fake_data_1L_lgl <- x@b_SpecificParameters@fake_1L_lgl
      paths_ls <- path_params_ls %>%
        ready4show::make_paths_ls(depth_1L_int = ifelse(is.null(transform_paths_ls),
                                                        1,
                                                        2))
      if(!is.null(transform_paths_ls)){
        paths_ls <- rlang::exec(transform_paths_ls$fn,
                                paths_ls,
                                !!!transform_paths_ls$args_ls)
      }
      paths_ls <- ready4show::write_all_outp_dirs(paths_ls = paths_ls)
      x@b_SpecificParameters@paths_ls <- paths_ls
    }
    if(what_1L_chr %in% c("descriptives","all")){
      ds_descvs_ls <- manufacture_SpecificProject(x, # Replace with manufacture when Exporting
                                  what_1L_chr = "ds_descvs_ls")
      descv_tbl_ls <- write_descv_tbls(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,#add youthvars:: when Exporting
                                                  ds_descvs_ls = ds_descvs_ls,
                                                  predictors_lup = x@b_SpecificParameters@predictors_lup,
                                                  descv_outp_dir_1L_chr = x@b_SpecificParameters@paths_ls$descv_outp_dir_1L_chr,
                                                  nbr_of_digits_1L_int = digits_1L_int,
                                                  participation_var_1L_chr = if(!series_1L_lgl){character(0)}else{x@a_YouthvarsProfile@participation_var_1L_chr})
      descv_plts_paths_ls <- youthvars::write_descv_plots(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                                          ds_descvs_ls = ds_descvs_ls,
                                                          descv_outp_dir_1L_chr = x@b_SpecificParameters@paths_ls$descv_outp_dir_1L_chr,
                                                          lbl_nms_chr = x@b_SpecificParameters@itm_labels_chr, # Should be domain labels
                                                          maui_domains_pfxs_1L_chr = hutils::longest_prefix(x@b_SpecificParameters@domain_labels_chr))
      
    }
    if(what_1L_chr %in% c("models","all")){
      x <- investigate(x)
      if(!is.null(prefd_mdl_types_chr))
        x <- renew(x,
                   new_val_xx = prefd_mdl_types_chr,
                   type_1L_chr = "results",
                   what_1L_chr = "prefd_mdls")
      x <- investigate(x)
      if(!is.null(prefd_covars_chr))
        x <- renew(x,
                   new_val_xx = prefd_covars_chr,
                   type_1L_chr = "results",
                   what_1L_chr = "prefd_covars")
      x <- investigate(x)
      x <- investigate(x)
      x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$session_data_ls <- session_data_ls
      # author(x,
      #        type_1L_chr = "purge_write")
    }
  }else{
    methods::callNextMethod()
  }
  return(x)
}
manufacture_SpecificProject <- function(x,
                                        what_1L_chr = "ds_descvs_ls",
                                        scndry_anlys_params_ls = NULL){
  series_1L_lgl <- x@a_YouthvarsProfile %>% inherits("YouthvarsSeries")
  if(what_1L_chr %in% c("ds_descvs_ls","ds_smry_ls","input_params_ls")){
    ds_descvs_ls <- make_ds_descvs_ls(candidate_predrs_chr =  x@b_SpecificParameters@candidate_predrs_chr,
                                      candidate_covar_nms_chr = x@b_SpecificParameters@candidate_covars_chr,
                                      cohort_descv_var_nms_chr = x@b_SpecificParameters@descv_var_nms_chr,
                                      dictionary_tb = x@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3,
                                      id_var_nm_1L_chr = x@a_YouthvarsProfile@id_var_nm_1L_chr,
                                      is_fake_1L_lgl = x@b_SpecificParameters@fake_1L_lgl,
                                      msrmnt_date_var_nm_1L_chr = if(!series_1L_lgl){character(0)}else{x@b_SpecificParameters@msrmnt_date_var_nm_1L_chr},
                                      round_var_nm_1L_chr = if(!series_1L_lgl){character(0)}else{x@a_YouthvarsProfile@timepoint_var_nm_1L_chr},
                                      round_vals_chr = if(!series_1L_lgl){"Overall"}else{x@a_YouthvarsProfile@timepoint_vals_chr},
                                      utl_wtd_var_nm_1L_chr = x@b_SpecificParameters@depnt_var_nm_1L_chr,
                                      maui_item_pfx_1L_chr = x@b_SpecificParameters@itm_prefix_1L_chr,
                                      utl_unwtd_var_nm_1L_chr = x@b_SpecificParameters@total_unwtd_var_nm_1L_chr)
    ds_descvs_ls$nbr_obs_in_raw_ds_1L_dbl <- nrow(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb)
    ds_descvs_ls$nbr_participants_1L_int <- length(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb %>%
                                                     dplyr::pull(ds_descvs_ls$id_var_nm_1L_chr) %>%
                                                     unique())
    object_xx <- ds_descvs_ls
  }
  if(what_1L_chr %in% c("ds_smry_ls","input_params_ls")){
    ds_smry_ls <- ds_descvs_ls %>%
      make_analysis_ds_smry_ls(candidate_covar_nms_chr = x@b_SpecificParameters@candidate_covars_chr,
                               predictors_lup = x@b_SpecificParameters@predictors_lup)
    object_xx <- ds_smry_ls
  }
  if(what_1L_chr %in% c("mdl_smry_ls","input_params_ls")){
    if(is.na(x@b_SpecificParameters@candidate_mdls_chr)){
      mdl_types_chr <- NULL
    }else{
      mdl_types_chr <- x@b_SpecificParameters@candidate_mdls_chr
    }
    if(is.na(x@b_SpecificParameters@candidate_mdl_pfcs_chr)){
      choose_from_pfx_chr <- NULL
    }else{
      choose_from_pfx_chr <- x@b_SpecificParameters@candidate_mdl_pfcs_chr
    }
    mdl_smry_ls <- make_mdl_smry_ls(mdl_types_lup = x@b_SpecificParameters@candidate_mdls_lup,
                                    mdl_types_chr = mdl_types_chr,
                                    choose_from_pfx_chr = choose_from_pfx_chr,
                                    folds_1L_int = x@b_SpecificParameters@folds_1L_int,
                                    max_nbr_of_boruta_mdl_runs_int = x@b_SpecificParameters@max_mdl_runs_1L_int)
    object_xx <- mdl_smry_ls
  }
  if(what_1L_chr == "input_params_ls"){
    y <- SpecificSynopsis()
    header_yaml_args_ls <- ready4show::make_header_yaml_args_ls(authors_tb = y@authors_r3,
                                                                institutes_tb = y@institutes_r3,
                                                                title_1L_chr = y@title_1L_chr,
                                                                keywords_chr = y@keywords_chr)
    maui_params_ls <- make_maui_params_ls(maui_domains_pfxs_1L_chr = x@b_SpecificParameters@itm_prefix_1L_chr,
                                          maui_itm_short_nms_chr = x@b_SpecificParameters@itm_labels_chr,
                                          maui_scoring_fn = NULL)
    output_format_ls <- ready4show::make_output_format_ls(manuscript_outp_1L_chr = y@outp_formats_chr[1],
                                                          manuscript_digits_1L_int = y@digits_int[1],
                                                          supplementary_outp_1L_chr = ifelse(length(y@outp_formats_chr)>1,y@outp_formats_chr[2],y@outp_formats_chr[1]),
                                                          supplementary_digits_1L_int = ifelse(length(y@digits_int)>1,y@digits_int[2],y@digits_int[1]))
    # scndry_anlys_params_ls <- make_scndry_anlys_params(candidate_predrs_chr = c("SOFAS"),
    #                                                    prefd_covars_chr = NA_character_)
    if(is.null(procure(x,
                       what = "prefd_covars"))){
      prefd_covars_chr <- NA_character_
    }else{
      prefd_covars_chr <- procure(x,
                                  what = "prefd_covars")
    }
    object_xx <- make_input_params(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                   control_ls = NULL,#
                                   ds_descvs_ls = ds_descvs_ls,
                                   dv_ds_nm_and_url_chr = NULL,#
                                   header_yaml_args_ls = header_yaml_args_ls,
                                   maui_params_ls = maui_params_ls,
                                   output_format_ls = output_format_ls,
                                   predictors_lup = x@b_SpecificParameters@predictors_lup,
                                   prefd_covars_chr = prefd_covars_chr,
                                   prefd_mdl_types_chr = procure(x,
                                                                 what = "prefd_mdls"),
                                   scndry_anlys_params_ls = scndry_anlys_params_ls,
                                   write_new_dir_1L_lgl = F)#
    
  }
  return(object_xx)
}
ratify_SpecificModels <- function(x,
                                  class_fn_1L_chr = "as.numeric",
                                  prototype_lup = NULL,
                                  scndry_anlys_params_ls = NULL){
  series_1L_lgl <- x@a_YouthvarsProfile %>% inherits("YouthvarsSeries")
  if(series_1L_lgl){
    x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb <- x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb %>%
      youthvars::add_interval_var(id_var_nm_1L_chr = x@a_YouthvarsProfile@id_var_nm_1L_chr,
                                  msrmnt_date_var_nm_1L_chr = ifelse(!is.na(x@b_SpecificParameters@msrmnt_date_var_nm_1L_chr),
                                                                     x@b_SpecificParameters@msrmnt_date_var_nm_1L_chr,
                                                                     "d_interview_date")) %>%
      youthvars::add_participation_var(id_var_nm_1L_chr = x@a_YouthvarsProfile@id_var_nm_1L_chr,
                                       fup_round_nbr_1L_int = length( x@a_YouthvarsProfile@timepoint_vals_chr))
    if(is.na(x@a_YouthvarsProfile@participation_var_1L_chr)){
      x@a_YouthvarsProfile@participation_var_1L_chr <- "participation"
    }
    x@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3 <- x@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3 %>%
      ready4::renew(var_nm_chr = c("bl_date_dtm",
                                   "interval_dbl",
                                   x@a_YouthvarsProfile@participation_var_1L_chr),
                    var_ctg_chr = c("Temporal",
                                    "Temporal",
                                    "Temporal"),
                    var_desc_chr = c("Date of baseline assessment",
                                     "Interval between baseline and follow-up assessments",
                                     "Rounds participated in"),
                    var_type_chr = c("date",
                                     "interval",
                                     "character"))
  }
  x@a_YouthvarsProfile@a_Ready4useDyad <- renew_Ready4useDyad(x@a_YouthvarsProfile@a_Ready4useDyad, type_1L_chr = "label") # DELETE WHEN EXPORTING
  # x <- renewSlot(x, # ADD BACK IN WHEN EXPORTING
  #                "a_YouthvarsProfile@a_Ready4useDyad",
  #                type_1L_chr = "label")
  x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb <- x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb %>%
    transform_mdl_vars_with_clss(predictors_lup = x@b_SpecificParameters@predictors_lup,
                                 depnt_var_nm_1L_chr = x@b_SpecificParameters@depnt_var_nm_1L_chr,
                                 prototype_lup = prototype_lup,
                                 class_fn_1L_chr = class_fn_1L_chr)
  ## Pick up here
  input_params_ls <- manufacture_SpecificProject(x, #### REPLACE WITH manufacture when EXPORTING
                                 what_1L_chr = "input_params_ls",
                                 scndry_anlys_params_ls = scndry_anlys_params_ls)
  x <- renewSlot(x,
                 "a_YouthvarsProfile@a_Ready4useDyad@ds_tb",
                 input_params_ls$params_ls$ds_tb) %>%
    renewSlot("a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3",
              input_params_ls$params_ls$ds_descvs_ls$dictionary_tb) %>%
    renewSlot("b_SpecificParameters@candidate_covars_chr",
              input_params_ls$params_ls$candidate_covar_nms_chr) %>%
    renewSlot("b_SpecificParameters@descv_var_nms_chr",
              input_params_ls$params_ls$ds_descvs_ls$cohort_descv_var_nms_chr) %>%
    renewSlot("b_SpecificParameters@candidate_predrs_chr",
              input_params_ls$params_ls$ds_descvs_ls$candidate_predrs_chr) %>%
    renewSlot("b_SpecificParameters@predictors_lup",
              input_params_ls$params_ls$predictors_lup) %>%
    renewSlot("b_SpecificParameters@depnt_var_nm_1L_chr",
              input_params_ls$params_ls$ds_descvs_ls$utl_wtd_var_nm_1L_chr) %>%
    renewSlot("c_SpecificResults@a_SpecificShareable@shareable_outp_ls",
              list(rename_lup = input_params_ls$rename_lup))
  #X@b_SpecificParameters@depnt_var_nm_1L_chr
  #X@b_SpecificParameters@predictors_lup
  return(x)
  
}
#### FROM THIS POINT ON IS TEMPORARY AND NOT FOR EXPORT
renew_ScorzAqol6Adol <- function(x, # DO NOT EXPORT TO SCORZ LIBRARY
                                 label_ds_1L_lgl = T,
                                 type_1L_chr = "score"){
  y <- x@a_YouthvarsProfile@a_Ready4useDyad
  y <- renew(y, type_1L_chr = "unlabel")
  if(type_1L_chr == "score"){
    if(identical(x@scrg_dss_ls,list(list()))){
      x@scrg_dss_ls <- get_aqol6d_scrg_dss()
    }
    select_chr <- setdiff(names(y@ds_tb),
                          x@instrument_dict_r3$var_nm_chr[!x@instrument_dict_r3$var_nm_chr %>%
                                                            startsWith(x@itm_prefix_1L_chr)] %>%
                            as.vector())
    y@ds_tb <- y@ds_tb %>%
      dplyr::select(tidyselect::all_of(select_chr))
    y@ds_tb <- add_adol6d_scores(y@ds_tb, # Make an enhanceSlot method - then generalise renew mthd to parent class
                                 aqol6d_scrg_dss_ls = x@scrg_dss_ls,
                                 prefix_1L_chr =  x@itm_prefix_1L_chr,
                                 id_var_nm_1L_chr = x@a_YouthvarsProfile@id_var_nm_1L_chr,
                                 total_aqol_var_nm_1L_chr = x@total_unwtd_var_nm_1L_chr,
                                 wtd_aqol_var_nm_1L_chr = x@total_wtd_var_nm_1L_chr)
    y@dictionary_r3 <- ready4::renew(y@dictionary_r3,
                                     new_cases_r3 = x@instrument_dict_r3)
  }
  if(label_ds_1L_lgl)
    y <- renew_Ready4useDyad(y)
  if(type_1L_chr == "score"){
    x@a_YouthvarsProfile@a_Ready4useDyad <- y
  }
  return(x)
}
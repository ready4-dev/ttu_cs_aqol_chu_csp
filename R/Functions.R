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
                              mdl_scaling_dbl = 0.1,
                              covariate_lgl = F) %>%
    dplyr::mutate(covariate_lgl = dplyr::case_when(short_name_chr == "SOFAS" ~ F,
                                                   T ~ covariate_lgl))
  return(predictors_r3)
}
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

#### FROM THIS POINT ON IS TEMPORARY
renew_ScorzAqol6Adol <- function(x, # DO NOT ADD TO SCORZ LIBRARY
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
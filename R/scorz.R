manufacture_ScorzProfile <- function(x,
                                     custom_args_ls = NULL,
                                     custom_fn = NULL,#identity,
                                     domain_pfxs_1L_chr = character(0),
                                     what_1L_chr = "domains",
                                     ...){
  object_xx <- NULL
  if(what_1L_chr == "domains"){
    if(is.null(custom_fn)){
      if(!identical(domain_pfxs_1L_chr, character(0))){
        object_xx <- purrr::map_chr(names(dplyr::select(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb, dplyr::starts_with(domain_pfxs_1L_chr))),
                                    ~ {
                                      domain_1L_chr <- eval(parse(text=paste0("attributes(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb$",.x,")$label")))
                                      domain_1L_chr <- ifelse(!stringr::str_detect(domain_1L_chr," - "),
                                                              domain_1L_chr,
                                                              stringr::str_sub(domain_1L_chr,
                                                                               start = stringi::stri_locate_last_fixed(domain_1L_chr," - ")[1,1] %>%
                                                                                 unname() + 2))
                                      domain_1L_chr
                                    })
        
        
      }
    }else{
      object_xx <-  rlang::exec(custom_fn, !!!custom_args_ls)
    }
  }
  return(object_xx)
}
manufacture_ScorzAqol6 <- function(x,
                                   what_1L_chr = "domains",
                                   ...){
  object_xx <- NULL
  if(what_1L_chr == "domains"){
    object_xx <- c("Independent Living", "Relationships", "Mental Health", "Coping", "Pain", "Senses")
  }else{
    object_xx <- methods::callNextMethod()
  }
  return(object_xx)
}
manufacture_ScorzEuroQol5 <- function(x,
                                      what_1L_chr = "domains",
                                      ...){
  object_xx <- NULL
  if(what_1L_chr == "domains"){
    object_xx <- x@itm_labels_chr
  }else{
    object_xx <- methods::callNextMethod()
  }
  return(object_xx)
}
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
    y <- renew_Ready4useDyad(y) # change to renew
  if(type_1L_chr == "score"){
    x@a_YouthvarsProfile@a_Ready4useDyad <- y
  }
  return(x)
}
renew_ScorzProfile <- function(x,
                               drop_msng_1L_lgl = F,
                               item_type_1L_chr = "numeric",
                               scoring_fn = identity,
                               scorz_args_ls = NULL,
                               label_ds_1L_lgl = T,
                               type_1L_chr = "score"){
  
  if(type_1L_chr %in% c("score","score-c","score-w")){
    if(type_1L_chr == "score"){
      x <-  renew_ScorzProfile(x,# Change to renew
                               scoring_fn = scoring_fn,
                               scorz_args_ls = scorz_args_ls,
                               label_ds_1L_lgl = label_ds_1L_lgl,
                               type_1L_chr = "score-w") %>% 
        renew_ScorzProfile(label_ds_1L_lgl = label_ds_1L_lgl,# Change to renew
                           type_1L_chr = "score-c") 
    }
    if(type_1L_chr %in% c("score-c","score-w")){
      y <- x@a_YouthvarsProfile@a_Ready4useDyad
      y <- renew(y, type_1L_chr = "unlabel")
    }
    if(type_1L_chr == "score-c"){
      y@ds_tb <- y@ds_tb %>%
        dplyr::mutate(`:=`(!!rlang::sym(x@total_unwtd_var_nm_1L_chr),
                           rowSums(dplyr::select(.,dplyr::starts_with(x@itm_prefix_1L_chr)))))
      if(drop_msng_1L_lgl)
        y@ds_tb <- y@ds_tb  %>%
          dplyr::filter(!is.na(!!rlang::sym(x@total_unwtd_var_nm_1L_chr)))
      if(!x@total_unwtd_var_nm_1L_chr %in% x@instrument_dict_r3$var_nm_chr){
        x@instrument_dict_r3 <- ready4use::renew.ready4use_dictionary(x@instrument_dict_r3,
                                                                      var_nm_chr = x@total_unwtd_var_nm_1L_chr,
                                                                      var_ctg_chr = "multi-attribute utility instrument unweighted total score",
                                                                      var_desc_chr = paste0(x@instrument_nm_1L_chr, " (unweighted total)"),
                                                                      var_type_chr = item_type_1L_chr)
      }
    }
    if(type_1L_chr == "score-w"){
      y@ds_tb <- rlang::exec(scoring_fn, y@ds_tb, !!!scorz_args_ls)
      if(!x@total_wtd_var_nm_1L_chr %in% x@instrument_dict_r3$var_nm_chr){
        x@instrument_dict_r3 <- ready4use::renew.ready4use_dictionary(x@instrument_dict_r3,
                                                                      var_nm_chr = x@total_wtd_var_nm_1L_chr,
                                                                      var_ctg_chr = "health utility",
                                                                      var_desc_chr = paste0(x@instrument_nm_1L_chr, " total score"),
                                                                      var_type_chr = item_type_1L_chr)
      }
    }
  }
  if(type_1L_chr %in% c("score-c","score-w")){
    y@dictionary_r3 <- ready4use::renew.ready4use_dictionary(y@dictionary_r3,
                                                             new_cases_r3 = x@instrument_dict_r3)
    if(label_ds_1L_lgl)
      y <- renew_Ready4useDyad(y) # change to renew
    x@a_YouthvarsProfile@a_Ready4useDyad <- y
  }
  return(x)
}
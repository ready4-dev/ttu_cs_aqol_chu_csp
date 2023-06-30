author_TTUProject <- function(x,
                              consent_1L_chr = "",
                              custom_args_ls = NULL,
                              custom_fn = NULL,
                              depnt_var_min_val_1L_dbl = numeric(0),
                              digits_1L_int = 2L,
                              download_tmpl_1L_lgl = T,
                              fl_nm_1L_chr = "TTUProject",
                              items_as_domains_1L_lgl = F,
                              supplement_fl_nm_1L_chr = "TA_PDF",
                              timepoint_new_nms_chr = NA_character_,
                              type_1L_chr = "auto",
                              what_1L_chr = "default",
                              ...){
  if(what_1L_chr %in% c("catalogue","Catalogue","dependencies", "Dependencies", "descriptives", "Descriptives", "manuscript", "Manuscript", "models", "Models",
                        "plots", "Plots", "purge", "Purge", "self", "Self",  "supplement", "Supplement")){
    if(what_1L_chr %in% c("self", "Self")){
      to_1L_chr <- paste0(x@c_SpecificProject@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                          "/",
                          fl_nm_1L_chr,
                          ".RDS")
      ready4::write_with_consent(consented_fn = saveRDS,
                                 prompt_1L_chr = paste0("Do you confirm that you want to write a copy of this TTUProject module to ",
                                                        to_1L_chr,
                                                        "?"),
                                 consent_1L_chr = consent_1L_chr,
                                 consented_args_ls = list(object = x,
                                                          file = to_1L_chr),
                                 consented_msg_1L_chr = paste0("A copy of this TTUProject module has been written to ",
                                                               to_1L_chr,
                                                               "."),
                                 declined_msg_1L_chr = "Write request cancelled - no new file has been written.")
    }
    
    if(what_1L_chr %in% c("catalogue","Catalogue")){
      authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, download_tmpl_1L_lgl = download_tmpl_1L_lgl, what_1L_chr = Hmisc::capitalize(what_1L_chr))
    }
    if(what_1L_chr %in% c("descriptives","Descriptives")){
      if(items_as_domains_1L_lgl == T){
        x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains",
                                    custom_args_ls = list(string = x@b_SpecificParameters@itm_labels_chr), custom_fn = Hmisc::capitalize)
      }else{
        x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains")
      }
      x <- renewSlot(x, "c_SpecificProject",
                     authorSlot(x, "c_SpecificProject",
                                consent_1L_chr = consent_1L_chr,
                                digits_1L_int = digits_1L_int,
                                what_1L_chr = tolower(what_1L_chr),
                                x_labels_chr = x_labels_chr))
    }
    if(what_1L_chr %in% c("manuscript", "Manuscript")){
      if(type_1L_chr=="auto"){
        x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0(Hmisc::capitalize(what_1L_chr),"_",Hmisc::capitalize(type_1L_chr))
        authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, download_tmpl_1L_lgl = download_tmpl_1L_lgl, type_1L_chr = "Report", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
      }
      if(type_1L_chr == "copy"){
        from_1L_chr <- paste0(x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr,
                              "/",
                              x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr,
                              "/Manuscript_Auto")
        to_1L_chr <- paste0(x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr,
                            "/",
                            x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr,
                            "/Manuscript_Submission")
        ready4::write_with_consent(consented_fn = R.utils::copyDirectory,
                                   prompt_1L_chr = paste0("Do you confirm that you want to copy the directory ",
                                                          from_1L_chr, #"packages.RDS",
                                                          " (and all its contents) to ",
                                                          to_1L_chr,
                                                          "?"),
                                   consent_1L_chr = consent_1L_chr,
                                   consented_args_ls = list(from = from_1L_chr,
                                                            to = to_1L_chr),
                                   consented_msg_1L_chr = paste0("The directory ",
                                                                 from_1L_chr,
                                                                 " has been copied to ",
                                                                 to_1L_chr,
                                                                 "."),
                                   declined_msg_1L_chr = "Write request cancelled - no new directory copy has been written.")
        
      }
      if(type_1L_chr %in% c("dependencies", "Dependencies")){
        author(x@d_TTUReports, consent_1L_chr = consent_1L_chr, type_1L_chr = Hmisc::capitalize(type_1L_chr), what_1L_chr = "Manuscript_Submission")
      }
      if(type_1L_chr %in% c("plots","Plots")){
        if(items_as_domains_1L_lgl == T){
          x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains",
                                      custom_args_ls = list(string = x@b_SpecificParameters@itm_labels_chr), custom_fn = Hmisc::capitalize)
        }else{
          x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains")
        }
        authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, depnt_var_desc_1L_chr = x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls$health_utl_nm_1L_chr,
                   depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, timepoint_new_nms_chr = timepoint_new_nms_chr, type_1L_chr = Hmisc::capitalize(type_1L_chr), what_1L_chr = "Manuscript_Submission", x_labels_chr = x_labels_chr)
      }
      if(type_1L_chr=="submission"){
        x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0(Hmisc::capitalize(what_1L_chr),"_",Hmisc::capitalize(type_1L_chr))
        authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr, download_tmpl_1L_lgl = download_tmpl_1L_lgl, type_1L_chr = "Report", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
      }
    }
    if(what_1L_chr %in% c("models", "Models")){
      x <- renewSlot(x, "c_SpecificProject",
                     authorData(procureSlot(x, "c_SpecificProject"), consent_1L_chr = consent_1L_chr))
    }
    if(what_1L_chr %in% c("purge")){
      authorSlot(x,"c_SpecificProject", type_1L_chr = "purge_write", consent_1L_chr = consent_1L_chr)
    }
    if(what_1L_chr %in% c("plots","Plots")){
      if(items_as_domains_1L_lgl == T){
        x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains",
                                    custom_args_ls = list(string = x@b_SpecificParameters@itm_labels_chr), custom_fn = Hmisc::capitalize)
      }else{
        x_labels_chr <- manufacture(x@a_ScorzProfile, what_1L_chr = "domains")
      }
      authorSlot(x, "d_TTUReports", consent_1L_chr = consent_1L_chr,
                 depnt_var_desc_1L_chr = x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls$health_utl_nm_1L_chr,
                 type_1L_chr = Hmisc::capitalize(what_1L_chr),
                 x_labels_chr = x_labels_chr)
      
      # authorSlot(x, "c_SpecificProject",
      #            consent_1L_chr = consent_1L_chr,
      #            digits_1L_int = digits_1L_int,
      #            what_1L_chr = tolower(what_1L_chr),
      #            x_labels_chr = x_labels_chr)
    }
    if(what_1L_chr %in% c("supplement","Supplement")){
      x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr <- paste0("Manuscript_",Hmisc::capitalize(type_1L_chr))
      authorReport(procureSlot(x, "d_TTUReports") %>%
                     renewSlot("a_TTUSynopsis@rmd_fl_nms_ls", ready4show::make_rmd_fl_nms_ls(pdf_fl_nm_1L_chr = supplement_fl_nm_1L_chr)) %>%
                     renewSlot("a_TTUSynopsis@outp_formats_chr", rep(x@d_TTUReports@a_TTUSynopsis@outp_formats_chr[2],2)) %>%
                     procureSlot("a_TTUSynopsis"),
                   consent_1L_chr = consent_1L_chr, fl_nm_1L_chr = "Supplement", what_1L_chr = x@d_TTUReports@a_TTUSynopsis@a_Ready4showPaths@ms_dir_1L_chr)
    }
  }else{
    x <- methods::callNextMethod()
  }
  return(x)
}
author_TTUReports <- function(x,
                              args_ls = NULL,
                              consent_1L_chr = "",
                              depnt_var_desc_1L_chr = NA_character_,
                              depnt_var_min_val_1L_dbl = numeric(0),
                              download_tmpl_1L_lgl = T,
                              fl_type_1L_chr = ".eps",
                              timepoint_new_nms_chr = NA_character_,
                              type_1L_chr = "Report",
                              what_1L_chr = NA_character_,
                              x_labels_chr = character(0),
                              ...){
  if(type_1L_chr == "Report"){
    if(download_tmpl_1L_lgl){
      authorData(x@a_TTUSynopsis,
                 consent_1L_chr = consent_1L_chr,
                 tmpl_url_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                          x@catalogue_tmpl_chr[1],
                                          x@manuscript_tmpl_chr[1]),
                 tmpl_version_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                              x@catalogue_tmpl_chr[2],
                                              x@manuscript_tmpl_chr[2]),
                 what_1L_chr = what_1L_chr)
    }
    if(what_1L_chr == "Catalogue"){
      x@a_TTUSynopsis@rmd_fl_nms_ls <- x@catalogue_fl_nms_ls
    }else{
      x@a_TTUSynopsis@rmd_fl_nms_ls <- x@manuscript_fl_nms_ls
    }
    if(what_1L_chr == "Catalogue"){
      author(x@a_TTUSynopsis,
             args_ls = args_ls,
             consent_1L_chr = consent_1L_chr,
             type_1L_chr = type_1L_chr,
             what_1L_chr = what_1L_chr)
    }else{
      authorReport(x@a_TTUSynopsis,
                   args_ls = args_ls,
                   consent_1L_chr = consent_1L_chr,
                   type_1L_chr = type_1L_chr,
                   what_1L_chr = what_1L_chr)
    }
  }else{
    dir_1L_chr <- paste0(x@a_TTUSynopsis@a_Ready4showPaths@outp_data_dir_1L_chr,
                         "/",
                         x@a_TTUSynopsis@a_Ready4showPaths@mkdn_data_dir_1L_chr,
                         "/",
                         what_1L_chr)
    if(type_1L_chr == "Dependencies"){
      df <- data.frame(Package = c("youthvars","scorz","specific","TTU") %>%
                         purrr::map(~ {
                           desc_ls <- utils::packageDescription(.x)
                           desc_ls[c("Depends", "Imports")] %>%
                             # `[`(c("Depends", "Imports")) %>%
                             purrr::map(~{
                               if(is.null(.x)){
                                 character(0)
                               }else{
                                 .x %>%
                                   strsplit(",\\n") %>%
                                   purrr::flatten_chr() %>%
                                   purrr::map(~strsplit(.x,", ") %>%
                                                purrr::flatten_chr()) %>%
                                   purrr::flatten_chr() %>% sort() %>%
                                   purrr::discard(~startsWith(.x,"R "))
                               }
                             }) %>%
                             purrr::flatten_chr() %>%
                             unique() %>%
                             sort()
                         }) %>%
                         purrr::reduce(~c(.x,.y)) %>%
                         purrr::map_chr(~{
                           updated_1L_chr <- stringr::str_replace_all(.x,"\\n"," ")
                           problem_idx_1L_chr <- stringr::str_locate(updated_1L_chr," ")[1,1] %>%
                             unname()
                           if(!is.na(problem_idx_1L_chr))
                             updated_1L_chr <- updated_1L_chr %>%
                             stringr::str_sub(end = problem_idx_1L_chr-1)
                           updated_1L_chr %>% trimws(which = "left")
                         }) %>% unique() %>% sort())
      df <- df %>%
        dplyr::mutate(Version = Package %>%
                        purrr::map_chr(~utils::packageDescription(.x) %>%
                                         purrr::pluck("Version")),
                      Citation = Package %>%
                        purrr::map_chr(~get_pkg_citation(.x)))
      ready4::write_with_consent(consented_fn = saveRDS,
                                 prompt_1L_chr = paste0("Do you confirm that you want to write the file ",
                                                        "packages.RDS",
                                                        " to ",
                                                        dir_1L_chr,
                                                        "?"),
                                 consent_1L_chr = consent_1L_chr,
                                 consented_args_ls = list(object = df,
                                                          file = paste0(dir_1L_chr, "/packages.RDS")),
                                 consented_msg_1L_chr = paste0("File ",
                                                               "packages.RDS",
                                                               " has been written to ",
                                                               dir_1L_chr,
                                                               "."),
                                 declined_msg_1L_chr = "Write request cancelled - no new files have been written.")
      
    }
    if(type_1L_chr == "Plots"){
      composite_1_plt <- depict(x@a_TTUSynopsis,#
                                consent_1L_chr = consent_1L_chr,
                                depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                timepoint_old_nms_chr = procureSlot(x,
                                                                    "a_TTUSynopsis@d_YouthvarsProfile@timepoint_vals_chr"),
                                timepoint_new_nms_chr = timepoint_new_nms_chr,
                                what_1L_chr = "composite_mdl",
                                write_1L_lgl = T)
      composite_2_plt <- depict(x@a_TTUSynopsis,#
                                consent_1L_chr = consent_1L_chr,
                                what_1L_chr = "composite_utl",
                                write_1L_lgl = T,
                                x_labels_chr = x_labels_chr)
      if(!is.na(what_1L_chr)){
        consented_fn <- function(composite_1_plt,
                                 composite_2_plt,
                                 dir_1L_chr,
                                 fl_type_1L_chr){
          ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig1", fl_type_1L_chr),
                          composite_2_plt)
          ggplot2::ggsave(file = paste0(dir_1L_chr, "/fig2", fl_type_1L_chr),
                          composite_1_plt)
        }
        ready4::write_with_consent(consented_fn = consented_fn,
                                   prompt_1L_chr = paste0("Do you confirm that you want to write the files ",
                                                          ready4::make_list_phrase(paste0("fig",1:2,fl_type_1L_chr)), #"packages.RDS",
                                                          " to ",
                                                          dir_1L_chr,
                                                          "?"),
                                   consent_1L_chr = consent_1L_chr,
                                   consented_args_ls = list(composite_1_plt = composite_1_plt,
                                                            composite_2_plt = composite_2_plt,
                                                            dir_1L_chr = dir_1L_chr,
                                                            fl_type_1L_chr = fl_type_1L_chr),
                                   consented_msg_1L_chr = paste0("Files ",
                                                                 ready4::make_list_phrase(paste0("fig",1:2,fl_type_1L_chr)),
                                                                 " have been written to ",
                                                                 dir_1L_chr,
                                                                 "."),
                                   declined_msg_1L_chr = "Write request cancelled - no new files have been written.")
        
      }
      
    }
  }
}
exhibit_TTUProject <- function(x,
                               captions_chr = NULL, 
                               display_1L_chr = "all",
                               header_1L_chr = "",
                               header_col_nms_chr = " ",
                               mkdn_tbl_refs_chr = NULL, 
                               profile_idx_int = NA_integer_,
                               output_type_1L_chr = "HTML",
                               type_1L_chr = "default",
                               use_lbls_as_col_nms_1L_lgl = T,
                               use_rdocx_1L_lgl = F,
                               variables_chr = character(0),
                               what_1L_chr = "predictors",
                               ...){
  if(what_1L_chr == "predictors"){
    exhibitSlot(x, "b_SpecificParameters@predictors_lup",
                ... = ...)
  }
  if(what_1L_chr == "profile"){
    exhibit(procure(x, variables_chr = variables_chr, type_1L_chr == "default", what_1L_chr = what_1L_chr, ...=...),
            captions_chr = captions_chr,
            header_1L_chr = header_1L_chr,
            header_col_nms_chr = header_col_nms_chr,
            mkdn_tbl_refs_chr = mkdn_tbl_refs_chr,
            profile_idx_int = profile_idx_int,
            output_type_1L_chr = output_type_1L_chr,
            what_1L_chr = type_1L_chr) #descriptives
  }
  if(what_1L_chr == "records"){
    if(type_1L_chr %in% c("ds","dict")){
      exhibit(x@a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad,
              caption_1L_chr = {if(is.null(captions_chr)){NA_character_}else{captions_chr[1]}},
              display_1L_chr = display_1L_chr,
              mkdn_tbl_ref_1L_chr = {if(is.null(mkdn_tbl_refs_chr)){""}else{mkdn_tbl_refs_chr[1]}},
              output_type_1L_chr = "HTML",
              type_1L_chr = type_1L_chr,
              use_lbls_as_col_nms_1L_lgl = T,
              use_rdocx_1L_lgl = F,
              ... = ...)
    }
  }
}
manufacture_TTUProject <- function(x,
                                   type_1L_chr = "dummys", 
                                   what_1L_chr = "factors",
                                   ...){
  object_xx <- manufacture(x@c_SpecificProject@a_YouthvarsProfile@a_Ready4useDyad, type_1L_chr = type_1L_chr, what_1L_chr = what_1L_chr, 
                           restrict_to_chr = x@c_SpecificProject@b_SpecificParameters@candidate_covars_chr) 
  return(object_xx)
  }
procure_TTUProject <- function(x,
                               type_1L_chr = "default",
                               variables_chr = character(0),
                               what_1L_chr = "records",
                               ...){
  if(what_1L_chr == "parameters"){
    if(type_1L_chr == "models_lup"){
      object_xx <- x@b_SpecificParameters@candidate_mdls_lup
    }
  }
  if(what_1L_chr %in% c("project")){
    if(type_1L_chr == "models"){
      object_xx <- procureSlot(x,"c_SpecificProject", use_procure_mthd_1L_lgl = T, what_1L_chr = "prefd_mdls")
    }
  }
  if(what_1L_chr %in% c("profile", "records")){
    object_xx <- x@a_ScorzProfile@a_YouthvarsProfile
    if(!identical(variables_chr, character(0))){
      object_xx@a_Ready4useDyad@ds_tb <- object_xx@a_Ready4useDyad@ds_tb %>% dplyr::select(variables_chr)
    }
    if(what_1L_chr == "records"){
      if(type_1L_chr %in% c("default")){
        object_xx <- object_xx@a_Ready4useDyad@ds_tb 
      }
    }
  }
  return(object_xx)
}
renew_TTUProject <- function(x,
                             new_val_xx = NULL,
                             consent_1L_chr = "",
                             depnt_var_min_val_1L_dbl = numeric(0),
                             fl_nm_1L_chr = character(0),
                             paths_chr = character(0),
                             type_1L_chr = "default",
                             y_Ready4useRepos = ready4use::Ready4useRepos(),
                             what_1L_chr = "utility",
                             ...){
  if(what_1L_chr == "parameters"){
    if(type_1L_chr=="default"){
      x <- renewSlot(x, "b_SpecificParameters", SpecificConverter(a_ScorzProfile = x@a_ScorzProfile) %>% metamorphose() %>% procureSlot("b_SpecificParameters"))
    }
    if(type_1L_chr == "range"){
      x <- renewSlot(x, "b_SpecificParameters@depnt_var_min_max_dbl", new_val_xx)
    }
    if(type_1L_chr=="predictors_lup"){
      if(identical(new_val_xx, "use_renew_mthd")){
        predictors_lup <- y_Ready4useRepos %>%
          ingest(fls_to_ingest_chr = c(fl_nm_1L_chr), metadata_1L_lgl = F)
      }else{
        predictors_lup <- new_val_xx
      }
      x <- renewSlot(x, "b_SpecificParameters@predictors_lup", predictors_lup)
    }
    if(type_1L_chr == "predictors_vars"){
      x <- renewSlot(x, "b_SpecificParameters@candidate_predrs_chr", new_val_xx)
    }
    if(type_1L_chr == "covariates"){
      x <- renewSlot(x, "b_SpecificParameters@candidate_covars_chr", new_val_xx)
    }
    if(type_1L_chr == "descriptives"){
      x <- renewSlot(x, "b_SpecificParameters@descv_var_nms_chr", new_val_xx)
    }
    if(type_1L_chr == "is_fake"){
      x <- renewSlot(x, "b_SpecificParameters@fake_1L_lgl", new_val_xx)
    }
    if(type_1L_chr == "temporal"){
      x <- renewSlot(x, "b_SpecificParameters@msrmnt_date_var_nm_1L_chr", new_val_xx)
    }
  }
  if(what_1L_chr == "project"){
    if(type_1L_chr == "default"){
      x <- renewSlot(x, "c_SpecificProject", SpecificModels(a_YouthvarsProfile = x@a_ScorzProfile@a_YouthvarsProfile,
                                                            b_SpecificParameters = x@b_SpecificParameters, paths_chr = paths_chr))
      x <- ratifySlot(x, "c_SpecificProject")
      x <- renewSlot(x, "c_SpecificProject", authorSlot(x, "c_SpecificProject", consent_1L_chr = consent_1L_chr, what_1L_chr = "workspace")) 
    }
    if(type_1L_chr == "dummys"){
      x <- renewSlot(x, "c_SpecificProject", renew(x@c_SpecificProject, new_val_xx, what_1L_chr = type_1L_chr))
    }
  }
  if(what_1L_chr == "records"){
    if(type_1L_chr == "ds"){
      x <- renewSlot(x, "a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@ds_tb", new_val_xx)
    }
    if(type_1L_chr == "dict"){
      x <- renewSlot(x, "a_ScorzProfile@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3", new_val_xx)
    }
    
  }
  if(what_1L_chr == "reporting"){
    if(type_1L_chr == "default"){
      Y <- metamorphoseSlot(x, "c_SpecificProject")
      Y <- TTUSynopsis(a_Ready4showPaths = Y@a_Ready4showPaths, b_SpecificResults = Y@b_SpecificResults, c_SpecificParameters = Y@c_SpecificParameters,
                       d_YouthvarsProfile = Y@d_YouthvarsProfile, rmd_fl_nms_ls = Y@rmd_fl_nms_ls)
      Y <- TTUReports(a_TTUSynopsis = Y)
      x <- renewSlot(x, "d_TTUReports", Y)
    }
    if(type_1L_chr == "abstract"){
      if(identical(new_val_xx,"use_renew_mthd")){
        descs_ls <- x@d_TTUReports@a_TTUSynopsis@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$results_ls$study_descs_ls
        x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis@abstract_args_ls",
                       manufactureSlot(x,"d_TTUReports@a_TTUSynopsis", what_1L_chr = "abstract_args_ls", depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                       depnt_var_nms_chr = c(descs_ls$health_utl_nm_1L_chr,descs_ls$health_utl_long_nm_1L_chr)))
        
      }else{
        x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("abstract_args_ls", new_val_xx))
      }
    }
    if(type_1L_chr == "authors"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("authors_r3", new_val_xx))
    }
    if(type_1L_chr == "background"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("background_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "changes"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("correspondences_r3", new_val_xx))
    }
    if(type_1L_chr == "conflicts"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("coi_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "conclusion"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("conclusion_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "digits"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("digits_int", new_val_xx))
    }
    if(type_1L_chr == "ethics"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("ethics_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "formats"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("outp_formats_chr", new_val_xx))
    }
    if(type_1L_chr == "figures-body"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("figures_in_body_lgl", new_val_xx))
    }
    if(type_1L_chr == "funding"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("funding_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "institutes"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("institutes_r3", new_val_xx))
    }
    if(type_1L_chr == "interval"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("interval_chr", new_val_xx))
    }
    if(type_1L_chr == "keywords"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("keywords_chr", new_val_xx))
    }
    if(type_1L_chr == "naming"){
      x <- enhanceSlot(x, "d_TTUReports@a_TTUSynopsis", depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, depnt_var_nms_chr = new_val_xx, with_1L_chr = "results_ls")
    }
    if(type_1L_chr == "repos"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("e_Ready4useRepos", new_val_xx))
    }
    if(type_1L_chr == "sample"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("sample_desc_1L_chr", new_val_xx))
    }
    if(type_1L_chr == "tables-body"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("tables_in_body_lgl", new_val_xx))
    }
    if(type_1L_chr == "template-catalaogue"){
      x <- renewSlot(x, "d_TTUReports", procureSlot(x, "d_TTUReports") %>% renewSlot("catalogue_tmpl_chr", new_val_xx))
    }
    if(type_1L_chr == "template-manuscript"){
      x <- renewSlot(x, "d_TTUReports", procureSlot(x, "d_TTUReports") %>% renewSlot("manuscript_tmpl_chr", new_val_xx))
    }
    if(type_1L_chr == "title"){
      x <- renewSlot(x, "d_TTUReports@a_TTUSynopsis", procureSlot(x, "d_TTUReports@a_TTUSynopsis") %>% renewSlot("title_1L_chr", new_val_xx))
    }
  }
  if(what_1L_chr == "results"){
    if(type_1L_chr =="covariates"){
      x <- renewSlot(x, "c_SpecificProject", renew(procureSlot(x, "c_SpecificProject"), new_val_xx = new_val_xx, type_1L_chr = "results", what_1L_chr = "prefd_covars"))
    }
    if(type_1L_chr == "models"){
      x <- renewSlot(x, "c_SpecificProject", renew(procureSlot(x, "c_SpecificProject"), new_val_xx = new_val_xx, type_1L_chr = "results", what_1L_chr = "prefd_mdls"))
    }
    
  }
  if(what_1L_chr == "utility"){
    x <- renewSlot(x, "a_ScorzProfile")
  }
  return(x)
}
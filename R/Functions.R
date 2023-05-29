## Rename _baseine to _scaled or no extension
## Editing: 
##    
## fit_ts_model_with_brm 
## add_utility_predn_to_ds make_mdl_coef_ratio_ls make_scaling_text  make_mdl_desc_lines  make_shareable_mdl make_smry_of_brm_mdl transform_nms_in_mdl_tbl transform_tb_to_mdl_inp transform_ts_mdl_data write_shareable_mdls
# Add to ctlg: make_predr_vars_nms_ls
## From ready4show
manufacture.ready4show_correspondences <- function(x,
                                                   data_ls,
                                                   flatten_1L_chr = F,
                                                   type_1L_chr = "new",
                                                   what_1L_chr = "names",
                                                   ...){
  if(what_1L_chr == "names"){
    object_xx <- data_ls %>% purrr::map(~.x %>%
                                          purrr::map_chr(~ifelse(.x %in% x$old_nms_chr,
                                                                 ready4::get_from_lup_obj(x,
                                                                                          match_value_xx = .x,
                                                                                          match_var_nm_1L_chr = ifelse(type_1L_chr == "old","new_nms_chr","old_nms_chr"),
                                                                                          target_var_nm_1L_chr = ifelse(type_1L_chr == "old","old_nms_chr","new_nms_chr")),
                                                                 .x)))
  }
  if(flatten_1L_lgl){
    object_xx <- object_xx %>% purrr::flatten_chr()
  }
  return(object_xx)
}
## From ready4use
manufacture_Ready4useDyad <- function(x,
                                      correspondences_r3 = ready4show::ready4show_correspondences(),
                                      dummys_ls = list(),
                                      indcs_int = 1L,
                                      flatten_1L_lgl = F,
                                      match_1L_chr = character(0),
                                      restrict_to_chr = character(0),
                                      remove_ctg_1L_lgl = F,
                                      target_1L_chr = "var_ctg_chr",
                                      type_1L_chr = "all",
                                      what_1L_chr = "factors",
                                      ...){
  if(what_1L_chr %in% c("factors","factors-c","factors-d")){
    factor_vars_chr <- x@dictionary_r3 %>% 
      ready4use::renew.ready4use_dictionary(filter_cdn_1L_chr = "var_type_chr == 'factor'") %>% 
      dplyr::pull(var_nm_chr) %>% as.character()
    if(!identical(restrict_to_chr, character(0))){
      factor_vars_chr <- intersect(factor_vars_chr, restrict_to_chr) 
    }
    if(!identical(dummys_ls, list())){
      old_ls <- manufacture_Ready4useDyad(x, flatten_1L_lgl = F, type_1L_chr = type_1L_chr, what_1L_chr = "factors", restrict_to_chr = names(dummys_ls)) ## Replace
      dummys_ls <- dummys_ls[names(old_ls)]
      dummys_ls <- purrr::map2(dummys_ls, names(dummys_ls), ~ paste0(.y,.x))
      object_xx <- manufacture_Ready4useDyad(x, 
                                             flatten_1L_lgl = F, 
                                             type_1L_chr = type_1L_chr, 
                                             what_1L_chr = what_1L_chr, 
                                             restrict_to_chr = names(dummys_ls), 
                                             new_nms_chr = dummys_ls %>% purrr::flatten_chr()) 
    }else{
      if(type_1L_chr == "vars" && what_1L_chr == "factors"){
        object_xx <- factor_vars_chr 
      }else{
        object_xx <- purrr::map(factor_vars_chr,
                                ~ paste0(.x,eval(parse(text=paste0("levels(x@ds_tb$",.x,")"))))) %>%
          stats::setNames(factor_vars_chr)
        if(what_1L_chr == "factors-d"){
          lookup_ls <- object_xx 
        }
        if(type_1L_chr == "base"){
          object_xx <- object_xx %>% purrr::map2(indcs_int, ~ .x[.y])
        }
        if(type_1L_chr == "dummys"){
          object_xx <- object_xx %>% purrr::map2(indcs_int, ~ .x[-.y])
        }
        if(remove_ctg_1L_lgl){
          object_xx <- object_xx %>% purrr::map2(names(object_xx), ~ stringr::str_remove(.x,.y))
        }
        if(!identical(correspondences_r3, ready4show::ready4show_correspondences())){
          object_xx <- manufacture.ready4show_correspondences(correspondences_r3, # update to ready4show::manufacture.ready4show_correspondences
                                                              data_ls = object_xx,
                                                              flatten_1L_chr = F, # Correct
                                                              type_1L_chr = "new",
                                                              what_1L_chr = "names",
                                                              ...) 
        }
        if(flatten_1L_lgl | what_1L_chr %in% c("factors-c","factors-d")){
          object_xx <- object_xx %>% purrr::flatten_chr()
        }
        if(what_1L_chr %in% c("factors-c", "factors-d")){
          object_xx <- ready4show::ready4show_correspondences(ready4show::make_pt_ready4show_correspondences(old_nms_chr = object_xx, ...))
          if(what_1L_chr == "factors-d"){
            object_xx <- ready4use_dictionary(make_pt_ready4use_dictionary(var_nm_chr = object_xx %>% purrr::pmap_chr(~ifelse(is.na(as.character(..2)),
                                                                                                                              as.character(..1),as.character(..2))),
                                                                           var_ctg_chr = object_xx$old_nms_chr %>%  purrr::map_chr(~{
                                                                             level_1L_chr <- .x
                                                                             names(lookup_ls)[purrr::map_lgl(lookup_ls,~level_1L_chr %in% .x)]
                                                                           }),
                                                                           var_desc_chr = lookup_ls %>% 
                                                                             purrr::map2(names(lookup_ls), 
                                                                                         ~ stringr::str_remove(.x,.y)) %>% 
                                                                             purrr::map2(indcs_int, 
                                                                                         ~ {
                                                                                           if(!type_1L_chr %in% c("base","dummys")){
                                                                                             description_chr <- .x
                                                                                           }
                                                                                           if(type_1L_chr == "base"){
                                                                                             description_chr <- .x[.y]
                                                                                           }
                                                                                           if(type_1L_chr == "dummys"){
                                                                                             description_chr <- .x[-.y]
                                                                                           }
                                                                                           description_chr
                                                                                         }) %>% purrr::flatten_chr(),
                                                                           var_type_chr = object_xx$old_nms_chr %>%  
                                                                             purrr::map(~{
                                                                               level_1L_chr <- .x
                                                                               purrr::map2(lookup_ls, indcs_int,
                                                                                           ~if(level_1L_chr %in% .x[.y]){
                                                                                             "base"
                                                                                           }else{
                                                                                             if(level_1L_chr %in% .x[-.y]){
                                                                                               "dummy"#rep("dummy", length(.x[-.y]))
                                                                                             }else{
                                                                                               NULL
                                                                                             }
                                                                                           }) %>% purrr::discard(is.null) %>%
                                                                                 purrr::flatten_chr()
                                                                             }) %>% purrr::flatten_chr()))
          }
          if(!identical(match_1L_chr, character(0))){
            object_xx <- ready4::get_from_lup_obj(object_xx,
                                                  match_value_xx = match_1L_chr,
                                                  match_var_nm_1L_chr = ifelse(what_1L_chr == "factors-d","var_nm_chr","old_nms_chr"),
                                                  target_var_nm_1L_chr = ifelse(what_1L_chr == "factors-d",target_1L_chr,"new_nms_chr"))
          }
          
        }
      }
    }


  }
  return(object_xx)
}
renew_Ready4useDyad <- function(x,
                                new_val_xx = NULL,
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
  if(type_1L_chr %in% c("base","dummys","levels")){
    dummys_dict_r3 <- manufacture_Ready4useDyad(x, # replace
                                                dummys_ls = new_val_xx,
                                                flatten_1L_lgl = F, 
                                                type_1L_chr = ifelse(type_1L_chr=="levels","all",type_1L_chr), 
                                                what_1L_chr = "factors-d") 
    x@dictionary_r3 <- renew.ready4use_dictionary(x@dictionary_r3,
                                                  new_cases_r3 = dummys_dict_r3)
    x@ds_tb <- purrr::reduce(dummys_dict_r3$var_ctg_chr %>% unique(),
                             .init = x@ds_tb,
                             ~{
                               var_nm_1L_chr <- .y
                               val_1_1L_chr <- if("base" %in% ready4::get_from_lup_obj(dummys_dict_r3, match_value_xx = .y, match_var_nm_1L_chr = "var_ctg_chr", target_var_nm_1L_chr = "var_type_chr")){character(0)}else{levels(.x %>% dplyr::pull(!!rlang::sym(.y)))[1]}
                               .x %>% dplyr::mutate(!!rlang::sym(.y) := factor(!!rlang::sym(.y), labels = c(val_1_1L_chr,
                                                                                                            ready4::get_from_lup_obj(dummys_dict_r3, match_value_xx = .y, match_var_nm_1L_chr = "var_ctg_chr", target_var_nm_1L_chr = "var_nm_chr") %>% 
                                                                                                              purrr::map_chr(~stringi::stri_replace_first_fixed(.x,var_nm_1L_chr,""))
                               )))
                             } )
     
  }
  return(x)
}
## From youthvars
make_csnl_example_dict <- function(ds_tb){ # New to youthvars
  dictionary_r3 <- Ready4useRepos(dv_nm_1L_chr = "TTU",
                                  dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                  dv_server_1L_chr = "dataverse.harvard.edu") %>% 
    ingest(fls_to_ingest_chr = c("dictionary_r3"), 
           metadata_1L_lgl = F) 
  dictionary_r3 <- dictionary_r3 %>% # Modify CALD to cald
    dplyr::filter(var_nm_chr %in% names(ds_tb)) %>%
    dplyr::filter(!var_nm_chr %in% c("CALD", "c_p_diag_s", "d_ATSI", "d_gender", "d_studying_working")) %>%
    renew.ready4use_dictionary(var_nm_chr = c(setdiff(names(ds_tb),dictionary_r3$var_nm_chr),"c_p_diag_s", "d_ATSI", "d_gender","d_studying_working") %>% sort(),
                               var_ctg_chr = c(rep("clinical symptom",5),
                                               rep("multi-attribute utility instrument question",9),
                                               "health_utility",
                                               rep("demographic",7),
                                               rep("difference",2),
                                               "psychological distress",
                                               "quality of life",
                                               rep("spatial",2),
                                               rep("validation",2)),
                               var_desc_chr = c("days unable to perform usual activities",
                                                "days out of role",
                                                "days cut back on usual activities",
                                                "primary diagnosis group",
                                                "primary diagnosis",
                                                paste0("Child Health Utility (9 Dimension) question ",1:9),
                                                "Child Health Utility (9 Dimension) total score",
                                                "Aboriginal and Torres Strait Islander",
                                                "culturally and linguistically diverse",
                                                "in employment",
                                                "employment type",
                                                "gender",
                                                "in education",
                                                "education and employment",
                                                "Difference between AQoL-6D and CHU-9D total scores",
                                                "Difference between AQoL-6D total scores and validation values",
                                                "Kessler Psychological Distress Scale (10 Item)",
                                                "My Life Tracker", 
                                                "area index of relative social disadvantage",
                                                "area remoteness",
                                                "validation unweighted aqol total",
                                                "validation weighted aqol total"
                               ),
                               var_type_chr = c(setdiff(names(ds_tb),dictionary_r3$var_nm_chr),"c_p_diag_s", "d_ATSI", "d_gender","d_studying_working") %>% 
                                 sort() %>% purrr::map_chr(~{
                                   classes_chr <- ds_tb[,.x][[1]] %>% class()
                                   ifelse("numeric" %in% classes_chr,
                                          ifelse(is.integer(ds_tb[,.x][[1]]),
                                                 "integer",
                                                 "double"),
                                          classes_chr[1])
                                 })) %>% 
    dplyr::arrange(var_ctg_chr, var_nm_chr)
  dictionary_r3
}
make_item_plt <- function(tfd_data_tb,
                          var_nm_1L_chr,
                          x_label_1L_chr,
                          fill_label_1L_chr = "Data collection",
                          legend_position_1L_chr = "none",
                          round_var_nm_1L_chr = "round",
                          sngl_round_lbl_1L_chr = "n",
                          sngl_round_var_nm_1L_chr = "n",
                          use_bw_theme_1L_lgl = F,
                          y_label_1L_chr = "Percentage",
                          y_scale_scl_fn = NULL
                          ){
  item_plt <- ggplot2::ggplot(tfd_data_tb %>%
                                dplyr::with_groups(NULL,ready4use::remove_labels_from_ds),#ready4use::remove_labels_from_ds(),
                              ggplot2::aes_string(var_nm_1L_chr)) +
    ggplot2::geom_bar(ggplot2::aes(y = y, fill = !!rlang::sym(ifelse(identical(round_var_nm_1L_chr,character(0)),
                                                                     sngl_round_var_nm_1L_chr,
                                                                     round_var_nm_1L_chr))
                                   ),
                      stat = "identity",
                      na.rm = TRUE,
                      position = "dodge",
                      colour = "white",
                      alpha= 0.7)
  if(!is.null(y_scale_scl_fn)){
    item_plt <- item_plt +
      ggplot2::scale_y_continuous(labels = y_scale_scl_fn)
  }
  item_plt <- item_plt +
    ggplot2::labs(x = x_label_1L_chr,
                  y = y_label_1L_chr,
                  fill = ifelse(identical(round_var_nm_1L_chr,character(0)),
                                sngl_round_lbl_1L_chr,
                                fill_label_1L_chr)
                  )
  if(use_bw_theme_1L_lgl){
    item_plt <- item_plt +
      ggplot2::theme_bw()
  }
  item_plt <- item_plt +
    ggplot2::theme(legend.position = legend_position_1L_chr)
  if(!identical(round_var_nm_1L_chr,character(0))){
    item_plt <- item_plt +
      ggplot2::scale_fill_manual(values = c("#de2d26",
                                            "#fc9272"
      ))
  }
  return(item_plt)
}
make_itm_resp_plts <- function(data_tb,
                               col_nms_chr,
                               lbl_nms_chr,
                               plot_rows_cols_pair_int,
                               heights_int,
                               round_var_nm_1L_chr = "round",
                               y_label_1L_chr = "Percentage"){
  plots_ls <- list()
  j=1
  for(i in col_nms_chr){
    tfd_data_tb <- data_tb %>%
      transform_ds_for_item_plt(var_nm_1L_chr = i,
                                round_var_nm_1L_chr = round_var_nm_1L_chr)
    labelx <- lbl_nms_chr[j]
    j = j+1
    plots_ls[[i]]<- make_item_plt(tfd_data_tb,
                                  var_nm_1L_chr = i,
                                  round_var_nm_1L_chr = round_var_nm_1L_chr,
                                  x_label_1L_chr = labelx,
                                  y_label_1L_chr = y_label_1L_chr,
                                  y_scale_scl_fn = scales::percent_format(),
                                  use_bw_theme_1L_lgl = T,
                                  legend_position_1L_chr = "none")
  }
  if(!identical(round_var_nm_1L_chr, character(0))){
    plot_plt <- make_item_plt(tfd_data_tb,
                              var_nm_1L_chr = i,
                              round_var_nm_1L_chr = round_var_nm_1L_chr,
                              x_label_1L_chr = labelx,
                              y_label_1L_chr = y_label_1L_chr,
                              y_scale_scl_fn = NULL,
                              use_bw_theme_1L_lgl = F,
                              legend_position_1L_chr = "bottom")
    legend_ls <- get_guide_box_lgd(plot_plt)
  }else{
    legend_ls <- NULL
  }
  composite_plt <- gridExtra::grid.arrange(ggpubr::ggarrange(plotlist=plots_ls,
                                                             nrow = plot_rows_cols_pair_int[1],
                                                             ncol = plot_rows_cols_pair_int[2]),
                                           legend_ls,
                                           nrow = length(heights_int),
                                           heights = heights_int)
  composite_plt$grobs <- composite_plt$grobs %>% purrr::discard(is.null)
  return(composite_plt)
}
make_predr_pars_and_cors_tbl <- function(data_tb,
                                         ds_descvs_ls,
                                         descv_tbl_ls,
                                         dictionary_tb,
                                         nbr_of_digits_1L_int = 2L,
                                         predictors_lup = NULL){
  predr_pars_and_cors_tb <- make_cors_with_utl_tbl(data_tb,
                                                   ds_descvs_ls = ds_descvs_ls,
                                                   dictionary_tb = dictionary_tb) %>%
    dplyr::mutate(label = paste0("Correlation with ",
                                 ready4::get_from_lup_obj(dictionary_tb,
                                                          match_var_nm_1L_chr = "var_nm_chr",
                                                          match_value_xx = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                          target_var_nm_1L_chr = "var_desc_chr",
                                                          evaluate_1L_lgl = F)))
  
  predr_pars_and_cors_tb <- purrr::map_dfr(1:nrow(predr_pars_and_cors_tb),
                                           ~ predr_pars_and_cors_tb %>%
                                             dplyr::slice(.x) %>%
                                             dplyr::mutate(dplyr::across(paste0(ds_descvs_ls$round_vals_chr,"_sig_dbl"),
                                                                         ~ format(round(.x, nbr_of_digits_1L_int),
                                                                                  nsmall = nbr_of_digits_1L_int))) %>%
                                             
                                             dplyr::mutate(p.value = paste0(c( !!!rlang::syms(paste0(ds_descvs_ls$round_vals_chr,"_sig_dbl"))),
                                                                            collapse = ", "))) %>%
    dplyr::mutate(dplyr::across(paste0(ds_descvs_ls$round_vals_chr,"_sig_dbl"), ~ "")) %>%
    dplyr::mutate(dplyr::across(paste0(ds_descvs_ls$round_vals_chr,"_cor_dbl"),
                                ~ format(round(.x, nbr_of_digits_1L_int),
                                         nsmall = nbr_of_digits_1L_int))) %>%
    dplyr::rename_with(~stringr::str_replace(.x,"_cor_dbl","_val_1_chr") %>%
                         stringr::str_replace("_sig_dbl","_val_2_chr")) %>%
    dplyr::select(variable_chr,
                  label,
                  dplyr::everything())
  
  
  main_outc_tbl_tb <- descv_tbl_ls$main_outc_tbl_tb %>%
    dplyr::filter(label %in% c("Mean (SD)", "Missing")) %>%
    dplyr::rename_with(~stringr::str_replace(.x,"_val_1_dbl","_val_1_chr") %>%
                         stringr::str_replace("_val_2_ls","_val_2_chr") %>%
                         stringr::str_replace("variable", "variable_chr")) %>%
    dplyr::filter(variable_chr %in% purrr::map_chr(ds_descvs_ls$candidate_predrs_chr,
                                                   ~ready4::get_from_lup_obj(dictionary_tb,
                                                                             match_var_nm_1L_chr = "var_nm_chr",
                                                                             match_value_xx = .x,
                                                                             target_var_nm_1L_chr = "var_desc_chr",
                                                                             evaluate_1L_lgl = F)))
  if("p.value" %in% names(predr_pars_and_cors_tb) & !"p.value" %in% names(main_outc_tbl_tb))
    main_outc_tbl_tb <- main_outc_tbl_tb %>%
    dplyr::mutate(p.value = "")
  predr_pars_and_cors_tb <- main_outc_tbl_tb$variable_chr %>% unique() %>%
    purrr::map_dfr(~tibble::add_case(main_outc_tbl_tb %>%
                                       dplyr::filter(variable_chr == .x),
                                     predr_pars_and_cors_tb %>%
                                       dplyr::filter(variable_chr == .x)))
  if(!is.null(predictors_lup)){
    predr_pars_and_cors_tb <- predr_pars_and_cors_tb %>%
      dplyr::mutate(variable_chr = purrr::map_chr(variable_chr,
                                                  ~ {
                                                    var_nm_1L_chr <- ready4::get_from_lup_obj(dictionary_tb,
                                                                                              match_var_nm_1L_chr = "var_desc_chr",
                                                                                              match_value_xx = .x,
                                                                                              target_var_nm_1L_chr = "var_nm_chr",
                                                                                              evaluate_1L_lgl = F)
                                                    paste0 (.x,
                                                            " (",
                                                            ready4::get_from_lup_obj(predictors_lup,
                                                                                     match_var_nm_1L_chr = "short_name_chr",
                                                                                     match_value_xx = var_nm_1L_chr,
                                                                                     target_var_nm_1L_chr = "min_val_dbl",
                                                                                     evaluate_1L_lgl = F),
                                                            "-",
                                                            ready4::get_from_lup_obj(predictors_lup,
                                                                                     match_var_nm_1L_chr = "short_name_chr",
                                                                                     match_value_xx = var_nm_1L_chr,
                                                                                     target_var_nm_1L_chr = "max_val_dbl",
                                                                                     evaluate_1L_lgl = F),
                                                            ")"
                                                    )
                                                  }
      ))
  }
  return(predr_pars_and_cors_tb)
}
make_subtotal_plt <- function(data_tb,
                              var_nm_1L_chr,
                              x_label_1L_chr,
                              legend_position_1L_chr = "none",
                              legend_sclg_1L_dbl = 1,
                              label_fill_1L_chr = NULL,
                              round_var_nm_1L_chr = "round",
                              axis_text_sclg_1L_dbl = 1,
                              axis_title_sclg_1L_dbl = 1,
                              use_bw_theme_1L_lgl = T,
                              y_label_1L_chr = "Percentage",
                              y_scale_scl_fn = scales::percent
){
  subtotal_plt <- ggplot2::ggplot(data_tb %>% ready4use::remove_labels_from_ds(),
                                  ggplot2::aes_string(var_nm_1L_chr)) +
    ggplot2::geom_histogram(bins=8,
                            color = "white",
                            if(identical(round_var_nm_1L_chr,character(0))){
                              ggplot2::aes(y = 2*(..density..)/sum(..density..))
                              }else{
                            ggplot2::aes(fill = !!rlang::sym(round_var_nm_1L_chr),
                                         y = 2*(..density..)/sum(..density..))},
                            position = 'dodge',
                            alpha=0.7)
  subtotal_plt <- subtotal_plt +
    ggplot2::labs(x = x_label_1L_chr,
                  y = y_label_1L_chr,                  
                  fill = if(identical(round_var_nm_1L_chr,character(0))){NULL}else{label_fill_1L_chr}
                  )
  if(use_bw_theme_1L_lgl){
    subtotal_plt <- subtotal_plt +
      ggplot2::theme_bw()
  }
  if(!is.null(y_scale_scl_fn)){
    subtotal_plt <- subtotal_plt +
      ggplot2::scale_y_continuous(labels = y_scale_scl_fn)
  }
  subtotal_plt <- subtotal_plt +
    ggplot2::theme(legend.position = legend_position_1L_chr,
                   legend.text = ggplot2::element_text(size=ggplot2::rel(legend_sclg_1L_dbl)),
                   legend.title = ggplot2::element_text(size=ggplot2::rel(legend_sclg_1L_dbl)),
                   axis.text = ggplot2::element_text(size=ggplot2::rel(axis_text_sclg_1L_dbl)),
                   axis.title = ggplot2::element_text(size=ggplot2::rel(axis_title_sclg_1L_dbl))) 
  if(!identical(character(0), round_var_nm_1L_chr)){
    subtotal_plt <- subtotal_plt +
      ggplot2::scale_fill_manual(values = c("#de2d26","#fc9272"))
  }
  return(subtotal_plt)
}
make_sub_tot_plts <- function(data_tb,
                              col_nms_chr,
                              heights_int,
                              plot_rows_cols_pair_int,
                              add_legend_1L_lgl = T,
                              axis_text_sclg_1L_dbl = 1,
                              axis_title_sclg_1L_dbl = 1,
                              legend_sclg_1L_dbl = 1,
                              make_log_log_tfmn_1L_lgl = F,
                              round_var_nm_1L_chr = "round",
                              x_labels_chr = character(0),
                              y_label_1L_chr = "Percentage"){
  if(!is.null(col_nms_chr)){
    plots_ls<-list()
    for(i in col_nms_chr){
      if(make_log_log_tfmn_1L_lgl){
        targetvar = paste0("tran_",i)
        data_tb <- dplyr::mutate(data_tb, !!targetvar := log(-log(1-!!as.name(i))))  %>%
          dplyr::mutate(!!targetvar :=ifelse(!!as.name(i)==1,log(-log(1-0.999)),!!as.name(targetvar)))
      }
      if(identical(x_labels_chr,character(0))){
        labelx <- eval(parse(text=paste0("attributes(data_tb$",i,")$label")))
        labelx <- stringr::str_sub(labelx,
                                   start = stringi::stri_locate_last_fixed(labelx," - ")[1,1] %>%
                                     unname() + 2)
      }else{
        labelx <- x_labels_chr[which(col_nms_chr==i)]
      }

      if(make_log_log_tfmn_1L_lgl){
        labelx<- paste0("log-log transformed ", labelx)
      }
      plots_ls[[i]]<- make_subtotal_plt(data_tb,
                                        legend_sclg_1L_dbl = legend_sclg_1L_dbl,
                                        round_var_nm_1L_chr = round_var_nm_1L_chr,
                                        axis_text_sclg_1L_dbl = axis_text_sclg_1L_dbl,
                                        axis_title_sclg_1L_dbl = axis_title_sclg_1L_dbl,
                                        var_nm_1L_chr = i,
                                        x_label_1L_chr = labelx,
                                        y_label_1L_chr = y_label_1L_chr)
    }
    if(add_legend_1L_lgl & !identical(character(0), round_var_nm_1L_chr)){
      plot_for_lgd_plt <- make_subtotal_plt(data_tb,
                                            legend_sclg_1L_dbl = legend_sclg_1L_dbl,
                                            round_var_nm_1L_chr = round_var_nm_1L_chr,
                                            var_nm_1L_chr = i,
                                            x_label_1L_chr = labelx,
                                            legend_position_1L_chr = "bottom",
                                            label_fill_1L_chr = "Data collection",
                                            axis_text_sclg_1L_dbl = axis_text_sclg_1L_dbl,
                                            axis_title_sclg_1L_dbl = axis_title_sclg_1L_dbl,
                                            y_label_1L_chr = y_label_1L_chr)
      legend_ls <- get_guide_box_lgd(plot_for_lgd_plt)
      composite_plt <- gridExtra::grid.arrange(ggpubr::ggarrange(plotlist = plots_ls,
                                                                 nrow = plot_rows_cols_pair_int[1],
                                                                 ncol = plot_rows_cols_pair_int[2]),
                                               legend_ls,
                                               nrow = length(heights_int),
                                               heights = heights_int)
    }else{
      legend_ls <- NULL
      heights_int <- min(heights_int[-length(heights_int)],length(plots_ls))
      composite_plt <- gridExtra::grid.arrange(ggpubr::ggarrange(plotlist=plots_ls,
                                                                 nrow = max(plot_rows_cols_pair_int[1], ceiling(length(plots_ls)/plot_rows_cols_pair_int[2])),
                                                                 ncol = plot_rows_cols_pair_int[2]),
                                               nrow = length(heights_int),
                                               heights = heights_int)
    }
  }else{
    composite_plt <- NULL
  }
  return(composite_plt)
}
make_var_by_round_plt <- function(data_tb,
                                  var_nm_1L_chr,
                                  x_label_1L_chr,
                                  label_fill_1L_chr = "Data collection",
                                  legend_sclg_1L_dbl = 1,
                                  axis_text_sclg_1L_dbl = 1,
                                  axis_title_sclg_1L_dbl = 1,
                                  round_var_nm_1L_chr = "round",
                                  y_label_1L_chr = "Percentage",
                                  y_scale_scl_fn = scales::percent){
  var_by_round_plt <- ggplot2::ggplot(data_tb %>% ready4use::remove_labels_from_ds(),
                                      if(identical(round_var_nm_1L_chr, character(0))){
                                        ggplot2::aes(x = !!rlang::sym(var_nm_1L_chr))
                                      }else{
                                        ggplot2::aes(x = !!rlang::sym(var_nm_1L_chr),
                                                     fill = !!rlang::sym(round_var_nm_1L_chr))
                                      }
                                      
                                      ) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(width*density)),
                            bins = 10,
                            position="dodge",
                            colour="white",
                            alpha=0.7)
  
  if(!is.null(y_scale_scl_fn)){
    var_by_round_plt <- var_by_round_plt +
      ggplot2::scale_y_continuous(labels = y_scale_scl_fn)
  }
  var_by_round_plt <- var_by_round_plt +
    #ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(y = y_label_1L_chr, x= x_label_1L_chr, 
                  fill = if(identical(round_var_nm_1L_chr,character(0))){NULL}else{label_fill_1L_chr}
                  )  +
    ggplot2::scale_fill_manual(values=c("#de2d26","#fc9272"))  +
    ggplot2::theme(legend.position="bottom",
                   legend.text = ggplot2::element_text(size=ggplot2::rel(legend_sclg_1L_dbl)),
                   legend.title = ggplot2::element_text(size=ggplot2::rel(legend_sclg_1L_dbl)),
                   axis.text = ggplot2::element_text(size=ggplot2::rel(axis_text_sclg_1L_dbl)),
                   axis.title = ggplot2::element_text(size=ggplot2::rel(axis_title_sclg_1L_dbl)))
  return(var_by_round_plt)
}
transform_csnl_example_ds <- function(ds_df){ # New to youthvars
  ds_tb <- ds_df %>% 
    tibble::as_tibble() %>%
    dplyr::mutate(c_p_diag_grouped = dplyr::case_when(as.character(DiagnosisPrimary) %in% c("Acute stress disorder",
                                                                                            "Adjustment disorder",
                                                                                            "Agoraphobia",
                                                                                            "Anxiety symptoms",
                                                                                            "Bipolar disorder",
                                                                                            "Cyclothymic disorder",
                                                                                            "Depressive disorder NOS",
                                                                                            "Depressive symptoms",
                                                                                            "Dysthymia",
                                                                                            "Generalised anxiety disorder",
                                                                                            "Major depressive disorder",
                                                                                            "Mixed anxiety and depressive symptoms",
                                                                                            "Obsessive-compulsive disorder",
                                                                                            "Other affective disorder",
                                                                                            "Other anxiety disorder",
                                                                                            "Panic disorder",
                                                                                            "Post-traumatic stress disorder",
                                                                                            "Separation anxiety disorder",
                                                                                            "Social phobia",
                                                                                            "Stress related") ~ "Anxiety and Depression",
                                                      as.character(DiagnosisPrimary) %in% c("Attention deficit hyperactivity disorder (ADHD)",
                                                                                            "Conduct disorder",
                                                                                            "Feeding and Eating Disorders",
                                                                                            "Gender Dysphoria",
                                                                                            "Neurocognitive Disorders",
                                                                                            "Neurodevelopmental Disorders",
                                                                                            "Oppositional defiant disorder",
                                                                                            "Other",
                                                                                            "Personality Disorders",
                                                                                            "Pervasive developmental disorder",
                                                                                            "Schizoaffective disorder",
                                                                                            "Schizophrenia",
                                                                                            "Sleep-Wake Disorders",
                                                                                            "Somatic Symptom and Related Disorders" ) ~ "Other Mental Disorder",
                                                      as.character(DiagnosisPrimary) %in% c("Alcohol dependence",
                                                                                            "Other drug dependence") ~ "Substance Use",
                                                      as.character(DiagnosisPrimary) %in% c("Not applicable (e.g. for non-Mental Health related services, or service provider not qualified to give diagnosis)",
                                                                                            "Diagnosis not yet assessed or requires further assessment",
                                                                                            "No diagnosis (and no sub-syndromal mental health problems)") ~ "Not applicable",
                                                      is.na(DiagnosisPrimary) ~ NA_character_,
                                                      T ~ "Uncategorised") %>%
                    as.factor()) %>%
    dplyr::rename(c_days_cut_back = K12_DaysCutDown,
                  c_days_unable = K11_DaysTotallyUnable,
                  CHU9D = chu9_total_w,
                  c_p_diag_s = DiagnosisPrimary,
                  d_age = Age,
                  d_ATSI = ATSI,
                  d_CALD = CALD, # Uncomment when dictionary is updated
                  d_employed = Working,
                  d_employment_type = EmploymentType,
                  d_gender = Gender,
                  d_studying = Studying,
                  K10 = K10_total,
                  MLT = MLT_mean, 
                  s_IRSD = IRSD,
                  s_remoteness = Remoteness,
                  validation_aqol_c = aqol6d_total_c,
                  validation_aqol_w = aqol6d_total_w) 
  ds_tb <- ds_tb %>%
    dplyr::mutate(dplyr::across(c(dplyr::starts_with("aqol6d_q"),
                                  dplyr::starts_with("chu9_q"),
                                  dplyr::starts_with("c_days_"),
                                  d_age,
                                  K10,
                                  s_IRSD,
                                  SOFAS), ~as.integer(.x))) %>%
  dplyr::mutate(dplyr::across(c(c_p_diag_s,
                                d_ATSI,
                                d_gender), ~as.factor(.x)))
  ds_tb <- ds_tb %>%
    dplyr::mutate(c_days_oor = c_days_cut_back + c_days_unable) %>%
    dplyr::mutate(d_studying_working = dplyr::case_when(purrr::map2_lgl(as.character(d_employed), as.character(d_studying), ~ is.na(.x) | is.na(.y)) ~ NA_character_,
                                                        purrr::map2_lgl(as.character(d_employed), as.character(d_studying), ~ .x == "No" && .y == "No") ~ "Not studying or working",
                                                        purrr::map2_lgl(as.character(d_employed), as.character(d_studying), ~ .x == "No" && .y == "Yes") ~ "Studying only",
                                                        purrr::map2_lgl(as.character(d_employed), as.character(d_studying), ~ .x == "Yes" && .y == "No") ~ "Working only",
                                                        purrr::map2_lgl(as.character(d_employed), as.character(d_studying), ~ .x == "Yes" && .y == "Yes") ~ "Studying and working",
                                                        T ~ "Uncategorised"
                                                        ) %>% as.factor()) %>%
    dplyr::mutate(difference_mauis = validation_aqol_w - CHU9D) %>%
    dplyr::mutate(difference_aqol_calcs = NA_real_)
  ds_tb <- youthvars::add_uids_to_tbs_ls(list(ds_tb),"Participant_") %>% purrr::pluck(1)  %>% ready4::remove_lbls_from_df()
  return(ds_tb)
}
transform_ds_for_item_plt <- function(data_tb,
                                      var_nm_1L_chr,
                                      round_var_nm_1L_chr = "round"){
  tfd_data_tb <- data_tb %>%
    dplyr::filter(!is.na(!!as.name(var_nm_1L_chr) ))
  if(!identical(round_var_nm_1L_chr, character(0))){
    tfd_data_tb <- tfd_data_tb %>%
      dplyr::group_by(!!rlang::sym(round_var_nm_1L_chr), !!as.name(var_nm_1L_chr) )
  }else{
    tfd_data_tb <- tfd_data_tb %>%
      dplyr::group_by(!!as.name(var_nm_1L_chr) )
  }
  tfd_data_tb <- tfd_data_tb %>%
    dplyr::summarise(n = dplyr::n()) 
  if(!identical(round_var_nm_1L_chr, character(0)))
    tfd_data_tb <- tfd_data_tb %>%
    dplyr::group_by(!!rlang::sym(round_var_nm_1L_chr)) 
  tfd_data_tb <- tfd_data_tb %>%
    dplyr::mutate(y = n/sum(n))
  return(tfd_data_tb)
}
transform_ds_for_tstng <- function (data_tb, depnt_var_nm_1L_chr = "aqol6d_total_w", depnt_var_max_val_1L_dbl = 0.999,
                                    candidate_predrs_chr = NA_character_, covar_var_nms_chr = NA_character_,
                                    round_var_nm_1L_chr = "round", round_val_1L_chr = "Baseline",
                                    remove_all_msng_1L_lgl = F)
{
  vars_to_keep_chr <- c(depnt_var_nm_1L_chr, candidate_predrs_chr,
                        covar_var_nms_chr) %>% purrr::discard(is.na)
  tfd_data_tb <- data_tb
  if(!identical(round_var_nm_1L_chr, character(0)) && (!is.na(round_var_nm_1L_chr) & !is.na(round_val_1L_chr)))
    tfd_data_tb <- tfd_data_tb %>%
    dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == round_val_1L_chr)
  tfd_data_tb <- tfd_data_tb %>%
    dplyr::select(!!!rlang::syms(vars_to_keep_chr)) %>%
    dplyr::mutate(`:=`(!!rlang::sym(depnt_var_nm_1L_chr), ifelse(!!rlang::sym(depnt_var_nm_1L_chr) >
                                                                   depnt_var_max_val_1L_dbl, depnt_var_max_val_1L_dbl, !!rlang::sym(depnt_var_nm_1L_chr))))
  if (remove_all_msng_1L_lgl)
    tfd_data_tb <- tfd_data_tb %>% stats::na.omit()
  return(tfd_data_tb)
}
write_descv_plots <- function(data_tb,
                              ds_descvs_ls,
                              descv_outp_dir_1L_chr,
                              lbl_nms_chr = c("Household tasks", "Getting around",
                                              "Morbility","Self care","Enjoy close rels",
                                              "Family rels", "Community involvement",
                                              "Despair","Worry", "Sad", "Agitated",
                                              "Energy level", "Control", "Coping",
                                              "Frequency of pain", "Degree of pain",
                                              "Pain interference","Vision", "Hearing",
                                              "Communication"),
                              maui_domains_pfxs_1L_chr = "vD",
                              item_plots_params_ls = list(plot_rows_cols_pair_int = c(5L,4L),
                                                          heights_int = c(10L, 1L),
                                                          width_1L_dbl = 9),
                              dim_plots_params_ls = list(plot_rows_cols_pair_int = c(3L,2L),
                                                         heights_int = c(10L, 1L),
                                                         width_1L_dbl = 8),
                              utl_by_rnd_plots_params_ls = list(width_1L_dbl = 6,
                                                                height_1L_dbl = 4),
                              combined_plot_params_ls = list(nrow_1L_int = 2L,
                                                             rel_heights_dbl = c(4,10),
                                                             scale_dbl = c(0.9,0.9),
                                                             base_height_dbl = 10),
                              x_labels_chr = character(0)
){
  if(is.null(maui_domains_pfxs_1L_chr)){
    maui_domains_col_nms_chr <- NULL
  }else{
    maui_domains_col_nms_chr <- names(dplyr::select(data_tb, dplyr::starts_with(maui_domains_pfxs_1L_chr)))
  }
  plots_params_ls <- list(qstn_rspns = list(plt_fn = make_itm_resp_plts,
                                            fn_args_ls = list(data_tb,
                                                              col_nms_chr = names(dplyr::select(data_tb,
                                                                                                starts_with(ds_descvs_ls$maui_item_pfx_1L_chr))),
                                                              lbl_nms_chr = lbl_nms_chr,
                                                              plot_rows_cols_pair_int = item_plots_params_ls$plot_rows_cols_pair_int,
                                                              heights_int = item_plots_params_ls$heights_int,
                                                              round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr),
                                            width_1L_dbl = item_plots_params_ls$width_1L_dbl,
                                            height_1L_dbl = sum(item_plots_params_ls$heights_int),
                                            path_to_write_to_1L_chr = descv_outp_dir_1L_chr,
                                            plt_nm_1L_chr = "qstn_rspns"),
                          wtd_sub_tots = list(plt_fn = make_sub_tot_plts,
                                              fn_args_ls = list(data_tb,
                                                                col_nms_chr = maui_domains_col_nms_chr,
                                                                plot_rows_cols_pair_int = dim_plots_params_ls$plot_rows_cols_pair_int,
                                                                round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                                heights_int = dim_plots_params_ls$heights_int,
                                                                x_labels_chr = x_labels_chr),
                                              width_1L_dbl = dim_plots_params_ls$width_1L_dbl,
                                              height_1L_dbl = sum(dim_plots_params_ls$heights_int),
                                              path_to_write_to_1L_chr = descv_outp_dir_1L_chr,
                                              plt_nm_1L_chr = "wtd_sub_tots"),
                          ll_sub_tot = list(plt_fn = make_sub_tot_plts,
                                            fn_args_ls = list(data_tb,
                                                              col_nms_chr = maui_domains_col_nms_chr,
                                                              plot_rows_cols_pair_int = dim_plots_params_ls$plot_rows_cols_pair_int,
                                                              round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                              heights_int = dim_plots_params_ls$heights_int,
                                                              make_log_log_tfmn_1L_lgl = T,
                                                              x_labels_chr = x_labels_chr),
                                            width_1L_dbl = dim_plots_params_ls$width_1L_dbl,
                                            height_1L_dbl = sum(dim_plots_params_ls$heights_int),
                                            path_to_write_to_1L_chr = descv_outp_dir_1L_chr,
                                            plt_nm_1L_chr = "ll_sub_tot"),
                          utl_by_rnd = list(plt_fn = make_var_by_round_plt,
                                            fn_args_ls = list(data_tb,
                                                              var_nm_1L_chr = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                              round_var_nm_1L_chr = ds_descvs_ls$round_var_nm_1L_chr,
                                                              x_label_1L_chr = ds_descvs_ls$dictionary_tb %>%
                                                                ready4::get_from_lup_obj(match_value_xx = ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                                                                                         match_var_nm_1L_chr = "var_nm_chr",
                                                                                         target_var_nm_1L_chr = "var_desc_chr",
                                                                                         evaluate_1L_lgl = F) %>% as.vector()),
                                            width_1L_dbl = utl_by_rnd_plots_params_ls$width_1L_dbl,
                                            height_1L_dbl = utl_by_rnd_plots_params_ls$height_1L_dbl,
                                            path_to_write_to_1L_chr = descv_outp_dir_1L_chr,
                                            plt_nm_1L_chr = "utl_by_rnd")
  )
  descv_plts_paths_ls <- purrr::map(plots_params_ls,
                                    ~ rlang::exec(ready4show::write_mdl_plt_fl,!!!.x)) %>%
    stats::setNames(names(plots_params_ls))
  combined_plt <- cowplot::plot_grid(rlang::exec(plots_params_ls$utl_by_rnd$plt_fn,!!!plots_params_ls$utl_by_rnd$fn_args_ls) + ggplot2::theme(legend.position = 'none'),
                                     rlang::exec(plots_params_ls$wtd_sub_tots$plt_fn,!!!plots_params_ls$wtd_sub_tots$fn_args_ls),
                                     nrow = combined_plot_params_ls$nrow_1L_int,
                                     rel_heights = combined_plot_params_ls$rel_heights_dbl,
                                     scale = combined_plot_params_ls$scale_dbl
  )
  descv_plts_paths_ls$combined_utl <- paste0(descv_outp_dir_1L_chr,
                                             "/combined_utl.png")
  cowplot::save_plot(descv_plts_paths_ls$combined_utl,
                     combined_plt,
                     base_height = combined_plot_params_ls$base_height_dbl)
  return(descv_plts_paths_ls)
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
  descv_tbl_ls$predr_pars_and_cors_tb <- make_predr_pars_and_cors_tbl(data_tb %>%
                                                                        dplyr::mutate(catch_all_round_chr = ds_descvs_ls$round_vals_chr[1]), ## ERROR HERE
                                                                      ds_descvs_ls = if(ds_descvs_ls$round_vals_chr == "Overall" & identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){ds_descvs_ls %>% purrr::list_modify(round_var_nm_1L_chr = "catch_all_round_chr")}else{ds_descvs_ls},
                                                                      descv_tbl_ls = descv_tbl_ls,
                                                                      dictionary_tb = ds_descvs_ls$dictionary_tb,
                                                                      nbr_of_digits_1L_int = nbr_of_digits_1L_int,
                                                                      predictors_lup = predictors_lup)
  saveRDS(descv_tbl_ls,paste0(descv_outp_dir_1L_chr,"/descv_tbls_ls.RDS"))
  return(descv_tbl_ls)
}
# From scorz
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
# From specific
add_utility_predn_to_ds <- function (data_tb,
                                     model_mdl,
                                     tfmn_1L_chr,
                                     depnt_var_nm_1L_chr,
                                     force_min_max_1L_lgl = T,
                                     force_new_data_1L_lgl = F,
                                     impute_1L_lgl = T,
                                     is_brms_mdl_1L_lgl = T,
                                     new_data_is_1L_chr = "Predicted",
                                     predn_type_1L_chr = NULL,
                                     predr_vars_nms_chr = NULL,
                                     rmv_tfd_depnt_var_1L_lgl = F,
                                     sd_dbl = NA_real_,
                                     utl_cls_fn = NULL,
                                     utl_min_val_1L_dbl = -1)
{
  depnt_vars_chr <- c(depnt_var_nm_1L_chr, transform_depnt_var_nm(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                  tfmn_1L_chr = tfmn_1L_chr)) %>% unique()
  data_tb <- purrr::reduce(depnt_vars_chr, .init = data_tb, ~dplyr::mutate(.x,
                                                                           !!rlang::sym(.y):= NA_real_))
  if(force_min_max_1L_lgl){
    min_max_vals_dbl <- c(utl_min_val_1L_dbl,1)
  }else{
    min_max_vals_dbl <- numeric(0)
  }
  predictions_dbl <- predict_vals(data_tb = data_tb, tfmn_1L_chr = tfmn_1L_chr,
                                  model_mdl = model_mdl,
                                  #force_min_max_1L_lgl = force_min_max_1L_lgl,
                                  #utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                  min_max_vals_dbl = min_max_vals_dbl,
                                  impute_1L_lgl = impute_1L_lgl,
                                  new_data_is_1L_chr = new_data_is_1L_chr,
                                  var_cls_fn = utl_cls_fn,
                                  is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl,
                                  force_new_data_1L_lgl = force_new_data_1L_lgl,
                                  predn_type_1L_chr = predn_type_1L_chr,
                                  sd_dbl = sd_dbl)
  data_tb <- data_tb %>% dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr):=predictions_dbl)
  if(!is.null(predr_vars_nms_chr)){
    derived_predrs_chr <- purrr::map(predr_vars_nms_chr, ~ paste0(.x,c("_baseline","_change","_scaled","_unscaled"))) %>%
      purrr::flatten_chr() %>% intersect(names(data_tb))
    data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(derived_predrs_chr))
  }
  if(rmv_tfd_depnt_var_1L_lgl){
    data_tb <- data_tb %>% dplyr::select(-tidyselect::all_of(depnt_vars_chr[depnt_vars_chr!=depnt_var_nm_1L_chr]))
  }
  return(data_tb)
}
add_utl_predn_to_new_ds <- function(data_tb,
                                    ingredients_ls,
                                    mdl_nm_1L_chr,
                                    analysis_1L_chr = NULL,
                                    deterministic_1L_lgl = T,
                                    force_min_max_1L_lgl = T,
                                    id_var_nm_1L_chr = NULL,
                                    model_mdl = NULL,
                                    #impute_1L_lgl = T,
                                    new_data_is_1L_chr = "Simulated",
                                    predr_vars_nms_chr = NULL,
                                    round_var_nm_1L_chr = "Timepoint",
                                    round_bl_val_1L_chr = "BL",
                                    utl_cls_fn = NULL,
                                    utl_var_nm_1L_chr = NULL){
  if(is.null(model_mdl))
    model_mdl <- get_table_predn_mdl(mdl_nm_1L_chr,
                                     ingredients_ls = ingredients_ls,
                                     analysis_1L_chr = analysis_1L_chr)
  mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr)
  tfmn_1L_chr <- ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                          match_value_xx = mdl_type_1L_chr,
                                          match_var_nm_1L_chr = "short_name_chr",
                                          target_var_nm_1L_chr = "tfmn_chr",
                                          evaluate_1L_lgl = F)
  predn_type_1L_chr <- ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                                match_value_xx = mdl_type_1L_chr,
                                                match_var_nm_1L_chr = "short_name_chr",
                                                target_var_nm_1L_chr = "predn_type_chr",
                                                evaluate_1L_lgl = F)
  if(is.na(predn_type_1L_chr))
    predn_type_1L_chr <- NULL
  id_var_nm_1L_chr <- ifelse(is.null(id_var_nm_1L_chr),
                             ingredients_ls$id_var_nm_1L_chr,
                             id_var_nm_1L_chr)
  if(!is.null(predr_vars_nms_chr)){
    data_tb <- rename_from_nmd_vec(data_tb,
                                   nmd_vec_chr = predr_vars_nms_chr,
                                   vec_nms_as_new_1L_lgl = T)
  }
  mdl_predr_terms_chr <- ingredients_ls$mdls_lup %>%
    dplyr::filter(mdl_nms_chr == mdl_nm_1L_chr) %>%
    dplyr::pull(predrs_ls) %>%
    purrr::flatten_chr()
  original_ds_vars_chr <- names(data_tb)[!names(data_tb) %in% c(mdl_predr_terms_chr,
                                                                ifelse(!is.null(utl_var_nm_1L_chr),
                                                                       utl_var_nm_1L_chr,
                                                                       ingredients_ls$depnt_var_nm_1L_chr))]
  updated_tb <- data_tb %>%
    transform_ds_to_predn_ds(predr_vars_nms_chr = mdl_predr_terms_chr,
                             tfmn_1L_chr = tfmn_1L_chr,
                             depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                             id_var_nm_1L_chr = id_var_nm_1L_chr,
                             round_var_nm_1L_chr = round_var_nm_1L_chr,
                             round_bl_val_1L_chr = round_bl_val_1L_chr,
                             predictors_lup = ingredients_ls$predictors_lup) %>%
    add_utility_predn_to_ds(model_mdl = model_mdl,
                            tfmn_1L_chr = tfmn_1L_chr,
                            depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                            predr_vars_nms_chr = mdl_predr_terms_chr,
                            force_min_max_1L_lgl = force_min_max_1L_lgl,
                            force_new_data_1L_lgl = T,
                            impute_1L_lgl = T, # Redundant?
                            is_brms_mdl_1L_lgl = inherits(model_mdl,"brmsfit"),
                            new_data_is_1L_chr = new_data_is_1L_chr,
                            predn_type_1L_chr = NULL,
                            rmv_tfd_depnt_var_1L_lgl = T,
                            utl_cls_fn = utl_cls_fn,
                            utl_min_val_1L_dbl = ingredients_ls$utl_min_val_1L_dbl,
                            sd_dbl = get_random_intercept(ingredients_ls$mdls_smry_tb,
                                                          mdl_nm_1L_chr = mdl_nm_1L_chr,
                                                          deterministic_1L_lgl = deterministic_1L_lgl))
  
  if(!is.null(utl_var_nm_1L_chr)){
    updated_tb <- updated_tb %>%
      dplyr::rename(!!rlang::sym(utl_var_nm_1L_chr):=tidyselect::all_of(ingredients_ls$depnt_var_nm_1L_chr))
  }
  if(!is.null(names(predr_vars_nms_chr))){
    updated_tb <- rename_from_nmd_vec(updated_tb,
                                      nmd_vec_chr = predr_vars_nms_chr,
                                      vec_nms_as_new_1L_lgl = F)
  }
  names_to_inc_chr <- c(names(updated_tb),
                        setdiff(names(data_tb),
                                names(updated_tb)))
  rename_tb <- make_uid_rename_lup(data_tb,
                                   id_var_nm_1L_chr = id_var_nm_1L_chr)
  updated_tb <- dplyr::left_join(data_tb %>%
                                   dplyr::select(tidyselect::all_of(original_ds_vars_chr)) %>%
                                   transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr,
                                                     rename_tb = rename_tb),
                                 updated_tb) %>%
    transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr,
                      rename_tb = rename_tb,
                      old_new_chr = c("new_id_int","old_id_xx"))
  return(updated_tb)
}
authorData_SpecificMixed <- function(x,
                                     depnt_var_min_val_1L_dbl = numeric(0),
                                     title_1L_chr = "An R model object",
                                     what_1L_chr = "Shareable"){
  if(what_1L_chr == "Shareable"){
    results_ls <- purrr::map(manufacture(x@c_SpecificResults,
                                         what_1L_chr = "indexed_shareable"),
                             ~{
                               outp_smry_ls <- append(procureSlot(x,
                                                                  "c_SpecificResults@b_SpecificPrivate@private_outp_ls"),
                                                      .x)
                               outp_smry_ls <- outp_smry_ls %>%
                                 write_shareable_mdls(new_dir_nm_1L_chr = "G_Shareable",
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
author_SpecificModels <- function(x,
                                  prefd_mdl_types_chr = NULL,
                                  what_1L_chr = "all",
                                  digits_1L_int = 3L,
                                  reference_1L_int = NULL,
                                  x_labels_chr = character(0)){
  series_1L_lgl <- x@a_YouthvarsProfile %>% inherits("YouthvarsSeries")
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
      descv_plts_paths_ls <- write_descv_plots(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb, #add youthvars:: when Exporting
                                               ds_descvs_ls = ds_descvs_ls,
                                               descv_outp_dir_1L_chr = x@b_SpecificParameters@paths_ls$descv_outp_dir_1L_chr,
                                               lbl_nms_chr = x@b_SpecificParameters@itm_labels_chr, # Should be domain labels
                                               x_labels_chr = x_labels_chr,
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
author_SpecificSynopsis <- function(x,
                                    reference_1L_int = NA_integer_,
                                    type_1L_chr = "Report",
                                    what_1L_chr = "Catalogue",
                                    ...){
  if(what_1L_chr == "Catalogue"){
    outp_smry_ls_ls <- manufacture(x@b_SpecificResults,
                                   what_1L_chr = "indexed_shareable")
    refs_int <- 1:length(outp_smry_ls_ls)
    if(!is.na(reference_1L_int)){
      outp_smry_ls_ls <- outp_smry_ls_ls[reference_1L_int]
      refs_int <- reference_1L_int
    }
    ctlg_nms_chr <- purrr::map2_chr(outp_smry_ls_ls,
                                    refs_int,
                                    ~ {
                                      fl_nm_1L_chr <- paste0("AAA_TTU_MDL_CTG",
                                                             ifelse(.y==1,
                                                                    "",
                                                                    paste0("-",(.y-1))))
                                      authorReport(x %>%
                                                     renewSlot("b_SpecificResults@a_SpecificShareable@shareable_outp_ls",
                                                               .x),
                                                   fl_nm_1L_chr = fl_nm_1L_chr,
                                                   what_1L_chr = "Catalogue",
                                                   ...)
                                      fl_nm_1L_chr
                                    }
    )
  }
}
depict_SpecificSynopsis <- function(x,
                                    axis_text_sclg_1L_dbl = 1.5,
                                    axis_title_sclg_1L_dbl = 2,
                                    base_height_1L_dbl = 13,
                                    base_size_1L_dbl = 30,
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
                                    y_label_1L_chr = " "){
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
    write_path_1L_chr <- paste0(outp_smry_ls$path_to_write_to_1L_chr,
                                "/dens_and_sctr.png")
  }
  if(what_1L_chr == "composite_utl"){
    ds_descvs_ls <- manufacture_SpecificSynopsis(x,#manufacture when exporting
                                                 depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                 what_1L_chr = "ds_descvs_ls")
    outp_smry_ls <- append(x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                           x@b_SpecificResults@b_SpecificPrivate@private_outp_ls)
    maui_domains_col_nms_chr <- x@c_SpecificParameters@domain_labels_chr
    first_plt <- rlang::exec(make_var_by_round_plt,#youthvars::make_var_by_round_plt, when exporting
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
    second_plt <- rlang::exec(make_sub_tot_plts,#youthvars::make_sub_tot_plts, when exporting
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
                                      y_label_1L_chr = y_label_1L_chr))
    legend_ls <- cowplot::get_legend(first_plt)
    
    plt <- cowplot::plot_grid(first_plt +
                                ggplot2::theme(legend.position="none"),
                              second_plt,
                              legend_ls,
                              nrow = 3L,
                              rel_heights = rel_heights_dbl,
                              scale = scale_dbl)
    write_path_1L_chr <- paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr,
                                "/Output/_Descriptives/combined_utl.png")
  }
  if(write_1L_lgl){
    cowplot::save_plot(write_path_1L_chr,
                       plt,
                       base_height = base_height_1L_dbl)
  }
  return(plt)
}
enhance_SpecificSynopsis <- function(x,
                                     depnt_var_nms_chr = NA_character_,
                                     depnt_var_min_val_1L_dbl = numeric(0),
                                     what_1L_chr = "shareable_outp_ls",
                                     with_1L_chr = "results_ls"){
  if(what_1L_chr == "shareable_outp_ls"){
    outp_smry_ls <- x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls
    if(with_1L_chr == "results_ls"){
      outp_smry_ls$results_ls <- manufacture_SpecificSynopsis(x,#manufacture when exporting
                                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                              depnt_var_nms_chr = depnt_var_nms_chr,
                                                              what_1L_chr = "results_ls")
      # outp_smry_ls$results_ls$abstract_args_ls <- x@abstract_args_ls
      
    }
    x <- renewSlot(x,
                   "b_SpecificResults@a_SpecificShareable@shareable_outp_ls",
                   outp_smry_ls)
  }
  return(x)
}
fit_ts_model_with_brm <- function (data_tb,# rename lngl ?
                                   depnt_var_nm_1L_chr, 
                                   family_fn_1L_chr,
                                   id_var_nm_1L_chr,
                                   predr_vars_nms_chr,
                                   backend_1L_chr = getOption("brms.backend", "rstan"), 
                                   control_ls = NULL,
                                   #link_1L_chr = "identity",
                                   is_csnl_1L_lgl = F,
                                   iters_1L_int = 4000L, seed_1L_int = 1000L, prior_ls = NULL)
{
  mdl_ls <- brms::brm(formula = stats::as.formula(paste0(depnt_var_nm_1L_chr,
                                                         " ~ ", purrr::map_chr(predr_vars_nms_chr, ~paste0(.x, 
                                                                                                           ifelse(is_csnl_1L_lgl,
                                                                                                                  ifelse(paste0(.x,"_scaled") %in% names(data_tb),
                                                                                                                         "_scaled + ",
                                                                                                                         " + "),
                                                                                                                  "_baseline + "), 
                                                                                                           ifelse(is_csnl_1L_lgl,"",paste0(.x, "_change + "))
                                                                                                           )) %>% paste0(collapse = ""),
                                                         "(1|", id_var_nm_1L_chr, ")")), backend = backend_1L_chr,
                      data = data_tb, family = eval(parse(text = family_fn_1L_chr)),
                      iter = iters_1L_int, seed = seed_1L_int, prior = prior_ls, control = control_ls)
  return(mdl_ls)
}
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr, what_1L_chr = "any", X_Ready4useDyad = NULL)
{

    signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>%
      purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>%
      purrr::flatten_chr() %>% unique()
    signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in%
                                            signif_vars_chr]
    if(what_1L_chr == "all"){
      signt_covars_chr <- signt_covars_chr[signt_covars_chr %>% purrr::map_lgl(~sum((mdls_with_covars_smry_tb$Significant %>%
                                                                                      purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>%
                                                                                      purrr::flatten_chr()) ==.x)==length(signt_covars_chr))]
    }
    if(!is.null(X_Ready4useDyad)){
      dummys_chr <- manufacture_Ready4useDyad(X_Ready4useDyad, flatten_1L_lgl = T, type_1L_chr = "dummys", what_1L_chr = "factors") 
      signt_dumys_ls <- mdls_with_covars_smry_tb$Significant %>%
        purrr::map(~{
          terms_1L_chr <- .x
          dummys_chr[dummys_chr %>%
            purrr::map_lgl(~stringr::str_detect(terms_1L_chr,.x))]
        })
      signt_dumys_chr <- signt_dumys_ls %>% purrr::flatten_chr() %>% unique()
      if(what_1L_chr == "all" && !identical(signt_dumys_chr, character(0))){
        signt_dumys_chr  <- signt_dumys_chr[signt_dumys_chr %>% purrr::map_lgl(~sum((signt_dumys_ls %>% purrr::flatten_chr())==.x)==length(signt_dumys_chr))]
      }
      signt_fctrs_chr <- signt_dumys_chr  %>% 
        purrr::map_chr(~manufacture_Ready4useDyad(X_Ready4useDyad, flatten_1L_lgl = T, type_1L_chr = "dummys", what_1L_chr = "factors-d", match_1L_chr = .x)) %>% # update when exported
        unique()
      signt_covars_chr <- c(signt_covars_chr, signt_fctrs_chr) %>% sort()
      }
    if(identical(signt_covars_chr, character(0))){
      signt_covars_chr <- NA_character_
    }
  return(signt_covars_chr)
}
get_mdls_with_signft_covars <- function(outp_smry_ls, 
                                        params_ls_ls){
  signft_covars_chr <- outp_smry_ls$mdls_with_covars_smry_tb %>%
    get_signft_covars(covar_var_nms_chr = params_ls_ls$params_ls$candidate_covar_nms_chr) # (Maybe) Needs editing to account for dummy variables - Need to check.  Would then need to update make_results_ls_spine
  signft_vars_ls <- outp_smry_ls[["mdls_with_covars_smry_tb"]]$Significant %>%
    purrr::map(~strsplit(.x, " ")) %>% purrr::flatten()
  mdls_with_signft_covars_ls <- signft_covars_chr %>%
    purrr::map(~{
      covar_nm_1L_chr <- .x
      mdls_chr <- outp_smry_ls$mdls_with_covars_smry_tb %>%
        dplyr::filter(purrr::map_lgl(signft_vars_ls,
                                     ~ any(.x == covar_nm_1L_chr))) %>%
        dplyr::pull(variable)
      mdls_chr
      
    }) %>%
    stats::setNames(signft_covars_chr)
  return(mdls_with_signft_covars_ls)
}
get_table_predn_mdl <- function(mdl_nm_1L_chr,
                                ingredients_ls,
                                analysis_1L_chr = NULL){
  mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr,
                                          mdl_types_lup = ingredients_ls$mdl_types_lup)
  tfmn_1L_chr <- ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                          match_value_xx = mdl_type_1L_chr,
                                          match_var_nm_1L_chr = "short_name_chr",
                                          target_var_nm_1L_chr = "tfmn_chr",
                                          evaluate_1L_lgl = F)
  if(is.null(analysis_1L_chr)){
    fake_ds_tb <- ingredients_ls$fake_ds_tb
  }else{
    reference_1L_chr <- ifelse(analysis_1L_chr == "Primary Analysis",
                               "Primary",
                               paste0("secondary_",which(LETTERS == stringr::str_sub(analysis_1L_chr,start=-1))))
    fake_ds_tb <- ingredients_ls %>% purrr::pluck(reference_1L_chr) %>% purrr::pluck("fake_ds_tb")
  }
  fake_ds_tb <- fake_ds_tb %>%
    add_tfd_var_to_ds(depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                      tfmn_1L_chr = tfmn_1L_chr)
  table_predn_mdl <- make_shareable_mdl(fake_ds_tb = fake_ds_tb,
                                        mdl_smry_tb = ingredients_ls$mdls_smry_tb %>% dplyr::filter(Model == mdl_nm_1L_chr),
                                        depnt_var_nm_1L_chr = ingredients_ls$depnt_var_nm_1L_chr,
                                        id_var_nm_1L_chr = ingredients_ls$id_var_nm_1L_chr,
                                        tfmn_1L_chr = tfmn_1L_chr,
                                        mdl_type_1L_chr = mdl_type_1L_chr,
                                        mdl_types_lup = ingredients_ls$mdl_types_lup,
                                        control_1L_chr = ready4::get_from_lup_obj(ingredients_ls$mdl_types_lup,
                                                                                  match_value_xx = mdl_type_1L_chr,
                                                                                  match_var_nm_1L_chr = "short_name_chr",
                                                                                  target_var_nm_1L_chr = "control_chr",
                                                                                  evaluate_1L_lgl = F),
                                        start_1L_chr = NA_character_,
                                        seed_1L_int = ingredients_ls$seed_1L_int)
  return(table_predn_mdl)
}
investigate_SpecificFixed <- function(x,
                                      combinations_1L_lgl = F, 
                                      depnt_var_min_val_1L_dbl = numeric(0),
                                      existing_predrs_ls = NULL, 
                                      max_nbr_of_covars_1L_int = integer(0)){
  results_ls <- write_mdls_with_covars_cmprsn(scored_data_tb = x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                              bl_tb = x@c_SpecificResults@b_SpecificPrivate@private_outp_ls$bl_tb,
                                              combinations_1L_lgl = combinations_1L_lgl, 
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              ds_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$ds_smry_ls,
                                              mdl_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$mdl_smry_ls,
                                              existing_predrs_ls = existing_predrs_ls, 
                                              max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                                              output_data_dir_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                                              seed_1L_int = x@b_SpecificParameters@seed_1L_int)
  rename_lup <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$rename_lup
  session_ls <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$session_ls
  x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(results_ls[-1],
                                                                      list(rename_lup = rename_lup,
                                                                           session_ls = session_ls)) # EDIT TO REMOVE INPUTS
  x@c_SpecificResults@b_SpecificPrivate@private_outp_ls <- results_ls[1]
  x_SpecificMixed <- SpecificMixed(a_YouthvarsProfile = x@a_YouthvarsProfile,
                                   b_SpecificParameters = x@b_SpecificParameters,
                                   c_SpecificResults = x@c_SpecificResults,
                                   paths_chr = x@paths_chr,
                                   dissemination_1L_chr = x@dissemination_1L_chr)
  return(x_SpecificMixed)
}
investigate_SpecificMixed <- function(x,
                                      backend_1L_chr = "cmdstanr",
                                      combinations_1L_lgl = F, 
                                      depnt_var_min_val_1L_dbl = numeric(0),
                                      existing_predrs_ls = NULL, 
                                      max_nbr_of_covars_1L_int = integer(0),
                                      new_dir_nm_1L_chr = "F_TS_Mdls",
                                      scndry_anlys_params_ls = NULL){
  if(identical(x@b_SpecificParameters@prior_ls,list(list()))){
    prior_ls <- NULL
  }else{
    prior_ls <- x@b_SpecificParameters@prior_ls
  }
  if(identical(x@b_SpecificParameters@control_ls,list(list()))){
    control_ls <- NULL
  }else{
    control_ls <- x@b_SpecificParameters@control_ls
  }
  if(is.null(scndry_anlys_params_ls)){
    results_ls <- write_ts_mdls_from_alg_outp(outp_smry_ls = append(x@c_SpecificResults@b_SpecificPrivate@private_outp_ls,
                                                                    x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls),
                                              combinations_1L_lgl = combinations_1L_lgl, 
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              existing_predrs_ls = existing_predrs_ls, 
                                              max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                                              predictors_lup = x@b_SpecificParameters@predictors_lup,
                                              utl_min_val_1L_dbl = x@b_SpecificParameters@depnt_var_min_max_dbl[1],# Change
                                              backend_1L_chr = backend_1L_chr,
                                              new_dir_nm_1L_chr = new_dir_nm_1L_chr, # Method Arg
                                              iters_1L_int = x@b_SpecificParameters@iters_1L_int,
                                              path_to_write_to_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,#x@a_Ready4showPaths@outp_data_dir_1L_chr,
                                              prior_ls = prior_ls,
                                              control_ls = control_ls)
    rename_lup <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$rename_lup
    session_ls <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$session_ls
    x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(results_ls[-1],
                                                                        list(rename_lup = rename_lup,
                                                                             session_ls = session_ls)) # EDIT TO REMOVE INPUTS
    x@c_SpecificResults@b_SpecificPrivate@private_outp_ls <- results_ls[1]
  }else{
    input_params_ls <- manufacture_SpecificProject(x, # change to manufacture
                                                   what_1L_chr = "input_params_ls")
    input_params_ls$rename_lup <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$rename_lup
    input_params_ls$scndry_anlys_params_ls <- scndry_anlys_params_ls
    input_params_ls$path_params_ls$paths_ls <- list(write_to_dir_nm_1L_chr = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$path_to_write_to_1L_chr %>%
                                                      stringr::str_sub(end=-8))
    # changed from x@paths_chr
    input_params_ls$outp_smry_ls <- append(x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                                           x@c_SpecificResults@b_SpecificPrivate@private_outp_ls)
    input_params_ls$params_ls$control_ls <- control_ls
    input_params_ls$params_ls$prior_ls <- prior_ls
    input_params_ls$params_ls$iters_1L_int <- x@b_SpecificParameters@iters_1L_int
    results_ls_ls <- write_secondary_analyses(input_params_ls,
                                              backend_1L_chr = backend_1L_chr,
                                              combinations_1L_lgl = combinations_1L_lgl, 
                                              existing_predrs_ls = existing_predrs_ls, 
                                              max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              new_dir_nm_1L_chr = new_dir_nm_1L_chr) %>%
      stats::setNames(names(scndry_anlys_params_ls))
    x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                                                                        results_ls_ls %>%
                                                                          purrr::map(~.x[-1]))
  }
  return(x)
}
investigate_SpecificModels <- function(x,
                                       depnt_var_max_val_1L_dbl = Inf,
                                       depnt_var_min_val_1L_dbl = 0.00001, 
                                       session_ls = NULL){
  results_ls <- write_mdl_cmprsn(scored_data_tb = x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                 depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                 depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl, 
                                 ds_smry_ls = manufacture_SpecificProject(x, # Replace with manufacture when Exporting
                                                          what_1L_chr = "ds_smry_ls"),
                                 mdl_smry_ls = manufacture_SpecificProject(x, # Replace with manufacture when Exporting
                                                           what_1L_chr = "mdl_smry_ls"),
                                 output_data_dir_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                                 seed_1L_int = x@b_SpecificParameters@seed_1L_int)
  rename_lup <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$rename_lup
  x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(results_ls[-1],
                                                                      list(rename_lup = rename_lup,
                                                                           session_ls = session_ls))
  x@c_SpecificResults@b_SpecificPrivate@private_outp_ls <- results_ls[1]
  
  x_SpecificPredictors <- SpecificPredictors(a_YouthvarsProfile = x@a_YouthvarsProfile,
                                             b_SpecificParameters = x@b_SpecificParameters,
                                             c_SpecificResults = x@c_SpecificResults,
                                             paths_chr = x@paths_chr,
                                             dissemination_1L_chr = x@dissemination_1L_chr)
  return(x_SpecificPredictors)
}
investigate_SpecificPredictors <- function(x,
                                           depnt_var_min_val_1L_dbl = numeric(0),
                                           signft_covars_cndn_1L_chr = "any"){
  results_ls <- write_predr_and_covars_cmprsn(scored_data_tb = x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                              bl_tb = x@c_SpecificResults@b_SpecificPrivate@private_outp_ls$bl_tb,
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              ds_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$ds_smry_ls,
                                              mdl_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$mdl_smry_ls,
                                              output_data_dir_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                                              seed_1L_int = x@b_SpecificParameters@seed_1L_int,
                                              signft_covars_cndn_1L_chr = signft_covars_cndn_1L_chr)
  rename_lup <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$rename_lup
  session_ls <- x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$session_ls
  x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(results_ls[-1],
                                                                      list(rename_lup = rename_lup,
                                                                           session_ls = session_ls))
  x@c_SpecificResults@b_SpecificPrivate@private_outp_ls <- results_ls[1]
  x_SpecificFixed <- SpecificFixed(a_YouthvarsProfile = x@a_YouthvarsProfile,
                                   b_SpecificParameters = x@b_SpecificParameters,
                                   c_SpecificResults = x@c_SpecificResults,
                                   paths_chr = x@paths_chr,
                                   dissemination_1L_chr = x@dissemination_1L_chr)
  return(x_SpecificFixed)
}
make_abstract_args_ls <- function(results_ls, # Rename - duplicate name from ready4show
                                  fl_nm_1L_chr = "abstract.txt"){
  mdl_cmprsns_ls <- get_mdl_cmprsns(results_ls, as_list_1L_lgl = T)
  abstract_args_ls <- list(abstract_ls = list(Background = get_background_text(results_ls),
                                              Objectives = paste0("We aimed to identify the best regression models to predict ",
                                                                  get_hlth_utl_nm(results_ls, short_nm_1L_lgl = F),
                                                                  " (",
                                                                  get_hlth_utl_nm(results_ls),
                                                                  ") utility and evaluate the predictive ability of ",
                                                                  get_nbr_of_predrs(results_ls),
                                                                  " candidate measure",
                                                                  ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1,"s",""),
                                                                  " of ",
                                                                  get_predr_ctgs(results_ls),
                                                                  "."),
                                              Methods = paste0(results_ls$study_descs_ls$sample_desc_1L_chr,
                                                               ifelse(is.na(results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr),
                                                                      "",
                                                                      paste0(" Follow-up measurements were ",
                                                                             results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr,
                                                                             " after baseline. ")),
                                                               paste0(length(mdl_cmprsns_ls$OLS) %>%
                                                                        xfun::numbers_to_words() %>%
                                                                        Hmisc::capitalize()),
                                                               " ordinary least squares (OLS) and ",
                                                               length(mdl_cmprsns_ls$GLM) %>%
                                                                 xfun::numbers_to_words(),
                                                               " generalised linear models (GLMs) were explored to identify the best algorithm. ",
                                                               " Predictive ability of ",
                                                               get_nbr_of_predrs(results_ls),
                                                               " candidate measure",
                                                               ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1,"s",""),
                                                               " of ",
                                                               get_predr_ctgs(results_ls),
                                                               " were assessed using ten fold cross validation",
                                                               ifelse(get_nbr_of_predrs(results_ls, as_words_1L_lgl = F)>1," and forest models",""),
                                                               ifelse(is.na(results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr),
                                                                      ". ",
                                                                      paste0(". Linear / generalised linear mixed effect models were then used to construct longitudinal predictive models for ",
                                                                             get_hlth_utl_nm(results_ls),
                                                                             " change."))),
                                              Results = paste0(make_ten_fold_text(results_ls, for_abstract_1L_lgl = T),
                                                               #make_random_forest_text(results_ls, for_abstract_1L_lgl = T),
                                                               ". ",
                                                               make_selected_mdl_text(results_ls, for_abstract_1L_lgl = T),
                                                               ifelse(is.na(results_ls$study_descs_ls$time_btwn_bl_and_fup_1L_chr),
                                                                      "",
                                                                      paste0(" The mean ratio between the within-person and between-person associated coefficients was ",
                                                                             make_within_between_ratios_text(results_ls,
                                                                                                             exclude_covars_1L_lgl = T),
                                                                             "."))
                                                               
                                                               ),
                                              Conclusions = get_conclusion_text(results_ls),
                                              Data = make_data_availability_text(results_ls)),
                           fl_nm_1L_chr = fl_nm_1L_chr)
  return(abstract_args_ls)
}
make_brms_mdl_plt <- function(outp_smry_ls,
                              depnt_var_min_val_1L_dbl = numeric(0),
                              depnt_var_desc_1L_chr,
                              mdl_nm_1L_chr,
                              type_1L_chr,
                              base_size_1L_dbl = 8,
                              brms_mdl = NULL,
                              correspondences_lup = NULL,
                              new_var_nm_1L_chr = "Predicted",
                              predn_type_1L_chr = NULL,
                              x_lbl_1L_chr = NA_character_,
                              y_lbl_1L_chr = NA_character_){
  sfx_1L_chr <- " from brmsfit"
  mdl_types_lup <- outp_smry_ls$mdl_types_lup
  if(is.null(brms_mdl)){
    incld_mdl_paths_chr <- make_incld_mdl_paths(outp_smry_ls)
    brms_mdl <- readRDS(paste0(outp_smry_ls$path_to_write_to_1L_chr,"/",
                               incld_mdl_paths_chr[incld_mdl_paths_chr %>% endsWith(paste0(mdl_nm_1L_chr,".RDS"))]))
    
  }
  mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr,
                                          mdl_types_lup = mdl_types_lup)
  tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                          match_value_xx = mdl_type_1L_chr,
                                          match_var_nm_1L_chr = "short_name_chr",
                                          target_var_nm_1L_chr = "tfmn_chr",
                                          evaluate_1L_lgl = F)
  plot_fn_and_args_ls <- make_plot_fn_and_args_ls(brms_mdl = brms_mdl,
                                                  tfd_data_tb = outp_smry_ls$scored_data_tb %>%
                                                    transform_tb_to_mdl_inp(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,##
                                                                            depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                                                                            predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique(),
                                                                            id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr,
                                                                            round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                                                                            round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr,
                                                                            scaling_fctr_dbl = make_scaling_fctr_dbl(outp_smry_ls)),
                                                  base_size_1L_dbl = base_size_1L_dbl,
                                                  correspondences_lup = correspondences_lup,
                                                  depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                                                  depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                                  new_var_nm_1L_chr = new_var_nm_1L_chr,
                                                  #mdl_nm_1L_chr = mdl_nm_1L_chr,
                                                  predn_type_1L_chr = predn_type_1L_chr,
                                                  round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                                                  sd_dbl = NA_real_,
                                                  sfx_1L_chr = sfx_1L_chr,
                                                  tfmn_1L_chr = tfmn_1L_chr,
                                                  type_1L_chr = type_1L_chr,
                                                  utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl),
                                                                              outp_smry_ls$utl_min_val_1L_dbl,
                                                                              -1),
                                                  x_lbl_1L_chr = x_lbl_1L_chr,
                                                  y_lbl_1L_chr = y_lbl_1L_chr)
  plt <- rlang::exec(plot_fn_and_args_ls$plt_fn, !!!plot_fn_and_args_ls$fn_args_ls)
  return(plt)
}
make_cmpst_sctr_and_dnst_plt <- function(outp_smry_ls,
                                         output_data_dir_1L_chr = NA_character_,
                                         predr_var_nms_chr = NA_character_,
                                         base_size_1L_dbl = 16,
                                         correspondences_lup = NULL,
                                         depnt_var_min_val_1L_dbl = numeric(0),
                                         depnt_var_desc_1L_chr = NA_character_,
                                         labels_chr = c("A","B","C","D"),
                                         label_x_1L_dbl = 0.1,
                                         label_y_1L_dbl = 0.9,
                                         label_size_1L_dbl = 22,
                                         mdl_indcs_int = 1:2,
                                         use_png_fls_1L_lgl = T){
  if(use_png_fls_1L_lgl){
    filtered_paths_chr <- outp_smry_ls$file_paths_chr %>% purrr::discard(~endsWith(.x,"_sim_sctr.png")|endsWith(.x,"_sim_dnst.png")|endsWith(.x,"_cnstrd_sctr_plt.png")|endsWith(.x,"_cnstrd_dnst.png"))
    filtered_paths_chr <- paste0(output_data_dir_1L_chr,"/",filtered_paths_chr[filtered_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,paste0(predr_var_nms_chr,"_1")) & (stringr::str_detect(.x,"_dnst.png") | stringr::str_detect(.x,"_sctr_plt.png")))])
    mdl_types_chr <- filtered_paths_chr %>% purrr::map_chr(~DescTools::SplitPath(.x)$filename %>%
                                                             stringr::str_remove("_dnst") %>%
                                                             stringr::str_remove("_sctr_plt") %>%
                                                             get_mdl_type_from_nm())
    ordered_paths_chr <- outp_smry_ls$prefd_mdl_types_chr %>%
      purrr::map(~filtered_paths_chr[which(mdl_types_chr==.x)]) %>%
      purrr::flatten_chr()
    plot_ls <- ordered_paths_chr %>% purrr::map(~cowplot::ggdraw() + cowplot::draw_image(.x))
  }else{
    plots_chr <- outp_smry_ls$mdl_nms_ls %>%
      purrr::flatten_chr()
    plots_chr <- plots_chr[mdl_indcs_int]
    # %>%
    #   `[`(mdl_indcs_int)
    plot_ls <- plots_chr %>%
      purrr::map(~{
        mdl_nm_1L_chr <- .x
        brms_mdl <- get_brms_mdl(outp_smry_ls,
                                 mdl_nm_1L_chr = mdl_nm_1L_chr)
        if(is.na(depnt_var_desc_1L_chr)){
          depnt_var_desc_1L_chr <- get_hlth_utl_nm(outp_smry_ls$results_ls,
                                                   short_nm_1L_lgl = T)
        }
        purrr::map(c("dnst", "sctr_plt"),
                   ~ {
                     make_brms_mdl_plt(outp_smry_ls,
                                       base_size_1L_dbl = base_size_1L_dbl,
                                       brms_mdl = brms_mdl,
                                       correspondences_lup = correspondences_lup,
                                       depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,#
                                       depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                       mdl_nm_1L_chr = mdl_nm_1L_chr,
                                       type_1L_chr = .x,
                                       predn_type_1L_chr = NULL,
                                       x_lbl_1L_chr = paste0("Observed ", depnt_var_desc_1L_chr),
                                       y_lbl_1L_chr = paste0("Predicted ", depnt_var_desc_1L_chr))
                   })
      }) %>%
      purrr::flatten()
    
  }
  composite_plt <- cowplot::plot_grid(plot_ls[[1]],plot_ls[[2]],plot_ls[[3]],plot_ls[[4]],
                                      nrow = 2,
                                      labels = labels_chr,
                                      label_x = label_x_1L_dbl,
                                      label_y = label_y_1L_dbl,
                                      label_size = label_size_1L_dbl)
  return(composite_plt)
}
make_csnl_example_predrs <- function(){ # New to specific
  predictors_r3 <- Ready4useRepos(dv_nm_1L_chr = "TTU", 
                                  dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0", 
                                  dv_server_1L_chr = "dataverse.harvard.edu") %>%
    ingest(fls_to_ingest_chr = c("predictors_r3"),
           metadata_1L_lgl = F)
  predictors_r3 <- renew.specific_predictors(predictors_r3,
                                             filter_cdn_1L_chr = "short_name_chr == 'SOFAS'") %>%
    renew.specific_predictors(short_name_chr = c("K10", "MLT", "CHU9D", "AQOL6D"),
                              long_name_chr = c("K10 total score", "MLT total score", "CHU9D health utility", "AQOL6D health utility"),
                              min_val_dbl = c(10,0,-0.1059,0.03),
                              max_val_dbl = c(50,100,1,1),
                              class_chr = c("integer","numeric","numeric","numeric"),
                              increment_dbl = 1,
                              class_fn_chr = c("as.integer","as.double","as.double","youthvars::youthvars_aqol6d_adol"), # update when new youthvars classes are created
                              mdl_scaling_dbl = c(0.01,0.01,1,1),
                              covariate_lgl = F) %>%
    dplyr::mutate(covariate_lgl = dplyr::case_when(short_name_chr == "SOFAS" ~ F,
                                                   T ~ covariate_lgl))
  return(predictors_r3)
}
make_fake_ts_data <- function (outp_smry_ls, # rename lngl
                               depnt_var_min_val_1L_dbl = numeric(0),
                               depnt_vars_are_NA_1L_lgl = T)
{
  data_tb <- outp_smry_ls$scored_data_tb %>% transform_tb_to_mdl_inp(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,#
                                                                     depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                                                                     predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique(),
                                                                     id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                                                                     round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
  if(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)) | ifelse(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)),T,is.na(outp_smry_ls$round_var_nm_1L_chr))){
    data_tb <- data_tb %>%
      dplyr::select(-(outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique() %>% paste0("_change")))
  }
  fk_data_ls <- synthpop::syn(data_tb, visit.sequence = names(data_tb)[names(data_tb) !=
                                                                         outp_smry_ls$id_var_nm_1L_chr], seed = outp_smry_ls$seed_1L_int)
  fk_data_tb <- fk_data_ls$syn
  if(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)) | ifelse(identical(outp_smry_ls$round_var_nm_1L_chr, character(0)),T,is.na(outp_smry_ls$round_var_nm_1L_chr))){
    fk_data_tb <- fk_data_tb %>%
      dplyr::ungroup()
  } else {
    fk_data_tb <- fk_data_tb %>%
      dplyr::mutate(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) := as.character(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr))) %>%
      dplyr::group_by(!!rlang::sym(outp_smry_ls$id_var_nm_1L_chr)) %>%
      dplyr::mutate(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) := !!rlang::sym(outp_smry_ls$round_var_nm_1L_chr) %>%
                      transform_timepoint_vals(timepoint_levels_chr = outp_smry_ls$scored_data_tb %>%
                                                 dplyr::pull(!!rlang::sym(outp_smry_ls$round_var_nm_1L_chr)) %>%
                                                 unique(),
                                               bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr)
      )%>%
      dplyr::ungroup()
  }
  if(depnt_vars_are_NA_1L_lgl){
    depnt_vars_chr <- names(fk_data_tb)[names(fk_data_tb) %>%
                                          purrr::map_lgl(~startsWith(.x, outp_smry_ls$depnt_var_nm_1L_chr))]
    fk_data_tb <- fk_data_tb %>% dplyr::mutate(dplyr::across(dplyr::all_of(depnt_vars_chr),
                                                             ~NA_real_))
  }
  return(fk_data_tb)
}
make_hlth_utl_and_predrs_ls <- function(outp_smry_ls, # Generalise from HU
                                        descv_tbls_ls,
                                        nbr_of_digits_1L_int = 2L,
                                        old_nms_chr = NULL,
                                        new_nms_chr = NULL){
  ranked_predrs_ls <- make_ranked_predrs_ls(descv_tbls_ls,
                                            old_nms_chr = old_nms_chr,
                                            new_nms_chr = new_nms_chr)
  var_nm_1L_chr <- descv_tbls_ls$ds_descvs_ls$dictionary_tb %>%
    ready4::get_from_lup_obj(match_var_nm_1L_chr = "var_nm_chr",
                             match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr,
                             target_var_nm_1L_chr = "var_desc_chr",
                             evaluate_1L_lgl = F) %>% as.vector()
  
  
  hlth_utl_and_predrs_ls = list(bl_hu_mean_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                           match_value_xx = var_nm_1L_chr,
                                                           target_var_nm_1L_chr = paste0(descv_tbls_ls$ds_descvs_ls$round_vals_chr[1],
                                                                                         "_val_1_dbl"),
                                                           evaluate_1L_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                bl_hu_sd_1L_dbl = descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                           match_value_xx = var_nm_1L_chr,
                                                           target_var_nm_1L_chr = paste0(descv_tbls_ls$ds_descvs_ls$round_vals_chr[1],
                                                                                         "_val_2_ls"),
                                                           evaluate_1L_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int),
                                fup_hu_mean_1L_dbl = ifelse(length(descv_tbls_ls$ds_descvs_ls$round_vals_chr)<2,
                                                            NA_real_,
                                                            descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                           match_value_xx = var_nm_1L_chr,
                                                           target_var_nm_1L_chr = paste0(descv_tbls_ls$ds_descvs_ls$round_vals_chr[2],
                                                                                         "_val_1_dbl"),
                                                           evaluate_1L_lgl = F) %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int)),
                                fup_hu_sd_1L_dbl = ifelse(length(descv_tbls_ls$ds_descvs_ls$round_vals_chr)<2,
                                                          NA_real_,
                                                          descv_tbls_ls$main_outc_tbl_tb %>%
                                  dplyr::filter(label == "Mean (SD)") %>%
                                  ready4::get_from_lup_obj(match_var_nm_1L_chr = "variable",
                                                           match_value_xx = var_nm_1L_chr,
                                                           target_var_nm_1L_chr = paste0(descv_tbls_ls$ds_descvs_ls$round_vals_chr[2],
                                                                                         "_val_2_ls"),
                                                           evaluate_1L_lgl = F) %>%
                                  stringr::str_remove("\\(") %>%
                                  stringr::str_remove("\\)") %>%
                                  as.numeric() %>%
                                  round(nbr_of_digits_1L_int)),
                                predrs_nartv_seq_chr = ranked_predrs_ls$unranked_predrs_chr,
                                cor_seq_dscdng_chr =  ranked_predrs_ls$ranked_predrs_chr)
  return(hlth_utl_and_predrs_ls)
}
make_input_params <- function(ds_tb, # Generalise MAUI
                              ds_descvs_ls,
                              header_yaml_args_ls,
                              maui_params_ls,
                              predictors_lup,
                              control_ls = NULL,
                              dv_ds_nm_and_url_chr = NULL,
                              iters_1L_int = 4000L,
                              mdl_smry_ls = make_mdl_smry_ls(),
                              output_format_ls = make_output_format_ls(),
                              path_params_ls = NULL,
                              prefd_covars_chr = NULL,
                              prefd_mdl_types_chr = NULL,
                              prior_ls = NULL,
                              seed_1L_int = 12345,
                              scndry_anlys_params_ls = NULL,
                              write_new_dir_1L_lgl = T){
  path_params_ls <- ready4show::make_path_params_ls(use_fake_data_1L_lgl = ds_descvs_ls$is_fake_1L_lgl,
                                                    dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr,
                                                    write_new_dir_1L_lgl = write_new_dir_1L_lgl)
  params_ls_ls <- make_analysis_core_params_ls(ds_descvs_ls = ds_descvs_ls,
                                               output_format_ls = output_format_ls,
                                               predictors_lup = predictors_lup,
                                               prefd_covars_chr = prefd_covars_chr,
                                               prefd_mdl_types_chr = prefd_mdl_types_chr,
                                               mdl_smry_ls = mdl_smry_ls,
                                               control_ls = control_ls,
                                               iters_1L_int = iters_1L_int,
                                               prior_ls = prior_ls,
                                               seed_1L_int = seed_1L_int) %>%
    make_valid_params_ls_ls(ds_tb = ds_tb,
                            maui_params_ls = maui_params_ls,
                            path_params_ls = path_params_ls)
  params_ls_ls$header_yaml_args_ls <- header_yaml_args_ls
  params_ls_ls$output_format_ls <- output_format_ls
  params_ls_ls$scndry_anlys_params_ls <- scndry_anlys_params_ls
  return(params_ls_ls)
}
make_mdl <- function (data_tb,
                      depnt_var_min_val_1L_dbl = numeric(0),
                      depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                      tfmn_1L_chr = "NTF",
                      predr_var_nm_1L_chr, covar_var_nms_chr = NA_character_, mdl_type_1L_chr = "OLS_NTF",
                      mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NULL)
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
  if (is.null(start_1L_chr)) {
    start_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                             match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                             target_var_nm_1L_chr = "start_chr", evaluate_1L_lgl = F)
  }
  if (!is.na(control_1L_chr)) {
    idx_1L_int <- 1 + stringi::stri_locate_last_fixed(mdl_type_1L_chr,
                                                      "_")[1, 1] %>% as.vector()
    link_1L_chr <- get_link_from_tfmn(stringr::str_sub(mdl_type_1L_chr, start = idx_1L_int))
  }
  mdl_1L_chr <- paste0(ready4::get_from_lup_obj(mdl_types_lup,
                                                match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                                target_var_nm_1L_chr = "fn_chr", evaluate_1L_lgl = F), "(",
                       transform_depnt_var_nm(depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr),
                       " ~ ", predr_var_nm_1L_chr, ifelse(is.na(covar_var_nms_chr[1]),
                                                          "", paste0(" + ", paste0(covar_var_nms_chr, collapse = " + "))),
                       ", data = data_tb", ifelse(!is.na(ready4::get_from_lup_obj(mdl_types_lup,
                                                                                  match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                                                                  target_var_nm_1L_chr = "family_chr", evaluate_1L_lgl = F)),
                                                  paste0(", family = ", ready4::get_from_lup_obj(mdl_types_lup,
                                                                                                 match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                                                                                 target_var_nm_1L_chr = "family_chr", evaluate_1L_lgl = F)),
                                                  ""), ifelse((!is.na(start_1L_chr)|(is.na(start_1L_chr) & !is.na(control_1L_chr))), ", ", ""), ifelse(!is.na(control_1L_chr),
                                                                                                                                                       paste0("link=\"", link_1L_chr, "\",control=", control_1L_chr,
                                                                                                                                                              "("), ""), ifelse(!is.na(start_1L_chr), paste0("start=c(",
                                                                                                                                                                                                             start_1L_chr, ")"), ""), ifelse(!is.na(control_1L_chr),
                                                                                                                                                                                                                                             ")", ""), ")")
  model_mdl <- eval(parse(text = mdl_1L_chr))
  return(model_mdl)
}
make_mdl_coef_ratio_ls <- function(mdl_ingredients_ls,
                                   predr_ctgs_ls = NULL){
  predrs_chr <- mdl_ingredients_ls$predictors_lup$short_name_chr
  mdl_type_chr <- mdl_ingredients_ls$mdls_smry_tb$Model %>%
    unique() %>%
    purrr::map_chr(~get_mdl_type_from_nm(.x,
                                         mdl_types_lup = mdl_ingredients_ls$mdl_types_lup)) %>%
    unique()
  main_mdls_ls <- predrs_chr %>% purrr::map(~paste0(paste0(.x,"_1_"),
                                                    mdl_type_chr))
  ratios_ls <- main_mdls_ls %>%
    purrr::map2(predrs_chr, ~{
      mdls_chr <- .x
      predr_1L_chr <- .y
      mdls_chr %>% purrr::map_dbl(~{
        coefs_dbl <- mdl_ingredients_ls$mdls_smry_tb %>% dplyr::filter(Model %in% .x) %>%
          dplyr::filter(Parameter %in% paste0(predr_1L_chr,(c(" baseline", " scaled", " unscaled",
                                                              " change")))) %>%
          dplyr::pull(Estimate)
        coefs_dbl[2]/coefs_dbl[1]
      })
    }
    ) %>% stats::setNames(predrs_chr)
  mdl_coef_ratios_ls = list(mean_ratios_dbl = ratios_ls %>% purrr::map_dbl(~mean(.x)))
  if(!is.null(predr_ctgs_ls)){
    append_ls <- purrr::map(predr_ctgs_ls,
                            ~ mdl_coef_ratios_ls$mean_ratios_dbl[predrs_chr %in% .x]) %>%
      stats::setNames(names(predr_ctgs_ls))
    mdl_coef_ratios_ls <- append(mdl_coef_ratios_ls, append_ls)
  }
  return(mdl_coef_ratios_ls)
}
make_mdl_desc_lines <- function(outp_smry_ls,
                                mdl_nm_1L_chr,
                                output_type_1L_chr = "PDF"){
  mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>%
    dplyr::filter(Model == mdl_nm_1L_chr)
  predictors_chr <- mdl_smry_tb$Parameter[!mdl_smry_tb$Parameter %in% c("SD (Intercept)","Intercept","R2","RMSE","Sigma")] %>%
    purrr::map_chr(~stringr::str_remove(.x," baseline") %>% stringr::str_remove(" change") %>% stringr::str_remove(.x," scaled") %>% stringr::str_remove(" unscaled") ) %>% unique()
  predictors_desc_chr <- predictors_chr %>%
    purrr::map_chr(~{
      scaling_1L_dbl <- ready4::get_from_lup_obj(outp_smry_ls$predictors_lup,
                                                 match_value_xx = .x,
                                                 match_var_nm_1L_chr = "short_name_chr",
                                                 target_var_nm_1L_chr = "mdl_scaling_dbl",
                                                 evaluate_1L_lgl = F)
      paste0(.x,
             " (",
             ready4::get_from_lup_obj(outp_smry_ls$dictionary_tb,
                                      match_value_xx = .x,
                                      match_var_nm_1L_chr = "var_nm_chr",
                                      target_var_nm_1L_chr = "var_desc_chr",
                                      evaluate_1L_lgl = F),
             ifelse(scaling_1L_dbl == 1,
                    "",
                    paste0(" (multiplied by ", scaling_1L_dbl,")")),
             ")")
    })
  if(length(predictors_desc_chr) > 1)
    predictors_desc_chr <- paste0(c(paste0("\n - ",
                                           predictors_desc_chr[-length(predictors_desc_chr)],
                                           collapse = ";"),
                                    paste0("\n - ",predictors_desc_chr[length(predictors_desc_chr)])),
                                  collapse = "; and")
  
  
  mdl_desc_lines_chr <- paste0(paste0("This model predicts values at two timepoints for ",
                                      ready4::get_from_lup_obj(outp_smry_ls$dictionary_tb,
                                                               match_value_xx = outp_smry_ls$depnt_var_nm_1L_chr,
                                                               match_var_nm_1L_chr = "var_nm_chr",
                                                               target_var_nm_1L_chr = "var_desc_chr",
                                                               evaluate_1L_lgl = F),
                                      ". The predictor variables are ",
                                      "baseline values and subsequent changes in ",
                                      collapse = ""), predictors_desc_chr,". ",
                               "The catalogue reference for this model is ",
                               ifelse(output_type_1L_chr == "PDF",
                                      paste0("\\texttt{\\detokenize{",mdl_nm_1L_chr,"}}"),
                                      mdl_nm_1L_chr),
                               ".")
  return(mdl_desc_lines_chr)
}
make_mdl_nms_ls <- function (predr_vars_nms_ls, mdl_types_chr)
{
  mdl_nms_ls <- purrr::map2(predr_vars_nms_ls, make_unique_ls_elmt_idx_int(predr_vars_nms_ls),
                            ~paste0(.x[1], "_", ifelse(is.na(.x[2]), "", paste0(.x[2],
                                                                                "_")), .y, "_", mdl_types_chr))
  return(mdl_nms_ls)
}
make_plot_fn_and_args_ls <- function(type_1L_chr,
                                     depnt_var_desc_1L_chr,
                                     args_ls = NULL,
                                     base_size_1L_dbl = 11,
                                     brms_mdl = NULL,
                                     correspondences_lup = NULL,
                                     depnt_var_nm_1L_chr = NULL,
                                     new_var_nm_1L_chr = NA_character_,
                                     predn_type_1L_chr = NULL,
                                     round_var_nm_1L_chr = NULL,
                                     sd_dbl = NA_real_,
                                     seed_1L_dbl = 23456,
                                     sfx_1L_chr = " from table",
                                     table_predn_mdl = NULL,
                                     tfmn_1L_chr = "NTF",
                                     tfd_data_tb = NULL,
                                     utl_min_val_1L_dbl = -1,
                                     x_lbl_1L_chr = NA_character_,
                                     y_lbl_1L_chr = NA_character_){
  if(!is.null(brms_mdl)){
    set.seed(seed_1L_dbl)
    tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb,
                                                    model_mdl = brms_mdl,
                                                    depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                    is_brms_mdl_1L_lgl = inherits(brms_mdl,"brmsfit"),
                                                    predn_type_1L_chr = predn_type_1L_chr,
                                                    sd_dbl = NA_real_,
                                                    sfx_1L_chr = ifelse(is.null(table_predn_mdl),
                                                                        " from brmsfit",
                                                                        sfx_1L_chr),
                                                    tfmn_1L_chr = tfmn_1L_chr,
                                                    utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    if(!is.null(table_predn_mdl)){
      tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb,
                                                      model_mdl = table_predn_mdl,
                                                      depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                      is_brms_mdl_1L_lgl = F,
                                                      predn_type_1L_chr = predn_type_1L_chr,
                                                      sd_dbl = sd_dbl,
                                                      sfx_1L_chr = ifelse(!is.null(brms_mdl),
                                                                          " from table",
                                                                          sfx_1L_chr),
                                                      tfmn_1L_chr = tfmn_1L_chr,
                                                      utl_min_val_1L_dbl = utl_min_val_1L_dbl)
    }
  }
  ref_idx_1L_int <- which(type_1L_chr == c("coefs", "hetg",
                                           "dnst", "sctr_plt",
                                           "sim_dnst", "sim_sctr",
                                           "cnstrd_dnst", "cnstrd_sctr_plt",
                                           "cnstrd_sim_dnst", "cnstrd_sim_sctr"))
  
  if (ref_idx_1L_int %in% c(3,5,7,9)) {
    plt_fn <- plot_obsd_predd_dnst
    fn_args_ls <- list(tfd_data_tb = tfd_data_tb,
                       base_size_1L_dbl = base_size_1L_dbl,
                       depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                       depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                       new_var_nm_1L_chr = new_var_nm_1L_chr,
                       predd_val_var_nm_1L_chr = ifelse(ref_idx_1L_int %in% c(3,7),
                                                        transform_predd_var_nm("Predicted",
                                                                               sfx_1L_chr = ifelse(!is.null(table_predn_mdl),
                                                                                                   " from brmsfit",
                                                                                                   sfx_1L_chr),
                                                                               utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 3, ####
                                                                                                           NA_real_,
                                                                                                           utl_min_val_1L_dbl)),
                                                        transform_predd_var_nm("Simulated",
                                                                               sfx_1L_chr = ifelse(!is.null(table_predn_mdl),
                                                                                                   " from brmsfit",
                                                                                                   sfx_1L_chr),
                                                                               utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 5,
                                                                                                           NA_real_,
                                                                                                           utl_min_val_1L_dbl))),
                       cmprsn_predd_var_nm_1L_chr = ifelse(is.null(table_predn_mdl),
                                                           NA_character_,
                                                           ifelse(ref_idx_1L_int %in% c(3,7), ##
                                                                  transform_predd_var_nm("Predicted",
                                                                                         sfx_1L_chr = " from table",
                                                                                         utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 3, ##
                                                                                                                     NA_real_,
                                                                                                                     utl_min_val_1L_dbl)),
                                                                  transform_predd_var_nm("Simulated",
                                                                                         sfx_1L_chr =" from table",
                                                                                         utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 5, ##
                                                                                                                     NA_real_,
                                                                                                                     utl_min_val_1L_dbl)))))
  }else{
    plt_fn <- plot_obsd_predd_sctr_cmprsn
    fn_args_ls <- list(tfd_data_tb = tfd_data_tb,
                       base_size_1L_dbl = base_size_1L_dbl,
                       correspondences_lup = correspondences_lup,
                       depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                       depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                       round_var_nm_1L_chr = round_var_nm_1L_chr,
                       predd_val_var_nm_1L_chr = ifelse(ref_idx_1L_int %in% c(4,8),
                                                        transform_predd_var_nm("Predicted",
                                                                               sfx_1L_chr = ifelse(!is.null(table_predn_mdl),
                                                                                                   " from brmsfit",
                                                                                                   sfx_1L_chr),
                                                                               utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 4,
                                                                                                           NA_real_,
                                                                                                           utl_min_val_1L_dbl)),
                                                        transform_predd_var_nm("Simulated",
                                                                               sfx_1L_chr = ifelse(!is.null(table_predn_mdl),
                                                                                                   " from brmsfit",
                                                                                                   sfx_1L_chr),
                                                                               utl_min_val_1L_dbl = ifelse(ref_idx_1L_int == 6,
                                                                                                           NA_real_,
                                                                                                           utl_min_val_1L_dbl))),
                       args_ls = args_ls,
                       x_lbl_1L_chr = x_lbl_1L_chr,
                       y_lbl_1L_chr = y_lbl_1L_chr)
  }
  plot_fn_and_args_ls <- list(plt_fn = plt_fn,
                              fn_args_ls = fn_args_ls)
  return(plot_fn_and_args_ls)
}
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls, combinations_1L_lgl = F, existing_predrs_ls = NULL, max_nbr_of_covars_1L_int = integer(0)) # May need copy in catalogue
{
  predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
    covars_chr <- .x
    purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr) %>% purrr::discard(is.na))) %>%
      purrr::flatten()
  }) %>% purrr::flatten() %>% unique()
  predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls,
                                                      length))]
  if(combinations_1L_lgl){
    main_predrs_ls <- predr_vars_nms_ls[predr_vars_nms_ls %>% purrr::map_lgl(~length(.x)==1)]
    combinations_from_chr <- setdiff(predr_vars_nms_ls[predr_vars_nms_ls %>% purrr::map_lgl(~length(.x)>1)] %>% purrr::flatten_chr(), main_predrs_chr)
    combinations_ls <- 1:(ifelse(identical(max_nbr_of_covars_1L_int,integer(0)),length(combinations_from_chr),min(max_nbr_of_covars_1L_int,length(combinations_from_chr)))) %>% 
      purrr::map(~gtools::combinations(length(combinations_from_chr),.x,combinations_from_chr) %>% t() %>% as.data.frame() %>% as.list() %>% unname()) %>% purrr::flatten()
    combinations_ls <- main_predrs_ls %>%
      purrr::map(~
                   {
                     main_1L_chr <- .x
                     combinations_ls %>%
                       purrr::map(~c(main_1L_chr,.x))
      }) %>% purrr::flatten()
    predr_vars_nms_ls <- append(main_predrs_ls,combinations_ls)
  }
  if(!is.null(existing_predrs_ls)){
    predr_vars_nms_ls <- predr_vars_nms_ls[predr_vars_nms_ls %>% purrr::map_lgl(~{
      test_chr <- .x
      !any(existing_predrs_ls %>% purrr::map_lgl(~identical(.x,test_chr))
      )})]
  }
  return(predr_vars_nms_ls)
}
make_results_ls <- function(spine_of_results_ls = NULL, # CORE OF S4 Classes - rename ts to lngl # and metamorphose methods
                            abstract_args_ls = NULL,
                            depnt_var_min_val_1L_dbl = numeric(0),
                            dv_ds_nm_and_url_chr = NULL,
                            output_format_ls = NULL,
                            params_ls_ls = NULL,
                            path_params_ls = NULL,
                            study_descs_ls = NULL,
                            fn_ls = NULL,
                            include_idx_int = NULL,
                            var_nm_change_lup = NULL,
                            ctgl_vars_regrouping_ls = NULL,
                            make_cmpst_plt_1L_lgl = T,
                            outp_smry_ls = NULL,
                            sig_covars_some_predrs_mdls_tb = NULL,
                            sig_thresh_covars_1L_chr = NULL,
                            version_1L_chr = NULL){
  if(is.null(spine_of_results_ls)){
    spine_of_results_ls <- make_results_ls_spine(output_format_ls = output_format_ls,
                                                 params_ls_ls = params_ls_ls,
                                                 path_params_ls = path_params_ls,
                                                 study_descs_ls = study_descs_ls,
                                                 fn_ls = fn_ls,
                                                 include_idx_int = include_idx_int,
                                                 outp_smry_ls = outp_smry_ls,
                                                 var_nm_change_lup = var_nm_change_lup)
  }
  mdls_smry_tbls_ls <- make_mdls_smry_tbls_ls(spine_of_results_ls$outp_smry_ls,
                                              nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int)
  covars_mdls_ls <- make_mdls_ls(spine_of_results_ls$outp_smry_ls,
                                 mdls_tb = mdls_smry_tbls_ls$covar_mdls_tb)
  descv_tbls_ls <- paste0(spine_of_results_ls$output_data_dir_1L_chr,"/",spine_of_results_ls$outp_smry_ls$file_paths_chr[spine_of_results_ls$outp_smry_ls$file_paths_chr %>% purrr::map_lgl(~stringr::str_detect(.x,"descv_tbls_ls.RDS"))]) %>% readRDS()
  if(make_cmpst_plt_1L_lgl){
    composite_plt <- make_cmpst_sctr_and_dnst_plt(spine_of_results_ls$outp_smry_ls,
                                                  depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                  output_data_dir_1L_chr = spine_of_results_ls$output_data_dir_1L_chr,
                                                  predr_var_nms_chr = spine_of_results_ls$outp_smry_ls$predr_vars_nms_ls[[1]])
    cowplot::save_plot(paste0(spine_of_results_ls$output_data_dir_1L_chr,"/dens_and_sctr.png"), composite_plt, base_height = 20)
    
  }
  ttu_cs_ls <- make_ttu_cs_ls(spine_of_results_ls$outp_smry_ls,
                              sig_covars_some_predrs_mdls_tb = sig_covars_some_predrs_mdls_tb,
                              sig_thresh_covars_1L_chr = sig_thresh_covars_1L_chr)
  mdl_types_chr <- mdls_smry_tbls_ls$prefd_predr_mdl_smry_tb$Model %>%
    purrr::map_chr(~get_mdl_type_from_nm(.x)) %>%
    unique()
  ttu_cs_ls$rf_seq_dscdng_chr <- ttu_cs_ls$rf_seq_dscdng_chr %>%
    purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                           .x %>%
                             ready4::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F),
                           .x))
  ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr <- ttu_cs_ls$cs_mdls_predrs_seq_dscdng_chr  %>%
    purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                           .x %>%
                             ready4::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F),
                           .x))
  ttu_lngl_ls = list(best_mdls_tb = tibble::tibble(model_type = mdl_types_chr %>%
                                                     purrr::map_chr(~ready4::get_from_lup_obj(spine_of_results_ls$outp_smry_ls$mdl_types_lup,
                                                                                              match_value_xx = .x,
                                                                                              match_var_nm_1L_chr = "short_name_chr",
                                                                                              target_var_nm_1L_chr = "mixed_acronym_chr",
                                                                                              evaluate_1L_lgl = F)),
                                                   link_and_tfmn_chr = mdl_types_chr %>%
                                                     purrr::map_chr(~ready4::get_from_lup_obj(spine_of_results_ls$outp_smry_ls$mdl_types_lup,
                                                                                              match_value_xx = .x,
                                                                                              match_var_nm_1L_chr = "short_name_chr",
                                                                                              target_var_nm_1L_chr = "with_chr",
                                                                                              evaluate_1L_lgl = F)),
                                                   name_chr = make_predrs_for_best_mdls(spine_of_results_ls$outp_smry_ls,
                                                                                        old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                                                                        new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr),
                                                   r2_dbl = mdls_smry_tbls_ls$prefd_predr_mdl_smry_tb %>%
                                                     dplyr::filter(Parameter == "R2") %>%
                                                     dplyr::pull(Estimate)),
                     cs_ts_ratios_tb = spine_of_results_ls$cs_ts_ratios_tb,
                     incld_covars_chr = spine_of_results_ls$outp_smry_ls$prefd_covars_chr)
  results_ls <- list(abstract_args_ls = abstract_args_ls,
                     candidate_covars_ls = spine_of_results_ls$candidate_covars_ls,
                     candidate_predrs_chr = spine_of_results_ls$candidate_predrs_chr,
                     cohort_ls = make_cohort_ls(descv_tbls_ls,
                                                ctgl_vars_regrouping_ls = ctgl_vars_regrouping_ls,
                                                nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int),
                     dv_ds_nm_and_url_chr = dv_ds_nm_and_url_chr,
                     header_yaml_args_ls = params_ls_ls$header_yaml_args_ls,
                     hlth_utl_and_predrs_ls = make_hlth_utl_and_predrs_ls(spine_of_results_ls$outp_smry_ls,
                                                                          descv_tbls_ls = descv_tbls_ls,
                                                                          nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int,
                                                                          old_nms_chr = spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                                                          new_nms_chr = spine_of_results_ls$var_nm_change_lup$new_nms_chr),
                     mdl_coef_ratios_ls = spine_of_results_ls$mdl_coef_ratios_ls,
                     mdl_ingredients_ls = spine_of_results_ls$mdl_ingredients_ls,
                     mdls_with_signft_covars_ls = spine_of_results_ls$mdls_with_signft_covars_ls,
                     output_format_ls = params_ls_ls$output_format_ls,
                     path_params_ls = params_ls_ls$path_params_ls,
                     paths_to_figs_ls = make_paths_to_ss_plts_ls(spine_of_results_ls$output_data_dir_1L_chr,
                                                                 outp_smry_ls = spine_of_results_ls$outp_smry_ls),
                     predr_var_nms_chr = spine_of_results_ls$outp_smry_ls$predr_vars_nms_ls[[1]] %>%
                       purrr::map_chr(~ifelse(.x %in% spine_of_results_ls$var_nm_change_lup$old_nms_chr,
                                              .x %>%
                                                ready4::get_from_lup_obj(data_lookup_tb = spine_of_results_ls$var_nm_change_lup,
                                                                         match_var_nm_1L_chr = "old_nms_chr",
                                                                         target_var_nm_1L_chr = "new_nms_chr",
                                                                         evaluate_1L_lgl = F),
                                              .x)),
                     r_version_1L_chr = paste0(spine_of_results_ls$outp_smry_ls$session_ls$R.version$major,
                                               ".",
                                               spine_of_results_ls$outp_smry_ls$session_ls$R.version$minor),
                     study_descs_ls = spine_of_results_ls$study_descs_ls,
                     tables_ls = make_ss_tbls_ls(spine_of_results_ls$outp_smry_ls,
                                                 mdls_smry_tbls_ls = mdls_smry_tbls_ls,
                                                 covars_mdls_ls = covars_mdls_ls,
                                                 descv_tbls_ls = descv_tbls_ls,
                                                 nbr_of_digits_1L_int = spine_of_results_ls$nbr_of_digits_1L_int),
                     ttu_cs_ls = ttu_cs_ls,
                     ttu_lngl_ls = ttu_lngl_ls,
                     ttu_version_1L_chr = spine_of_results_ls$outp_smry_ls$session_ls$otherPkgs$TTU$Version,
                     var_nm_change_lup = spine_of_results_ls$var_nm_change_lup,
                     version_1L_chr = version_1L_chr)
  results_ls <- transform_tbls_for_covar_nms(results_ls) %>% transform_tbls_for_csnl_mdls()
  return(results_ls)
}
make_scaling_text <- function(results_ls,
                              table_1L_chr = "cfscl"){
  if(startsWith(table_1L_chr,"cfscl")){
    table_df <- results_ls$tables_ls$ind_preds_coefs_tbl
  }else{
    if(startsWith(table_1L_chr,"coefscovarstype")){
      table_df <- results_ls$tables_ls %>%
        purrr::pluck(paste0("mdl_type_",
                            table_1L_chr %>% stringr::str_remove("coefscovarstype"),
                            "_covar_mdls_tb"))
    }
  }
  predrs_chr <- table_df$Parameter %>% setdiff(c("SD (Intercept)","Intercept")) %>% stringr::str_replace_all(" model","") %>% 
    stringr::str_replace_all(" baseline","") %>% stringr::str_replace_all(" change","") %>% 
    stringr::str_replace_all(" scaled","") %>% stringr::str_replace_all(" unscaled","") %>%unique()
  predrs_lup <- results_ls$mdl_ingredients_ls$predictors_lup %>%
    dplyr::filter(short_name_chr %in% predrs_chr)
  scaling_dbl <- predrs_lup$mdl_scaling_dbl %>% unique()
  text_1L_chr <- ifelse(all(scaling_dbl==1),
                        "",
                        paste0("Note: ",
                               scaling_dbl %>%
                                 purrr::map_chr(~ {
                                   scaled_predrs_chr <- predrs_lup %>%
                                     dplyr::filter(mdl_scaling_dbl == .x) %>%
                                     dplyr::pull(short_name_chr) %>%
                                     sort() %>%
                                     transform_names(rename_lup = results_ls$var_nm_change_lup)
                                   ifelse(.x ==1,
                                          "",
                                          paste0("The ",
                                                 scaled_predrs_chr %>%
                                                   paste0(collapse = ", ") %>%
                                                   stringi::stri_replace_last(fixed = ",", " and"),
                                                 " parameter",
                                                 ifelse(length(scaled_predrs_chr) == 1,
                                                        " was",
                                                        "s were"),
                                                 " first multiplied by ",
                                                 .x,
                                                 "."))
                                 }
                                 ) %>%
                                 paste0(collapse = " ")))
  return(text_1L_chr)
}
make_shareable_mdl <- function (fake_ds_tb, mdl_smry_tb,
                                depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                                id_var_nm_1L_chr = "fkClientID", tfmn_1L_chr = "CLL", mdl_type_1L_chr = "OLS_CLL",
                                mdl_types_lup = NULL, control_1L_chr = NA_character_, start_1L_chr = NA_character_,
                                seed_1L_int = 12345L)
{
  if (is.null(mdl_types_lup))
    utils::data(mdl_types_lup, envir = environment())
  if (is.na(tfmn_1L_chr))
    tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                            match_value_xx = mdl_type_1L_chr,
                                            match_var_nm_1L_chr = "short_name_chr",
                                            target_var_nm_1L_chr = "tfmn_chr",
                                            evaluate_1L_lgl = F)
  predr_var_nms_chr <- mdl_smry_tb$Parameter[!mdl_smry_tb$Parameter %in% c("SD (Intercept)", "Intercept",
                                                                           "R2", "RMSE","Sigma")] %>%
    stringi::stri_replace_last_fixed(" baseline","_baseline") %>%
    stringi::stri_replace_last_fixed(" change","_change") %>%
    stringi::stri_replace_last_fixed(" scaled","_scaled") %>%
    stringi::stri_replace_last_fixed(" unscaled","_unscaled") 
  tfd_depnt_var_nm_1L_chr <- transform_depnt_var_nm(depnt_var_nm_1L_chr,
                                                    tfmn_1L_chr = tfmn_1L_chr)
  if (length(predr_var_nms_chr) > 1) {
    covar_var_nms_chr <- predr_var_nms_chr[2:length(predr_var_nms_chr)]
  }else{
    covar_var_nms_chr <- NA_character_
  }
  model_mdl <- make_mdl(fake_ds_tb %>%
                          dplyr::select(tidyselect::all_of(c(id_var_nm_1L_chr,
                                                             tfd_depnt_var_nm_1L_chr,
                                                             predr_var_nms_chr))),
                        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                        predr_var_nm_1L_chr = predr_var_nms_chr[1], covar_var_nms_chr = covar_var_nms_chr,
                        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr,
                        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr,
                        start_1L_chr = start_1L_chr)
  if(ready4::get_from_lup_obj(mdl_types_lup,
                              match_value_xx = mdl_type_1L_chr,
                              match_var_nm_1L_chr = "short_name_chr",
                              target_var_nm_1L_chr = "fn_chr",
                              evaluate_1L_lgl = F) == "betareg::betareg"){
    model_coeffs_dbl <- model_mdl$coefficients$mean
    
  }else{
    model_coeffs_dbl <- model_mdl$coefficients
  }
  param_nms_chr <- model_coeffs_dbl %>% names()
  mdl_smry_tb <- mdl_smry_tb %>% dplyr::mutate(Parameter = dplyr::case_when(Parameter ==
                                                                              "Intercept" ~ "(Intercept)", TRUE ~ purrr::map_chr(Parameter,
                                                                                                                                 ~stringr::str_replace_all(.x, " ", "_")))) %>% dplyr::filter(Parameter %in%
                                                                                                                                                                                                param_nms_chr) %>% dplyr::slice(match(param_nms_chr, Parameter))
  assertthat::assert_that(all(param_nms_chr == mdl_smry_tb$Parameter),
                          msg = "Parameter names mismatch between data and model summary table")
  model_coeffs_dbl <- mdl_smry_tb$Estimate
  names(model_coeffs_dbl) <- param_nms_chr
  if(ready4::get_from_lup_obj(mdl_types_lup,
                              match_value_xx = mdl_type_1L_chr,
                              match_var_nm_1L_chr = "short_name_chr",
                              target_var_nm_1L_chr = "fn_chr",
                              evaluate_1L_lgl = F) == "betareg::betareg"){
    model_mdl$coefficients$mean <- model_coeffs_dbl
    
  }else{
    model_mdl$coefficients <- model_coeffs_dbl
  }
  return(model_mdl)
}
make_smry_of_brm_mdl <- function (mdl_ls,
                                  data_tb,
                                  depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                                  predr_vars_nms_chr,
                                  mdl_nm_1L_chr = NA_character_,
                                  seed_1L_dbl = 23456,
                                  tfmn_1L_chr) {
  if (is.na(mdl_nm_1L_chr))
    mdl_nm_1L_chr <- predr_vars_nms_chr[1]
  set.seed(seed_1L_dbl)
  predictions <- stats::predict(mdl_ls, summary = F) %>%
    calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                             tfmn_is_outp_1L_lgl = T)
  sd_intcpt_df <- summary(mdl_ls, digits = 4)$random[[1]]
  sd_intcpt_df <- sd_intcpt_df[1:nrow(sd_intcpt_df), 1:4]  %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
  coef <- summary(mdl_ls, digits = 4)$fixed #
  coef <- coef[1:nrow(coef), 1:4] %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
  R2 <- brms::bayes_R2(mdl_ls) %>%
    as.vector()#
  RMSE <- psych::describe(apply(predictions, 1, calculate_rmse, y_dbl = data_tb %>%
                                  dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))), quant = c(0.25,
                                                                                             0.75), skew = F, ranges = F)
  RMSE <- cbind(RMSE$mean, RMSE$sd, RMSE$Q0.25, RMSE$Q0.75) %>%
    as.vector()
  Sigma <- summary(mdl_ls, digits = 4)$spec_par[1:4]
  smry_of_brm_mdl_tb <- data.frame(round(rbind(sd_intcpt_df,
                                               coef,
                                               R2,
                                               RMSE,
                                               Sigma), 3)) %>%
    dplyr::mutate(Parameter = c("SD (Intercept)","Intercept",
                                purrr::map(predr_vars_nms_chr,
                                           ~paste0(.x, 
                                                   c(" baseline",
                                                     " change",
                                                     " scaled",
                                                     " unscaled")
                                                   )) %>%
                                  purrr::flatten_chr() %>%
                                  intersect(purrr::map_chr(names(mdl_ls$data),
                                            ~ stringi::stri_replace_last_fixed(.x,"_baseline"," baseline") %>%
                                              stringi::stri_replace_last_fixed("_change", " change") %>%
                                              stringi::stri_replace_last_fixed("_scaled", " scaled") %>%
                                              stringi::stri_replace_last_fixed("_unscaled", " unscaled"))),
                                "R2", "RMSE", "Sigma"),
                  Model = mdl_nm_1L_chr) %>%
    dplyr::mutate(`95% CI` = paste(l.95..CI,
                                   ",",
                                   u.95..CI)) %>%
    dplyr::rename(SE = Est.Error) %>%
    dplyr::select(Model, Parameter, Estimate, SE, `95% CI`)
  rownames(smry_of_brm_mdl_tb)<-NULL
  return(smry_of_brm_mdl_tb)
}
make_smry_of_mdl_outp <- function (data_tb,
                                   #model_mdl = NULL, # TEST REMOVE
                                   folds_1L_int = 10,
                                   depnt_var_min_val_1L_dbl = numeric(0),
                                   depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                                   start_1L_chr = NULL,
                                   tfmn_1L_chr = "NTF",
                                   predr_var_nm_1L_chr,
                                   covar_var_nms_chr = NA_character_,
                                   mdl_type_1L_chr = "OLS_NTF",
                                   mdl_types_lup = NULL,
                                   predn_type_1L_chr = NULL)
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  data_tb <- data_tb %>% dplyr::filter(!is.na(!!rlang::sym(predr_var_nm_1L_chr)))
  data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
  mdl_desc_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                              match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                              target_var_nm_1L_chr = "long_name_chr", evaluate_1L_lgl = F)
  folds_ls <- make_folds_ls(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                            folds_1L_int = folds_1L_int)
  control_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                             match_var_nm_1L_chr = "short_name_chr",
                                             match_value_xx = mdl_type_1L_chr,
                                             target_var_nm_1L_chr = "control_chr",
                                             evaluate_1L_lgl = F)
  smry_of_one_predr_mdl_tb <- purrr::map_dfr(folds_ls, ~{
    model_mdl <- make_mdl(data_tb[-.x,], depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                          depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, # Unnecessary
                          start_1L_chr = start_1L_chr, tfmn_1L_chr = tfmn_1L_chr,
                          predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
                          mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr)
    predd_old_dbl <- stats::predict(model_mdl, type = predn_type_1L_chr) %>%
      calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr, ## TEST
                               tfmn_is_outp_1L_lgl = T)
    predd_new_dbl <- stats::predict(model_mdl, newdata = data_tb[.x, ],
                                    type = predn_type_1L_chr) %>% calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                                                                                           tfmn_is_outp_1L_lgl = T)
    tibble::tibble(Rsquared = caret::R2(predd_old_dbl, data_tb[-.x,
    ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
    form = "traditional"), RMSE = caret::RMSE(predd_old_dbl,
                                              data_tb[-.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
    MAE = caret::MAE(predd_old_dbl, data_tb[-.x, ] %>%
                       dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
    RsquaredP = caret::R2(predd_new_dbl, data_tb[.x, ] %>%
                            dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                          form = "traditional"), RMSEP = caret::RMSE(predd_new_dbl,
                                                                     data_tb[.x, ] %>% dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))),
    MAEP = caret::MAE(predd_new_dbl, data_tb[.x, ] %>%
                        dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr))))
  }) %>% dplyr::summarise_all(mean) %>% dplyr::mutate(Model = mdl_desc_1L_chr) %>%
    dplyr::select(Model, dplyr::everything())
  return(smry_of_one_predr_mdl_tb)
}
make_smry_of_ts_mdl_outp <- function (data_tb, # rename ts to lngl
                                      predr_vars_nms_chr, mdl_nm_1L_chr, path_to_write_to_1L_chr = NA_character_,
                                      depnt_var_min_val_1L_dbl = numeric(0),
                                      depnt_var_nm_1L_chr = "utl_total_w", # Remove Default
                                      id_var_nm_1L_chr = "fkClientID",
                                      round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", predictors_lup, utl_min_val_1L_dbl = -1,
                                      backend_1L_chr = getOption("brms.backend", "rstan"), iters_1L_int = 4000L, mdl_types_lup,
                                      seed_1L_int = 1000L, prior_ls = NULL, control_ls = NULL)
{
  scaling_fctr_dbl <- predr_vars_nms_chr %>% purrr::map_dbl(~
                                                              ifelse(.x %in% predictors_lup$short_name_chr,
                                                                     ready4::get_from_lup_obj(predictors_lup,
                                                                                              target_var_nm_1L_chr = "mdl_scaling_dbl",
                                                                                              match_value_xx = .x,
                                                                                              match_var_nm_1L_chr = "short_name_chr",
                                                                                              evaluate_1L_lgl = F),
                                                                     1))
  mdl_type_1L_chr <- mdl_nm_1L_chr %>%
    stringr::str_remove(paste0(predr_vars_nms_chr[1], "_", ifelse(is.na(predr_vars_nms_chr[2]), "", paste0(predr_vars_nms_chr[2],
                                                                                                           "_"))))
  mdl_type_1L_chr <- stringr::str_sub(mdl_type_1L_chr, start = 1 + (mdl_type_1L_chr %>% stringi::stri_locate_first_fixed("_"))[1,2] %>% as.vector())
  tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                          target_var_nm_1L_chr = "tfmn_chr",
                                          match_value_xx = mdl_type_1L_chr,
                                          match_var_nm_1L_chr = "short_name_chr",
                                          evaluate_1L_lgl = F)
  tfd_data_tb <- transform_tb_to_mdl_inp(data_tb, 
                                         depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                         depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, 
                                         predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
                                         round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr,
                                         scaling_fctr_dbl = scaling_fctr_dbl,
                                         tfmn_1L_chr = tfmn_1L_chr) ## PICK UP FROM HERE
  tfd_depnt_var_nm_1L_chr <- transform_depnt_var_nm(depnt_var_nm_1L_chr,
                                                    tfmn_1L_chr = tfmn_1L_chr)
  family_fn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                               match_var_nm_1L_chr = "short_name_chr",
                                               match_value_xx = mdl_type_1L_chr,
                                               target_var_nm_1L_chr = "family_chr",
                                               evaluate_1L_lgl = F)
  family_fn_1L_chr <- ifelse(is.na(family_fn_1L_chr),
                             ifelse(startsWith(mdl_type_1L_chr,"BET"),
                                    paste0("brms::Beta(link = \"",
                                           get_link_from_tfmn(stringr::str_sub(mdl_type_1L_chr,
                                                                               start = -3)),
                                           "\")"),
                                    "gaussian(identity)"),
                             family_fn_1L_chr)
  args_ls <- list(data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
                  predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
                  is_csnl_1L_lgl = !(!identical(round_var_nm_1L_chr, character(0)) && ifelse(identical(round_var_nm_1L_chr, character(0)),T,!is.na(round_var_nm_1L_chr))),
                  iters_1L_int = iters_1L_int,
                  backend_1L_chr = backend_1L_chr,
                  family_fn_1L_chr = family_fn_1L_chr,
                  seed_1L_int = seed_1L_int,
                  prior_ls = prior_ls, control_ls = control_ls)
  # if(startsWith(mdl_type_1L_chr, "GLM_BNL")){
  # WRITE FN
  # }else{
  mdl_ls <- rlang::exec(fit_ts_model_with_brm, !!!args_ls)
  # }
  smry_of_ts_mdl_ls <- list(smry_of_ts_mdl_tb = make_smry_of_brm_mdl(mdl_ls,
                                                                     data_tb = tfd_data_tb, depnt_var_nm_1L_chr = tfd_depnt_var_nm_1L_chr,
                                                                     predr_vars_nms_chr = predr_vars_nms_chr,
                                                                     mdl_nm_1L_chr = mdl_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr))
  if (!is.na(path_to_write_to_1L_chr)) {
    smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr <- paste0(path_to_write_to_1L_chr,
                                                      "/", mdl_nm_1L_chr, ".RDS")
    if (file.exists(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr))
      file.remove(smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
    saveRDS(mdl_ls, smry_of_ts_mdl_ls$path_to_mdl_ls_1L_chr)
    smry_of_ts_mdl_ls$paths_to_mdl_plts_chr <- write_ts_mdl_plts(mdl_ls,
                                                                 tfd_data_tb = tfd_data_tb,
                                                                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                 mdl_nm_1L_chr = mdl_nm_1L_chr,
                                                                 path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                                 round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                                 tfmn_1L_chr = tfmn_1L_chr,
                                                                 utl_min_val_1L_dbl = utl_min_val_1L_dbl)
  }
  return(smry_of_ts_mdl_ls)
}
make_unique_ls_elmt_idx_int <- function (data_ls)
{
  combos_tb <- tibble::tibble(names_chr = data_ls %>% purrr::map_chr(~paste0(.x[1],ifelse(length(.x)>1,.x[2],""))),
                              indices_int = NA_integer_)
  combos_tb <- combos_tb %>% dplyr::group_by(names_chr) %>%
    dplyr::mutate(indices_int = dplyr::row_number())
  unique_ls_elmt_idx_int <- combos_tb$indices_int
  return(unique_ls_elmt_idx_int)
}
make_valid_params_ls_ls <- function(analysis_core_params_ls, # Generalise utility component
                                    ds_tb,
                                    path_params_ls,
                                    maui_params_ls,
                                    candidate_covar_nms_chr = NA_character_,
                                    prefd_covars_chr = NULL,
                                    prefd_mdl_types_chr = NULL,
                                    raw_ds_tfmn_fn = NULL,
                                    scndry_analysis_extra_vars_chr = NA_character_,
                                    subtitle_1L_chr = "Methods Report 1: Analysis Program (Primary Analysis)",
                                    utl_class_fn_1L_chr = "as.numeric"){
  valid_params_ls_ls <- make_prmry_analysis_params_ls(analysis_core_params_ls = analysis_core_params_ls,
                                                      candidate_covar_nms_chr = candidate_covar_nms_chr,
                                                      ds_tb = ds_tb,
                                                      path_params_ls = path_params_ls,
                                                      maui_params_ls = maui_params_ls,
                                                      prefd_covars_chr = prefd_covars_chr,
                                                      prefd_mdl_types_chr = prefd_mdl_types_chr,
                                                      raw_ds_tfmn_fn = raw_ds_tfmn_fn,
                                                      subtitle_1L_chr = subtitle_1L_chr,
                                                      utl_class_fn_1L_chr = utl_class_fn_1L_chr) %>%
    transform_params_ls_to_valid(scndry_analysis_extra_vars_chr = scndry_analysis_extra_vars_chr)
  valid_params_ls_ls$params_ls$short_and_long_nm <- NULL
  valid_params_ls_ls$short_and_long_nm <- maui_params_ls$short_and_long_nm
  valid_params_ls_ls$path_params_ls <- path_params_ls
  return(valid_params_ls_ls)
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
manufacture_SpecificSynopsis <- function(x,
                                         depnt_var_min_val_1L_dbl = numeric(0),
                                         depnt_var_nms_chr = NA_character_,
                                         make_cmpst_plt_1L_lgl = F,
                                         scndry_anlys_params_ls = NULL,
                                         version_1L_chr = "",
                                         what_1L_chr = "input_params_ls"){
  if(what_1L_chr %in% c("abstract_args_ls","ds_descvs_ls","ds_smry_ls","input_params_ls","results_ls","mdl_smry_ls")){
    y_SpecificMixed <- SpecificMixed(a_YouthvarsProfile = x@d_YouthvarsProfile,
                                     b_SpecificParameters = x@c_SpecificParameters,
                                     c_SpecificResults = x@b_SpecificResults,
                                     paths_chr = x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls$path_to_write_to_1L_chr)
    
    if(what_1L_chr %in% c("ds_descvs_ls","ds_smry_ls","mdl_smry_ls")){ # Could add input_params_ls to this logic once corresponding SpecificProject methd is updated
      # Would then need to pass scndry_anlys_params_ls to mthd
      object_xx <- manufacture_SpecificProject(y_SpecificMixed, # manufacture when exported
                                               depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                               what_1L_chr = what_1L_chr)
    }
    if(what_1L_chr %in% c("abstract_args_ls","input_params_ls","results_ls")){
      header_yaml_args_ls <- ready4show::make_header_yaml_args_ls(authors_tb = x@authors_r3,
                                                                  institutes_tb = x@institutes_r3,
                                                                  title_1L_chr = x@title_1L_chr,
                                                                  keywords_chr = x@keywords_chr)
      maui_params_ls <- make_maui_params_ls(maui_domains_pfxs_1L_chr = y_SpecificMixed@b_SpecificParameters@itm_prefix_1L_chr,
                                            maui_itm_short_nms_chr = y_SpecificMixed@b_SpecificParameters@itm_labels_chr,
                                            maui_scoring_fn = NULL)
      output_format_ls <- ready4show::make_output_format_ls(manuscript_outp_1L_chr = x@outp_formats_chr[1],
                                                            manuscript_digits_1L_int = x@digits_int[1],
                                                            supplementary_outp_1L_chr = ifelse(length(x@outp_formats_chr)>1,x@outp_formats_chr[2],x@outp_formats_chr[1]),
                                                            supplementary_digits_1L_int = ifelse(length(x@digits_int)>1,x@digits_int[2],x@digits_int[1]))
      # scndry_anlys_params_ls <- make_scndry_anlys_params(candidate_predrs_chr = c("SOFAS"),
      #                                                    prefd_covars_chr = NA_character_)
      object_xx <- make_input_params(y_SpecificMixed@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                     control_ls = y_SpecificMixed@b_SpecificParameters@control_ls,
                                     ds_descvs_ls = manufacture_SpecificProject(y_SpecificMixed, #manufacture when expporting
                                                                what_1L_chr = "ds_descvs_ls"),
                                     dv_ds_nm_and_url_chr = c(x@e_Ready4useRepos@dv_nm_1L_chr,
                                                              x@e_Ready4useRepos@dv_ds_nm_1L_chr),
                                     header_yaml_args_ls = header_yaml_args_ls,
                                     maui_params_ls = maui_params_ls,
                                     output_format_ls = output_format_ls,
                                     predictors_lup = y_SpecificMixed@b_SpecificParameters@predictors_lup,
                                     prefd_covars_chr = ifelse(is.null(procure(y_SpecificMixed,
                                                                               what = "prefd_covars")),
                                                               NA_character_,
                                                               procure(y_SpecificMixed,
                                                                       what = "prefd_covars")),
                                     prefd_mdl_types_chr = procure(y_SpecificMixed,
                                                                   what = "prefd_mdls"),
                                     scndry_anlys_params_ls = scndry_anlys_params_ls,
                                     write_new_dir_1L_lgl = F)
      if(is.na(depnt_var_nms_chr[1]))
        depnt_var_nms_chr <- c(y_SpecificMixed@b_SpecificParameters@depnt_var_nm_1L_chr,
                               y_SpecificMixed@a_YouthvarsProfile@a_Ready4useDyad@dictionary_r3 %>%
                                 ready4::get_from_lup_obj(match_value_xx = y_SpecificMixed@b_SpecificParameters@depnt_var_nm_1L_chr,
                                                          match_var_nm_1L_chr = "var_nm_chr",
                                                          target_var_nm_1L_chr = "var_desc_chr"))
      object_xx$short_and_long_nm <- depnt_var_nms_chr
      object_xx <- object_xx %>%
        make_study_descs_ls(time_btwn_bl_and_fup_1L_chr = x@interval_chr,
                            background_1L_chr = x@background_1L_chr,
                            coi_1L_chr = x@coi_1L_chr,
                            conclusion_1L_chr = x@conclusion_1L_chr,
                            ethics_1L_chr = x@ethics_1L_chr,
                            funding_1L_chr = x@funding_1L_chr,
                            #predr_ctgs_ls = make_predr_ctgs_ls(x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls),
                            sample_desc_1L_chr = x@sample_desc_1L_chr,
                            var_nm_change_lup = x@correspondences_r3)
      if(what_1L_chr %in% c("abstract_args_ls","results_ls")){
        object_xx$study_descs_ls$predr_ctgs_ls <- make_predr_ctgs_ls(x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls)
        object_xx <- make_results_ls(dv_ds_nm_and_url_chr = object_xx$path_params_ls$dv_ds_nm_and_url_chr,
                                     depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                     make_cmpst_plt_1L_lgl = make_cmpst_plt_1L_lgl,
                                     outp_smry_ls = x@b_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                                     output_format_ls = object_xx$output_format_ls,
                                     params_ls_ls = object_xx,
                                     path_params_ls = list(paths_ls = list(output_data_dir_1L_chr = paste0(x@a_Ready4showPaths@outp_data_dir_1L_chr,"/Output"))),
                                     study_descs_ls = object_xx$study_descs_ls,
                                     var_nm_change_lup = object_xx$study_descs_ls$var_nm_change_lup,
                                     version_1L_chr = version_1L_chr)
        # object_xx$abstract_args_ls <- x@abstract_args_ls
        # object_xx$abstract_args_ls$abstract_ls$Background <- x@background_1L_chr
        object_xx$abstract_args_ls$abstract_ls$Conclusions <- x@conclusion_1L_chr
        object_xx$study_descs_ls$background_1L_chr <- x@background_1L_chr
        object_xx$study_descs_ls$conclusion_1L_chr <- x@conclusion_1L_chr
      }
      if(what_1L_chr == "abstract_args_ls"){
        object_xx <- make_abstract_args_ls(object_xx)
      }
    }
  }else{
    object_xx <- methods::callNextMethod()
  }
  return(object_xx)
}
plot_obsd_predd_sctr_cmprsn <- function (tfd_data_tb, depnt_var_nm_1L_chr = "utl_total_w",# Remove defaults
                                         depnt_var_desc_1L_chr = "Total weighted utility score", # Remove defaults
                                         round_var_nm_1L_chr = "round",
                                         args_ls = NULL, base_size_1L_dbl = 11, correspondences_lup = NULL, predd_val_var_nm_1L_chr = "Predicted", x_lbl_1L_chr = NA_character_, y_lbl_1L_chr = NA_character_)
{
  if(is.na(x_lbl_1L_chr))
    x_lbl_1L_chr <- paste0("Observed ", depnt_var_desc_1L_chr)
  
  if(is.na(y_lbl_1L_chr))
    y_lbl_1L_chr <- paste0(predd_val_var_nm_1L_chr, " ", depnt_var_desc_1L_chr)
  if(!is.null(correspondences_lup) && !identical(round_var_nm_1L_chr, character(0)) && ifelse(identical(round_var_nm_1L_chr, character(0)),T,!is.na(round_var_nm_1L_chr)))
    tfd_data_tb <- tfd_data_tb %>%
      dplyr::mutate(!!rlang::sym(round_var_nm_1L_chr) := !!rlang::sym(round_var_nm_1L_chr) %>%
                      purrr::map_chr(~ready4::get_from_lup_obj(correspondences_lup,
                                                               match_value_xx = .x,
                                                               match_var_nm_1L_chr = "old_nms_chr",
                                                               target_var_nm_1L_chr = "new_nms_chr")))
  if(!identical(round_var_nm_1L_chr, character(0)) && ifelse(identical(round_var_nm_1L_chr, character(0)),T,!is.na(round_var_nm_1L_chr))){
    mapping_aes <- ggplot2::aes(x = !!rlang::sym(depnt_var_nm_1L_chr), y = !!rlang::sym(predd_val_var_nm_1L_chr),
                                col = !!rlang::sym(round_var_nm_1L_chr))
    }else{
      mapping_aes <- ggplot2::aes(x = !!rlang::sym(depnt_var_nm_1L_chr), y = !!rlang::sym(predd_val_var_nm_1L_chr))
      }
  ggplot2::ggplot(tfd_data_tb) + 
    rlang::exec(ggplot2::geom_point,
                mapping_aes,
                size = 1,
                !!!args_ls) + 
    ggplot2::theme_bw(base_size = base_size_1L_dbl) +
    ggplot2::xlim(0,1) + 
    ggplot2::ylim(0, 1) + 
    ggplot2::scale_color_manual(values = c("#D55E00", "#56B4E9")) + 
    ggplot2::labs(x = x_lbl_1L_chr,
                  y = y_lbl_1L_chr,
                  col = "") + 
    ggplot2::theme(legend.position = "bottom")
}
print_covar_ttu_tbls <- function(params_ls, # Rename and generalise from TTU
                                 caption_1L_chr,
                                 table_1L_chr,
                                 ref_1L_int = 1){
  results_ls <- params_ls$results_ls
  df <- results_ls$tables_ls %>%
    purrr::pluck(paste0("mdl_type_",ref_1L_int,"_covar_mdls_tb"))
  if(!is.null(df)){
    df <- df %>%
      transform_nms_in_mdl_tbl(col_nm_1L_chr = "Parameter",
                               var_nm_change_lup = results_ls$var_nm_change_lup) %>%
      dplyr::mutate(Parameter = Parameter %>% purrr::map_chr(~stringr::str_replace_all(.x,"_"," ")))
    df %>%
      print_lngl_ttu_tbls(caption_1L_chr = caption_1L_chr,
                          params_ls = params_ls,
                          ref_1L_int = ref_1L_int,
                          table_1L_chr = table_1L_chr)
  }
}
print_indpnt_predrs_coefs_tbl <- function(params_ls,
                                          caption_1L_chr,
                                          mkdn_tbl_ref_1L_chr){
  results_ls <- params_ls$results_ls
  tb <- results_ls$tables_ls$ind_preds_coefs_tbl %>%
    transform_nms_in_mdl_tbl(col_nm_1L_chr = "Parameter",
                             var_nm_change_lup = results_ls$var_nm_change_lup)
  if(params_ls$output_type_1L_chr =="Word"){
    tb$Parameter <- stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(tb$Parameter, '\\\\textbf', ''), '\\{', ''), '\\}', '')
  }
  if(params_ls$output_type_1L_chr == "PDF"){
    tb %>%
      dplyr::mutate(dplyr::across(.cols = everything(),
                                  ~ dplyr::case_when(is.na(.x) ~ "",
                                                     T ~ .x))) %>%
      kableExtra::kbl(booktabs = T,
                      caption = knitr::opts_current$get("tab.cap"),
                      escape = F,
                      longtable = T,
                      col.names = c("Parameter", "Estimate",	"SE",	"CI (95\\%)", "R2", "Sigma", "Estimate", "SE", "CI (95\\%)", "R2", "Sigma")) %>%
      kableExtra::kable_styling(latex_options = c("repeat_header")) %>%
      kableExtra::column_spec(3:6, width = "3em") %>%
      kableExtra::row_spec(which(!is.na(tb[,5])), bold = T) %>%
      #kableExtra::column_spec(1, bold = T, width = "14em") %>%
      kableExtra::add_header_above(parse(text=paste0("c(",
                                                     "\" \"",
                                                     ", ",
                                                     paste0("\"",
                                                            results_ls$ttu_lngl_ls$best_mdls_tb[[1,"model_type"]],
                                                            " - ",
                                                            results_ls$ttu_lngl_ls$best_mdls_tb[[1,"link_and_tfmn_chr"]],
                                                            "\" = 5"),
                                                     ", ",
                                                     paste0("\"",
                                                            results_ls$ttu_lngl_ls$best_mdls_tb[[2,"model_type"]],
                                                            " - ",
                                                            results_ls$ttu_lngl_ls$best_mdls_tb[[2,"link_and_tfmn_chr"]],
                                                            "\" = 5"),
                                                     ")"
      )) %>% eval()
      ) %>%
      kableExtra::collapse_rows(columns = 1) %>%
      kableExtra::footnote(general = make_scaling_text(results_ls),
                           general_title = " ")
  }else{
    add_to_row_ls <- NULL
    tb %>%
      ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                              caption_1L_chr = caption_1L_chr,
                              mkdn_tbl_ref_1L_chr = mkdn_tbl_ref_1L_chr,
                              use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                              add_to_row_ls = add_to_row_ls,
                              footnotes_chr = make_scaling_text(results_ls),
                              sanitize_fn = force)
  }
}
print_lngl_ttu_tbls <- function(table_df, # Rename and generalise from TTU
                                params_ls,
                                caption_1L_chr,
                                table_1L_chr,
                                column_1_width_1L_chr = "25em",
                                ref_1L_int = 1){
  results_ls <- params_ls$results_ls
  if(params_ls$output_type_1L_chr == "PDF"){
    add_to_row_ls <- list()
    add_to_row_ls$pos <- list(0, nrow(table_df))
    add_to_row_ls$command <- c("Parameter & Estimate	& SE	& 95CI & R2	& Sigma\\\\\n",
                               paste0("\\hline\n",
                                      "{\\footnotesize ",
                                      make_scaling_text(results_ls,
                                                        table_1L_chr = knitr::opts_current$get("tab.id")),
                                      "}\n"))
    table_df$Parameter <- table_df$Parameter %>% purrr::map_chr(~ifelse(endsWith(.x,
                                                                                 " model"),
                                                                        paste0("\\textbf{",.x,"}"),.x))
  }else{
    add_to_row_ls <- NULL
  }
  if(params_ls$output_type_1L_chr == "PDF"){
    table_df %>%
      dplyr::mutate(dplyr::across(.cols = everything(),
                                  ~ dplyr::case_when(is.na(.x) ~ "",
                                                     T ~ .x))) %>%
      kableExtra::kbl(booktabs = T,
                      caption = knitr::opts_current$get("tab.cap"),
                      escape = F,
                      longtable = T,
                      col.names = c("Parameter", "Estimate",	"SE",	"CI (95\\%)", "R2", "Sigma")) %>%
      kableExtra::kable_styling(latex_options = c("repeat_header"),
                                full_width = F) %>%
      kableExtra::column_spec(1, width = column_1_width_1L_chr) %>%
      kableExtra::row_spec(which(!is.na(table_df[,5])), bold = T) %>%
      kableExtra::collapse_rows(columns = 1) %>%
      kableExtra::footnote(general = make_scaling_text(results_ls,
                                                       table_1L_chr = table_1L_chr),
                           general_title = " ")
  }else{
    table_df %>%
      ready4show::print_table(output_type_1L_chr = params_ls$output_type_1L_chr,
                              caption_1L_chr = caption_1L_chr,
                              mkdn_tbl_ref_1L_chr = paste0("tab:",table_1L_chr),
                              use_rdocx_1L_lgl = ifelse(params_ls$output_type_1L_chr=="Word",T,F),
                              add_to_row_ls = add_to_row_ls,
                              footnotes_chr = make_scaling_text(results_ls,
                                                                table_1L_chr = table_1L_chr),
                              hline_after_ls = c(-1,0),
                              sanitize_fn = force)
  }
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
renew_SpecificProject <- function(x,
                                  new_val_xx,
                                  type_1L_chr = "results",
                                  what_1L_chr = "prefd_mdls"
){
  if(type_1L_chr == "results"){
    if(what_1L_chr == "prefd_mdls") # Make version for SpecificFixed and above
      x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$mdl_smry_ls$prefd_mdl_types_chr <- new_val_xx
  }
  if(what_1L_chr == "prefd_covars"){
    x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$mdl_smry_ls$prefd_covars_chr <- new_val_xx
  }
  if(what_1L_chr %in% c("base","dummys","levels")){
    x@a_YouthvarsProfile@a_Ready4useDyad <- renew_Ready4useDyad(x@a_YouthvarsProfile@a_Ready4useDyad,
                                                                new_val_xx = new_val_xx,
                                                                type_1L_chr = what_1L_chr)
  }

  return(x)
}
transform_ds_for_mdlng <- function (data_tb,
                                    depnt_var_nm_1L_chr = "utl_total_w", # Undo defaults
                                    depnt_var_min_val_1L_dbl = numeric(0),
                                    predr_var_nm_1L_chr,
                                    covar_var_nms_chr = NA_character_)
{
  mdl_vars_chr <- c(names(data_tb)[names(data_tb) %>% startsWith(depnt_var_nm_1L_chr)],
                    predr_var_nm_1L_chr, covar_var_nms_chr) %>% purrr::discard(is.na)
  tfd_data_tb <- data_tb %>% tidyr::drop_na(!!!rlang::syms(mdl_vars_chr)) %>%
    dplyr::select(!!!rlang::syms(mdl_vars_chr))
  if(!identical(depnt_var_min_val_1L_dbl, numeric(0)))
    tfd_data_tb <- tfd_data_tb %>% dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr) := !!rlang::sym(depnt_var_nm_1L_chr) %>%
                                                   purrr::map_dbl(~max(.x,depnt_var_min_val_1L_dbl)))
  return(tfd_data_tb)
}
transform_ds_to_predn_ds <- function(data_tb,
                                     predr_vars_nms_chr,
                                     tfmn_1L_chr,
                                     depnt_var_nm_1L_chr,
                                     id_var_nm_1L_chr,
                                     round_var_nm_1L_chr,
                                     round_bl_val_1L_chr,
                                     predictors_lup){
  data_tb <- data_tb %>%
    dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr):= NA_real_)
  data_tb <- purrr::reduce(predr_vars_nms_chr,
                           .init = data_tb,
                           ~ {
                             predr_cls_fn <- eval(parse(text=ready4::get_from_lup_obj(predictors_lup,
                                                                                      match_var_nm_1L_chr = "short_name_chr",
                                                                                      match_value_xx = .y,
                                                                                      target_var_nm_1L_chr = "class_fn_chr",
                                                                                      evaluate_1L_lgl = F)))
                             dplyr::mutate(.x,
                                           !!rlang::sym(.y) := !!rlang::sym(.y) %>% rlang::exec(.fn = predr_cls_fn))
                           })
  data_tb <- data_tb %>% transform_tb_to_mdl_inp(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                 predr_vars_nms_chr = predr_vars_nms_chr,
                                                 id_var_nm_1L_chr = id_var_nm_1L_chr,
                                                 round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                 round_bl_val_1L_chr = round_bl_val_1L_chr,
                                                 drop_all_msng_1L_lgl = F,
                                                 scaling_fctr_dbl = purrr::map_dbl(predr_vars_nms_chr,
                                                                                   ~ ready4::get_from_lup_obj(predictors_lup,
                                                                                                              target_var_nm_1L_chr = "mdl_scaling_dbl",
                                                                                                              match_var_nm_1L_chr = "short_name_chr",
                                                                                                              match_value_xx = .x,
                                                                                                              evaluate_1L_lgl = F)),
                                                 ungroup_1L_lgl = T,
                                                 tfmn_1L_chr = tfmn_1L_chr)
  return(data_tb)
}
transform_nms_in_mdl_tbl <- function(mdl_tbl_tb,
                                     col_nm_1L_chr = "Parameter",
                                     var_nm_change_lup = NULL){
  if(is.null(var_nm_change_lup)){
    tfd_mdl_tbl_tb <-  mdl_tbl_tb
  }else{
    tfd_mdl_tbl_tb <-  mdl_tbl_tb %>%
      dplyr::mutate(!!rlang::sym(col_nm_1L_chr) := dplyr::case_when(!!rlang::sym(col_nm_1L_chr) %>%
                                                                      purrr::map_lgl(~(endsWith(.x," model") | endsWith(.x," baseline") | endsWith(.x," change") | endsWith(.x," scaled") | endsWith(.x," unscaled"))) ~ !!rlang::sym(col_nm_1L_chr) %>% purrr::map_chr(~{
                                                                        sfx_starts_1L_int <- stringi::stri_locate_first_fixed(.x," ")[[1,1]]
                                                                        paste0(stringr::str_sub(.x,end=(sfx_starts_1L_int-1)) %>%
                                                                                 strsplit("_") %>%
                                                                                 purrr::pluck(1) %>%
                                                                                 transform_names(rename_lup = var_nm_change_lup) %>%
                                                                                 paste0(collapse = "_"),
                                                                               stringr::str_sub(.x,start=sfx_starts_1L_int))}),
                                                                    T ~ !!rlang::sym(col_nm_1L_chr)))
  }
  return(tfd_mdl_tbl_tb)
}
transform_params_ls_from_lup <- function(params_ls,
                                         rename_lup){
  if(!is.null(params_ls$ds_descvs_ls)){
    params_ls$ds_descvs_ls$candidate_predrs_chr <- params_ls$ds_descvs_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
    params_ls$ds_descvs_ls$cohort_descv_var_nms_chr <- params_ls$ds_descvs_ls$cohort_descv_var_nms_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
    params_ls$ds_descvs_ls$utl_wtd_var_nm_1L_chr <- params_ls$ds_descvs_ls$utl_wtd_var_nm_1L_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  if(!is.null(params_ls$predictors_lup)){
    params_ls$predictors_lup$short_name_chr <-  params_ls$predictors_lup$short_name_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  params_ls$candidate_covar_nms_chr <- params_ls$candidate_covar_nms_chr %>%
    purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                           .x,
                           ready4::get_from_lup_obj(rename_lup,
                                                    match_value_xx = .x,
                                                    match_var_nm_1L_chr = "old_nms_chr",
                                                    target_var_nm_1L_chr = "new_nms_chr",
                                                    evaluate_1L_lgl = F)))
  if(!is.na(params_ls$prefd_covars_chr[1])){
    params_ls$prefd_covars_chr <- params_ls$prefd_covars_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  if(!is.null(params_ls$candidate_predrs_chr)){
    params_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  return(params_ls)
}
transform_params_ls_to_valid <- function(params_ls,
                                         scndry_analysis_extra_vars_chr = NA_character_){
  target_var_nms_chr <- c(params_ls$ds_descvs_ls$utl_wtd_var_nm_1L_chr,
                          params_ls$ds_descvs_ls$candidate_predrs_chr,
                          params_ls$candidate_covar_nms_chr,
                          scndry_analysis_extra_vars_chr) %>%
    purrr::discard(is.na) %>%
    unique()
  valid_var_nms_chr <- target_var_nms_chr %>%
    stringi::stri_replace_last_fixed("_dbl","") %>%
    stringi::stri_replace_last_fixed("_int","") %>%
    stringi::stri_replace_all_fixed("_","")
  unchanged_var_nms_chr <- setdiff(params_ls$ds_descvs_ls$dictionary_tb$var_nm_chr,
                                   target_var_nms_chr)
  rename_lup <- tibble::tibble(old_nms_chr = c(unchanged_var_nms_chr,target_var_nms_chr),
                               new_nms_chr = make.unique(c(unchanged_var_nms_chr,
                                                           valid_var_nms_chr), sep="V")) %>%
    dplyr::filter(!old_nms_chr %in% unchanged_var_nms_chr)
  params_ls$ds_tb <- youthvars::transform_ds_with_rename_lup(params_ls$ds_tb,
                                                             rename_lup = rename_lup,
                                                             target_var_nms_chr = target_var_nms_chr)
  params_ls$ds_descvs_ls$dictionary_tb <- params_ls$ds_descvs_ls$dictionary_tb %>%
    transform_dict_with_rename_lup(rename_lup = rename_lup)
  rename_lup <- rename_lup %>%
    dplyr::filter(old_nms_chr != new_nms_chr)
  valid_params_ls_ls <- list(params_ls = params_ls %>%
                               transform_params_ls_from_lup(rename_lup = rename_lup),
                             rename_lup = rename_lup)
  return(valid_params_ls_ls)
}
transform_tb_to_mdl_inp <- function (data_tb,
                                     depnt_var_min_val_1L_dbl = numeric(0),
                                     depnt_var_max_val_1L_dbl = 0.99999,
                                     depnt_var_nm_1L_chr = "utl_total_w", # remove default
                                     predr_vars_nms_chr,
                                     id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
                                     round_bl_val_1L_chr = "Baseline", drop_all_msng_1L_lgl = T, scaling_fctr_dbl = 1,
                                     tfmn_1L_chr = "NTF", ungroup_1L_lgl = F)
{
  if(length(scaling_fctr_dbl)!=length(predr_vars_nms_chr)){
    scaling_fctr_dbl <- rep(scaling_fctr_dbl[1],length(predr_vars_nms_chr))
  }
  data_tb <- data.frame(data_tb) %>%
    ready4use::remove_labels_from_ds()
  tfd_for_mdl_inp_tb <- data_tb %>% dplyr::select(dplyr::all_of(id_var_nm_1L_chr), dplyr::all_of(round_var_nm_1L_chr),
                                                  dplyr::all_of(predr_vars_nms_chr),
                                                  dplyr::all_of(depnt_var_nm_1L_chr) # Moved from first var in tb
  ) %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) 
  
  tfd_for_mdl_inp_tb <- if(!identical(round_var_nm_1L_chr, character(0)) && ifelse(identical(round_var_nm_1L_chr, character(0)),T,!is.na(round_var_nm_1L_chr))){
    tfd_for_mdl_inp_tb %>% dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr), 
                                          !!rlang::sym(round_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr),
                  .init = tfd_for_mdl_inp_tb,
                  ~ {
                    idx_1L_int <- as.integer(.y)
                    .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]),
                                                       .fns = list(baseline = ~dplyr::first(.)*scaling_fctr_dbl[idx_1L_int],
                                                                   change = ~ifelse(!!rlang::sym(round_var_nm_1L_chr) == round_bl_val_1L_chr,
                                                                                    0,
                                                                                    (. - dplyr::lag(.))*scaling_fctr_dbl[idx_1L_int]))))
                  })
  }else{
    tfd_for_mdl_inp_tb %>% dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr))
    tfd_for_mdl_inp_tb <- purrr::reduce(1:length(predr_vars_nms_chr),
                                        .init = tfd_for_mdl_inp_tb,
                                        ~ {
                                          idx_1L_int <- as.integer(.y)
                                          table_tb <- .x %>% dplyr::mutate(dplyr::across(dplyr::all_of(predr_vars_nms_chr[idx_1L_int]),
                                                                             .fns = list(baseline = ~dplyr::first(.)*scaling_fctr_dbl[idx_1L_int],
                                                                                         change = ~ 0)))
                                          old_name_1L_chr <- paste0(predr_vars_nms_chr[idx_1L_int],"_baseline")
                                          new_name_1L_chr <- paste0(predr_vars_nms_chr[idx_1L_int],ifelse(scaling_fctr_dbl[idx_1L_int]==1,"_unscaled","_scaled"))
                                          table_tb <- table_tb %>% dplyr::rename(!!rlang::sym(new_name_1L_chr) := !!rlang::sym(old_name_1L_chr))
                                        })
  }

  if(!identical(depnt_var_min_val_1L_dbl, numeric(0))){
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>% dplyr::mutate(!!rlang::sym(depnt_var_nm_1L_chr) := !!rlang::sym(depnt_var_nm_1L_chr) %>%
                                                                 purrr::map_dbl(~max(.x,depnt_var_min_val_1L_dbl)))
  }
  tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
    add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                      tfmn_1L_chr = tfmn_1L_chr,
                      depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl)
  if(drop_all_msng_1L_lgl){
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
      stats::na.omit()
  }
  if(ungroup_1L_lgl){
    tfd_for_mdl_inp_tb <- tfd_for_mdl_inp_tb %>%
      dplyr::ungroup()
  }
  tfd_for_mdl_inp_tb  <- tfd_for_mdl_inp_tb %>%
    transform_uid_var(id_var_nm_1L_chr = id_var_nm_1L_chr)
  

  return(tfd_for_mdl_inp_tb)
}
transform_tbls_for_csnl_mdls <- function(results_ls){
  if(is.na(results_ls$cohort_ls$n_fup_1L_dbl)){
    results_ls$tables_ls <- results_ls$tables_ls %>% 
      purrr::map(~{
        column_nm_1L_chr <- names(.x)[1]
        .x %>%
          dplyr::mutate(!!rlang::sym(column_nm_1L_chr) := !!rlang::sym(column_nm_1L_chr) %>%
                          purrr::map_chr(~ifelse(endsWith(.x, " baseline"), stringi::stri_replace_last_fixed(.x, " baseline","") ,.x)))
      })
  }
  return(results_ls)
}
transform_tbls_for_covar_nms <- function(results_ls){
  results_ls$tables_ls <- results_ls$tables_ls %>% 
    purrr::map(~{
      column_nm_1L_chr <- names(.x)[1]
      predr_vars_nms_chr <- get_predrs_by_ctg(results_ls,collapse_1L_lgl = T) %>% purrr::flatten_chr()
      .x %>%
        dplyr::mutate(!!rlang::sym(column_nm_1L_chr) := !!rlang::sym(column_nm_1L_chr) %>%
                        purrr::map_chr(~
                                         {
                                           var_nm_1L_chr <- .x
                                           purrr::reduce(c(" baseline"," change"),
                                                         .init = var_nm_1L_chr,
                                                         ~ ifelse(endsWith(.x, .y) &&  !(stringi::stri_replace_last_fixed(.x, .y,"") %in% predr_vars_nms_chr),
                                                                  ready4::get_from_lup_obj(results_ls$mdl_ingredients_ls$dictionary_tb,
                                                                                           match_value_xx = stringi::stri_replace_last_fixed(.x, .y,""),
                                                                                           match_var_nm_1L_chr = "var_nm_chr",
                                                                                           target_var_nm_1L_chr = "var_desc_chr") %>%
                                                                    Hmisc::capitalize(),
                                                                  .x))
                                         }
                        ))
    })
  return(results_ls)
}
# transform_ts_mdl_data <- function (mdl_ls, data_tb,
#                                    depnt_var_nm_1L_chr = "utl_total_w", # remove default
#                                    predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", mdl_nm_1L_chr)
# {
#   old_data_tb <- data_tb %>% dplyr::select(c(dplyr::all_of(id_var_nm_1L_chr),
#                                              dplyr::all_of(depnt_var_nm_1L_chr), predr_vars_nms_chr %>%
#                                                purrr::map(~paste0(.x, c("", "_baseline", "_change"))) %>% # Condidtionally replace "_baseline" with "_scaled" or "_unscaled"
#                                                purrr::flatten_chr()))
#   cnfdl_mdl_ls <- mdl_ls
#   cnfdl_mdl_ls$data <- old_data_tb %>% as.data.frame() %>%
#     dplyr::summarise(dplyr::across(dplyr::everything(), ~sample(.x,
#                                                                 1)))
#   return(cnfdl_mdl_ls)
# }
write_analyses <- function(input_params_ls,
                           abstract_args_ls = NULL,
                           combinations_1L_lgl = F, 
                           depnt_var_min_val_1L_dbl = numeric(0),
                           existing_predrs_ls = NULL, 
                           max_nbr_of_covars_1L_int = integer(0),
                           start_at_int = c(2,1)){
  ready4show::write_report(params_ls = input_params_ls$params_ls,
                           paths_ls = input_params_ls$path_params_ls$paths_ls,
                           rprt_nm_1L_chr = "AAA_PMRY_ANLYS_MTH",
                           abstract_args_ls = abstract_args_ls,
                           header_yaml_args_ls = input_params_ls$header_yaml_args_ls)
  if(!is.null(input_params_ls$scndry_anlys_params_ls)){
    write_secondary_analyses(input_params_ls,
                             combinations_1L_lgl = combinations_1L_lgl, 
                             existing_predrs_ls = existing_predrs_ls, 
                             max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                             depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,)
    
  }
}
write_mdl_cmprsn <- function(scored_data_tb,
                             ds_smry_ls,
                             mdl_smry_ls,
                             output_data_dir_1L_chr,
                             depnt_var_max_val_1L_dbl = 0.99999, # MODIFIED
                             depnt_var_min_val_1L_dbl = 0.00001, #NEW
                             seed_1L_int = 1234){
  bl_tb <- transform_ds_for_tstng(scored_data_tb,# youthvars::
                                  depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                  candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                                  depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl,
                                  # round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                                  # round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr,
                                  round_var_nm_1L_chr = if(identical(ds_smry_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_smry_ls$round_var_nm_1L_chr},
                                  round_val_1L_chr = if(identical(ds_smry_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_smry_ls$round_bl_val_1L_chr}
  )
  ds_smry_ls$candidate_predrs_chr <- reorder_cndt_predrs_chr(ds_smry_ls$candidate_predrs_chr,
                                                             data_tb = bl_tb,
                                                             depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr)
  mdl_smry_ls <- add_prefd_predr_var_to_mdl_smry_ls(mdl_smry_ls,
                                                    ds_smry_ls = ds_smry_ls)
  mdl_smry_ls$smry_of_sngl_predr_mdls_tb <- write_sngl_predr_multi_mdls_outps(data_tb = bl_tb,
                                                                              folds_1L_int = mdl_smry_ls$folds_1L_int,
                                                                              mdl_types_chr = mdl_smry_ls$mdl_types_chr,
                                                                              depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                              predr_var_nm_1L_chr = mdl_smry_ls$predr_var_nm_1L_chr,
                                                                              predr_var_desc_1L_chr = mdl_smry_ls$predr_var_desc_1L_chr,
                                                                              predr_vals_dbl = mdl_smry_ls$predr_vals_dbl,
                                                                              path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                                                              new_dir_nm_1L_chr =  "A_Candidate_Mdls_Cmprsn",
                                                                              #start_1L_chr = NA_character_,
                                                                              mdl_types_lup = mdl_smry_ls$mdl_types_lup,
                                                                              dictionary_tb = ds_smry_ls$dictionary_tb)
  mdl_smry_ls$prefd_mdl_types_chr <- make_prefd_mdls_vec(mdl_smry_ls$smry_of_sngl_predr_mdls_tb,
                                                         choose_from_pfx_chr = mdl_smry_ls$choose_from_pfx_chr,
                                                         mdl_types_lup = mdl_smry_ls$mdl_types_lup)
  mdl_cmprsn_ls <- list(bl_tb = bl_tb,
                        ds_smry_ls = ds_smry_ls,
                        mdl_smry_ls = mdl_smry_ls)
  return(mdl_cmprsn_ls)
}
write_mdl_plts <- function (data_tb, model_mdl, mdl_fl_nm_1L_chr = "OLS_NTF", 
                            depnt_var_min_val_1L_dbl = numeric(0),
                            depnt_var_nm_1L_chr = "utl_total_w",
                            depnt_var_desc_1L_chr = "Utility score", # Remove defaults
                            tfmn_1L_chr = "NTF", predr_var_nm_1L_chr, predr_var_desc_1L_chr,
                            predr_vals_dbl, covar_var_nms_chr = NA_character_, path_to_write_to_1L_chr,
                            predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_,
                            plt_indcs_int = 1:5)
{
  data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                    depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
  tfd_data_tb <- transform_data_tb_for_cmprsn(data_tb, model_mdl = model_mdl, 
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                              predn_type_1L_chr = predn_type_1L_chr, tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                              family_1L_chr = family_1L_chr, tfmn_1L_chr = tfmn_1L_chr)
  if (1 %in% plt_indcs_int) {
    predn_ds_tb <- make_predn_ds_with_one_predr(model_mdl, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                                predr_vals_dbl = predr_vals_dbl, predn_type_1L_chr = predn_type_1L_chr)
  }else{
    predn_ds_tb <- NULL
  }
  purrr::pwalk(list(plt_fn_ls = list(plot_lnr_cmprsn,
                                     plot_auto_lm,
                                     plot_obsd_predd_dnst,
                                     plot_obsd_predd_dnst,
                                     plot_sctr_plt_cmprsn)[plt_indcs_int],
                    fn_args_ls_ls = list(list(data_tb = data_tb,
                                              predn_ds_tb = predn_ds_tb,
                                              predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                              predr_var_desc_1L_chr = predr_var_desc_1L_chr,
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                              depnt_var_desc_1L_chr = depnt_var_desc_1L_chr),
                                         list(model_mdl, which_dbl = 1:6, ncol_1L_int = 3L, label_size_1L_int = 3),
                                         list(tfd_data_tb = tfd_data_tb,
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                              depnt_var_desc_1L_chr = depnt_var_desc_1L_chr),
                                         list(tfd_data_tb = transform_data_tb_for_cmprsn(data_tb,
                                                                                         model_mdl = model_mdl,
                                                                                         depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                                         new_data_is_1L_chr = ifelse(!4 %in% plt_indcs_int,
                                                                                                                     "Predicted", "Simulated"),
                                                                                         predn_type_1L_chr = NULL,
                                                                                         tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl,
                                                                                         family_1L_chr = family_1L_chr,
                                                                                         tfmn_1L_chr = tfmn_1L_chr),
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                              depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                              predd_val_var_nm_1L_chr = "Simulated"),
                                         list(tfd_data_tb = tfd_data_tb,
                                              depnt_var_nm_1L_chr = depnt_var_nm_1L_chr))[plt_indcs_int],
                    plt_nm_sfx_chr = c("_LNR_CMPRSN",
                                       "_AUTOPLT", "_PRED_DNSTY", "_SIM_DNSTY", "_PRED_SCTR")[plt_indcs_int],
                    size_ls = list(c(6, 6), c(4, 7), c(6, 6), c(6, 6), c(6,
                                                                         6))[plt_indcs_int]), ~ ready4show::write_mdl_plt_fl(plt_fn = ..1,
                                                                                                                             fn_args_ls = ..2, path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                                                                                                             plt_nm_1L_chr = paste0(mdl_fl_nm_1L_chr, ifelse(!is.na(covar_var_nms_chr[1]),
                                                                                                                                                                             paste("_", paste0(covar_var_nms_chr[1:min(length(covar_var_nms_chr),
                                                                                                                                                                                                                       3)], collapse = "")), ""), ..3), height_1L_dbl = ..4[1],
                                                                                                                             width_1L_dbl = ..4[2]))
}
write_mdl_type_covars_mdls <- function (data_tb,
                                        depnt_var_min_val_1L_dbl = numeric(0),
                                        depnt_var_nm_1L_chr = "utl_total_w",# Remove default
                                        predrs_var_nms_chr,
                                        covar_var_nms_chr, mdl_type_1L_chr, path_to_write_to_1L_chr, new_dir_nm_1L_chr = "D_Covars_Selection",
                                        fl_nm_pfx_1L_chr = "D_CT", mdl_types_lup = NULL, start_1L_chr = NA_character_)
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  arg_vals_chr <- c("control_chr", "predn_type_chr","tfmn_chr") %>%
    purrr::map_chr(~ready4::get_from_lup_obj(mdl_types_lup,
                                             match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                             target_var_nm_1L_chr = .x, evaluate_1L_lgl = F))
  control_1L_chr <- arg_vals_chr[1]
  predn_type_1L_chr <- arg_vals_chr[2]
  tfmn_1L_chr <- arg_vals_chr[3]
  if (is.na(predn_type_1L_chr))
    predn_type_1L_chr <- NULL
  data_tb <- data_tb %>%
    add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                      tfmn_1L_chr = tfmn_1L_chr)
  output_dir_1L_chr <- output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr,
                                                               new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  smry_of_mdls_with_covars_tb <- purrr::map_dfr(predrs_var_nms_chr,
                                                ~{
                                                  model_mdl <- make_mdl(data_tb, 
                                                                        depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                        depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                        predr_var_nm_1L_chr = .x, covar_var_nms_chr = covar_var_nms_chr,
                                                                        tfmn_1L_chr = tfmn_1L_chr, mdl_type_1L_chr = mdl_type_1L_chr,
                                                                        control_1L_chr = control_1L_chr, mdl_types_lup = mdl_types_lup,
                                                                        start_1L_chr = start_1L_chr)
                                                  mdl_fl_nm_1L_chr <- paste0(fl_nm_pfx_1L_chr, "_",
                                                                             .x, "_", mdl_type_1L_chr)
                                                  saveRDS(model_mdl, paste0(output_dir_1L_chr, "/",
                                                                            mdl_fl_nm_1L_chr, ".RDS"))
                                                  if("summary.betareg" %in% class(summary(model_mdl))){
                                                    coefficients_mat <- summary(model_mdl)$coefficients$mean
                                                  }else{
                                                    coefficients_mat <- summary(model_mdl)$coefficients
                                                  }
                                                  tibble::tibble(variable = .x,
                                                                 Rsquare = caret::R2(stats::predict(model_mdl, type = predn_type_1L_chr) %>% calculate_depnt_var_tfmn(tfmn_1L_chr = tfmn_1L_chr,
                                                                                                                                                                      tfmn_is_outp_1L_lgl = T),
                                                                                     data_tb %>%
                                                                                       dplyr::pull(!!rlang::sym(depnt_var_nm_1L_chr)),
                                                                                     form = "traditional"),
                                                                 AIC = stats::AIC(model_mdl),
                                                                 BIC = stats::BIC(model_mdl),
                                                                 Significant = paste(names(which(coefficients_mat[,4] < 0.01)), collapse = " "))
                                                })
  smry_of_mdls_with_covars_tb <- smry_of_mdls_with_covars_tb %>%
    dplyr::arrange(dplyr::desc(AIC))
  saveRDS(smry_of_mdls_with_covars_tb, paste0(output_dir_1L_chr,
                                              "/", paste0(fl_nm_pfx_1L_chr, "_", "SMRY", "_", mdl_type_1L_chr),
                                              ".RDS"))
  return(smry_of_mdls_with_covars_tb)
}
write_mdl_type_multi_outps <- function (data_tb, folds_1L_int = 10, predrs_var_nms_chr, covar_var_nms_chr = NA_character_,
                                        start_1L_chr = NULL, mdl_type_1L_chr,
                                        depnt_var_min_val_1L_dbl = numeric(0),
                                        depnt_var_nm_1L_chr = "utl_total_w",# remove default
                                        path_to_write_to_1L_chr, new_dir_nm_1L_chr, mdl_types_lup = NULL, fl_nm_pfx_1L_chr = "C_PREDR",
                                        plt_indcs_int = c(3, 5))
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  smry_of_mdl_sngl_predrs_tb <- purrr::map_dfr(predrs_var_nms_chr,
                                               ~{
                                                 tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                                                                         match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                                                                         target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
                                                 mdl_smry_tb <- write_mdl_type_sngl_outps(data_tb,
                                                                                          folds_1L_int = folds_1L_int, 
                                                                                          depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                                          depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                                                          tfmn_1L_chr = tfmn_1L_chr, start_1L_chr = start_1L_chr,
                                                                                          predr_var_nm_1L_chr = .x, predr_var_desc_1L_chr = NA_character_,
                                                                                          predr_vals_dbl = NA_real_, covar_var_nms_chr = covar_var_nms_chr,
                                                                                          mdl_type_1L_chr = mdl_type_1L_chr, path_to_write_to_1L_chr = output_dir_1L_chr,
                                                                                          mdl_types_lup = mdl_types_lup, mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr,
                                                                                                                                                   "_", .x, "_", mdl_type_1L_chr), plt_indcs_int = plt_indcs_int)
                                                 if (!is.null(folds_1L_int)) {
                                                   mdl_smry_tb <- mdl_smry_tb %>% dplyr::select((-Model)) %>%
                                                     dplyr::mutate(Predictor = .x) %>% dplyr::select(Predictor,
                                                                                                     dplyr::everything())
                                                 }
                                                 mdl_smry_tb
                                               })
  if (!is.null(folds_1L_int)) {
    smry_of_mdl_sngl_predrs_tb <- smry_of_mdl_sngl_predrs_tb %>%
      dplyr::arrange(dplyr::desc(RsquaredP))
  }
  return(smry_of_mdl_sngl_predrs_tb)
}
write_mdl_type_sngl_outps <- function (data_tb, folds_1L_int = 10, depnt_var_min_val_1L_dbl = numeric(0),
                                       depnt_var_nm_1L_chr = "utl_total_w",
                                       start_1L_chr = NULL, tfmn_1L_chr = "NTF", predr_var_nm_1L_chr,
                                       predr_var_desc_1L_chr, predr_vals_dbl, covar_var_nms_chr = NA_character_,
                                       mdl_type_1L_chr = "OLS_NTF", mdl_types_lup = NULL, path_to_write_to_1L_chr,
                                       mdl_fl_nm_1L_chr, plt_indcs_int = NA_integer_)
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  arg_vals_chr <- c("control_chr", "family_chr", "predn_type_chr") %>%
    purrr::map_chr(~ready4::get_from_lup_obj(mdl_types_lup,
                                             match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                             target_var_nm_1L_chr = .x, evaluate_1L_lgl = F))
  control_1L_chr <- arg_vals_chr[1]
  family_1L_chr <- arg_vals_chr[2]
  predn_type_1L_chr <- arg_vals_chr[3]
  if (is.na(predn_type_1L_chr))
    predn_type_1L_chr <- NULL
  if (is.na(plt_indcs_int[1])) {
    plt_indcs_int <- 1:5
    if (!is.na(control_1L_chr)) {
      if (control_1L_chr %>% startsWith("betareg"))
        plt_indcs_int <- c(1, 3, 4, 5)
    }
  }
  tfmn_for_bnml_1L_lgl <- ready4::get_from_lup_obj(mdl_types_lup,
                                                   match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_type_1L_chr,
                                                   target_var_nm_1L_chr = "tfmn_for_bnml_lgl", evaluate_1L_lgl = F)
  data_tb <- data_tb %>%
    add_tfd_var_to_ds(depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                      tfmn_1L_chr = tfmn_1L_chr)
  model_mdl <- make_mdl(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                        depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                        tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                        covar_var_nms_chr = covar_var_nms_chr, mdl_type_1L_chr = mdl_type_1L_chr,
                        mdl_types_lup = mdl_types_lup, control_1L_chr = control_1L_chr,
                        start_1L_chr = start_1L_chr)
  write_mdl_plts(data_tb, model_mdl = model_mdl, mdl_fl_nm_1L_chr = mdl_fl_nm_1L_chr,
                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, tfmn_1L_chr = tfmn_1L_chr,
                 predr_var_nm_1L_chr = predr_var_nm_1L_chr, predr_var_desc_1L_chr = predr_var_desc_1L_chr,
                 predr_vals_dbl = predr_vals_dbl, covar_var_nms_chr = covar_var_nms_chr,
                 path_to_write_to_1L_chr = path_to_write_to_1L_chr, predn_type_1L_chr = predn_type_1L_chr,
                 tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, family_1L_chr = family_1L_chr,
                 plt_indcs_int = plt_indcs_int)
  if (!is.null(folds_1L_int)) {
    smry_of_one_predr_mdl_tb <- make_smry_of_mdl_outp(data_tb,
                                                      folds_1L_int = folds_1L_int, 
                                                      depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                      depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                      tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
                                                      mdl_type_1L_chr = mdl_type_1L_chr, mdl_types_lup = mdl_types_lup, start_1L_chr = start_1L_chr,
                                                      predn_type_1L_chr = predn_type_1L_chr)
  }
  else {
    smry_of_one_predr_mdl_tb <- tibble::tibble()
  }
  saveRDS(model_mdl, paste0(path_to_write_to_1L_chr, "/", mdl_fl_nm_1L_chr,
                            ".RDS"))
  return(smry_of_one_predr_mdl_tb)
}
write_mdls_to_dv <- function(outp_smry_ls,
                             new_dir_nm_1L_chr = "G_Shareable",
                             shareable_title_detail_1L_chr = "",
                             output_dir_chr = NA_character_){
  if(is.na(output_dir_chr[1]))
    output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  if (!is.null(outp_smry_ls$dv_ls)) {
    write_shareable_mdls_to_dv(outp_smry_ls,
                               new_dir_nm_1L_chr = new_dir_nm_1L_chr,
                               share_ingredients_1L_lgl = T,
                               output_dir_chr = output_dir_chr)
    #if(write_mdls_to_dv_1L_lgl){
    outp_smry_ls$shareable_mdls_tb <- write_shareable_mdls_to_dv(outp_smry_ls,
                                                                 new_dir_nm_1L_chr = new_dir_nm_1L_chr,
                                                                 shareable_title_detail_1L_chr = shareable_title_detail_1L_chr,
                                                                 share_ingredients_1L_lgl = F,
                                                                 output_dir_chr = output_dir_chr)
    #}
  }
  return(outp_smry_ls)
}
write_mdls_with_covars_cmprsn <- function(scored_data_tb,
                                          bl_tb,
                                          combinations_1L_lgl = F, 
                                          depnt_var_min_val_1L_dbl = numeric(0),
                                          ds_smry_ls,
                                          existing_predrs_ls = NULL, 
                                          max_nbr_of_covars_1L_int = integer(0),
                                          mdl_smry_ls,
                                          output_data_dir_1L_chr,
                                          seed_1L_int = 1234,
                                          session_data_ls = NULL){
  empty_tb <- write_mdl_type_multi_outps(data_tb = bl_tb,
                                         folds_1L_int = NULL,
                                         start_1L_chr = NA_character_,
                                         predrs_var_nms_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr,
                                         covar_var_nms_chr = mdl_smry_ls$prefd_covars_chr,
                                         mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1],
                                         depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                         depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                         path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                         new_dir_nm_1L_chr = "E_Predrs_W_Covars_Sngl_Mdl_Cmprsn",
                                         mdl_types_lup = mdl_smry_ls$mdl_types_lup,
                                         fl_nm_pfx_1L_chr = "E_CK_CV")
  mdl_smry_ls$predr_vars_nms_ls <- make_predr_vars_nms_ls(main_predrs_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr,
                                                          covars_ls = list(mdl_smry_ls$prefd_covars_chr),
                                                          combinations_1L_lgl = combinations_1L_lgl,
                                                          existing_predrs_ls = existing_predrs_ls,
                                                          max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int)
  mdl_smry_ls$mdl_nms_ls <- make_mdl_nms_ls(mdl_smry_ls$predr_vars_nms_ls,
                                            mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr)
  outp_smry_ls <- list(scored_data_tb = scored_data_tb,
                       dictionary_tb = ds_smry_ls$dictionary_tb,
                       predictors_lup = ds_smry_ls$predictors_lup,
                       smry_of_sngl_predr_mdls_tb = mdl_smry_ls$smry_of_sngl_predr_mdls_tb,
                       prefd_mdl_types_chr = mdl_smry_ls$prefd_mdl_types_chr,
                       predr_cmprsn_tb = mdl_smry_ls$predr_cmprsn_tb,
                       smry_of_mdl_sngl_predrs_tb = mdl_smry_ls$smry_of_mdl_sngl_predrs_tb,
                       mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb,
                       signt_covars_chr = mdl_smry_ls$signt_covars_chr,
                       prefd_covars_chr = mdl_smry_ls$prefd_covars_chr,
                       depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                       predr_vars_nms_ls = mdl_smry_ls$predr_vars_nms_ls,
                       mdl_nms_ls = mdl_smry_ls$mdl_nms_ls,
                       id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr,
                       round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                       round_bl_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr,
                       path_to_write_to_1L_chr = output_data_dir_1L_chr,
                       seed_1L_int = seed_1L_int,
                       folds_1L_int = mdl_smry_ls$folds_1L_int,
                       max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int,
                       mdl_types_lup = mdl_smry_ls$mdl_types_lup,
                       file_paths_chr = list.files(output_data_dir_1L_chr, recursive = T),
                       session_data_ls = session_data_ls)
  saveRDS(outp_smry_ls,
          paste0(outp_smry_ls$path_to_write_to_1L_chr,"/I_ALL_OUTPUT_.RDS"))
  return(outp_smry_ls)
}
write_predr_and_covars_cmprsn <- function(scored_data_tb,
                                          bl_tb,
                                          depnt_var_min_val_1L_dbl = numeric(0),
                                          ds_smry_ls,
                                          mdl_smry_ls,
                                          output_data_dir_1L_chr,
                                          seed_1L_int = 1234,
                                          signft_covars_cndn_1L_chr = "any"){
  mdl_smry_ls$predr_cmprsn_tb <- write_predr_cmprsn_outps(data_tb = bl_tb,
                                                          depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                          path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                                          new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsn",
                                                          depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                                          candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                                                          max_nbr_of_boruta_mdl_runs_int = mdl_smry_ls$max_nbr_of_boruta_mdl_runs_int)
  if(identical(mdl_smry_ls$predr_cmprsn_tb$predr_chr, character(0))){
    stop("No important predictors identified - execution aborted. Try specifying other predictors.")
  }
  mdl_smry_ls$smry_of_mdl_sngl_predrs_tb <- write_mdl_type_multi_outps(data_tb = bl_tb,
                                                                       folds_1L_int = mdl_smry_ls$folds_1L_int,
                                                                       predrs_var_nms_chr = mdl_smry_ls$predr_cmprsn_tb$predr_chr,
                                                                       mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1],
                                                                       depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                       depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                                                       path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                                                       new_dir_nm_1L_chr = "C_Predrs_Sngl_Mdl_Cmprsn",
                                                                       fl_nm_pfx_1L_chr = "C_PREDR",
                                                                       start_1L_chr = NA_character_,
                                                                       mdl_types_lup = mdl_smry_ls$mdl_types_lup)
  bl_tb <- scored_data_tb %>%
    transform_ds_for_tstng(depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,#youthvars::
                           candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                           covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr,
                           remove_all_msng_1L_lgl = T,
                           round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                           round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr)
  mdl_smry_ls$mdls_with_covars_smry_tb <- write_mdl_type_covars_mdls(bl_tb,
                                                                     depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                     depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                                                     predrs_var_nms_chr = ds_smry_ls$candidate_predrs_chr,
                                                                     covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr,
                                                                     mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1],
                                                                     path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                                                     new_dir_nm_1L_chr = "D_Predr_Covars_Cmprsn",
                                                                     fl_nm_pfx_1L_chr = "D_CT",
                                                                     mdl_types_lup = mdl_smry_ls$mdl_types_lup)
  mdl_smry_ls$signt_covars_chr <- get_signft_covars(mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb,
                                                    covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr,
                                                    X_Ready4useDyad = Ready4useDyad(ds_tb = scored_data_tb,
                                                                                    dictionary_r3 = ds_smry_ls$dictionary_tb),
                                                    what_1L_chr = signft_covars_cndn_1L_chr)
  predr_and_covars_cmprsn_ls <- list(bl_tb = bl_tb,
                                     ds_smry_ls = ds_smry_ls,
                                     mdl_smry_ls = mdl_smry_ls)
  return(predr_and_covars_cmprsn_ls)
}
write_predr_and_mdl_tstng_results <- function(scored_data_tb,
                                              combinations_1L_lgl = F,
                                              depnt_var_max_val_1L_dbl = 0.99999, 
                                              depnt_var_min_val_1L_dbl = 0.00001, 
                                              ds_smry_ls,
                                              existing_predrs_ls = NULL,
                                              max_nbr_of_covars_1L_int = integer(0),
                                              mdl_smry_ls,
                                              session_data_ls,
                                              output_data_dir_1L_chr,
                                              seed_1L_int = 1234,
                                              signft_covars_cndn_1L_chr = "any"){
  
  cmprsn_ls <- write_mdl_cmprsn(scored_data_tb = scored_data_tb,
                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl, 
                                ds_smry_ls = ds_smry_ls,
                                mdl_smry_ls = mdl_smry_ls,
                                output_data_dir_1L_chr = output_data_dir_1L_chr,
                                seed_1L_int = seed_1L_int)
  cmprsn_ls <- write_predr_and_covars_cmprsn(scored_data_tb = scored_data_tb,
                                             bl_tb = cmprsn_ls$bl_tb,
                                             depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                             ds_smry_ls = cmprsn_ls$ds_smry_ls,
                                             mdl_smry_ls  = cmprsn_ls$mdl_smry_ls,
                                             output_data_dir_1L_chr = output_data_dir_1L_chr,
                                             seed_1L_int = seed_1L_int,
                                             signft_covars_cndn_1L_chr  = signft_covars_cndn_1L_chr)
  if(ifelse(is.null(cmprsn_ls$mdl_smry_ls$prefd_covars_chr),
            T,
            is.na(cmprsn_ls$mdl_smry_ls$prefd_covars_chr[1]))){
    cmprsn_ls$mdl_smry_ls$prefd_covars_chr <- cmprsn_ls$mdl_smry_ls$signt_covars_chr
  }
  outp_smry_ls <- write_mdls_with_covars_cmprsn(scored_data_tb = scored_data_tb,
                                                bl_tb = cmprsn_ls$bl_tb,
                                                combinations_1L_lgl = combinations_1L_lgl,
                                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                ds_smry_ls = cmprsn_ls$ds_smry_ls,
                                                existing_predrs_ls = existing_predrs_ls,
                                                max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                                                mdl_smry_ls = cmprsn_ls$mdl_smry_ls,
                                                output_data_dir_1L_chr = output_data_dir_1L_chr,
                                                seed_1L_int = seed_1L_int,
                                                session_data_ls = session_data_ls)
  return(outp_smry_ls)
}
write_predr_cmprsn_outps <- function (data_tb, path_to_write_to_1L_chr, new_dir_nm_1L_chr = "B_Candidate_Predrs_Cmprsn", 
                                      depnt_var_min_val_1L_dbl = numeric(0),
                                      depnt_var_nm_1L_chr = "utl_total_w",
                                      candidate_predrs_chr, max_nbr_of_boruta_mdl_runs_int = 300L)
{
  if (length(candidate_predrs_chr) > 1) {
    covar_var_nms_chr <- candidate_predrs_chr[2:length(candidate_predrs_chr)]
  } else {
    covar_var_nms_chr <- NA_character_
  }
  data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                    depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                    predr_var_nm_1L_chr = candidate_predrs_chr[1], covar_var_nms_chr = covar_var_nms_chr)
  rf_mdl <- randomForest::randomForest(stats::as.formula(paste0(depnt_var_nm_1L_chr,
                                                                " ~ .")), data = data_tb, importance = TRUE)
  boruta_mdl <- Boruta::Boruta(stats::as.formula(paste0(depnt_var_nm_1L_chr,
                                                        " ~ .")), data = data_tb, maxRuns = max_nbr_of_boruta_mdl_runs_int)
  output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  purrr::pwalk(list(fn_ls = list(randomForest::varImpPlot,
                                 plot), fn_args_ls_ls = list(list(rf_mdl, main = ""),
                                                             list(boruta_mdl, cex = 1.5, cex.axis = 0.8, las = 2,
                                                                  xlab = "", main = "")), plt_nm_sfx_chr = c("_RF_VAR_IMP",
                                                                                                             "_BORUTA_VAR_IMP"), size_ls = list(c(6, 6), c(4, 6))),
               ~ ready4show::write_mdl_plt_fl(plt_fn = ..1, fn_args_ls = ..2,
                                              path_to_write_to_1L_chr = output_dir_1L_chr,
                                              plt_nm_1L_chr = paste0("B_PRED_CMPRSN", ..3), height_1L_dbl = ..4[1],
                                              width_1L_dbl = ..4[2]))
  confirmed_predrs_chr <- names(boruta_mdl$finalDecision)[boruta_mdl$finalDecision ==
                                                            "Confirmed"]
  confirmed_predrs_tb <- rf_mdl$importance %>% tibble::as_tibble(rownames = "predr_chr") %>%
    dplyr::arrange(dplyr::desc(`%IncMSE`)) %>% dplyr::filter(predr_chr %in%
                                                               confirmed_predrs_chr)
  return(confirmed_predrs_tb)
}
write_scndry_analysis <- function(valid_params_ls_ls,
                                  candidate_covar_nms_chr,
                                  path_params_ls,
                                  reference_1L_int,
                                  combinations_1L_lgl = F, 
                                  depnt_var_max_val_1L_dbl = 0.99999, 
                                  depnt_var_min_val_1L_dbl = 0.00001, 
                                  existing_predrs_ls = NULL, 
                                  max_nbr_of_covars_1L_int = integer(0),
                                  backend_1L_chr = "cmdstanr",
                                  candidate_predrs_chr = NULL,
                                  new_dir_nm_1L_chr = "F_TS_Mdls",
                                  predictors_lup = NULL,
                                  prefd_covars_chr = NA_character_,
                                  signft_covars_cndn_1L_chr = "any"){
  analysis_params_ls <- valid_params_ls_ls$params_ls %>%
    append(path_params_ls[1:2])
  rename_lup <- valid_params_ls_ls$rename_lup
  if(!is.null(predictors_lup)){
    predictors_lup$short_name_chr <- predictors_lup$short_name_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
    analysis_params_ls$predictors_lup <- predictors_lup
  }
  if(!is.null(candidate_predrs_chr)){
    candidate_predrs_chr <- candidate_predrs_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
    analysis_params_ls$ds_descvs_ls$candidate_predrs_chr <- candidate_predrs_chr
  }
  if(!is.null(candidate_covar_nms_chr)){
    candidate_covar_nms_chr <- candidate_covar_nms_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  if(ifelse(is.null(prefd_covars_chr),F,!is.na(prefd_covars_chr))){
    prefd_covars_chr <- prefd_covars_chr %>%
      purrr::map_chr(~ifelse(!.x %in% rename_lup$old_nms_chr,
                             .x,
                             ready4::get_from_lup_obj(rename_lup,
                                                      match_value_xx = .x,
                                                      match_var_nm_1L_chr = "old_nms_chr",
                                                      target_var_nm_1L_chr = "new_nms_chr",
                                                      evaluate_1L_lgl = F)))
  }
  analysis_params_ls$prefd_covars_chr <- prefd_covars_chr
  analysis_params_ls$candidate_covar_nms_chr <- candidate_covar_nms_chr
  path_params_ls$paths_ls <- write_scndry_analysis_dir(path_params_ls$paths_ls,
                                                       reference_1L_int = reference_1L_int)
  params_ls <- list(candidate_predrs_chr = candidate_predrs_chr,
                    transform_paths_ls = list(fn = transform_paths_ls_for_scndry,
                                              args_ls = list(reference_1L_int = reference_1L_int))) %>%
    append(analysis_params_ls)
  params_ls$utl_class_fn_1L_chr <- params_ls$raw_ds_tfmn_fn <- NULL
  params_ls_ls <- transform_params_ls_to_valid(params_ls)
  params_ls <- params_ls_ls %>%
    purrr::pluck("params_ls") %>%
    append(list(rename_lup = params_ls_ls$rename_lup))
  outp_smry_ls <- valid_params_ls_ls$outp_smry_ls
  mdl_smry_ls <- params_ls$mdl_smry_ls
  data_tb <- outp_smry_ls$scored_data_tb
  ds_smry_ls <- params_ls$ds_descvs_ls %>%
    make_analysis_ds_smry_ls(candidate_covar_nms_chr = params_ls$candidate_covar_nms_chr,
                             predictors_lup = params_ls$predictors_lup)
  ds_smry_ls$candidate_predrs_chr <- params_ls$candidate_predrs_chr
  existing_mdls_chr <- outp_smry_ls[["mdl_nms_ls"]] %>% purrr::flatten_chr()
  existing_predrs_ls <- outp_smry_ls$predr_vars_nms_ls
  ### WRITE CNDT MDL TESTS
  cmprsns_ls <- write_mdl_cmprsn(scored_data_tb = data_tb,
                                 depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                 depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl, 
                                 ds_smry_ls = ds_smry_ls,
                                 mdl_smry_ls = mdl_smry_ls,
                                 output_data_dir_1L_chr = path_params_ls$paths_ls$write_to_dir_nm_1L_chr,
                                 seed_1L_int = params_ls$seed_1L_int)
  if(!is.null(params_ls$prefd_mdl_types_chr)){
    cmprsns_ls$mdl_smry_ls$prefd_mdl_types_chr <- params_ls$prefd_mdl_types_chr
  }
  cmprsns_ls <- write_predr_and_covars_cmprsn(scored_data_tb = data_tb,
                                              bl_tb = cmprsns_ls$bl_tb,
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              ds_smry_ls = cmprsns_ls$ds_smry_ls,
                                              mdl_smry_ls  = cmprsns_ls$mdl_smry_ls,
                                              output_data_dir_1L_chr = path_params_ls$paths_ls$write_to_dir_nm_1L_chr,
                                              seed_1L_int = params_ls$seed_1L_int,
                                              signft_covars_cndn_1L_chr = signft_covars_cndn_1L_chr)
  if(!is.null(params_ls$prefd_covars_chr)){
    cmprsns_ls$mdl_smry_ls$prefd_covars_chr <- params_ls$prefd_covars_chr
  }
  outp_smry_ls <- write_mdls_with_covars_cmprsn(scored_data_tb = data_tb,
                                                bl_tb = cmprsns_ls$bl_tb,
                                                combinations_1L_lgl = F, # Correct
                                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                ds_smry_ls = cmprsns_ls$ds_smry_ls,
                                                existing_predrs_ls = NULL, # Correct
                                                max_nbr_of_covars_1L_int = integer(0), # Correct
                                                mdl_smry_ls = cmprsns_ls$mdl_smry_ls,
                                                output_data_dir_1L_chr = path_params_ls$paths_ls$write_to_dir_nm_1L_chr,
                                                seed_1L_int = params_ls$seed_1L_int,
                                                session_data_ls = sessionInfo())
  outp_smry_ls$mdl_nms_ls <- outp_smry_ls$mdl_nms_ls %>%
    purrr::map(~.x[!.x %in% existing_mdls_chr]) %>%
    purrr::compact()
  outp_smry_ls$predr_vars_nms_ls <- outp_smry_ls$predr_vars_nms_ls[outp_smry_ls$predr_vars_nms_ls %>% purrr::map_lgl(~{
    test_chr <- .x
    !any(existing_predrs_ls %>% purrr::map_lgl(~identical(.x,test_chr))
    )})]
  outp_smry_ls <- write_ts_mdls_from_alg_outp(outp_smry_ls = outp_smry_ls,
                                              backend_1L_chr = backend_1L_chr,
                                              combinations_1L_lgl = combinations_1L_lgl, 
                                              control_ls = params_ls$control_ls,
                                              depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                              existing_predrs_ls = existing_predrs_ls, 
                                              iters_1L_int = params_ls$iters_1L_int,
                                              max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                                              new_dir_nm_1L_chr = new_dir_nm_1L_chr,
                                              path_to_write_to_1L_chr = outp_smry_ls$path_to_write_to_1L_chr,
                                              predictors_lup = params_ls$predictors_lup,
                                              prior_ls = params_ls$prior_ls,
                                              utl_min_val_1L_dbl = params_ls$utl_min_val_1L_dbl,)
  return(outp_smry_ls)
}
write_scndry_analysis_dir <- function(paths_ls,
                                      reference_1L_int = 1){
  paths_ls <- transform_paths_ls_for_scndry(paths_ls,
                                            reference_1L_int = reference_1L_int)
  ready4::write_new_dirs(paths_ls$write_to_dir_nm_1L_chr)
  return(paths_ls)
}
write_secondary_analyses <- function(input_params_ls,
                                     backend_1L_chr = "cmdstanr",
                                     combinations_1L_lgl = F, 
                                     existing_predrs_ls = NULL, 
                                     max_nbr_of_covars_1L_int = integer(0),
                                     depnt_var_min_val_1L_dbl = numeric(0),
                                     new_dir_nm_1L_chr = "F_TS_Mdls"){
  references_int <- 1:length(input_params_ls$scndry_anlys_params_ls)
  results_ls <- references_int %>%
    purrr::map(~{
      changes_ls <- input_params_ls$scndry_anlys_params_ls %>%
        purrr::pluck(.x)
      if(is.null(changes_ls$candidate_covar_nms_chr))
        changes_ls$candidate_covar_nms_chr <- input_params_ls$params_ls$candidate_covar_nms_chr %>% transform_names(input_params_ls$rename_lup, invert_1L_lgl = T)
      if(is.null(changes_ls$candidate_predrs_chr)){
        changes_ls$candidate_covar_nms_chr <- changes_ls$candidate_covar_nms_chr[!changes_ls$candidate_covar_nms_chr %in% changes_ls$candidate_predrs_chr]
      }
      write_scndry_analysis(valid_params_ls_ls = input_params_ls,
                            backend_1L_chr = backend_1L_chr,
                            candidate_covar_nms_chr = changes_ls$candidate_covar_nms_chr,
                            candidate_predrs_chr = changes_ls$candidate_predrs_chr,
                            combinations_1L_lgl = combinations_1L_lgl, 
                            depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                            existing_predrs_ls = existing_predrs_ls, 
                            max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int,
                            new_dir_nm_1L_chr = new_dir_nm_1L_chr,
                            path_params_ls = input_params_ls$path_params_ls,
                            predictors_lup = changes_ls$predictors_lup,
                            prefd_covars_chr = changes_ls$prefd_covars_chr,
                            reference_1L_int = .x)})
  return(results_ls)
}
write_shareable_mdls <- function (outp_smry_ls,
                                  depnt_var_min_val_1L_dbl = numeric(0),
                                  new_dir_nm_1L_chr = "G_Shareable",
                                  shareable_title_detail_1L_chr = "",
                                  write_mdls_to_dv_1L_lgl = F)
{
  output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls,
                                        new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  incld_mdl_paths_chr <- make_incld_mdl_paths(outp_smry_ls)
  fake_ds_tb <- make_fake_ts_data(outp_smry_ls, depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, depnt_vars_are_NA_1L_lgl = F)
  mdl_types_lup <- outp_smry_ls$mdl_types_lup
  shareable_mdls_ls <- outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr() %>%
    purrr::map2(incld_mdl_paths_chr, ~{
      model_mdl <- readRDS(paste0(outp_smry_ls$path_to_write_to_1L_chr,"/",.y))
      
      mdl_smry_tb <- outp_smry_ls$mdls_smry_tb %>% dplyr::filter(Model ==
                                                                   .x)
      mdl_nm_1L_chr <- .x
      mdl_type_1L_chr <- get_mdl_type_from_nm(mdl_nm_1L_chr,
                                              mdl_types_lup = mdl_types_lup)
      tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                              match_value_xx = mdl_type_1L_chr,
                                              match_var_nm_1L_chr = "short_name_chr",
                                              target_var_nm_1L_chr = "tfmn_chr",
                                              evaluate_1L_lgl = F)
      predn_type_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                                    match_value_xx = mdl_type_1L_chr,
                                                    match_var_nm_1L_chr = "short_name_chr",
                                                    target_var_nm_1L_chr = "predn_type_chr",
                                                    evaluate_1L_lgl = F)
      if (is.na(predn_type_1L_chr))
        predn_type_1L_chr <- NULL
      control_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                                 match_value_xx = mdl_type_1L_chr,
                                                 match_var_nm_1L_chr = "short_name_chr",
                                                 target_var_nm_1L_chr = "control_chr",
                                                 evaluate_1L_lgl = F)
      sd_dbl <- mdl_smry_tb %>%
        dplyr::filter(Parameter == "SD (Intercept)") %>%
        dplyr::select(Estimate, SE) %>%
        t() %>%
        as.vector()
      mdl_fake_ds_tb <- fake_ds_tb %>%
        add_tfd_var_to_ds(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                          tfmn_1L_chr = tfmn_1L_chr,
                          depnt_var_max_val_1L_dbl = 0.999) %>%
        dplyr::select(names(model_mdl$data))
      model_mdl$data <- mdl_fake_ds_tb
      table_predn_mdl <- make_shareable_mdl(fake_ds_tb = mdl_fake_ds_tb,
                                            mdl_smry_tb = mdl_smry_tb, depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                                            id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr,
                                            tfmn_1L_chr = tfmn_1L_chr,
                                            mdl_type_1L_chr = mdl_type_1L_chr,
                                            mdl_types_lup = mdl_types_lup,
                                            control_1L_chr = control_1L_chr,
                                            start_1L_chr = NA_character_,
                                            seed_1L_int = outp_smry_ls$seed_1L_int)
      saveRDS(table_predn_mdl, paste0(output_dir_chr[4], "/", .x,
                                      ".RDS"))
      saveRDS(model_mdl, paste0(output_dir_chr[3], "/", .x,
                                ".RDS"))
      scaling_fctr_dbl <- make_scaling_fctr_dbl(outp_smry_ls)
      write_ts_mdl_plts(brms_mdl = model_mdl,
                        table_predn_mdl = table_predn_mdl,
                        tfd_data_tb = outp_smry_ls$scored_data_tb %>%
                          transform_tb_to_mdl_inp(depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,##
                                                  depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                                                  predr_vars_nms_chr = outp_smry_ls$predr_vars_nms_ls %>% purrr::flatten_chr() %>% unique(),
                                                  id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr,
                                                  round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                                                  round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr,
                                                  scaling_fctr_dbl = scaling_fctr_dbl),
                        depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                        mdl_nm_1L_chr = mdl_nm_1L_chr,
                        path_to_write_to_1L_chr = output_dir_chr[3],
                        predn_type_1L_chr = predn_type_1L_chr,
                        round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                        sd_dbl = sd_dbl,
                        sfx_1L_chr = " from table",
                        tfmn_1L_chr = tfmn_1L_chr,
                        utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl),
                                                    outp_smry_ls$utl_min_val_1L_dbl,
                                                    -1))
      table_predn_mdl
    }) %>% stats::setNames(outp_smry_ls$mdl_nms_ls %>% purrr::flatten_chr())
  outp_smry_ls$shareable_mdls_ls <- shareable_mdls_ls
  outp_smry_ls$shareable_mdls_tb <-  NULL
  ingredients_ls <- list(depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr,
                         dictionary_tb = outp_smry_ls$dictionary_tb %>%
                           dplyr::filter(var_nm_chr %in% names(fake_ds_tb)),
                         id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr,
                         fake_ds_tb = fake_ds_tb,
                         mdls_lup = outp_smry_ls$shareable_mdls_ls %>%
                           purrr::map2_dfr(names(outp_smry_ls$shareable_mdls_ls),
                                           ~{
                                             if(inherits(.x,"betareg")){
                                               coeffs_dbl <- .x$coefficients$mean
                                             }else{
                                               coeffs_dbl <- .x$coefficients
                                             }
                                             mdl_type_1L_chr = get_mdl_type_from_nm(.y,
                                                                                    mdl_types_lup = outp_smry_ls$mdl_types_lup)
                                             tibble::tibble(mdl_nms_chr = .y) %>%
                                               dplyr::mutate(predrs_ls = list(coeffs_dbl %>%
                                                                                names() %>%
                                                                                stringr::str_remove_all("_change") %>%
                                                                                stringr::str_remove_all("_baseline") %>%
                                                                                stringr::str_remove_all("_scaled") %>%
                                                                                stringr::str_remove_all("_unscaled") %>%
                                                                                unique() %>%
                                                                                purrr::discard(~ .x== "(Intercept)")),
                                                             mdl_type_chr = mdl_type_1L_chr,
                                                             tfmn_chr = ready4::get_from_lup_obj(outp_smry_ls$mdl_types_lup,
                                                                                                 match_value_xx = mdl_type_1L_chr,
                                                                                                 match_var_nm_1L_chr = "short_name_chr",
                                                                                                 target_var_nm_1L_chr = "tfmn_chr",
                                                                                                 evaluate_1L_lgl = F))
                                           }), #
                         mdls_smry_tb = outp_smry_ls$mdls_smry_tb,#
                         mdl_types_lup = mdl_types_lup,
                         predictors_lup = outp_smry_ls$predictors_lup,#
                         round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                         seed_1L_int = outp_smry_ls$seed_1L_int,
                         utl_min_val_1L_dbl = ifelse(!is.null(outp_smry_ls$utl_min_val_1L_dbl),
                                                     outp_smry_ls$utl_min_val_1L_dbl,
                                                     -1))
  saveRDS(ingredients_ls, paste0(output_dir_chr[2], "/mdl_ingredients",
                                 ".RDS"))
  outp_smry_ls <- write_mdls_to_dv(outp_smry_ls,
                                   new_dir_nm_1L_chr = new_dir_nm_1L_chr,
                                   shareable_title_detail_1L_chr = shareable_title_detail_1L_chr,
                                   output_dir_chr = output_dir_chr)
  return(outp_smry_ls)
}
write_shareable_mdls_to_dv <- function (outp_smry_ls,
                                        new_dir_nm_1L_chr = "G_Shareable",
                                        shareable_title_detail_1L_chr = "",
                                        share_ingredients_1L_lgl = T,
                                        output_dir_chr = NA_character_){
  if(is.na(output_dir_chr[1]))
    output_dir_chr <- write_shareable_dir(outp_smry_ls = outp_smry_ls,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  if(share_ingredients_1L_lgl){
    shareable_mdls_tb <- tibble::tibble(ds_obj_nm_chr = "mdl_ingredients",
                                        title_chr = "An R object that can be used to construct model objects from tables of coefficients. Contains a synthetic dataset.")
  }else{
    shareable_mdls_tb <- tibble::tibble(ds_obj_nm_chr = names(outp_smry_ls$shareable_mdls_ls),
                                        title_chr = paste0("A shareable (contains no confidential data) statistical model, ",
                                                           names(outp_smry_ls$shareable_mdls_ls), ".",
                                                           shareable_title_detail_1L_chr))
  }
  ready4::write_to_dv_with_wait(shareable_mdls_tb, dv_nm_1L_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr,
                                ds_url_1L_chr = outp_smry_ls$dv_ls$ds_url_1L_chr, parent_dv_dir_1L_chr = outp_smry_ls$dv_ls$parent_dv_dir_1L_chr,
                                paths_to_dirs_chr = output_dir_chr[ifelse(share_ingredients_1L_lgl,
                                                                          2,
                                                                          3)], paths_are_rltv_1L_lgl = F,inc_fl_types_chr = ".RDS")
  if(!share_ingredients_1L_lgl){
    ds_ls <- dataverse::get_dataset(outp_smry_ls$dv_ls$ds_url_1L_chr)
    shareable_mdls_tb <- shareable_mdls_tb %>%
      dplyr::mutate(dv_nm_chr = outp_smry_ls$dv_ls$dv_nm_1L_chr,
                    fl_ids_int = ds_obj_nm_chr %>%
                      purrr::map_int(~ready4::get_fl_id_from_dv_ls(ds_ls,
                                                                   fl_nm_1L_chr = paste0(.x, ".RDS")) %>%
                                       as.integer()))
    
  }
  return(shareable_mdls_tb)
}
write_sngl_predr_multi_mdls_outps <- function (data_tb, mdl_types_chr, predr_var_nm_1L_chr, predr_var_desc_1L_chr,
                                               predr_vals_dbl, path_to_write_to_1L_chr, new_dir_nm_1L_chr =  "A_Candidate_Mdls_Cmprsn", start_1L_chr = NULL,
                                               covar_var_nms_chr = NA_character_,
                                               depnt_var_min_val_1L_dbl = numeric(0), # NEW
                                               depnt_var_nm_1L_chr = "utl_total_w", # Remove default
                                               folds_1L_int = 10, mdl_types_lup = NULL, fl_nm_pfx_1L_chr = "A_RT_",
                                               plt_indcs_int = NA_integer_, dictionary_tb)
{
  if (is.null(mdl_types_lup))
    utils::data("mdl_types_lup", envir = environment())
  data_tb <- transform_ds_for_mdlng(data_tb, depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                    depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl, 
                                    predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr)
  output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  ready4show::write_mdl_plt_fl(plt_fn = make_tfmn_cmprsn_plt,
                               fn_args_ls = list(data_tb = data_tb,
                                                 depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                 dictionary_tb = dictionary_tb),
                               path_to_write_to_1L_chr = output_dir_1L_chr,
                               plt_nm_1L_chr = "A_TFMN_CMPRSN_DNSTY",
                               height_1L_dbl = 6,
                               width_1L_dbl = 10)
  smry_of_sngl_predr_mdls_tb <- purrr::map_dfr(mdl_types_chr,
                                               ~{
                                                 tfmn_1L_chr <- ready4::get_from_lup_obj(mdl_types_lup,
                                                                                         match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x,
                                                                                         target_var_nm_1L_chr = "tfmn_chr", evaluate_1L_lgl = F)
                                                 write_mdl_type_sngl_outps(data_tb, folds_1L_int = folds_1L_int,
                                                                           depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                                                           depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, start_1L_chr = start_1L_chr,
                                                                           tfmn_1L_chr = tfmn_1L_chr, predr_var_nm_1L_chr = predr_var_nm_1L_chr,
                                                                           predr_var_desc_1L_chr = predr_var_desc_1L_chr,
                                                                           predr_vals_dbl = predr_vals_dbl, covar_var_nms_chr = covar_var_nms_chr,
                                                                           mdl_type_1L_chr = .x, path_to_write_to_1L_chr = output_dir_1L_chr,
                                                                           mdl_types_lup = mdl_types_lup, mdl_fl_nm_1L_chr = paste0(fl_nm_pfx_1L_chr,
                                                                                                                                    predr_var_nm_1L_chr, "_", .x), plt_indcs_int = plt_indcs_int)
                                               })
  if (!is.null(folds_1L_int))
    smry_of_sngl_predr_mdls_tb <- smry_of_sngl_predr_mdls_tb %>%
    dplyr::arrange(dplyr::desc(RsquaredP))
  return(smry_of_sngl_predr_mdls_tb)
}
write_ts_mdls <- function (data_tb, # Rename lngl
                           depnt_var_min_val_1L_dbl = numeric(0),
                           depnt_var_nm_1L_chr = "utl_total_w", #Remove default
                           predr_vars_nms_ls,
                           id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round",
                           round_bl_val_1L_chr = "Baseline", utl_min_val_1L_dbl = -1, backend_1L_chr = getOption("brms.backend",
                                                                                                                 "rstan"),
                           mdl_nms_ls, mdl_smry_dir_1L_chr, predictors_lup, iters_1L_int = 4000L,
                           mdl_types_lup, seed_1L_int = 1000L, prior_ls = NULL, control_ls = NULL){
  if (!dir.exists(mdl_smry_dir_1L_chr))
    dir.create(mdl_smry_dir_1L_chr)
  mdls_smry_tb <- purrr::map_dfr(1:length(mdl_nms_ls), ~{
    idx_1L_int <- .x
    purrr::map_dfr(mdl_nms_ls[[idx_1L_int]], ~{
      smry_ls <- make_smry_of_ts_mdl_outp(data_tb = data_tb,
                                          depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                          predr_vars_nms_chr = predr_vars_nms_ls[[idx_1L_int]],
                                          mdl_nm_1L_chr = .x,
                                          path_to_write_to_1L_chr = mdl_smry_dir_1L_chr,
                                          depnt_var_nm_1L_chr = depnt_var_nm_1L_chr, id_var_nm_1L_chr = id_var_nm_1L_chr,
                                          round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr,
                                          predictors_lup = predictors_lup, utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                          backend_1L_chr = backend_1L_chr, iters_1L_int = iters_1L_int,
                                          mdl_types_lup = mdl_types_lup, seed_1L_int = seed_1L_int, prior_ls = prior_ls, control_ls = control_ls)
      Sys.sleep(5)
      smry_ls$smry_of_ts_mdl_tb
    })
  })
  saveRDS(mdls_smry_tb, paste0(mdl_smry_dir_1L_chr, "/mdls_smry_tb.RDS"))
  return(mdls_smry_tb)
}
write_ts_mdls_from_alg_outp <- function (outp_smry_ls, # rename lngl
                                         predictors_lup,
                                         depnt_var_min_val_1L_dbl = numeric(0),
                                         utl_min_val_1L_dbl = -1,
                                         backend_1L_chr = getOption("brms.backend", "rstan"),
                                         combinations_1L_lgl = F, 
                                         existing_predrs_ls = NULL, 
                                         iters_1L_int = 4000L,
                                         max_nbr_of_covars_1L_int = integer(0),
                                         new_dir_nm_1L_chr = "F_TS_Mdls",
                                         path_to_write_to_1L_chr = NA_character_,
                                         prior_ls = NULL,
                                         control_ls = NULL){
  if(is.na(path_to_write_to_1L_chr)) # BIN THIS AFTER TESTING
    path_to_write_to_1L_chr <- outp_smry_ls$path_to_write_to_1L_chr %>%
      stringr::str_sub(end=-8)
  output_dir_1L_chr <- write_new_outp_dir(path_to_write_to_1L_chr,
                                          new_dir_nm_1L_chr = new_dir_nm_1L_chr)
  outp_smry_ls$predr_vars_nms_ls <- make_predr_vars_nms_ls(main_predrs_chr = outp_smry_ls$predr_cmprsn_tb$predr_chr,
                                                           covars_ls = list(outp_smry_ls$prefd_covars_chr),
                                                           combinations_1L_lgl = combinations_1L_lgl,
                                                           existing_predrs_ls = existing_predrs_ls,
                                                           max_nbr_of_covars_1L_int = max_nbr_of_covars_1L_int)
  outp_smry_ls$mdl_nms_ls <- make_mdl_nms_ls(outp_smry_ls$predr_vars_nms_ls,
                                             mdl_types_chr = outp_smry_ls$prefd_mdl_types_chr)
  mdls_smry_tb <- write_ts_mdls(data_tb = outp_smry_ls$scored_data_tb,
                                depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                depnt_var_nm_1L_chr = outp_smry_ls$depnt_var_nm_1L_chr, predr_vars_nms_ls = outp_smry_ls$predr_vars_nms_ls,
                                id_var_nm_1L_chr = outp_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = outp_smry_ls$round_var_nm_1L_chr,
                                round_bl_val_1L_chr = outp_smry_ls$round_bl_val_1L_chr,
                                mdl_nms_ls = outp_smry_ls$mdl_nms_ls,
                                mdl_smry_dir_1L_chr = output_dir_1L_chr,
                                predictors_lup = predictors_lup,
                                utl_min_val_1L_dbl = utl_min_val_1L_dbl,
                                backend_1L_chr = backend_1L_chr, iters_1L_int = iters_1L_int, mdl_types_lup = outp_smry_ls$mdl_types_lup,
                                seed_1L_int = outp_smry_ls$seed_1L_int, prior_ls = prior_ls, control_ls = control_ls)
  outp_smry_ls$mdls_smry_tb <- mdls_smry_tb
  outp_smry_ls$utl_min_val_1L_dbl <- utl_min_val_1L_dbl
  outp_smry_ls$file_paths_chr <- list.files(outp_smry_ls$path_to_write_to_1L_chr, recursive = T)
  outp_smry_ls$combinations_1L_lgl <- combinations_1L_lgl
  outp_smry_ls$max_nbr_of_covars_1L_int <- max_nbr_of_covars_1L_int
  return(outp_smry_ls)
}
write_ts_mdl_plts <- function (brms_mdl, # Rename lngl
                               table_predn_mdl = NULL, tfd_data_tb, mdl_nm_1L_chr, path_to_write_to_1L_chr,
                               #depnt_var_min_val_1L_dbl = numeric(0),
                               depnt_var_nm_1L_chr = "utl_total_w", depnt_var_desc_1L_chr = "Utility score", 
                               predn_type_1L_chr = NULL, round_var_nm_1L_chr = "round", sd_dbl = NA_real_, sfx_1L_chr = " from table", tfmn_1L_chr = "NTF",
                               units_1L_chr = "in", height_dbl = c(rep(6, 2), rep(5,8)), width_dbl = c(rep(6, 2), rep(6, 8)),
                               rsl_dbl = rep(300,10), args_ls = NULL, seed_1L_dbl = 23456, utl_min_val_1L_dbl = -1)
{
  set.seed(seed_1L_dbl)
  tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb,
                                                  model_mdl = brms_mdl,
                                                  depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                  is_brms_mdl_1L_lgl = inherits(brms_mdl,"brmsfit"),
                                                  predn_type_1L_chr = predn_type_1L_chr,
                                                  sd_dbl = NA_real_,
                                                  sfx_1L_chr = ifelse(inherits(brms_mdl,"brmsfit"),
                                                                      " from brmsfit",
                                                                      sfx_1L_chr),
                                                  tfmn_1L_chr = tfmn_1L_chr,
                                                  utl_min_val_1L_dbl = utl_min_val_1L_dbl)
  if(!is.null(table_predn_mdl)){
    tfd_data_tb <- transform_ds_for_all_cmprsn_plts(tfd_data_tb = tfd_data_tb,
                                                    model_mdl = table_predn_mdl,
                                                    depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                    is_brms_mdl_1L_lgl = F,
                                                    predn_type_1L_chr = predn_type_1L_chr,
                                                    sd_dbl = sd_dbl,
                                                    sfx_1L_chr = ifelse(is.null(brms_mdl),
                                                                        " from table",
                                                                        sfx_1L_chr),
                                                    tfmn_1L_chr = tfmn_1L_chr,
                                                    utl_min_val_1L_dbl = utl_min_val_1L_dbl)
  }
  plt_nms_chr <- paste0(mdl_nm_1L_chr, "_",
                        c("coefs", "hetg",
                          "dnst", "sctr_plt",
                          "sim_dnst", "sim_sctr",
                          "cnstrd_dnst","cnstrd_sctr_plt",
                          "cnstrd_sim_dnst", "cnstrd_sim_sctr"))
  mdl_plts_paths_ls <- purrr::map(ifelse(inherits(brms_mdl,"brmsfit"),1,3):10, ~{
    plt_fn <- fn_args_ls <- NULL
    if (.x %in% c(1, 2)) {
      plt <- plot(brms_mdl, ask = F, plot = F)
      if (length(plt) >= .x) {
        fn_args_ls <- list(brms_mdl = brms_mdl, idx_1L_int = as.integer(.x))
        plt_fn <- function(brms_mdl, idx_1L_int) {
          plot(brms_mdl, ask = F, plot = F)[idx_1L_int]
        }
      }
    }  else {
      plot_fn_and_args_ls <- make_plot_fn_and_args_ls(tfd_data_tb = tfd_data_tb,
                                                      args_ls = args_ls,
                                                      depnt_var_nm_1L_chr = depnt_var_nm_1L_chr,
                                                      depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                                      round_var_nm_1L_chr = round_var_nm_1L_chr,
                                                      sfx_1L_chr = ifelse(is.null(table_predn_mdl),
                                                                          ifelse(inherits(brms_mdl,"brmsfit"),
                                                                                 " from brmsfit",
                                                                                 sfx_1L_chr),
                                                                          ifelse(is.null(brms_mdl),
                                                                                 " from table",
                                                                                 sfx_1L_chr)),
                                                      table_predn_mdl = table_predn_mdl,
                                                      tfmn_1L_chr = tfmn_1L_chr,
                                                      type_1L_chr = c("coefs", "hetg",
                                                                      "dnst", "sctr_plt",
                                                                      "sim_dnst", "sim_sctr",
                                                                      "cnstrd_dnst","cnstrd_sctr_plt",
                                                                      "cnstrd_sim_dnst", "cnstrd_sim_sctr")[.x],
                                                      brms_mdl = NULL, # This is correct
                                                      predn_type_1L_chr =  predn_type_1L_chr,
                                                      sd_dbl = sd_dbl,
                                                      seed_1L_dbl = seed_1L_dbl)
      plt_fn <- plot_fn_and_args_ls$plt_fn
      fn_args_ls <- plot_fn_and_args_ls$fn_args_ls
    }
    ready4show::write_mdl_plt_fl(plt_fn,
                                 fn_args_ls = fn_args_ls,
                                 path_to_write_to_1L_chr = path_to_write_to_1L_chr,
                                 plt_nm_1L_chr = plt_nms_chr[.x],
                                 units_1L_chr = units_1L_chr,
                                 width_1L_dbl = width_dbl[.x],
                                 height_1L_dbl = height_dbl[.x],
                                 rsl_1L_dbl = rsl_dbl[.x])
  }) %>%
    stats::setNames(plt_nms_chr[ifelse(inherits(brms_mdl,"brmsfit"),1,3):10]) %>%
    purrr::discard(is.na)
  return(mdl_plts_paths_ls)
}
## For TTU
author_TTUReports <- function(x,
                              depnt_var_desc_1L_chr = NA_character_,
                              depnt_var_min_val_1L_dbl = numeric(0),
                              download_tmpl_1L_lgl = T,
                              fl_type_1L_chr = ".eps",
                              timepoint_new_nms_chr = NA_character_,
                              type_1L_chr = "Report",
                              what_1L_chr = NA_character_,
                              ...){
  if(type_1L_chr == "Report"){
    if(download_tmpl_1L_lgl){
      authorData(x@a_TTUSynopsis,
                 tmpl_url_1L_chr = ifelse(what_1L_chr == "Catalogue",
                                          x@catalogue_tmpl_chr[1],
                                          x@manuscript_tmpl_chr[1]),
                 tmpl_version_1_L_chr = ifelse(what_1L_chr == "Catalogue",
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
             type_1L_chr = "Report",
             what_1L_chr = what_1L_chr)
    }else{
      author_SpecificSynopsis(x@a_TTUSynopsis, # authorReport when Exported
                   what_1L_chr = what_1L_chr,
                   ...)
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
                           utils::packageDescription(.x) %>%
                             `[`(c("Depends", "Imports")) %>%
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
      saveRDS(df,
              paste0(dir_1L_chr,
                     "/packages.RDS"))
      
    }
    if(type_1L_chr == "Plots"){
      composite_1_plt <- depict_SpecificSynopsis(x@a_TTUSynopsis,#depictSlot(x,"a_TTUSynopsis", when exported
                                    depnt_var_desc_1L_chr = depnt_var_desc_1L_chr,
                                    depnt_var_min_val_1L_dbl = depnt_var_min_val_1L_dbl,
                                    timepoint_old_nms_chr = procureSlot(x,
                                                                        "a_TTUSynopsis@d_YouthvarsProfile@timepoint_vals_chr"),
                                    timepoint_new_nms_chr = timepoint_new_nms_chr,
                                    what_1L_chr = "composite_mdl",
                                    write_1L_lgl = T)
      composite_2_plt <- depict_SpecificSynopsis(x@a_TTUSynopsis,#depictSlot(x,"a_TTUSynopsis", when exported
                                                 what_1L_chr = "composite_utl",
                                    write_1L_lgl = T)
      if(!is.na(what_1L_chr)){
        ggplot2::ggsave(file = paste0(dir_1L_chr,
                                      "/fig1",
                                      fl_type_1L_chr),
                        composite_2_plt)
        ggplot2::ggsave(file = paste0(dir_1L_chr,
                                      "/fig2",
                                      fl_type_1L_chr),
                        composite_1_plt)
        
      }
      
    }
  }
}
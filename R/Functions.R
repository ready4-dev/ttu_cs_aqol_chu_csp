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
                              y_label_1L_chr = "Percentage"){
  if(!is.null(col_nms_chr)){
    plots_ls<-list()
    for(i in col_nms_chr){
      if(make_log_log_tfmn_1L_lgl){
        targetvar = paste0("tran_",i)
        data_tb <- dplyr::mutate(data_tb, !!targetvar := log(-log(1-!!as.name(i))))  %>%
          dplyr::mutate(!!targetvar :=ifelse(!!as.name(i)==1,log(-log(1-0.999)),!!as.name(targetvar)))
      }
      labelx <- eval(parse(text=paste0("attributes(data_tb$",i,")$label")))
      labelx <- stringr::str_sub(labelx,
                                 start = stringi::stri_locate_last_fixed(labelx," - ")[1,1] %>%
                                   unname() + 2)
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
      heights_int <- heights_int[-length(heights_int)]
      composite_plt <- gridExtra::grid.arrange(ggpubr::ggarrange(plotlist=plots_ls,
                                                                 nrow = plot_rows_cols_pair_int[1],
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
                                                             base_height_dbl = 10)
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
                                                                heights_int = dim_plots_params_ls$heights_int),
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
                                                              make_log_log_tfmn_1L_lgl = T),
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
      descv_plts_paths_ls <- write_descv_plots(x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb, #add youthvars:: when Exporting
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
investigate_SpecificModels <- function(x,
                                       depnt_var_max_val_1L_dbl = Inf,
                                       session_ls = NULL){
  results_ls <- write_mdl_cmprsn(scored_data_tb = x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                 ds_smry_ls = manufacture_SpecificProject(x, # Replace with manufacture when Exporting
                                                          what_1L_chr = "ds_smry_ls"),
                                 mdl_smry_ls = manufacture_SpecificProject(x, # Replace with manufacture when Exporting
                                                           what_1L_chr = "mdl_smry_ls"),
                                 output_data_dir_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                                 depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl,
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
  # }
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
write_mdl_cmprsn <- function(scored_data_tb,
                             ds_smry_ls,
                             mdl_smry_ls,
                             output_data_dir_1L_chr,
                             depnt_var_max_val_1L_dbl = 0.99,
                             seed_1L_int = 1234){
  bl_tb <- transform_ds_for_tstng(scored_data_tb,# youthvars::
                                  depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                  candidate_predrs_chr = ds_smry_ls$candidate_predrs_chr,
                                  depnt_var_max_val_1L_dbl = depnt_var_max_val_1L_dbl,
                                  # round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                                  # round_val_1L_chr = ds_smry_ls$round_bl_val_1L_chr,
                                  round_var_nm_1L_chr = if(identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_descvs_ls$round_var_nm_1L_chr},
                                  round_val_1L_chr = if(identical(ds_descvs_ls$round_var_nm_1L_chr,character(0))){NA_character_}else{ds_smry_ls$round_bl_val_1L_chr}
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
investigate_SpecificMixed <- function(x,
                                      backend_1L_chr = "cmdstanr",
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
    input_params_ls <- manufacture(x,
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
                                              new_dir_nm_1L_chr = new_dir_nm_1L_chr) %>%
      stats::setNames(names(scndry_anlys_params_ls))
    x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls <- append(x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls,
                                                                        results_ls_ls %>%
                                                                          purrr::map(~.x[-1]))
  }
  return(x)
}
investigate_SpecificPredictors <- function(x){
  results_ls <- write_predr_and_covars_cmprsn(scored_data_tb = x@a_YouthvarsProfile@a_Ready4useDyad@ds_tb,
                                              bl_tb = x@c_SpecificResults@b_SpecificPrivate@private_outp_ls$bl_tb,
                                              ds_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$ds_smry_ls,
                                              mdl_smry_ls = x@c_SpecificResults@a_SpecificShareable@shareable_outp_ls$mdl_smry_ls,
                                              output_data_dir_1L_chr = x@b_SpecificParameters@paths_ls$output_data_dir_1L_chr,
                                              seed_1L_int = x@b_SpecificParameters@seed_1L_int)
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
write_predr_and_covars_cmprsn <- function(scored_data_tb,
                                          bl_tb,
                                          ds_smry_ls,
                                          mdl_smry_ls,
                                          output_data_dir_1L_chr,
                                          seed_1L_int = 1234){
  mdl_smry_ls$predr_cmprsn_tb <- write_predr_cmprsn_outps(data_tb = bl_tb,
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
                                                                     depnt_var_nm_1L_chr = ds_smry_ls$depnt_var_nm_1L_chr,
                                                                     predrs_var_nms_chr = ds_smry_ls$candidate_predrs_chr,
                                                                     covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr,
                                                                     mdl_type_1L_chr = mdl_smry_ls$prefd_mdl_types_chr[1],
                                                                     path_to_write_to_1L_chr = output_data_dir_1L_chr,
                                                                     new_dir_nm_1L_chr = "D_Predr_Covars_Cmprsn",
                                                                     fl_nm_pfx_1L_chr = "D_CT",
                                                                     mdl_types_lup = mdl_smry_ls$mdl_types_lup)
  mdl_smry_ls$signt_covars_chr <- get_signft_covars(mdls_with_covars_smry_tb = mdl_smry_ls$mdls_with_covars_smry_tb,
                                                    covar_var_nms_chr = ds_smry_ls$candidate_covar_nms_chr)
  predr_and_covars_cmprsn_ls <- list(bl_tb = bl_tb,
                                     ds_smry_ls = ds_smry_ls,
                                     mdl_smry_ls = mdl_smry_ls)
  return(predr_and_covars_cmprsn_ls)
}
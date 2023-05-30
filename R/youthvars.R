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
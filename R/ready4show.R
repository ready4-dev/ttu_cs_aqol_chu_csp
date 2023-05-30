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
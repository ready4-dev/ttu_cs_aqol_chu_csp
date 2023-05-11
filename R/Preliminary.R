source("R/Functions.R")
ymh_clinical_aqol_chu_dict_r3 <- readRDS(path_to_data_1L_chr) %>% # Path value needs to be assigned to path_to_data_1L_chr
  transform_csnl_example_ds() %>%
  make_csnl_example_dict()
predictors_r3 <- make_csnl_example_predrs()
Y <- Ready4useRepos(dv_nm_1L_chr = "TTU", # Replace with values for a dataverse & dataset for which
                    dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/FDRUXH", #  you have write permissions.
                    dv_server_1L_chr = "dataverse.harvard.edu")

Y <- share(Y,
           obj_to_share_xx = ymh_clinical_aqol_chu_dict_r3,
           fl_nm_1L_chr = "ymh_clinical_aqol_chu_dict_r3",
           description_1L_chr = "Data dictionary for the pre-analysis study dataset.")
Y <- share(Y,
           obj_to_share_xx = predictors_r3,
           fl_nm_1L_chr = "predictors_r3",
           description_1L_chr = "Candiate predictors for utility mapping models.")
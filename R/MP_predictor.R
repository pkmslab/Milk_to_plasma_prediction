#' MP predictor function
#'
#' Predict the MP using information from CHEMBL
#'
#' @description The function calculates the milk-to-plasma ratio (MP) using information pulled from the CHEMBL database API. It will acquire
#' the data required to predict the phase distribution equation, the equation proposed by Koshimichi et al., the logarithmically transformed
#' phase distribution equation, Meskin and Lien's equation, and finally the model trained using the XGBoost ML method. The relevant equations
#' can be found online or by viewing the equations nested within the function. Drug names or "custom" are all that is required to use the function
#' and the user is able to modify any of the input parameters as they see fit. The goal of this function is to facilitate MP predictions
#' and to provide the user with the choice of MP to use in PBPK models.

#' @param drug_names The name of the drug of interest or input custom to build the dataset from scratch
#' @param ID_url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search" Non-editable URL, required to pull CHEMBL IDs
#' @param base_url = "https://www.ebi.ac.uk/chembl/api/data/activity" Non-editable URL, required to pull unbound fractions and activity data
#' @param limit = 10000 Required to scan through the CHEMBL API pages

#' @return A data frame with MP predictions from the phase distribution equation, Koshimichi et al. equation, logarithmically transformed PD equation, Meskin and Lien equation, and
#' XGBoost model
#'
#' @details Enter the name of the drug of interest. The first step is to go through CHEMBL
#' and pull data from databases called resources. An entry may or may not have a record in a
#' particular resource.
#'
#' @details If you are unsure about which name to use, you could go onto the CHEMBL website: https://www.ebi.ac.uk/chembl/
#' and find the drug of interest. Then enter the name into the function and the function should work as intended.
#'
#' @import tidyverse
#' @importFrom stats predict
#' @importFrom dplyr %>%
#'
#' @note
#' Links to the webservice documentation:
#' \itemize{
#'      \item \url{https://www.ebi.ac.uk/chembl/api/data/docs}
#' }
#' @export

#### MP predictions ####

# Function to get ChEMBL data for multiple drug names
MP_predictor <- function(drug_names, ID_url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search",
                                   base_url = "https://www.ebi.ac.uk/chembl/api/data/activity", limit = 10000) {

  for (drug_name in drug_names) {
    if(!is.character(drug_name)) {
      return(paste0("Sorry, you need to ensure your input is a character. Your input is ", class(drug_name)))
      }
    }

  # Initialize data frames to store results
  ID_all_results <- data.frame()
  Fu_all_results <- data.frame()

  # Cycle through the drugs in the list
  for (drug_name in drug_names) {
    if(drug_name == stringr::regex("custom", ignore_case = TRUE)) {

      # Define variables with arbitrary inputs
      Drugs <- "Input drug name"
      Type <- "Input drug type"
      MW <- 0
      pka1 <- 0
      pka2 <- 0
      PSA <- 0
      HBD <- 0
      LogP <- 0
      LogD7.4 <- 0
      fup <- 0

      # Create a data frame
      drug_info <- data.frame(Drugs, Type, MW, pka1, pka2, PSA, HBD, LogP, LogD7.4, fup)
      drug_names_chembl <- NULL
      chembl_id <- NULL
    }
    else{
      drug_names_chembl <- drug_names
    }
  }

  tictoc::tic("Elapsed time")
  # Loop through each drug name in the list
  for(drug_name in drug_names_chembl) {

    # For pagination
    offset <- 0
    repeat {
      # Construct the molecule search URL
      url <- paste0(ID_url ,"?q=", utils::URLencode(drug_name), "&limit=", limit, "&offset=", offset)

      req <- httr2::request(url)
      resp <- req |> httr2::req_headers(Accept = "application/json")
      response <- httr2::req_perform(resp)
      content <- httr2::resp_body_string(response)
      parsed_result <- jsonlite::fromJSON(content, flatten = TRUE)

      # If no data is returned, break the loop
      if (is.null(parsed_result$molecules) || length(parsed_result$molecules) == 0) {
        message(paste0("Drug name search: ", drug_name))
        break
      }

      # Convert to data frame if necessary
      parsed_data <- as.data.frame(parsed_result$molecules)

      # Extract molecule data
      molecule_data <- parsed_data %>%
        dplyr::select(pref_name, molecule_chembl_id, molecule_properties.full_mwt,
               molecule_properties.cx_most_apka, molecule_properties.cx_most_bpka,
               molecule_properties.psa, molecule_properties.hbd, molecule_properties.cx_logp,
               molecule_properties.cx_logd, molecule_properties.molecular_species)

      # Add Drug Name
      molecule_data <- dplyr::mutate(molecule_data, Drug_Name = drug_name)

      # Combine results
      ID_all_results <- dplyr::bind_rows(ID_all_results, molecule_data)

      # Update the offset for pagination
      offset <- offset + limit
    }

    # Filter the results case-insensitively
    ID_all_results <- ID_all_results %>%
      dplyr::filter(tolower(Drug_Name) == tolower(pref_name))
  }

  # Loop through each drug name in the list
  for(chembl_id in ID_all_results$molecule_chembl_id) {

    # For pagination
    offset <- 0
    repeat {

      # Construct the activity search URL
      url <- paste0(base_url, "?molecule_chembl_id=", chembl_id, "&limit=", limit, "&offset=", offset)

      req <- httr2::request(url)
      resp <- req |> httr2::req_headers(Accept = "application/json")
      response <- httr2::req_perform(resp)
      content <- httr2::resp_body_string(response)
      parsed_result <- jsonlite::fromJSON(content, flatten = TRUE)

      # If no data is returned, break the loop
      if (is.null(parsed_result$activities) || length(parsed_result$activities) == 0) {
        message(paste0("Scanning the database for: ", chembl_id))
        break
      }

      # Convert to data frame if necessary
      parsed_data <- as.data.frame(parsed_result$activities)

      # Extract activity data
      activity_data <- parsed_data %>%
        dplyr::select(molecule_chembl_id, standard_type, standard_value, units, assay_description)

      # Add the ChEMBL ID for reference
      activity_data <- dplyr::mutate(activity_data, Entry = chembl_id)

      # Combine results
      Fu_all_results <- dplyr::bind_rows(Fu_all_results, activity_data)

      # Update the offset for pagination
      offset <- offset + limit
    }

    Fu_all_results <- Fu_all_results %>%
      dplyr::filter(tolower(Entry) == tolower(molecule_chembl_id), standard_type == "Fu") %>%
      dplyr::mutate(
        standard_value = ifelse(is.na(standard_value), NA,
                                ifelse(standard_value >= 0, standard_value, NA))
      )


    # Combine molecule and activity data
    drug_info <- dplyr::left_join(ID_all_results, Fu_all_results, by = c("molecule_chembl_id" = "Entry")) %>%
      dplyr::mutate(Drugs = pref_name,
                   CHEMBL = molecule_chembl_id,
                   Type = molecule_properties.molecular_species,
                   MW = as.numeric(molecule_properties.full_mwt),
                   pka1 = as.numeric(molecule_properties.cx_most_apka),
                   pka2 = as.numeric(molecule_properties.cx_most_bpka),
                   PSA = as.numeric(molecule_properties.psa),
                   HBD = as.numeric(molecule_properties.hbd),
                   LogP = as.numeric(molecule_properties.cx_logp),
                   LogD7.4 = as.numeric(molecule_properties.cx_logd)) %>%
      dplyr::mutate(pka1 = dplyr::coalesce(pka1, pka2),
             pka2 = ifelse(pka2 == pka1, NA, pka2)) %>%
      dplyr::distinct(pka1, pka2, .keep_all = TRUE) %>%
      dplyr::mutate(fup = as.numeric(standard_value))

    # Edit the data before running calculations/predictions
    drug_info <- drug_info %>%
      dplyr::select(Drugs, Type, MW, pka1, pka2, PSA, HBD, LogP, LogD7.4, fup)
  }

  message(paste0("Confirm/edit the data. When complete, close the editor."))

  # Edit the data
  prelim_drug_data <- utils::edit(drug_info, editor = "xedit")

  # Required model parameters
  drug_data <- prelim_drug_data %>%
    dplyr::mutate(
        funp = case_when(
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 1,
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(7.4-pka1)),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(pka1-7.4)),
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 1,
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.4 - pka1)+10^(7.4 - pka2)+10^(2*7.4 - (pka1+pka2))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(pka1-7.4)+10^(pka2-7.4)+10^((pka1+pka2)-2*7.4)),
          stringr::str_detect(Type, stringr::regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.4-pka1)+10^(pka2-7.4)+10^(pka2-pka1)),
          TRUE ~ NA_real_  # Default case for unmatched Types
        ),
        funsm = case_when(
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 1,
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(7.2-pka1)),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(pka1-7.2)),
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 1,
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.2 - pka1)+10^(7.2 - pka2)+10^(2*7.2 - (pka1+pka2))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(pka1-7.2)+10^(pka2-7.2)+10^((pka1+pka2)-2*7.2)),
          stringr::str_detect(Type, stringr::regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.2-pka1)+10^(pka2-7.2)+10^(pka2-pka1)),
          TRUE ~ NA_real_
        ),
        MuPu = funp/funsm
    )

  ## Phase distribution
  PD <- drug_data %>%
    dplyr::mutate(
        fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
        Dmilk7.2 = case_when(
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          TRUE ~ NA_real_
        ),
        Phase_Distribution_MP = (fup*MuPu*(0.955/fusm + 0.045*Dmilk7.2))) %>%
    dplyr::select(Phase_Distribution_MP)

  ## Koshimichi et al.
  Koshi <- drug_data %>%
    dplyr::mutate(
        fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
        Dmilk7.2 = case_when(
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
          stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          stringr::str_detect(Type, stringr::regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
          TRUE ~ NA_real_
        ),
        fm_total = (1/(0.955/fusm + 0.045*Dmilk7.2)),
        CLre = 10^(2.793 + 0.179*LogP - 0.132*HBD),
        CLsec = 10^(-3.912 - 0.015*PSA + 3.367*log10(MW) - 0.164*log10((10^LogP)/(10^LogD7.4))),
        Koshi_MP = ((CLsec/CLre)*(fup/fm_total))) %>%
    dplyr::select(Koshi_MP)

  ## Log Phase distribution
  log_PD <- drug_data %>%
    dplyr::mutate(
          fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
          Dmilk7.2 = case_when(
            stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
            stringr::str_detect(Type, stringr::regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
            stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
            stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
            stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
            stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
            stringr::str_detect(Type, stringr::regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
            TRUE ~ NA_real_
          ),
          log_Phase_Distribution_MP = case_when(
            stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) ~ exp(-0.405 + 9.36*log(MuPu) - 0.69*log(fup) - 1.54*log(0.955/fusm + 0.045*Dmilk7.2)),
            stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) ~ exp(0.03 + 2.28*log(MuPu) + 0.89*log(fup) + 0.51*log(0.955/fusm + 0.045*Dmilk7.2)),
            TRUE ~ NA_real_
          )) %>%
    dplyr::select(log_Phase_Distribution_MP)

  ## Meskin and Lien
  M_and_L <- drug_data %>%
    dplyr::mutate(
      M_and_L_MP = case_when(
        stringr::str_detect(Type, stringr::regex("ACID", ignore_case = TRUE)) ~ 10^(2.068-0.162*sqrt(MW)-0.185*LogP),
        stringr::str_detect(Type, stringr::regex("BASE", ignore_case = TRUE)) ~ 10^(0.265-0.153*LogP-0.128*(7.4-pka1)),
        TRUE ~ NA_real_
      )) %>%
    dplyr::select(M_and_L_MP)

  ## XGBoost
  # Select variables
  XGB_data <- drug_data %>%
    dplyr::mutate(pKa1 = pka1) %>%
    dplyr::select(pKa1, fup, PSA, LogP, LogD7.4, funp, MuPu)

  # Store XGBoost variables as a matrix
  test_matrix <- xgboost::xgb.DMatrix(data = as.matrix(XGB_data))

  ### NEED TO FIX XGBOOST
  xgbmodel <- system.file("extdata", "xgboost_model.rds", package = "MP.prediction")
  bstDMatrix <- readRDS(xgbmodel)

  # Make prediction using existing model
  pred_xgboost <- predict(bstDMatrix, test_matrix)
  xgb_pred <- data.frame(MP_ratio = (10^pred_xgboost))

tictoc::toc()

  # Report values and corresponding drugs
  MP_values <- data.frame(Drug = drug_data$Drugs,
                          Type = drug_data$Type,
                          Phase_Distribution_MP = PD$Phase_Distribution_MP,
                          Koshimichi_et_al_MP = Koshi$Koshi_MP,
                          Log_Phase_Distribution_MP = log_PD$log_Phase_Distribution_MP,
                          Meskin_and_Lien_MP = M_and_L$M_and_L_MP,
                          XGBoost_MP = xgb_pred$MP_ratio)

  return(MP_values)
}





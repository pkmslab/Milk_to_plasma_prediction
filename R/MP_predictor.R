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
#' @details Enter the name of the drug of interest. The first step is to go through the ChEMBL
#' and pull data from databases called resources. An entry may or may not have a record in a
#' particular resource.
#'
#' @details If you are unsurea about which name to use, you could go onto the CHEMBL website and find the drug of
#' interest. Then enter the name into the function and the function should work as intended.
#'
#' @note
#' Links to the webservice documentation:
#' \itemize{
#'      \item \url{https://www.ebi.ac.uk/chembl/api/data/docs}
#' }
#' @export

# Function to predict the MP ratio of drugs using a variety of methods
MP_prediction_function <- function(drug_names, ID_url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search",
                                   base_url = "https://www.ebi.ac.uk/chembl/api/data/activity", limit = 10000) {

  # Initialize data frames to store results
  ID_all_results <- data.frame()
  Fu_all_results <- data.frame()

  # Cycle through the drugs in the list
  for (drug_name in drug_names) {
    if(drug_name == regex("custom", ignore_case = TRUE)) {

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

  tic("Molecular Properties")

  # Loop through each drug name in the list
  for(drug_name in drug_names_chembl) {

    # For pagination
    offset <- 0
    repeat {

      # Construct the molecule search URL
      url <- paste0(ID_url, "?q=", URLencode(drug_name), "&limit=", limit, "&offset=", offset)

      # Make the GET request
      response <- GET(url, add_headers(Accept = "application/json"))

      if (status_code(response) == 200) {
        json_data <- content(response, "text", encoding = "UTF-8")
        parsed_data <- fromJSON(json_data, flatten = TRUE)

        # If no data is returned, break the loop
        if (is.null(parsed_data$molecules) || length(parsed_data$molecules) == 0) {
          message(drug_name)
          break
        }

        # Extract molecule data
        molecule_data <- parsed_data$molecules %>%
          select(pref_name, molecule_chembl_id, molecule_properties.full_mwt,
                 molecule_properties.cx_most_apka, molecule_properties.cx_most_bpka,
                 molecule_properties.psa, molecule_properties.hbd, molecule_properties.cx_logp,
                 molecule_properties.cx_logd, molecule_properties.molecular_species)

        # Add Drug Name
        molecule_data <- mutate(molecule_data, Drug_Name = drug_name)

        # Combine results
        ID_all_results <- bind_rows(ID_all_results, molecule_data)

        # Update the offset for pagination
        offset <- offset + limit
      }
      else {
        message("Failed to retrieve data for ", drug_name, ". Status Code: ", status_code(response))
        message("Response: ", content(response, "text"))
        break
      }

      # Filter the results case-insensitively
      ID_all_results <- ID_all_results %>%
        filter(tolower(Drug_Name) == tolower(pref_name))
    }
  }
  toc()

  tic("Activity data")

  for (chembl_id in ID_all_results$molecule_chembl_id) {

    # For pagination
    offset <- 0
    repeat{

      # Construct the activity search URL
      url <- paste0(base_url, "?molecule_chembl_id=", chembl_id, "&limit=", limit, "&offset=", offset)

      # Make the GET request
      response <- GET(url, add_headers(Accept = "application/json"))

      if (status_code(response) == 200) {
        json_data <- content(response, "text", encoding = "UTF-8")
        parsed_data <- fromJSON(json_data, flatten = TRUE)

        # If no data is returned, break the loop
        if (is.null(parsed_data$activities) || length(parsed_data$activities) == 0) {
          message(chembl_id)
          break
        }

        # Extract activity data
        activity_data <- parsed_data$activities %>%
          select(molecule_chembl_id, standard_type, standard_value, units, assay_description)

        # Add the ChEMBL ID for reference
        activity_data <- mutate(activity_data, Entry = chembl_id)

        # Combine results
        Fu_all_results <- bind_rows(Fu_all_results, activity_data)

        # Update the offset for pagination
        offset <- offset + limit
      }
      else {
        message("Failed to retrieve data for ", chembl_id, ". Status Code: ", status_code(response))
        message("Response: ", content(response, "text"))
        break
      }
    }

    # Filter the activity results for specific types (e.g., "Fu")
    Fu_all_results <- Fu_all_results %>%
      filter(tolower(Entry) == tolower(molecule_chembl_id), standard_type == "Fu")

    # Combine molecule and activity data
    drug_info <- left_join(ID_all_results, Fu_all_results, by = c("molecule_chembl_id" = "Entry")) %>%
      mutate(Drugs = pref_name,
             CHEMBL = molecule_chembl_id,
             Type = molecule_properties.molecular_species,
             MW = as.numeric(molecule_properties.full_mwt),
             pka1 = as.numeric(molecule_properties.cx_most_apka),
             pka2 = as.numeric(molecule_properties.cx_most_bpka),
             PSA = as.numeric(molecule_properties.psa),
             HBD = as.numeric(molecule_properties.hbd),
             LogP = as.numeric(molecule_properties.cx_logp),
             LogD7.4 = as.numeric(molecule_properties.cx_logd)) %>%
      filter(str_detect(assay_description, regex("human|plasma", ignore_case = TRUE))) %>%
      mutate(pka1 = coalesce(pka1, pka2),
             pka2 = ifelse(pka2 == pka1, NA, pka2)) %>%
      distinct(pka1, pka2, .keep_all = TRUE) %>%
      mutate(fup = as.numeric(standard_value))

    # Edit the data before running calculations/predictions
    drug_info <- drug_info %>%
      select(Drugs, Type, MW, pka1, pka2, PSA, HBD, LogP, LogD7.4, fup)
  }
  toc()

  # Edit the data
  prelim_drug_data <- edit(drug_info, editor = "xedit")

  tic("MP Predictions")

  # Required model parameters
  drug_data <- prelim_drug_data %>%
    mutate(
      funp = case_when(
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 1,
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(7.4-pka1)),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(pka1-7.4)),
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 1,
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.4 - pka1)+10^(7.4 - pka2)+10^(2*7.4 - (pka1+pka2))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(pka1-7.4)+10^(pka2-7.4)+10^((pka1+pka2)-2*7.4)),
        str_detect(Type, regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.4-pka1)+10^(pka2-7.4)+10^(pka2-pka1)),
        TRUE ~ NA_real_  # Default case for unmatched Types
      ),
      funsm = case_when(
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 1,
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(7.2-pka1)),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 1/(1+10^(pka1-7.2)),
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 1,
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.2 - pka1)+10^(7.2 - pka2)+10^(2*7.2 - (pka1+pka2))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(pka1-7.2)+10^(pka2-7.2)+10^((pka1+pka2)-2*7.2)),
        str_detect(Type, regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 1/(1+10^(7.2-pka1)+10^(pka2-7.2)+10^(pka2-pka1)),
        TRUE ~ NA_real_
      ),
      MuPu = funp/funsm
    )

  ## Phase distribution
  PD <- drug_data %>%
    mutate(
      fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
      Dmilk7.2 = case_when(
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        TRUE ~ NA_real_
      ),
      Phase_Distribution_MP = (fup*MuPu*(0.955/fusm + 0.045*Dmilk7.2))) %>%
    select(Phase_Distribution_MP)

  ## Koshimichi et al.
  Koshi <- drug_data %>%
    mutate(
      fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
      Dmilk7.2 = case_when(
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        TRUE ~ NA_real_
      ),
      fm_total = (1/(0.955/fusm + 0.045*Dmilk7.2)),
      CLre = 10^(2.793 + 0.179*LogP - 0.132*HBD),
      CLsec = 10^(-3.912 - 0.015*PSA + 3.367*log10(MW) - 0.164*log10((10^LogP)/(10^LogD7.4))),
      Koshi_MP = ((CLsec/CLre)*(fup/fm_total))) %>%
    select(Koshi_MP)

  ## Log Phase distribution
  log_PD <- drug_data %>%
    mutate(
      fusm = ((fup^0.448)/(0.000694^0.448 + fup^0.448)),
      Dmilk7.2 = case_when(
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("NEUTRAL", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*LogD7.4),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("ACID", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        str_detect(Type, regex("AMPHOLYTE|ZWITTERION", ignore_case = TRUE)) & !is.na(pka2) ~ 10^(-0.88 + 1.29*(log10((10^LogD7.4)*funsm/funp))),
        TRUE ~ NA_real_
      ),
      log_Phase_Distribution_MP = case_when(
        str_detect(Type, regex("ACID", ignore_case = TRUE)) ~ exp(-0.405 + 9.36*log(MuPu) - 0.69*log(fup) - 1.54*log(0.955/fusm + 0.045*Dmilk7.2)),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) ~ exp(0.03 + 2.28*log(MuPu) + 0.89*log(fup) + 0.51*log(0.955/fusm + 0.045*Dmilk7.2)),
        TRUE ~ NA_real_
      )) %>%
    select(log_Phase_Distribution_MP)

  ## Meskin and Lien
  M_and_L <- drug_data %>%
    mutate(
      M_and_L_MP = case_when(
        str_detect(Type, regex("ACID", ignore_case = TRUE)) ~ 10^(2.068-0.162*sqrt(MW)-0.185*LogP),
        str_detect(Type, regex("BASE", ignore_case = TRUE)) ~ 10^(0.265-0.153*LogP-0.128*(7.4-pka1)),
        TRUE ~ NA_real_
      )) %>%
    select(M_and_L_MP)

  ## XGBoost
  # Select variables
  XGB_data <- drug_data %>%
    mutate(pKa1 = pka1) %>%
    select(pKa1, fup, PSA, LogP, LogD7.4, funp, MuPu)

  # Store XGBoost variables as a matrix
  test_matrix <- xgb.DMatrix(data = as.matrix(XGB_data))

  # Make prediction using existing model
  pred_xgboost <- predict(bstDMatrix, test_matrix)
  xgb_pred <- data.frame(MP_ratio = (10^pred_xgboost))

  toc()

  # Report values and corresponding drugs
  MP_values <- data.frame(Drug = drug_data$Drugs,
                          Type = drug_data$Type,
                          Phase_Distribution_MP = PD,
                          Koshimichi_et_al_MP = Koshi,
                          log_Phase_Distribution_MP = log_PD,
                          M_and_L_MP = M_and_L,
                          XGBoost_MP = xgb_pred$MP_ratio)

  return(MP_values)
}



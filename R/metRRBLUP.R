#' @title Multi-Trial Analysis using RR-BLUP (Ridge Regression BLUP)
#'
#' @description Fits genomic prediction models that estimate marker effects directly
#' using Ridge Regression BLUP via the LMMsolver package. Supports four model types:
#' Main Effects (A), Main Effects (A+D), GCA, and SCA/GCA. Operates in two modes:
#' fitting models from phenotypic data ("fit") or predicting unphenotyped individuals
#' using pre-estimated marker effects ("predict").
#'
#' @param phenoDTfile The Bioflow data object containing predictions, metrics, modeling,
#'   status, and data slots.
#' @param analysisId Character or numeric. The STA analysisId(s) to use as input data
#'   (fit mode) or ignored in predict mode.
#' @param analysisIdGeno Character or numeric. The analysisId for genotypic QA modifications.
#'   If NULL or empty, raw genotype data is used.
#' @param trait Character vector. Trait(s) to analyze.
#' @param useWeights Logical. Whether to use reliability-based weights in the model.
#'   Default TRUE.
#' @param heritLB Numeric. Lower bound for heritability filter (inclusive). Default 0.
#' @param heritUB Numeric. Upper bound for heritability filter (inclusive). Default 1.
#' @param meanLB Numeric. Lower bound for trait mean filter (inclusive). Default -Inf.
#' @param meanUB Numeric. Upper bound for trait mean filter (inclusive). Default Inf.
#' @param maxIters Numeric. Maximum number of iterations for LMMsolver. Default 50.
#' @param verbose Logical. Whether to emit progress and warning messages. Default TRUE.
#' @param model_type Character. One of "Main Effects (A)", "Main Effects (A+D)",
#'   "GCA", or "SCA/GCA". Default "Main Effects (A)".
#' @param mode Character. One of "fit" (fit a model from phenotypic data) or
#'   "predict" (predict unphenotyped individuals using stored marker effects).
#'   Default "fit".
#' @param markerEffectsId Character or numeric. The MTA_Stamp analysisId containing
#'   marker effect estimates to use in predict mode. Required when mode = "predict".
#'
#' @return The modified phenoDTfile object with new rows appended to the predictions,
#'   metrics, modeling, and status tables.
#'
#' @export
met_rrblup <- function(
    phenoDTfile = NULL,
    analysisId = NULL,
    analysisIdGeno = NULL,
    trait = NULL,
    useWeights = TRUE,
    heritLB = 0,
    heritUB = 1,
    meanLB = -Inf,
    meanUB = Inf,
    maxIters = 50,
    verbose = TRUE,
    model_type = c("Main Effects (A)", "Main Effects (A+D)", "GCA", "SCA/GCA"),
    mode = c("fit", "predict"),
    markerEffectsId = NULL
) {

  ## ───────────────────────────────────────────────────────────────────────────

## INPUT VALIDATION
## ───────────────────────────────────────────────────────────────────────────

  if (is.null(phenoDTfile)) {
    stop("Please provide the phenotype file (phenoDTfile is NULL).", call. = FALSE)
  }
  if (is.null(analysisId)) {
    stop("Please provide the STA analysisId to be used as input (analysisId is NULL).", call. = FALSE)
  }
  if (is.null(trait)) {
    stop("Please provide traits to be analyzed (trait is NULL).", call. = FALSE)
  }


  # Validate that specified traits exist in the predictions table for the given analysisId
  baseData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId), ]
  availableTraits <- unique(baseData[, "trait"])
  # --- TPP: resolve pheno_traits for validation (TPP trait names won't be in data directly) ---
  tpp_config <- phenoDTfile$metadata$tpp_analysis_config
  traitsForValidation <- trait
  if (!is.null(tpp_config) && !is.null(tpp_config$trait_map)) {
    traitsForValidation <- vapply(trait, function(t) {
      if (t %in% names(tpp_config$trait_map)) tpp_config$trait_map[[t]] else t
    }, character(1), USE.NAMES = FALSE)
  }
  missingTraits <- setdiff(traitsForValidation, availableTraits)
  if (length(intersect(traitsForValidation, availableTraits)) == 0) {
    stop(
      paste0(
        "The traits you have specified are not present in the analysisId provided. ",
        "Missing traits: ", paste(missingTraits, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  # Match arguments
  model_type <- match.arg(model_type)
  mode <- match.arg(mode)

  ## ───────────────────────────────────────────────────────────────────────────
## GENERATE UNIQUE ANALYSIS ID FOR THIS RUN
## ───────────────────────────────────────────────────────────────────────────

  mtaAnalysisId <- as.numeric(Sys.time())

  ## ───────────────────────────────────────────────────────────────────────────
## ENVIRONMENT FILTERING (Task 1.2)
## Filter environments by heritability and mean bounds from the metrics table
## ───────────────────────────────────────────────────────────────────────────

  # Expand scalar filter bounds into per-trait named vectors
  heritLB <- rep(heritLB, length.out = length(trait)); names(heritLB) <- trait
  heritUB <- rep(heritUB, length.out = length(trait)); names(heritUB) <- trait
  meanLB  <- rep(meanLB,  length.out = length(trait)); names(meanLB)  <- trait
  meanUB  <- rep(meanUB,  length.out = length(trait)); names(meanUB)  <- trait

  # Load metrics from the STA analysisId(s)
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId), ]

  # Load predictions for the given STA analysisId(s)
  mydata <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId), ]

  # Per-trait environment filtering
  traitsToAnalyze <- character(0)
  filteredDataList <- list()
  filterWarnings <- list()

  # Heritability parameter names (same pattern as metLMMsolver)
  heritParams <- c("plotH2", "H2", "meanR2", "r2",
                   apply(expand.grid(c("plotH2", "H2", "meanR2", "r2"),
                                     c("designation", "mother", "father")),
                         1, function(f) { paste(f, collapse = "_") }))

  # Mean parameter names (same pattern as metLMMsolver)
  meanParams <- c("mean",
                  apply(expand.grid(c("mean"),
                                    c("designation", "mother", "father")),
                        1, function(f) { paste(f, collapse = "_") }))

  for (iTrait in trait) {
    # --- NEW: TPP trait resolution ---
    tpp_config <- phenoDTfile$metadata$tpp_analysis_config
    actual_trait <- iTrait
    tpp_env_filter <- NULL

    if (!is.null(tpp_config) && iTrait %in% names(tpp_config$trait_map)) {
      actual_trait <- tpp_config$trait_map[[iTrait]]
      tpp_env_filter <- tpp_config$env_map[[iTrait]]
    }
    # --- END TPP trait resolution ---

    # Subset data for this trait (use actual_trait for data lookup)
    provIdx <- which(mydata[, "trait"] == actual_trait)
    if (length(provIdx) == 0) next
    prov <- mydata[provIdx, ]

    # --- NEW: Apply TPP environment filter ---
    if (!is.null(tpp_env_filter)) {
      available_envs <- unique(prov$environment)
      valid_envs <- intersect(tpp_env_filter, available_envs)
      if (length(valid_envs) < 1) {
        warning(paste("TPP trait", iTrait, "skipped: no environments remained after filtering."))
        next
      }
      # Log info if some env_filter environments not found in data
      missing_envs <- setdiff(tpp_env_filter, available_envs)
      if (length(missing_envs) > 0 && verbose) {
        message(paste("TPP trait", iTrait, ": env_filter environments not found in data:",
                      paste(missing_envs, collapse = ", ")))
      }
      # Restrict to only environments in the TPP filter
      prov <- prov[which(prov$environment %in% valid_envs), ]
    }
    # --- END TPP environment filter ---

    # --- Filter by heritability bounds [heritLB, heritUB] ---
    metricsSub_h2 <- metrics[which(metrics$trait == actual_trait &
                                     metrics$parameter %in% heritParams), ]
    goodFields_h2 <- unique(metricsSub_h2[which(
      metricsSub_h2$value >= heritLB[iTrait] &
        metricsSub_h2$value <= heritUB[iTrait]
    ), "environment"])
    prov <- prov[which(prov$environment %in% goodFields_h2), ]

    # --- Filter by mean bounds [meanLB, meanUB] ---
    metricsSub_mean <- metrics[which(metrics$trait == actual_trait &
                                       metrics$parameter %in% meanParams), ]
    goodFields_mean <- unique(metricsSub_mean[which(
      metricsSub_mean$value >= meanLB[iTrait] &
        metricsSub_mean$value <= meanUB[iTrait]
    ), "environment"])
    prov <- prov[which(prov$environment %in% goodFields_mean), ]

    # --- Check minimum environment count ---
    # metRRBLUP is inherently genomic: a single environment is valid for GEBV estimation
    remainingEnvs <- unique(prov$environment)
    if (length(remainingEnvs) < 1) {
      # Record warning in modeling table and skip this trait
      warningEntry <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = "warningFilterEnv",
        value = paste0("Trait '", iTrait, "' skipped: no environments remained after heritability/mean filtering."),
        stringsAsFactors = FALSE
      )
      filterWarnings[[iTrait]] <- warningEntry
      if (verbose) {
        message(paste0("Warning: Trait '", iTrait, "' skipped - no environments remained after filtering."))
      }
      next
    }

    # Trait passed filtering; store for downstream processing
    traitsToAnalyze <- c(traitsToAnalyze, iTrait)
    filteredDataList[[iTrait]] <- prov

    if (verbose) {
      message(paste0("Trait '", iTrait, "': ", length(remainingEnvs),
                     " environments passed filtering."))
    }
  }

  # Append any filter warning entries to the modeling table
  if (length(filterWarnings) > 0) {
    warningDf <- do.call(rbind, filterWarnings)
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, warningDf)
  }

  ## ───────────────────────────────────────────────────────────────────────────
## MARKER MATRIX CONSTRUCTION (Task 1.3)
## Construct incidence matrices based on model_type
## ───────────────────────────────────────────────────────────────────────────

  # --- Load genotype matrix ---
  Markers <- NULL
  if (!is.null(phenoDTfile$data$geno)) {
    Markers <- as.matrix(phenoDTfile$data$geno)
  }
  if (is.null(Markers) || nrow(Markers) == 0 || ncol(Markers) == 0) {
    stop("Genotype matrix is not available. Please upload marker data and run Genotype QA/QC.", call. = FALSE)
  }

  # Apply genotype modifications if analysisIdGeno is provided

  if (!is.null(analysisIdGeno) && analysisIdGeno != "") {
    if (inherits(phenoDTfile$data$geno, "genlight")) {
      theresMatch <- which(as.character(analysisIdGeno) %in% names(phenoDTfile$data$geno_imp))
      if (length(theresMatch) > 0) {
        Markers <- as.matrix(phenoDTfile$data$geno_imp[[as.character(analysisIdGeno)]])
      }
    } else {
      modificationsMarkers <- phenoDTfile$modifications$geno
      theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdGeno)
      if (length(theresMatch) > 0) {
        modificationsMarkers <- modificationsMarkers[theresMatch, ]
        Markers <- cgiarBase::applyGenoModifications(M = Markers, modifications = modificationsMarkers)
      }
    }
  }

  # Check for missing marker data
  if (length(which(is.na(Markers))) > 0) {
    stop("Markers have missing data. Please run Genotype QA/QC to impute missing values before using this module.", call. = FALSE)
  }

  # --- Validate biallelic markers for A+D and SCA/GCA models ---
  if (model_type %in% c("Main Effects (A+D)", "SCA/GCA")) {
    unique_vals <- sort(unique(as.vector(Markers)))
    if (!all(unique_vals %in% c(0L, 1L, 2L))) {
      stop(
        "Dominance models require biallelic marker data coded as 0, 1, 2. ",
        "Found values: ", paste(head(unique_vals[!unique_vals %in% c(0, 1, 2)], 5), collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --- Validate pedigree data for GCA/SCA models ---
  if (model_type %in% c("GCA", "SCA/GCA")) {
    # Check mother/father columns exist in the predictions table
    if (!all(c("mother", "father") %in% colnames(baseData))) {
      stop("Pedigree information (mother and father columns) is required for the selected model type.", call. = FALSE)
    }
    # Get records with non-missing mother and father
    pedData <- baseData[!is.na(baseData$mother) & !is.na(baseData$father) &
                          baseData$mother != "" & baseData$father != "", ]
    if (nrow(pedData) == 0) {
      stop("Pedigree information (mother and father) is required for the selected model but no records have valid mother/father entries.", call. = FALSE)
    }
    unique_mothers <- unique(pedData$mother)
    unique_fathers <- unique(pedData$father)
    if (length(unique_mothers) < 2 || length(unique_fathers) < 2) {
      stop(
        "The SCA/GCA or GCA model requires at least 2 unique mothers and 2 unique fathers. ",
        "Found ", length(unique_mothers), " unique mother(s) and ",
        length(unique_fathers), " unique father(s).",
        call. = FALSE
      )
    }
    # Validate that parental genotypes exist in the marker matrix
    mothers_in_geno <- unique_mothers[unique_mothers %in% rownames(Markers)]
    fathers_in_geno <- unique_fathers[unique_fathers %in% rownames(Markers)]
    if (length(mothers_in_geno) == 0 || length(fathers_in_geno) == 0) {
      stop(
        "Parental genotypes not found in the marker matrix. ",
        "Found ", length(mothers_in_geno), " mother(s) and ",
        length(fathers_in_geno), " father(s) with genotype data.",
        call. = FALSE
      )
    }
  }

  # --- Construct marker matrices based on model_type ---
  M_a <- NULL
  M_d <- NULL
  M_mother <- NULL
  M_father <- NULL
  M_sca <- NULL
  het_covariate <- NULL

  if (model_type == "Main Effects (A)") {
    # Center the marker dosage matrix: M_a = M - 2p (where p = column means / 2)
    # Subset to individuals present in the analysis data
    indiv_ids <- unique(baseData$designation)
    common_ids <- intersect(indiv_ids, rownames(Markers))
    M_sub <- Markers[common_ids, , drop = FALSE]
    p <- colMeans(M_sub) / 2  # allele frequencies
    M_a <- sweep(M_sub, 2, 2 * p, FUN = "-")

  } else if (model_type == "Main Effects (A+D)") {
    # Additive component: M_a = M - 2p
    indiv_ids <- unique(baseData$designation)
    common_ids <- intersect(indiv_ids, rownames(Markers))
    M_sub <- Markers[common_ids, , drop = FALSE]
    p <- colMeans(M_sub) / 2  # allele frequencies
    M_a <- sweep(M_sub, 2, 2 * p, FUN = "-")

    # Dominance component: M_d[i,j] = 1 if M[i,j] == 1, else 0, minus column mean
    M_het <- ifelse(M_sub == 1, 1, 0)
    het_col_means <- colMeans(M_het)
    M_d <- sweep(M_het, 2, het_col_means, FUN = "-")

    # Mean heterozygosity per individual as fixed covariate (directional dominance)
    het_covariate <- rowMeans(M_het)

  } else if (model_type == "GCA") {
    # Construct separate maternal and paternal marker matrices
    # For each individual (hybrid), mother contributes half the alleles, father the other half
    # Use mother/father columns from predictions to look up rows in genotype matrix

    # Keep only records with valid pedigree and parental genotypes
    pedData <- baseData[!is.na(baseData$mother) & !is.na(baseData$father) &
                          baseData$mother != "" & baseData$father != "", ]
    pedData <- pedData[pedData$mother %in% rownames(Markers) &
                         pedData$father %in% rownames(Markers), ]

    # Get unique hybrid-mother-father mappings
    ped_map <- unique(pedData[, c("designation", "mother", "father")])

    # Construct maternal marker matrix (rows = hybrids, columns = markers)
    M_mother_raw <- Markers[ped_map$mother, , drop = FALSE]
    rownames(M_mother_raw) <- ped_map$designation

    # Construct paternal marker matrix
    M_father_raw <- Markers[ped_map$father, , drop = FALSE]
    rownames(M_father_raw) <- ped_map$designation

    # Center each by column means
    M_mother <- sweep(M_mother_raw, 2, colMeans(M_mother_raw), FUN = "-")
    M_father <- sweep(M_father_raw, 2, colMeans(M_father_raw), FUN = "-")

  } else if (model_type == "SCA/GCA") {
    # Construct maternal, paternal, and SCA interaction matrices
    # Following De Jong model: GCA (maternal + paternal) + SCA (interaction) + directional dominance

    # Keep only records with valid pedigree and parental genotypes
    pedData <- baseData[!is.na(baseData$mother) & !is.na(baseData$father) &
                          baseData$mother != "" & baseData$father != "", ]
    pedData <- pedData[pedData$mother %in% rownames(Markers) &
                         pedData$father %in% rownames(Markers), ]

    # Get unique hybrid-mother-father mappings
    ped_map <- unique(pedData[, c("designation", "mother", "father")])

    # Construct maternal marker matrix
    M_mother_raw <- Markers[ped_map$mother, , drop = FALSE]
    rownames(M_mother_raw) <- ped_map$designation

    # Construct paternal marker matrix
    M_father_raw <- Markers[ped_map$father, , drop = FALSE]
    rownames(M_father_raw) <- ped_map$designation

    # Center each by column means
    M_mother <- sweep(M_mother_raw, 2, colMeans(M_mother_raw), FUN = "-")
    M_father <- sweep(M_father_raw, 2, colMeans(M_father_raw), FUN = "-")

    # SCA interaction matrix: element-wise product of centered maternal and paternal matrices
    M_sca <- M_mother * M_father

    # Mean heterozygosity per individual as fixed covariate (directional dominance)
    # For hybrids: heterozygosity at locus j = 1 if mother allele != father allele
    # Using raw (non-centered) parental dosages to compute heterozygosity
    M_het_hybrid <- ifelse(M_mother_raw != M_father_raw, 1, 0)
    het_covariate <- rowMeans(M_het_hybrid)
  }

  if (verbose) {
    if (!is.null(M_a)) message("Additive marker matrix: ", nrow(M_a), " individuals x ", ncol(M_a), " markers")
    if (!is.null(M_d)) message("Dominance marker matrix: ", nrow(M_d), " individuals x ", ncol(M_d), " markers")
    if (!is.null(M_mother)) message("Maternal marker matrix: ", nrow(M_mother), " individuals x ", ncol(M_mother), " markers")
    if (!is.null(M_father)) message("Paternal marker matrix: ", nrow(M_father), " individuals x ", ncol(M_father), " markers")
    if (!is.null(M_sca)) message("SCA interaction matrix: ", nrow(M_sca), " individuals x ", ncol(M_sca), " markers")
    if (!is.null(het_covariate)) message("Heterozygosity covariate computed for ", length(het_covariate), " individuals")
  }

  ## ───────────────────────────────────────────────────────────────────────────
## MODEL FITTING (Task 2.1) - mode == "fit"
## Fit LMMsolver models for each trait x environment combination
## ───────────────────────────────────────────────────────────────────────────

  fittingResults <- list()  # stores marker effects per trait x environment

  if (mode == "fit" && length(traitsToAnalyze) > 0) {

    if (verbose) message("Fitting RR-BLUP models (mode = 'fit').")

    for (iTrait in traitsToAnalyze) {

      prov <- filteredDataList[[iTrait]]
      environments <- unique(prov$environment)

      if (verbose) message(paste("  Analyzing trait:", iTrait))

      for (iEnv in environments) {

        # Subset data for this trait x environment
        envData <- prov[which(prov$environment == iEnv), ]
        if (nrow(envData) == 0) next

        # --- Determine which individuals to include ---
        # Only individuals present in BOTH genotype matrix rows and predictions data
        if (model_type %in% c("Main Effects (A)", "Main Effects (A+D)")) {
          available_ids <- intersect(unique(envData$designation), rownames(M_a))
        } else {
          # GCA / SCA/GCA: individuals must be in M_mother rows
          available_ids <- intersect(unique(envData$designation), rownames(M_mother))
        }

        if (length(available_ids) < 2) {
          if (verbose) {
            message(paste("    Skipping environment", iEnv,
                          "- fewer than 2 individuals with genotype data."))
          }
          next
        }

        # Subset to common individuals
        envData <- envData[which(envData$designation %in% available_ids), ]
        # Remove duplicates (keep first occurrence per designation within the environment)
        envData <- envData[!duplicated(envData$designation), ]

        # --- Build the model data frame ---
        modelData <- data.frame(
          predictedValue = envData$predictedValue,
          stringsAsFactors = FALSE
        )
        rownames(modelData) <- envData$designation

        # --- Construct fixed formula and add fixed-effect columns ---
        if (model_type %in% c("Main Effects (A+D)", "SCA/GCA")) {
          # Directional dominance: mean heterozygosity as fixed covariate
          modelData$het_covariate <- het_covariate[envData$designation]
          fixFormula <- "predictedValue ~ het_covariate"
        } else {
          # Main Effects (A) and GCA: intercept only
          fixFormula <- "predictedValue ~ 1"
        }

        # --- Construct random effects (incidence matrices via group argument) ---
        grouping <- list()
        startCol <- ncol(modelData) + 1

        if (model_type == "Main Effects (A)") {
          # Additive marker matrix
          M_sub <- M_a[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_sub))
          endCol <- ncol(modelData)
          grouping[["marker"]] <- startCol:endCol

        } else if (model_type == "Main Effects (A+D)") {
          # Additive marker matrix
          M_a_sub <- M_a[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_a_sub))
          endCol <- ncol(modelData)
          grouping[["marker_a"]] <- startCol:endCol

          # Dominance marker matrix
          startCol2 <- ncol(modelData) + 1
          M_d_sub <- M_d[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_d_sub))
          endCol2 <- ncol(modelData)
          grouping[["marker_d"]] <- startCol2:endCol2

        } else if (model_type == "GCA") {
          # Maternal marker matrix
          M_m_sub <- M_mother[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_m_sub))
          endCol <- ncol(modelData)
          grouping[["gca_mother"]] <- startCol:endCol

          # Paternal marker matrix
          startCol2 <- ncol(modelData) + 1
          M_f_sub <- M_father[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_f_sub))
          endCol2 <- ncol(modelData)
          grouping[["gca_father"]] <- startCol2:endCol2

        } else if (model_type == "SCA/GCA") {
          # Maternal marker matrix
          M_m_sub <- M_mother[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_m_sub))
          endCol <- ncol(modelData)
          grouping[["gca_mother"]] <- startCol:endCol

          # Paternal marker matrix
          startCol2 <- ncol(modelData) + 1
          M_f_sub <- M_father[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_f_sub))
          endCol2 <- ncol(modelData)
          grouping[["gca_father"]] <- startCol2:endCol2

          # SCA interaction matrix
          startCol3 <- ncol(modelData) + 1
          M_sca_sub <- M_sca[envData$designation, , drop = FALSE]
          modelData <- cbind(modelData, as.data.frame(M_sca_sub))
          endCol3 <- ncol(modelData)
          grouping[["sca"]] <- startCol3:endCol3
        }

        # --- Construct random formula using grp() ---
        ranFormula <- paste("~", paste(paste0("grp(", names(grouping), ")"), collapse = " + "))

        # --- Weights ---
        weightsArg <- NULL
        if (useWeights && "reliability" %in% colnames(envData)) {
          rel_vals <- envData$reliability
          if (!all(is.na(rel_vals)) && any(rel_vals > 0, na.rm = TRUE)) {
            # Use reliability as weights; replace NA/zero with small value
            rel_vals[is.na(rel_vals) | rel_vals <= 0] <- min(rel_vals[rel_vals > 0], na.rm = TRUE) * 0.1
            modelData$w <- rel_vals
            weightsArg <- "w"
          }
        }

        # --- Fit the model with tryCatch for non-convergence ---
        mix <- tryCatch({
          LMMsolver::LMMsolve(
            fixed = as.formula(fixFormula),
            random = as.formula(ranFormula),
            group = grouping,
            weights = weightsArg,
            data = modelData,
            maxit = maxIters
          )
        }, error = function(e) {
          if (verbose) {
            message(paste("    WARNING: Model failed for trait", iTrait,
                          "environment", iEnv, ":", conditionMessage(e)))
          }
          return(NULL)
        })

        # Skip if model failed
        if (is.null(mix) || inherits(mix, "try-error")) {
          if (verbose) {
            message(paste("    Skipping trait", iTrait, "environment", iEnv,
                          "due to non-convergence."))
          }
          next
        }

        # --- Extract marker effects as BLUPs ---
        markerEffects <- list()

        if (model_type == "Main Effects (A)") {
          pick <- mix$ndxCoefficients[["marker"]]
          effects_a <- as.vector(mix$coefMME[pick])
          names(effects_a) <- colnames(M_a)
          markerEffects[["marker_a"]] <- effects_a

        } else if (model_type == "Main Effects (A+D)") {
          pick_a <- mix$ndxCoefficients[["marker_a"]]
          effects_a <- as.vector(mix$coefMME[pick_a])
          names(effects_a) <- colnames(M_a)
          markerEffects[["marker_a"]] <- effects_a

          pick_d <- mix$ndxCoefficients[["marker_d"]]
          effects_d <- as.vector(mix$coefMME[pick_d])
          names(effects_d) <- colnames(M_d)
          markerEffects[["marker_d"]] <- effects_d

        } else if (model_type == "GCA") {
          pick_m <- mix$ndxCoefficients[["gca_mother"]]
          effects_m <- as.vector(mix$coefMME[pick_m])
          names(effects_m) <- colnames(M_mother)
          markerEffects[["marker_gca_mother"]] <- effects_m

          pick_f <- mix$ndxCoefficients[["gca_father"]]
          effects_f <- as.vector(mix$coefMME[pick_f])
          names(effects_f) <- colnames(M_father)
          markerEffects[["marker_gca_father"]] <- effects_f

        } else if (model_type == "SCA/GCA") {
          pick_m <- mix$ndxCoefficients[["gca_mother"]]
          effects_m <- as.vector(mix$coefMME[pick_m])
          names(effects_m) <- colnames(M_mother)
          markerEffects[["marker_gca_mother"]] <- effects_m

          pick_f <- mix$ndxCoefficients[["gca_father"]]
          effects_f <- as.vector(mix$coefMME[pick_f])
          names(effects_f) <- colnames(M_father)
          markerEffects[["marker_gca_father"]] <- effects_f

          pick_sca <- mix$ndxCoefficients[["sca"]]
          effects_sca <- as.vector(mix$coefMME[pick_sca])
          names(effects_sca) <- colnames(M_sca)
          markerEffects[["marker_sca"]] <- effects_sca
        }

        # Store results keyed by trait and environment
        resultKey <- paste(iTrait, iEnv, sep = "||")
        fittingResults[[resultKey]] <- list(
          trait = iTrait,
          environment = iEnv,
          markerEffects = markerEffects,
          model = mix,
          individuals = envData$designation
        )

        if (verbose) {
          message(paste("    Successfully fitted model for environment:", iEnv))
        }

      } # end environment loop
    } # end trait loop

    if (verbose) {
      message(paste("  Model fitting complete.", length(fittingResults),
                    "trait x environment combinations fitted successfully."))
    }

  } # end if mode == "fit"

  ## ───────────────────────────────────────────────────────────────────────────
## INDIVIDUAL-LEVEL PREDICTION COMPUTATION (Task 2.2)
## Compute individual genomic values: M × β̂ for ALL individuals in genotype matrix
## ───────────────────────────────────────────────────────────────────────────

  individualPredictions <- list()

  if (mode == "fit" && length(fittingResults) > 0) {

    if (verbose) message("Computing individual-level predictions (M x beta_hat).")

    for (resultKey in names(fittingResults)) {

      res <- fittingResults[[resultKey]]
      iTrait <- res$trait
      iEnv <- res$environment
      markerEffects <- res$markerEffects

      # Initialize a list to collect data frames for this key
      predRows <- list()

      if (model_type == "Main Effects (A)") {
        # Compute for ALL individuals in M_a (full genotype matrix, already centered)
        marker_a <- markerEffects[["marker_a"]]
        # Genomic value = M_a %*% marker_a
        gv <- as.vector(M_a %*% marker_a)
        names(gv) <- rownames(M_a)

        predRows[["designation"]] <- data.frame(
          designation = rownames(M_a),
          predictedValue = gv,
          effectType = "designation",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "Main Effects (A+D)") {
        # Additive component
        marker_a <- markerEffects[["marker_a"]]
        gv_a <- as.vector(M_a %*% marker_a)
        names(gv_a) <- rownames(M_a)

        # Dominance component
        marker_d <- markerEffects[["marker_d"]]
        gv_d <- as.vector(M_d %*% marker_d)
        names(gv_d) <- rownames(M_d)

        # Total = additive + dominance
        gv_total <- gv_a + gv_d

        predRows[["designationA"]] <- data.frame(
          designation = rownames(M_a),
          predictedValue = gv_a,
          effectType = "designationA",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

        predRows[["designationD"]] <- data.frame(
          designation = rownames(M_d),
          predictedValue = gv_d,
          effectType = "designationD",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

        predRows[["designation"]] <- data.frame(
          designation = rownames(M_a),
          predictedValue = gv_total,
          effectType = "designation",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "GCA") {
        # GCA = maternal contribution + paternal contribution
        marker_gca_mother <- markerEffects[["marker_gca_mother"]]
        marker_gca_father <- markerEffects[["marker_gca_father"]]

        gv_gca <- as.vector(M_mother %*% marker_gca_mother) +
                   as.vector(M_father %*% marker_gca_father)
        names(gv_gca) <- rownames(M_mother)

        predRows[["designationGCA"]] <- data.frame(
          designation = rownames(M_mother),
          predictedValue = gv_gca,
          effectType = "designationGCA",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

        # For GCA model, total = GCA
        predRows[["designation"]] <- data.frame(
          designation = rownames(M_mother),
          predictedValue = gv_gca,
          effectType = "designation",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "SCA/GCA") {
        # GCA component = maternal + paternal
        marker_gca_mother <- markerEffects[["marker_gca_mother"]]
        marker_gca_father <- markerEffects[["marker_gca_father"]]

        gv_gca <- as.vector(M_mother %*% marker_gca_mother) +
                   as.vector(M_father %*% marker_gca_father)
        names(gv_gca) <- rownames(M_mother)

        # SCA component
        marker_sca <- markerEffects[["marker_sca"]]
        gv_sca <- as.vector(M_sca %*% marker_sca)
        names(gv_sca) <- rownames(M_sca)

        # Total = GCA + SCA
        gv_total <- gv_gca + gv_sca

        predRows[["designationGCA"]] <- data.frame(
          designation = rownames(M_mother),
          predictedValue = gv_gca,
          effectType = "designationGCA",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

        predRows[["designationSCA"]] <- data.frame(
          designation = rownames(M_sca),
          predictedValue = gv_sca,
          effectType = "designationSCA",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )

        predRows[["designation"]] <- data.frame(
          designation = rownames(M_mother),
          predictedValue = gv_total,
          effectType = "designation",
          trait = iTrait,
          environment = iEnv,
          stringsAsFactors = FALSE
        )
      }

      # Combine all component data frames for this key
      if (length(predRows) > 0) {
        individualPredictions[[resultKey]] <- do.call(rbind, predRows)
        rownames(individualPredictions[[resultKey]]) <- NULL
      }

    } # end loop over fittingResults

    if (verbose) {
      totalIndivPreds <- sum(vapply(individualPredictions, nrow, integer(1)))
      message(paste("  Individual-level predictions computed:", totalIndivPreds, "rows across",
                    length(individualPredictions), "trait x environment combinations."))
    }

  } # end if mode == "fit" for individual predictions

  ## ───────────────────────────────────────────────────────────────────────────
## PREDICT MODE (Task 3.1) - mode == "predict"
## Predict unphenotyped individuals using stored marker effects
## ───────────────────────────────────────────────────────────────────────────

  predictResults <- list()
  predictMarkerCounts <- list()

  if (mode == "predict") {

    if (is.null(markerEffectsId) || markerEffectsId == "") {
      stop("markerEffectsId is required when mode = 'predict'.", call. = FALSE)
    }

    if (verbose) message("Predict mode: retrieving marker effects for analysisId = ", markerEffectsId)

    # --- 1. Retrieve marker effect predictions from the predictions table ---
    allPreds <- phenoDTfile$predictions
    markerEffectRows <- allPreds[
      which(allPreds$analysisId == markerEffectsId &
              grepl("^marker_", allPreds$effectType)), ]

    if (nrow(markerEffectRows) == 0) {
      stop("No marker effects found in predictions table for markerEffectsId = '",
           markerEffectsId, "'. Ensure the selected MTA stamp contains marker-level predictions.",
           call. = FALSE)
    }

    # --- 2. The genotype matrix (Markers) is already loaded from earlier in the function ---
    # Get marker names from the genotype matrix columns
    genoMarkerNames <- colnames(Markers)

    # --- 3. Intersect marker names between effects and genotype matrix columns ---
    effectMarkerNames <- unique(markerEffectRows$designation)
    commonMarkers <- intersect(effectMarkerNames, genoMarkerNames)

    if (length(commonMarkers) == 0) {
      stop("Zero marker name overlap between the stored marker effects and the genotype matrix. ",
           "Cannot compute predictions.", call. = FALSE)
    }

    nMarkersMatched <- length(commonMarkers)
    nMarkersTotal <- length(effectMarkerNames)

    if (verbose) {
      message("  Marker overlap: ", nMarkersMatched, " out of ", nMarkersTotal,
              " effect markers found in genotype matrix (",
              ncol(Markers), " total genotype columns).")
    }

    # --- 4. Center the genotype matrix columns ---
    # Use column means of the full genotype matrix: p = colMeans(M) / 2, M_centered = M - 2p
    # This matches the centering applied during fitting
    M_predict <- Markers[, commonMarkers, drop = FALSE]
    p <- colMeans(M_predict) / 2
    M_centered <- sweep(M_predict, 2, 2 * p, FUN = "-")

    # --- 5. Compute genomic predictions for each trait × environment × effectType ---
    # Get unique trait × environment combinations
    traitEnvCombos <- unique(markerEffectRows[, c("trait", "environment")])

    for (i in seq_len(nrow(traitEnvCombos))) {
      iTrait <- traitEnvCombos$trait[i]
      iEnv   <- traitEnvCombos$environment[i]
      resultKey <- paste(iTrait, iEnv, sep = "||")

      # Get all marker effect rows for this trait × environment
      teRows <- markerEffectRows[
        which(markerEffectRows$trait == iTrait &
                markerEffectRows$environment == iEnv), ]

      # Get unique effectTypes for this combination
      uniqueEffectTypes <- unique(teRows$effectType)

      # Compute predictions for each effectType
      predsByEffectType <- list()

      for (effType in uniqueEffectTypes) {
        # Get effects for this effectType
        effRows <- teRows[which(teRows$effectType == effType), ]

        # Build named effect vector (marker name -> effect value)
        effectVec <- setNames(effRows$predictedValue, effRows$designation)

        # Subset to common markers (maintain order of commonMarkers)
        effectsCommon <- effectVec[commonMarkers]

        # Handle any NAs (markers in commonMarkers but missing from this effectType)
        effectsCommon[is.na(effectsCommon)] <- 0

        # Compute prediction: M_centered %*% beta_hat
        gv <- as.vector(M_centered %*% effectsCommon)
        names(gv) <- rownames(M_centered)

        predsByEffectType[[effType]] <- gv
      }

      # Store results keyed by "trait||environment"
      predictResults[[resultKey]] <- list(
        trait = iTrait,
        environment = iEnv,
        markerEffectTypes = uniqueEffectTypes,
        predictions = predsByEffectType,
        nMarkersMatched = nMarkersMatched,
        nMarkersTotal = nMarkersTotal
      )

    } # end loop over trait × environment combinations

    # Store marker counts for later use by Task 3.2
    predictMarkerCounts <- list(
      nMarkersMatched = nMarkersMatched,
      nMarkersTotal = nMarkersTotal
    )

    if (verbose) {
      message("  Predict mode complete: ", length(predictResults),
              " trait x environment combinations predicted for ",
              nrow(M_centered), " individuals.")
    }

  } # end if mode == "predict"

  ## ───────────────────────────────────────────────────────────────────────────
## PREDICT MODE OUTPUT ASSEMBLY (Task 3.2)
## Map effectTypes and assemble output tables for predict mode
## ───────────────────────────────────────────────────────────────────────────

  if (mode == "predict" && length(predictResults) > 0) {

    if (verbose) message("Assembling predict mode output tables.")

    predictPredRows <- list()

    for (resultKey in names(predictResults)) {

      res <- predictResults[[resultKey]]
      iTrait <- res$trait
      iEnv <- res$environment
      markerEffectTypes <- res$markerEffectTypes
      predictions <- res$predictions

      # --- effectType mapping logic ---
      if ("marker_a" %in% markerEffectTypes && "marker_d" %in% markerEffectTypes) {
        # A+D model: produce designation_A, designation_D, and designation (combined)
        gv_a <- predictions[["marker_a"]]
        gv_d <- predictions[["marker_d"]]
        gv_combined <- gv_a + gv_d

        predictPredRows[[paste0(resultKey, "_designation_A")]] <- data.frame(
          designation = names(gv_a),
          predictedValue = as.numeric(gv_a),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation_A",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

        predictPredRows[[paste0(resultKey, "_designation_D")]] <- data.frame(
          designation = names(gv_d),
          predictedValue = as.numeric(gv_d),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation_D",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

        predictPredRows[[paste0(resultKey, "_designation")]] <- data.frame(
          designation = names(gv_combined),
          predictedValue = as.numeric(gv_combined),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if ("marker_gca" %in% markerEffectTypes && "marker_sca" %in% markerEffectTypes) {
        # SCA/GCA model: produce designation_GCA, designation_SCA, and designation (combined)
        gv_gca <- predictions[["marker_gca"]]
        gv_sca <- predictions[["marker_sca"]]
        gv_combined <- gv_gca + gv_sca

        predictPredRows[[paste0(resultKey, "_designation_GCA")]] <- data.frame(
          designation = names(gv_gca),
          predictedValue = as.numeric(gv_gca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation_GCA",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

        predictPredRows[[paste0(resultKey, "_designation_SCA")]] <- data.frame(
          designation = names(gv_sca),
          predictedValue = as.numeric(gv_sca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation_SCA",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

        predictPredRows[[paste0(resultKey, "_designation")]] <- data.frame(
          designation = names(gv_combined),
          predictedValue = as.numeric(gv_combined),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if ("marker_a" %in% markerEffectTypes) {
        # Additive-only model: produce designation (combined = just additive)
        gv_a <- predictions[["marker_a"]]

        predictPredRows[[paste0(resultKey, "_designation")]] <- data.frame(
          designation = names(gv_a),
          predictedValue = as.numeric(gv_a),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if ("marker_gca" %in% markerEffectTypes) {
        # GCA-only model: produce designation (combined = just GCA)
        gv_gca <- predictions[["marker_gca"]]

        predictPredRows[[paste0(resultKey, "_designation")]] <- data.frame(
          designation = names(gv_gca),
          predictedValue = as.numeric(gv_gca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "designation",
          environment = iEnv,
          entryType = "test",
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )
      }

    } # end loop over predictResults

    ## ─── APPEND PREDICTIONS ──────────────────────────────────────────────────
    if (length(predictPredRows) > 0) {
      newPredictions <- do.call(rbind, predictPredRows)
      rownames(newPredictions) <- NULL

      # Ensure effectType column exists in existing predictions
      if (!is.null(phenoDTfile$predictions)) {
        if (!"effectType" %in% colnames(phenoDTfile$predictions)) {
          phenoDTfile$predictions$effectType <- NA
        }
      }

      # Append using rbind, subsetting to match existing columns
      phenoDTfile$predictions <- rbind(
        phenoDTfile$predictions,
        newPredictions[, colnames(phenoDTfile$predictions)]
      )
    }

    ## ─── STATUS TABLE ────────────────────────────────────────────────────────
    newStatus <- data.frame(
      module = "mta_rrblup",
      analysisId = mtaAnalysisId,
      analysisIdName = NA_character_,
      stringsAsFactors = FALSE
    )
    phenoDTfile$status <- rbind(
      phenoDTfile$status,
      newStatus[, colnames(phenoDTfile$status)]
    )

    ## ─── MODELING TABLE (marker count info + input analysisId) ────────────────
    predictModelingRows <- list()

    # Get unique traits from the predictions
    predictTraits <- unique(vapply(predictResults, function(r) r$trait, character(1)))

    for (iTrait in predictTraits) {
      # Find the first result for this trait to get marker counts
      traitKeys <- grep(paste0("^", iTrait, "\\|\\|"), names(predictResults), value = TRUE)
      if (length(traitKeys) > 0) {
        traitRes <- predictResults[[traitKeys[1]]]

        predictModelingRows[[paste0(iTrait, "_nMarkersMatched")]] <- data.frame(
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          trait = iTrait,
          environment = "general",
          parameter = "nMarkersMatched",
          value = as.character(traitRes$nMarkersMatched),
          stringsAsFactors = FALSE
        )

        predictModelingRows[[paste0(iTrait, "_nMarkersTotal")]] <- data.frame(
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          trait = iTrait,
          environment = "general",
          parameter = "nMarkersTotal",
          value = as.character(traitRes$nMarkersTotal),
          stringsAsFactors = FALSE
        )
      }
    }

    # Record input marker effects analysisId
    predictModelingRows[["input_markerEffects"]] <- data.frame(
      module = "mta_rrblup",
      analysisId = mtaAnalysisId,
      trait = "inputObject",
      environment = "general",
      parameter = "analysisId",
      value = as.character(markerEffectsId),
      stringsAsFactors = FALSE
    )

    if (length(predictModelingRows) > 0) {
      newModeling <- do.call(rbind, predictModelingRows)
      rownames(newModeling) <- NULL
      phenoDTfile$modeling <- rbind(
        phenoDTfile$modeling,
        newModeling[, colnames(phenoDTfile$modeling)]
      )
    }

    if (verbose) message("Predict mode output assembly complete. analysisId: ", mtaAnalysisId)

  } # end if mode == "predict" && length(predictResults) > 0

  ## ───────────────────────────────────────────────────────────────────────────
## OUTPUT ASSEMBLY (Task 2.3)
## Assemble predictions, metrics, modeling, and status tables
## ───────────────────────────────────────────────────────────────────────────

  if (mode == "fit" && length(fittingResults) > 0) {

    if (verbose) message("Assembling output tables.")

    ## ─── 1. MARKER-LEVEL PREDICTIONS ─────────────────────────────────────────
    markerPredRows <- list()

    for (resultKey in names(fittingResults)) {
      res <- fittingResults[[resultKey]]
      iTrait <- res$trait
      iEnv <- res$environment
      markerEffects <- res$markerEffects

      if (model_type == "Main Effects (A)") {
        # marker_a rows
        effects_a <- markerEffects[["marker_a"]]
        markerPredRows[[paste0(resultKey, "_marker_a")]] <- data.frame(
          designation = names(effects_a),
          predictedValue = as.numeric(effects_a),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_a",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "Main Effects (A+D)") {
        # marker_a rows
        effects_a <- markerEffects[["marker_a"]]
        markerPredRows[[paste0(resultKey, "_marker_a")]] <- data.frame(
          designation = names(effects_a),
          predictedValue = as.numeric(effects_a),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_a",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )
        # marker_d rows
        effects_d <- markerEffects[["marker_d"]]
        markerPredRows[[paste0(resultKey, "_marker_d")]] <- data.frame(
          designation = names(effects_d),
          predictedValue = as.numeric(effects_d),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_d",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "GCA") {
        # Combined GCA marker effects: (maternal + paternal) / 2
        effects_m <- markerEffects[["marker_gca_mother"]]
        effects_f <- markerEffects[["marker_gca_father"]]
        effects_gca <- (effects_m + effects_f) / 2
        markerPredRows[[paste0(resultKey, "_marker_gca")]] <- data.frame(
          designation = names(effects_gca),
          predictedValue = as.numeric(effects_gca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_gca",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )

      } else if (model_type == "SCA/GCA") {
        # Combined GCA marker effects: (maternal + paternal) / 2
        effects_m <- markerEffects[["marker_gca_mother"]]
        effects_f <- markerEffects[["marker_gca_father"]]
        effects_gca <- (effects_m + effects_f) / 2
        markerPredRows[[paste0(resultKey, "_marker_gca")]] <- data.frame(
          designation = names(effects_gca),
          predictedValue = as.numeric(effects_gca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_gca",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )
        # marker_sca rows
        effects_sca <- markerEffects[["marker_sca"]]
        markerPredRows[[paste0(resultKey, "_marker_sca")]] <- data.frame(
          designation = names(effects_sca),
          predictedValue = as.numeric(effects_sca),
          stdError = NA_real_,
          reliability = NA_real_,
          trait = iTrait,
          effectType = "marker_sca",
          environment = iEnv,
          entryType = NA_character_,
          gid = NA_character_,
          mother = NA_character_,
          father = NA_character_,
          pipeline = NA_character_,
          module = "mta_rrblup",
          analysisId = mtaAnalysisId,
          stringsAsFactors = FALSE
        )
      }
    }

    ## ─── 2. INDIVIDUAL-LEVEL PREDICTIONS ─────────────────────────────────────
    ## Add missing columns to the individual prediction data frames

    # Lookup table for mother/father/gid/pipeline from baseData
    .first_sorted <- function(x) (sort(x, decreasing = FALSE, na.last = TRUE))[1]
    boGid    <- aggregate(baseData[, "gid",    drop = FALSE], by = list(designation = baseData[, "designation"]), FUN = .first_sorted)
    boMother <- aggregate(baseData[, "mother", drop = FALSE], by = list(designation = baseData[, "designation"]), FUN = .first_sorted)
    boFather <- aggregate(baseData[, "father", drop = FALSE], by = list(designation = baseData[, "designation"]), FUN = .first_sorted)
    if ("pipeline" %in% colnames(baseData)) {
      boPipe <- aggregate(baseData[, "pipeline", drop = FALSE], by = list(designation = baseData[, "designation"]),
                          FUN = function(x) paste(unique(sort(x, decreasing = FALSE)), collapse = ", "))
      baseOrigin <- data.frame(designation = boGid$designation, gid = boGid$gid,
                               mother = boMother$mother, father = boFather$father,
                               pipeline = boPipe$pipeline, stringsAsFactors = FALSE)
    } else {
      baseOrigin <- data.frame(designation = boGid$designation, gid = boGid$gid,
                               mother = boMother$mother, father = boFather$father,
                               pipeline = NA_character_, stringsAsFactors = FALSE)
    }

    # Combine all individual prediction data frames
    if (length(individualPredictions) > 0) {
      indivPredBind <- do.call(rbind, individualPredictions)
      rownames(indivPredBind) <- NULL

      # Add remaining columns
      indivPredBind$stdError    <- NA_real_
      indivPredBind$reliability <- NA_real_
      indivPredBind$entryType   <- "test"
      indivPredBind$module      <- "mta_rrblup"
      indivPredBind$analysisId  <- mtaAnalysisId

      # Merge with baseOrigin to get gid, mother, father, pipeline
      indivPredBind <- merge(indivPredBind, baseOrigin, by = "designation", all.x = TRUE)
    } else {
      indivPredBind <- NULL
    }

    ## ─── 3. COMBINE ALL PREDICTIONS AND APPEND ───────────────────────────────
    # Combine marker-level and individual-level predictions
    allNewPreds <- list()
    if (length(markerPredRows) > 0) {
      allNewPreds[["markers"]] <- do.call(rbind, markerPredRows)
    }
    if (!is.null(indivPredBind)) {
      allNewPreds[["individuals"]] <- indivPredBind
    }

    if (length(allNewPreds) > 0) {
      newPredictions <- do.call(rbind, allNewPreds)
      rownames(newPredictions) <- NULL

      # Ensure effectType column exists in existing predictions
      if (!is.null(phenoDTfile$predictions)) {
        if (!"effectType" %in% colnames(phenoDTfile$predictions)) {
          phenoDTfile$predictions$effectType <- NA
        }
      }

      # Append using rbind, subsetting to match existing columns
      phenoDTfile$predictions <- rbind(
        phenoDTfile$predictions,
        newPredictions[, colnames(phenoDTfile$predictions)]
      )
    }

    ## ─── 4. METRICS TABLE ────────────────────────────────────────────────────
    # Add basic metrics per trait (nEnv, nEntries)
    metricsRows <- list()
    for (iTrait in traitsToAnalyze) {
      # Count how many environments were successfully fitted for this trait
      traitKeys <- grep(paste0("^", iTrait, "\\|\\|"), names(fittingResults), value = TRUE)
      nEnvFitted <- length(traitKeys)
      # Count unique individuals
      traitIndivPreds <- individualPredictions[traitKeys]
      nEntries <- 0
      if (length(traitIndivPreds) > 0) {
        allDesig <- unique(unlist(lapply(traitIndivPreds, function(df) {
          df$designation[df$effectType == "designation"]
        })))
        nEntries <- length(allDesig)
      }

      metricsRows[[iTrait]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = c("nEnv", "nEntries"),
        method = c("n", "n"),
        value = c(nEnvFitted, nEntries),
        stdError = c(NA_real_, NA_real_),
        stringsAsFactors = FALSE
      )
    }

    if (length(metricsRows) > 0) {
      newMetrics <- do.call(rbind, metricsRows)
      rownames(newMetrics) <- NULL
      phenoDTfile$metrics <- rbind(
        phenoDTfile$metrics,
        newMetrics[, colnames(phenoDTfile$metrics)]
      )
    }

    ## ─── 5. MODELING TABLE ───────────────────────────────────────────────────
    # Entries per trait: modelType, fixedFormula, randomFormula, selectedModel
    modelingRows <- list()
    for (iTrait in traitsToAnalyze) {
      # Determine the formulas used (they are the same for all envs within a trait)
      if (model_type %in% c("Main Effects (A+D)", "SCA/GCA")) {
        fixFormulaStr <- "predictedValue ~ het_covariate"
      } else {
        fixFormulaStr <- "predictedValue ~ 1"
      }

      if (model_type == "Main Effects (A)") {
        ranFormulaStr <- "~ grp(marker)"
      } else if (model_type == "Main Effects (A+D)") {
        ranFormulaStr <- "~ grp(marker_a) + grp(marker_d)"
      } else if (model_type == "GCA") {
        ranFormulaStr <- "~ grp(gca_mother) + grp(gca_father)"
      } else if (model_type == "SCA/GCA") {
        ranFormulaStr <- "~ grp(gca_mother) + grp(gca_father) + grp(sca)"
      }

      modelingRows[[paste0(iTrait, "_modelType")]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = "modelType",
        value = "rrblup",
        stringsAsFactors = FALSE
      )
      modelingRows[[paste0(iTrait, "_fixedFormula")]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = "fixedFormula",
        value = fixFormulaStr,
        stringsAsFactors = FALSE
      )
      modelingRows[[paste0(iTrait, "_randomFormula")]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = "randomFormula",
        value = ranFormulaStr,
        stringsAsFactors = FALSE
      )
      modelingRows[[paste0(iTrait, "_selectedModel")]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = "selectedModel",
        value = model_type,
        stringsAsFactors = FALSE
      )
    }

    # One entry per input STA stamp: parameter="analysisId", trait="inputObject", environment="general"
    for (inputId in analysisId) {
      modelingRows[[paste0("input_", inputId)]] <- data.frame(
        module = "mta_rrblup",
        analysisId = mtaAnalysisId,
        trait = "inputObject",
        environment = "general",
        parameter = "analysisId",
        value = as.character(inputId),
        stringsAsFactors = FALSE
      )
    }

    if (length(modelingRows) > 0) {
      newModeling <- do.call(rbind, modelingRows)
      rownames(newModeling) <- NULL
      phenoDTfile$modeling <- rbind(
        phenoDTfile$modeling,
        newModeling[, colnames(phenoDTfile$modeling)]
      )
    }

    ## ─── 6. STATUS TABLE ─────────────────────────────────────────────────────
    # Append exactly one status entry
    newStatus <- data.frame(
      module = "mta_rrblup",
      analysisId = mtaAnalysisId,
      analysisIdName = NA_character_,
      stringsAsFactors = FALSE
    )
    phenoDTfile$status <- rbind(
      phenoDTfile$status,
      newStatus[, colnames(phenoDTfile$status)]
    )

    if (verbose) message("Output assembly complete. analysisId: ", mtaAnalysisId)

  } # end if mode == "fit" && length(fittingResults) > 0

  ## ───────────────────────────────────────────────────────────────────────────
## RETURN
## ───────────────────────────────────────────────────────────────────────────

  return(phenoDTfile)
}

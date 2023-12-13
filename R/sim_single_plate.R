#' Simulate a single plate
#'
#' @param plate_layout data frame similar to platemap
#' @param hit_genes number of hit genes
#' @param hit_genes_with_interaction number of hit genes
#'                               with interaction with compound
#' @param desired_cells_per_well numeric value desired number of cells per well
#' @param base_level numeric value baseline expression of the target gene
#' @param mu_bg average gene effect on expression, percentage
#' @param sigma_bg between gene variation in effect on target expression
#' @param mu_btg average gene treatment interaction, percentage
#' @param sigma_btg between gene variation in interaction
#' @param cv_cell between cell coefficient of variation
#' @param measure_error measurement error log scale
#' @param spatial_effect data frame wells with three columns, Row, Column
#'                       and spatial bias size
#' @param sigma_bg0 between gene variation among the true negative gene editing
#'
#' @return data frame
#'
#' @importFrom rlang .data
#' @export
sim_single_plate <- function(plate_layout,
                             hit_genes,
                             hit_genes_with_interaction,
                             desired_cells_per_well,
                             base_level,
                             mu_bg,
                             sigma_bg,
                             mu_btg,
                             sigma_btg,
                             cv_cell,
                             measure_error,
                             spatial_effect,
                             sigma_bg0) {

  if (length(hit_genes_with_interaction) > length(hit_genes)) {
    stop("It is impossible to have more hit genes
         with interactions than the number of hit genes")
  }

  # plate information ----
  ## some wells have no gene perturbation, gene = NA
  my_genes <- unique(plate_layout$gene[!is.na(plate_layout$gene)])
  n_genes <- length(my_genes)
  n_compounds <- length(unique(plate_layout$compound))
  n_wells <- nrow(plate_layout)
  ## simulate the number of cells in each cell
  n_cells_per_well <- stats::rpois(n_wells, lambda = desired_cells_per_well)

  # true gene effect ----
  mu_bg <- mu_bg
  sigma_bg <- sigma_bg

  mu_btg <- mu_btg
  sigma_btg <- sigma_btg

  # mean endpoint per well: base.level + gKO.effect in log scale
  bg <- vector("numeric", length = n_genes)
  btg <- vector("numeric", length = n_genes)
  for (i in 1:n_genes) {
    if (my_genes[i] %in% hit_genes_with_interaction) {
      bg[i] <- stats::rnorm(1, mean = mu_bg, sd = sigma_bg)
      btg[i] <- stats::rnorm(1, mean = mu_btg, sd = sigma_btg)
    } else if (!(my_genes[i] %in% hit_genes_with_interaction) &&
                 (my_genes[i] %in% hit_genes)) {
      bg[i] <- stats::rnorm(1, mean = mu_bg, sd = sigma_bg)
      btg[i] <- 0
    } else {
      # gene editing effect of the true negative
      # there is a between gene variation
      bg[i] <- stats::rnorm(1, mean = 0, sd = sigma_bg0)
      btg[i] <- 0
    }
  }

  gene_effect <- plate_layout %>%
    dplyr::filter(!is.na(.data$gene)) %>%
    dplyr::select(.data$gene) %>%
    dplyr::distinct() %>%
    dplyr::mutate(bg = bg, btg = btg) %>%
    dplyr::bind_rows(data.frame(gene = NA, bg = 0, btg = 0)) %>%
    dplyr::mutate(hit = ifelse(.data$gene %in% hit_genes, 1L, 0L),
                  hit_with_interaction =
                    ifelse(.data$gene %in% hit_genes_with_interaction, 1L, 0L))

  # true compound effect ----
  bt <- stats::rnorm(n_compounds, mean = 0, sd = 0.0)

  compound_effect <- plate_layout %>%
    dplyr::select(.data$compound) %>%
    dplyr::distinct() %>%
    dplyr::mutate(bt = bt,
                  bt = ifelse(.data$compound == "DMSO", 0, bt))

  # true signal ----
  my_truth <- plate_layout %>%
    dplyr::mutate(a = base_level,
                  G = ifelse(is.na(.data$gene), 0L, 1L),
                  T = ifelse(.data$compound == "DMSO", 0L, 1L)) %>%
    dplyr::left_join(gene_effect, by = "gene") %>%
    dplyr::left_join(compound_effect, by = "compound") %>%
    dplyr::mutate(mu_per_cell = exp(.data$a + .data$bg * .data$G +
                                      (.data$bt +
                                         .data$btg * .data$G) * .data$T)) %>%
    dplyr::mutate(n_cells = n_cells_per_well) %>%
    dplyr::mutate(sigma_cell = .data$mu_per_cell * cv_cell,
                  x_per_cell = purrr::pmap(list(.data$n_cells,
                                                .data$mu_per_cell,
                                                .data$sigma_cell),
                                           ~stats::rnorm(..1,
                                                         mean = ..2,
                                                         sd = ..3)))

  my_sim <- my_truth %>%
    tidyr::unnest(.data$x_per_cell) %>%
    dplyr::full_join(spatial_effect, by = c("Row", "Column")) %>%
    dplyr::mutate(spatial_bias = ifelse(is.na(.data$spatial_bias),
                                        0,
                                        .data$spatial_bias),
                  # the signal will be measured in each cell and
                  # the systematic bias will influence all of them
                  # it is additive at original scale,
                  # x_per_cell is on original scale
                  z_per_cell = log(.data$x_per_cell + .data$spatial_bias),
                  y_per_cell =
                    purrr::map_dbl(.data$z_per_cell,
                                   ~stats::rlnorm(1,
                                                  meanlog = .x,
                                                  sdlog = measure_error))) %>%
    # the mean intensity value is derived from sum up all cells' intensity and
    # then divide the number of cells
    dplyr::group_by(.data$Row, .data$Column, .data$gene, .data$compound,
                    .data$hit, .data$hit_with_interaction, .data$a, .data$bg,
                    .data$btg, .data$bt, .data$n_cells, .data$spatial_bias) %>%
    dplyr::summarize(y_all_cells = sum(.data$y_per_cell)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(y = .data$y_all_cells / .data$n_cells)

  return(my_sim)
}

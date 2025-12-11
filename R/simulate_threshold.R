#' @title Simulating Dynamic Treatment Intervention (With Duration Control)
#' @description \code{simulate_network_threshold} simulates an intervention where the thresholds of targeted symptoms are directly manipulated. It allows defining how long the intervention lasts.
#'
#' @param W_init A square matrix representing the initial weighted connections.
#' @param b_init A numeric vector representing the initial thresholds.
#' @param target A numeric vector matching the number of symptoms. Values indicate the magnitude of intervention added to each symptom's threshold.
#' @param connectivity A numeric value (default: 1) controlling overall connection strength.
#' @param trial An integer (default: 10) specifying the number of simulation trials.
#' @param baseline_iteration An integer (default: 10) steps before intervention.
#' @param intervention_duration An integer (default: 20). How many steps the intervention lasts before being removed.
#' @param num_updates An integer (default: 5). How many times the network should be re-estimated after intervention.
#' @param update_interval An integer (default: 10). The number of simulation steps between each network re-estimation.
#' @param follow_up_iteration An integer (default: 10) steps after the dynamic update phase.
#' @param symptom_name A character vector for symptom names.
#'
#' @return A list containing plot_data, raw_data, final_network, network_history_full, result_plot, and result_text.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef sd rnorm
#' @importFrom qgraph qgraph
#' @importFrom tibble tibble
#' @importFrom dplyr mutate left_join
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_fill_manual scale_color_manual scale_y_continuous theme_classic theme element_text labs ggplotGrob geom_vline
#' @importFrom cowplot ggdraw draw_image
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off tiff
#' @importFrom scales squish
#' @export

simulate_network_threshold <- function(W_init,
                                       b_init,
                                       target,
                                       connectivity = 1,
                                       trial = 10,
                                       baseline_iteration = 10,
                                       intervention_duration = 20, #介入の持続期間
                                       num_updates = 5,
                                       update_interval = 10,
                                       follow_up_iteration = 10,
                                       symptom_name = NULL) {

  # --- 入力チェック ---
  if (nrow(W_init) != ncol(W_init)) {
    stop("Error: W_init must be a square matrix.")
  }
  if (length(b_init) != nrow(W_init)) {
    stop(paste0("Error: Length of b_init (", length(b_init), ") must match the number of symptoms in W_init (", nrow(W_init), ")."))
  }
  if (length(target) != nrow(W_init)) {
    stop(paste0("Error: Length of target (", length(target), ") must match the number of symptoms in W_init (", nrow(W_init), ")."))
  }

  # --- 内部関数: ロジスティック関数 ---
  sigmoid <- function(x) {
    x <- pmax(pmin(x, 700), -700)
    1 / (1 + exp(x))
  }

  # --- 内部関数: ネットワーク推定 ---
  estimate_symptom_network <- function(X_series) {
    P <- ncol(X_series)
    T_time <- nrow(X_series)

    Y_data <- as.matrix(X_series[2:T_time, ])
    Z_data <- as.matrix(X_series[1:(T_time - 1), ])

    B_temporal <- matrix(0, nrow = P, ncol = P)
    Alpha_threshold <- numeric(P)

    for (i in 1:P) {
      Y_i <- Y_data[, i]
      if (stats::sd(Y_i) == 0) {
        Alpha_threshold[i] <- 0
        next
      }
      tryCatch({
        n_folds <- min(10, length(Y_i))
        if(n_folds < 3) n_folds <- length(Y_i)

        if(nrow(Z_data) < 3) {
          Alpha_threshold[i] <- 0
          next
        }
        fit <- suppressWarnings(glmnet::cv.glmnet(x = Z_data, y = Y_i, family = "binomial", alpha = 1, nfolds = n_folds))
        coefs <- stats::coef(fit, s = "lambda.min")
        coefs_vec <- as.vector(coefs)
        Alpha_threshold[i] <- coefs_vec[1]
        B_temporal[i, ] <- coefs_vec[2:(P + 1)]
      }, error = function(e) {})
    }
    return(list(B = B_temporal, Alpha = Alpha_threshold))
  }

  # --- 内部関数: 1試行分のシミュレーション ---
  simulate_network <- function(W_init, b_init, target, connectivity,
                               baseline_iteration, intervention_duration,
                               num_updates, update_interval, follow_up_iteration, symptom_name) {

    W <- W_init
    b <- b_init
    num_symptom <- nrow(W)

    history_list <- list()

    if (is.null(symptom_name)) {
      symptom_name <- letters[1:num_symptom]
    }
    colnames(W) <- symptom_name
    rownames(W) <- symptom_name
    names(b) <- symptom_name

    history_list[[1]] <- list(step_index = 0, phase = "Initial", W = W, b = b)

    X <- matrix(numeric(num_symptom), 1, num_symptom)
    total_steps <- baseline_iteration + (num_updates * update_interval) + follow_up_iteration
    X_matrix_full <- matrix(0, nrow = total_steps, ncol = num_symptom)

    # 介入終了のタイミングを計算
    intervention_end_step <- baseline_iteration + intervention_duration

    safe_sample <- function(prob) {
      if (is.na(prob)) prob <- 0.5
      if (prob < 0) prob <- 0
      if (prob > 1) prob <- 1
      sample(x = c(1, 0), size = 1, replace = T, prob = c(prob, 1 - prob))
    }

    current_step <- 0

    # --- Phase 1: Baseline ---
    for (s in 1:baseline_iteration) {
      current_step <- current_step + 1
      A <- numeric(num_symptom)
      P_prob <- numeric(num_symptom)
      for (j in 1:num_symptom) {
        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      X_matrix_full[current_step, ] <- X
    }

    # --- Start Intervention ---
    # ここで閾値を上げる
    b <- b + target
    history_list[[2]] <- list(step_index = current_step, phase = "Intervention_Start", W = W, b = b)

    # --- Phase 2: Dynamic Updates ---
    for (i in 1:num_updates) {
      X_history_current <- matrix(0, nrow = update_interval, ncol = num_symptom)
      for (k in 1:update_interval) {
        current_step <- current_step + 1

        # === 介入終了チェック ===
        if (current_step == intervention_end_step) {
          b <- b - target # 閾値を下げる（介入除去）
          # 履歴に記録
          history_list[[length(history_list) + 1]] <- list(step_index = current_step, phase = "Intervention_End", W = W, b = b)
        }
        # ========================

        A <- numeric(num_symptom)
        P_prob <- numeric(num_symptom)
        for (j in 1:num_symptom) {
          input_val <- sum(connectivity * W[, j] * X)
          A[j] <- input_val
          exponent_val <- abs(b[j]) - A[j]
          exponent_val <- max(min(exponent_val, 100), -100)
          P_prob[j] <- 1 / (1 + exp(exponent_val))
          X[j] <- safe_sample(P_prob[j])
        }
        X_history_current[k, ] <- X
        X_matrix_full[current_step, ] <- X
      }

      # ネットワーク再推定
      symptom_data <- X_history_current
      if (stats::sd(as.vector(symptom_data)) > 0) {
        estimated <- estimate_symptom_network(symptom_data)
        B_new <- estimated$B
        Alpha_new <- estimated$Alpha
        B_new[is.na(B_new)] <- 0
        Alpha_new[is.na(Alpha_new)] <- 0

        W <- t(B_new)
        b <- Alpha_new
        colnames(W) <- symptom_name
        rownames(W) <- symptom_name
        names(b) <- symptom_name

        history_list[[length(history_list) + 1]] <- list(step_index = current_step, phase = paste0("Update_", i), W = W, b = b)
      } else {
        history_list[[length(history_list) + 1]] <- list(step_index = current_step, phase = paste0("Update_", i, "_NoChange"), W = W, b = b)
      }
    }

    # --- Phase 3: Follow-up ---
    for (k in 1:follow_up_iteration) {
      current_step <- current_step + 1

      # === 介入終了チェック (期間が長くFollow-up中に終了する場合用) ===
      if (current_step == intervention_end_step) {
        b <- b - target
        history_list[[length(history_list) + 1]] <- list(step_index = current_step, phase = "Intervention_End", W = W, b = b)
      }
      # ========================

      A <- numeric(num_symptom)
      P_prob <- numeric(num_symptom)
      for (j in 1:num_symptom) {
        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      X_matrix_full[current_step, ] <- X
    }
    return(list(X_full = X_matrix_full, W = W, b = b, history = history_list))
  }

  # --- メイン処理 ---
  target <- unlist(target)
  total_time <- baseline_iteration + num_updates * update_interval + follow_up_iteration
  initial_num_symptoms <- nrow(W_init)

  if (is.null(symptom_name)) {
    if (!is.null(colnames(W_init))) {
      symptom_name <- colnames(W_init)
    } else {
      symptom_name <- paste0("Symptom_", 1:initial_num_symptoms)
    }
  }

  raw_data_array <- array(0, dim = c(trial, total_time, initial_num_symptoms))
  final_W_array <- array(0, dim = c(trial, initial_num_symptoms, initial_num_symptoms))
  final_b_matrix <- matrix(0, nrow = trial, ncol = initial_num_symptoms)
  network_history_full <- vector("list", trial)
  names(network_history_full) <- paste0("Trial_", 1:trial)

  for (l in 1:trial) {
    result <- simulate_network(W_init, b_init, target, connectivity,
                               baseline_iteration, intervention_duration,
                               num_updates, update_interval, follow_up_iteration, symptom_name)
    raw_data_array[l, , ] <- result$X_full
    final_W_array[l, , ] <- result$W
    final_b_matrix[l, ] <- result$b
    network_history_full[[l]] <- result$history
  }

  # --- Plot Data Preparation ---
  mean_traj <- apply(raw_data_array, c(2, 3), mean, na.rm = TRUE)
  sd_traj <- apply(raw_data_array, c(2, 3), stats::sd, na.rm = TRUE)
  colnames(mean_traj) <- symptom_name
  colnames(sd_traj) <- symptom_name

  df_mean <- as.data.frame(mean_traj) |>
    dplyr::mutate(time = 1:total_time) |>
    tidyr::pivot_longer(cols = -time, names_to = "symptom", values_to = "mean_val")
  df_sd <- as.data.frame(sd_traj) |>
    dplyr::mutate(time = 1:total_time) |>
    tidyr::pivot_longer(cols = -time, names_to = "symptom", values_to = "sd_val")
  plot_data <- dplyr::left_join(df_mean, df_sd, by = c("time", "symptom"))

  mean_final_W <- apply(final_W_array, c(2, 3), mean, na.rm=TRUE)
  mean_final_b <- colMeans(final_b_matrix, na.rm=TRUE)
  colnames(mean_final_W) <- symptom_name
  rownames(mean_final_W) <- symptom_name

  # --- 最終試行のデータ抽出 ---
  last_trial_W <- final_W_array[trial, , ]
  colnames(last_trial_W) <- symptom_name
  rownames(last_trial_W) <- symptom_name

  last_trial_data <- raw_data_array[trial, , ]
  last_trial_mean_traj <- rowMeans(last_trial_data)
  last_trial_sd_traj <- apply(last_trial_data, 1, stats::sd)

  # --- プロット作成 ---
  v_fixed <- 8
  target_indices <- which(target != 0)

  if (length(target_indices) > 0) {
    targeted_names <- symptom_name[target_indices]
    plot_title <- paste("Final Trial Structure\nIntervention on:", paste(targeted_names, collapse = ", "))
  } else {
    plot_title <- "Final Trial Structure\nNo Intervention"
  }

  color_node <- rep("lightblue", initial_num_symptoms)
  color_node[target_indices] <- "orange"

  temp_file_path <- file.path(tempdir(), "tmp_net.tiff")
  grDevices::tiff(temp_file_path, width = 1200, height = 800, res = 400)
  tryCatch({
    qgraph::qgraph(last_trial_W, layout = "spring", theme = "colorblind", vsize = v_fixed,
                   color = color_node, posCol = "blue", negCol = "red", negDashed = F, title = plot_title)
  }, error = function(e) {
    plot(1, type="n", axes=F, xlab="", ylab="", main="Network Plot Error")
  })
  grDevices::dev.off()

  # 介入終了タイミング
  intervention_end_x <- baseline_iteration + intervention_duration

  p2 <- tibble::tibble(time = 1:total_time, mean = last_trial_mean_traj, sd = last_trial_sd_traj) |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = mean)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, fill = "blue") +
    # 介入開始ライン
    ggplot2::geom_vline(xintercept = baseline_iteration + 0.5, linetype="dashed", color="gray50") +
    # 介入終了ライン (もし全期間内であれば表示)
    {if(intervention_end_x < total_time) ggplot2::geom_vline(xintercept = intervention_end_x + 0.5, linetype="dotted", color="red")} +
    ggplot2::scale_y_continuous(limits = c(0.0, 1.0), oob = scales::squish) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Time Step", y = "Activation", title = "Global Response (Final Trial)",
                  caption = "Dashed: Intervention Start | Dotted: Intervention End")

  qgraph_image <- cowplot::ggdraw() + cowplot::draw_image(temp_file_path)
  qgraph_image_grob <- ggplot2::ggplotGrob(qgraph_image)
  result_plot <- gridExtra::grid.arrange(qgraph_image_grob, p2, ncol = 2)

  result_text <- paste0(sprintf("Global activation at final step (Last Trial): %.3f", last_trial_mean_traj[total_time]))

  return(list(plot_data = plot_data, raw_data = raw_data_array,
              final_network = list(W = mean_final_W, b = mean_final_b),
              network_history_full = network_history_full,
              result_plot = result_plot, result_text = result_text))
}

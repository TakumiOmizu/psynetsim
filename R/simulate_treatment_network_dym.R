#' @title Simulating Dynamic Treatment Intervention on Symptom Network
#' @description \code{simulate_treatment_network_dym} simulates treatment intervention on a symptom network, where **symptom network connections dynamically update** based on recent patient state data using regularized logistic regression (Lasso/Elastic Net), capturing how the system responds to treatment over time. This functionality is an extension of the base `simulate_treatment_network` function.
#'
#' @param W_init A square matrix representing the initial weighted connections (edges) between symptoms. The number of rows/columns determines the number of symptoms in the network.
#' @param b_init A numeric vector representing the initial thresholds (biases) for each symptom. Its length must match the number of symptoms in `W_init`.
#' @param target A list or vector indicating which symptoms are targeted by the treatment. Use 1 for targeted symptoms and 0 for non-targeted symptoms. Its length must match the number of symptoms.
#' @param connectivity A numeric value (default: 1) controlling the overall strength of connections within the network. Higher values increase the influence of connected nodes.
#' @param edge_between_TC A numeric value (default: 1) controlling the strength of connections between different treatment components (TCs).
#' @param weight_bias A numeric value (default: 1) influencing the weight of connections between treatment components and targeted symptoms. A value greater than 1 increases the negative (inhibitory) effect.
#' @param TB A numeric value (default: 1) influencing the threshold (bias) of treatment components. Higher values make TCs less likely to activate spontaneously.
#' @param trial An integer (default: 10) specifying the number of simulation trials to run. Results are averaged across trials.
#' @param baseline_iteration An integer (default: 10) specifying the number of iterations for the baseline (pre-treatment) simulation phase.
#' @param num_TC An integer (default: 5) specifying the total number of treatment components to be introduced.
#' @param TC_iteration_per_component An integer (default: 10) specifying the number of iterations after each TC introduction for the network to adapt, before re-estimation.
#' @param follow_up_iteration An integer (default: 10) specifying the number of iterations for the follow-up phase.
#' @param symptom_name A character vector (default: NULL) providing custom names for the symptoms. If NULL, symptoms will be named alphabetically.
#'
#' @return A list containing:
#'   \item{result_plot}{A `cowplot` object displaying two plots: a network visualization of the *initial* symptom connections (which are dynamically updated during simulation) and a time-series plot showing the mean and standard deviation of symptom and treatment component activation.}
#'   \item{result_text}{A character string summarizing the mean and standard deviation of symptom and treatment component activation at the final simulation step.}
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @importFrom qgraph qgraph centrality
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_fill_manual scale_color_manual scale_y_continuous theme_classic theme element_text labs ggplotGrob
#' @importFrom cowplot ggdraw draw_image
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off tiff
#' @importFrom stats rnorm sd
#' @importFrom scales squish
#' @export

simulate_treatment_network_dym <- function(W_init,
                                           b_init,
                                           target,
                                           connectivity = 1,
                                           edge_between_TC = 1,
                                           weight_bias = 1,
                                           TB = 1,
                                           trial = 10,
                                           baseline_iteration = 10,
                                           num_TC = 5,
                                           TC_iteration_per_component = 10,
                                           follow_up_iteration = 10,
                                           symptom_name = NULL) {



  # --- 内部関数: ロジスティック関数 ---
  sigmoid <- function(x) {
    x <- pmax(pmin(x, 700), -700)
    1 / (1 + exp(x))
  }

  # --- 内部関数: 症状ネットワークの推定 (正則化ロジスティック回帰) ---
  estimate_symptom_network <- function(X_series) {
    P <- ncol(X_series)
    T_time <- nrow(X_series)

    Y_data <- as.matrix(X_series[2:T_time, ])
    Z_data <- as.matrix(X_series[1:(T_time - 1), ])

    B_temporal <- matrix(0, nrow = P, ncol = P)
    Alpha_threshold <- numeric(P)

    for (i in 1:P) {
      Y_i <- Y_data[, i]

      if (sd(Y_i) == 0) {
        Alpha_threshold[i] <- 0
        next
      }

      tryCatch({
        n_folds <- min(10, length(Y_i))
        if(n_folds < 3) n_folds <- length(Y_i)

        fit <- suppressWarnings(cv.glmnet(x = Z_data, y = Y_i, family = "binomial", alpha = 1, nfolds = n_folds))

        coefs <- coef(fit, s = "lambda.min")
        coefs_vec <- as.vector(coefs)

        Alpha_threshold[i] <- coefs_vec[1]
        B_temporal[i, ] <- coefs_vec[2:(P + 1)]

      }, error = function(e) {
      })
    }

    return(list(B = B_temporal, Alpha = Alpha_threshold))
  }

  # --- 内部関数: 1試行分のシミュレーション ---
  simulate_network <- function(W_init, b_init, target, connectivity, edge_between_TC,
                               weight_bias, TB, baseline_iteration, num_TC,
                               TC_iteration_per_component, follow_up_iteration, symptom_name) {

    W <- W_init
    b <- b_init

    # TCの重み生成のための平均・分散計算用
    if (isSymmetric(as.matrix(W))) {
      # 無向グラフ: 上三角成分のみ取得（重複回避）
      edges <- W[row(W) < col(W)]
    } else {
      # 有向グラフ: 対角成分以外の非ゼロ要素を全て取得
      edges <- W[row(W) != col(W) & W != 0]
    }

    if(length(edges) == 0) {
      mean_weight <- 0
      sd_weight <- 0.1
    } else {
      mean_weight <- mean(edges, na.rm = TRUE)
      sd_weight <- sd(edges, na.rm = TRUE)
      if(is.na(sd_weight)) sd_weight <- 0.1
    }

    mean_threthold <- mean(b, na.rm = TRUE)
    sd_threthold <- sd(b, na.rm = TRUE)
    if(is.na(sd_threthold)) sd_threthold <- 0.1

    num_symptom <- nrow(W)
    X <- matrix(rep(0, length(b)), 1, length(b))

    D <- numeric(0)
    TC <- numeric(0)

    if (is.null(symptom_name)) {
      symptom_name <- letters[1:num_symptom]
    }

    safe_sample <- function(prob) {
      if (is.na(prob)) prob <- 0.5
      if (prob < 0) prob <- 0
      if (prob > 1) prob <- 1
      sample(x = c(1, 0), size = 1, replace = T, prob = c(prob, 1 - prob))
    }

    # --- Baseline ---
    X_history_baseline <- matrix(0, nrow = baseline_iteration, ncol = length(b))

    for (s in 1:baseline_iteration) {
      A <- numeric(length(b))
      P_prob <- numeric(length(b))

      for (j in 1:length(b)) {

        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      X_history_baseline[s, ] <- X
      D[s] <- sum(X[1:num_symptom]) / num_symptom
      TC[s] <- 0
    }

    last_phase_data <- X_history_baseline

    # --- TC Loop ---
    for (i in 1:num_TC) {

      #  次のTCが入る直前までの各症状の時系列データからネットワーク更新
      if (i > 1 && any(target == 1)) {
        symptom_data <- last_phase_data[, 1:num_symptom]

        if (sd(as.vector(symptom_data)) > 0) {
          estimated <- estimate_symptom_network(symptom_data)
          B_new <- estimated$B
          Alpha_new <- estimated$Alpha

          B_new[is.na(B_new)] <- 0
          Alpha_new[is.na(Alpha_new)] <- 0

          W[1:num_symptom, 1:num_symptom] <- t(B_new)
          b[1:num_symptom] <- Alpha_new
        }
      }

      # Wの拡張
      current_nodes_count <- nrow(W)
      new_dim <- current_nodes_count + 1
      W_expanded <- matrix(0, nrow = new_dim, ncol = new_dim)
      if (current_nodes_count > 0) {
        W_expanded[1:current_nodes_count, 1:current_nodes_count] <- W
      }
      W <- W_expanded

      # TCを入れるとtargetと負の相関を持つようになる
      for (j in 1:(num_symptom + i)) {
        if (j <= num_symptom) {
          if (target[j] == 1) {
            mw <- if(is.na(mean_weight)) 0 else mean_weight
            sw <- if(is.na(sd_weight)) 0.1 else sd_weight

            temp_weight <- -abs(rnorm(1, mean_weight, sw)) * weight_bias

            W[num_symptom + i, j] <- temp_weight
            W[j, num_symptom + i] <- temp_weight
          }
        } else if (j > num_symptom && j < (num_symptom + i)) {
          mw <- if(is.na(mean_weight)) 0 else mean_weight
          sw <- if(is.na(sd_weight)) 0.1 else sd_weight
          temp_weight <- rnorm(1, edge_between_TC * mw, sw)
          W[num_symptom + i, j] <- temp_weight
          W[j, num_symptom + i] <- temp_weight
        }
      }

      mt <- if(is.na(mean_threthold)) 0 else mean_threthold
      st <- if(is.na(sd_threthold)) 0.1 else sd_threthold
      b <- c(b, rnorm(1, mt * TB, st))

      TC_name <- paste0("TC", 1:i)
      colnames(W) <- c(symptom_name, TC_name)
      rownames(W) <- c(symptom_name, TC_name)

      X <- cbind(X, c(0))

      X_history_current <- matrix(0, nrow = TC_iteration_per_component, ncol = length(b))

      for (k in 1:TC_iteration_per_component) {
        A <- numeric(length(b))
        P_prob <- numeric(length(b))
        for (j in 1:length(b)) {
          input_val <- sum(connectivity * W[, j] * X)
          A[j] <- input_val
          exponent_val <- abs(b[j]) - A[j]
          exponent_val <- max(min(exponent_val, 100), -100)
          P_prob[j] <- 1 / (1 + exp(exponent_val))
          X[j] <- safe_sample(P_prob[j])
        }
        X_history_current[k, ] <- X
        tmp_D <- sum(X[1:num_symptom]) / num_symptom
        tmp_TC <- sum(X[(num_symptom + 1):(num_symptom + i)]) / i
        time_no <- baseline_iteration + (i - 1) * TC_iteration_per_component + k
        D[time_no] <- tmp_D
        TC[time_no] <- tmp_TC
      }
      last_phase_data <- X_history_current
    }

    # --- Follow-up ---
    for (k in 1:follow_up_iteration) {
      A <- numeric(length(b))
      P_prob <- numeric(length(b))
      for (j in 1:length(b)) {
        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      tmp_D <- sum(X[1:num_symptom]) / num_symptom
      tmp_TC <- sum(X[(num_symptom + 1):(num_symptom + i)]) / i
      time_no <- length(D) + 1
      D[time_no] <- tmp_D
      TC[time_no] <- tmp_TC
    }

    return(list(D = D, TC = TC, W = W))
  }

  # --- メイン処理 ---
  target <- unlist(target)
  total_time <- baseline_iteration + num_TC * TC_iteration_per_component + follow_up_iteration
  D_iteration <- matrix(0, trial, total_time)
  TC_iteration <- matrix(0, trial, total_time)

  initial_num_symptoms <- nrow(W_init)
  final_nodes <- initial_num_symptoms + num_TC
  W_iteration <- matrix(0, trial, final_nodes^2)

  # 症状名の取得
  if (is.null(symptom_name)) {
    if (!is.null(colnames(W_init))) {
      symptom_name <- colnames(W_init)
    } else {
      symptom_name <- letters[1:initial_num_symptoms]
    }
  }

  for (l in 1:trial) {
    result <- simulate_network(W_init, b_init, target, connectivity, edge_between_TC,
                               weight_bias, TB, baseline_iteration, num_TC,
                               TC_iteration_per_component, follow_up_iteration,
                               symptom_name)
    D_iteration[l, ] <- result$D
    TC_iteration[l, ] <- result$TC
    W_iteration[l, ] <- as.vector(result$W)
  }

  W_names_full <- colnames(result$W)
  D_plot <- colMeans(D_iteration, na.rm=TRUE)
  TC_plot <- colMeans(TC_iteration, na.rm=TRUE)
  D_sd <- apply(D_iteration, 2, sd, na.rm=TRUE)
  TC_sd <- apply(TC_iteration, 2, sd, na.rm=TRUE)

  # --- プロット用データの整形（症状のみ抽出） ---
  W_full_mean <- matrix(colMeans(W_iteration, na.rm=TRUE), final_nodes, final_nodes)
  colnames(W_full_mean) <- W_names_full
  rownames(W_full_mean) <- W_names_full

  W_symptoms_only <- W_full_mean[1:initial_num_symptoms, 1:initial_num_symptoms]

  color_node <- rep("lightblue", initial_num_symptoms)

  # 固定サイズ設定
  v_fixed <- 8

  # タイトル設定
  target_indices <- which(target == 1)
  if (length(target_indices) > 0) {
    targeted_names <- symptom_name[target_indices]
    plot_title <- paste("Intervention on:", paste(targeted_names, collapse = ", "))
  } else {
    plot_title <- "No Intervention"
  }

  temp_file_path <- file.path(tempdir(), "tmp.tiff")
  tiff(temp_file_path, width = 1200, height = 800, res = 400)
  tryCatch({
    qgraph(W_symptoms_only, layout = "spring", theme = "colorblind",
           vsize = v_fixed,
           color = color_node, posCol = "blue",
           negCol = "red", negDashed = F,
           title = plot_title)
  }, error = function(e) {
    plot(1, type="n", axes=F, xlab="", ylab="", main="Network Plot Error")
  })
  dev.off()

  p2 <- tibble(time = 1:total_time,
               D_plot = D_plot,
               TC_plot = TC_plot,
               D_sd = D_sd,
               TC_sd = TC_sd) |>
    pivot_longer(cols = c(D_plot, TC_plot), names_to = "variable", values_to = "mean") |>
    mutate(sd = ifelse(variable == "D_plot", D_sd, TC_sd)) |>
    ggplot(aes(x = time, y = mean, color = variable, fill = variable)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    scale_fill_manual(name = "Variable", values = c("D_plot" = "blue", "TC_plot" = "orange")) +
    scale_color_manual(name = "Variable", values = c("D_plot" = "blue", "TC_plot" = "orange")) +
    scale_y_continuous(limits = c(0.0, 1.0), oob = scales::squish) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = "", y = "", title = "")

  qgraph_image <- ggdraw() + draw_image(temp_file_path)
  qgraph_image_grob <- ggplotGrob(qgraph_image)
  result_plot <- grid.arrange(qgraph_image_grob, p2, ncol = 2)

  result_text <- paste0(
    sprintf("The mean value of symptom at the final step(t=%d). = %.3f", total_time, D_plot[total_time]), "\n",
    sprintf("The mean value of treatment component at the final step(t=%d). = %.3f", total_time, TC_plot[total_time]), "\n",
    sprintf("The SD value of symptom at the final step(t=%d). = %.3f", total_time, D_sd[total_time]), "\n",
    sprintf("The SD value of treatment component at the final step(t=%d). = %.3f", total_time, TC_sd[total_time])
  )
  return(list(result_plot = result_plot, result_text = result_text))
}

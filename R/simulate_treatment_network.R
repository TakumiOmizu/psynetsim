#' @title Simulating treatment intervention to symptom network
#' @description \code{simulate_treatment_network} simulates the effect of treatment interventions on a symptom network. This package, **ROXGAN**, provides tools for understanding and visualizing the dynamic impact of interventions within complex symptom interdependencies.
#'
#' @param W_init A square matrix representing the initial weighted connections (edges) between symptoms. The number of rows/columns determines the number of symptoms in the network.
#' @param b_init A numeric vector representing the initial thresholds (biases) for each symptom. Its length must match the number of symptoms in `W_init`.
#' @param target A list or vector indicating which symptoms are targeted by the treatment. Use 1 for targeted symptoms and 0 for non-targeted symptoms. Its length must match the number of symptoms.
#' @param connectivity A numeric value (default: 1) controlling the overall strength of connections within the network. Higher values increase the influence of connected nodes.
#' @param edge_between_TC A numeric value (default: 1) controlling the strength of connections between different treatment components (TCs).
#' @param weight_bias A numeric value (default: 1) influencing the weight of connections between treatment components and targeted symptoms. A value greater than 1 increases the negative (inhibitory) effect, while less than 1 decreases it.
#' @param TB A numeric value (default: 1) influencing the threshold (bias) of treatment components. Higher values make TCs less likely to activate spontaneously.
#' @param trial An integer (default: 10) specifying the number of simulation trials to run. The final results are averaged across these trials.
#' @param baseline_iteration An integer (default: 10) specifying the number of iterations for the baseline (pre-treatment) simulation phase.
#' @param num_TC An integer (default: 5) specifying the total number of treatment components to be introduced during the simulation.
#' @param TC_iteration_per_component An integer (default: 10) specifying the number of iterations for each treatment component introduction phase.
#' @param follow_up_iteration An integer (default: 10) specifying the number of iterations for the follow-up (post-treatment) simulation phase.
#' @param symptom_name A character vector (default: NULL) providing custom names for the symptoms. If NULL, symptoms will be named alphabetically (e.g., "a", "b", "c", ...). Its length must match the number of symptoms.
#'
#' @return A list containing:
#'   \item{result_plot}{A `cowplot` object displaying two plots: a network visualization of the average weights (symptom and treatment components) and a time-series plot showing the mean and standard deviation of symptom and treatment component activation over the simulation period.}
#'   \item{result_text}{A character string summarizing the mean and standard deviation of symptom and treatment component activation at the final simulation step.}
#' @examples
#' # Install and load necessary packages if not already installed
#' # install.packages(c("qgraph", "tidyverse", "cowplot", "gridExtra"))
#' library(qgraph)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(ggplot2)
#' library(cowplot)
#' library(gridExtra)
#'
#' # Example data for a 6-symptom network
#' set.seed(456)
#' weight_6 <- matrix(rnorm(6*6, mean = 0.2, sd = 0.08), nrow = 6, ncol = 6)
#' diag(weight_6) <- 0
#' weight_6[upper.tri(weight_6)] <- t(weight_6)[upper.tri(weight_6)]
#' threshold_6 <- data.frame(threshold = rnorm(6, mean = 0.3, sd = 0.05))
#' target_list_6 <- list(symptom1 = 1, symptom2 = 0, symptom3 = 1,
#'                       symptom4 = 0, symptom5 = 0, symptom6 = 1)
#' custom_symptom_names_6 <- c("Anxiety", "Sadness", "Fatigue",
#'                             "Insomnia", "Irritability", "Pain")
#'
#' # Run the simulation with custom parameters
#' sim_results <- simulate_treatment_network(
#'   W_init = weight_6,
#'   b_init = threshold_6$threshold,
#'   target = target_list_6,
#'   connectivity = 1.2,
#'   edge_between_TC = 0.8,
#'   weight_bias = 1.2,
#'   TB = 0.8,
#'   trial = 5, # Example: Overriding default 10
#'   baseline_iteration = 15, # Example: Overriding default 10
#'   num_TC = 4, # Example: Overriding default 5
#'   TC_iteration_per_component = 12, # Example: Overriding default 10
#'   follow_up_iteration = 15, # Example: Overriding default 10
#'   symptom_name = custom_symptom_names_6
#' )
#'
#' # Print summary text and display plots
#' print(sim_results$result_text)
#' print(sim_results$result_plot)
#'
#' # Run the simulation with only required parameters (using all default iteration/trial values)
#' sim_results_default <- simulate_treatment_network(
#'   W_init = weight_6,
#'   b_init = threshold_6$threshold,
#'   target = target_list_6
#' )
#' print(sim_results_default$result_text)
#' print(sim_results_default$result_plot)
#'
#' @importFrom qgraph qgraph centrality
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot ggplotGrob aes geom_line geom_ribbon scale_fill_manual scale_color_manual scale_y_continuous theme_classic theme element_text labs
#' @importFrom cowplot ggdraw draw_image
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off tiff
#' @importFrom stats rnorm sd
#' @importFrom scales squish
#' @export
#' @export

simulate_treatment_network<-function(W_init,
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
                                     symptom_name = NULL){

  simulate_network <- function(W_init,
                               b_init,
                               target,
                               connectivity = 1,
                               edge_between_TC = 1,
                               weight_bias = 1,
                               TB = 1,
                               baseline_iteration = 1000,
                               num_TC = 20,
                               TC_iteration_per_component = 100,
                               follow_up_iteration = 1000,
                               symptom_name){

    W <- W_init
    b <- b_init

    edges <- W[row(W) < col(W)]

    mean_weight <- mean(edges)
    sd_weight <- sd(edges)
    mean_threthold <- mean(b)
    sd_threthold <- sd(b)
    # baseline simulation
    num_symptom <- nrow(W)
    X <- matrix(rep(0, length(b)), 1, length(b))
    D <- NULL
    TC <- NULL

    # default symptom_name
    if (is.null(symptom_name)) {
      symptom_name <- letters[1:num_symptom]
    } else if (length(symptom_name) != num_symptom) {
      stop("The length of symptom_name must match the number of symptoms (nrow(W_init)).")
    }

    for (s in 1:baseline_iteration) {
      A <- NULL
      P <- NULL
      for (j in 1:length(b)) {
        A[j] <- sum(connectivity*W[,j]*X)
        P[j] <- 1/(1+exp(abs(b[j])-A[j]))
        X[j] <- sample(x=c(1,0),size=1,replace = T,prob=c(P[j],1-P[j]))
      }
      tmp_D <- sum(X[1:num_symptom])/num_symptom
      tmp_TC <- 0
      D[s] <- tmp_D
      TC[s] <- tmp_TC
      time <- 1:s
    }

    for (i in 1:num_TC) {

      current_nodes_count <- nrow(W)
      new_dim <- current_nodes_count + 1
      W_expanded <- matrix(0, nrow = new_dim, ncol = new_dim)

      if (current_nodes_count > 0) {
        W_expanded[1:current_nodes_count, 1:current_nodes_count] <- W
      }
      W <- W_expanded

      for (j in 1:(num_symptom+i)) {
        if(j <= num_symptom){

          if(target[j] == 1){

            temp_weight <- rnorm(1, weight_bias*-mean_weight, sd_weight)
            W[num_symptom+i,j] <- temp_weight
            W[j, num_symptom+i] <- temp_weight
          }else{

            W[num_symptom+i,j] <- 0
            W[j, num_symptom+i] <- 0
          }
        }else if(j > num_symptom && j < (num_symptom+i)){

          temp_weight <- rnorm(1, edge_between_TC*mean_weight, sd_weight)
          W[num_symptom+i,j] <- temp_weight
          W[j, num_symptom+i] <- temp_weight
        }else{

          W[j, j] <- 0
        }
      }
      b <- c(b, rnorm(1, mean_threthold*TB, sd_threthold))
      TC_name <- NULL
      for (s_tc in 1:i) {
        TC_name <- c(TC_name, paste0("TC",s_tc))
      }
      colnames(W) <- c(symptom_name, TC_name)
      rownames(W) <- c(symptom_name, TC_name)
      X <- cbind(X,c(0))
      for (k in 1:TC_iteration_per_component) {
        A <- NULL
        P <- NULL
        for (j in 1:length(b)) {
          A[j] <- sum(connectivity*W[,j]*X)
          P[j] <- 1/(1+exp(abs(b[j])-A[j]))
          X[j] <- sample(x=c(1,0),size=1,replace = T,prob=c(P[j],1-P[j]))
        }

        tmp_D <- sum(X[1:num_symptom])/num_symptom
        tmp_TC <- sum(X[(num_symptom+1):(num_symptom+i)])/i
        time_no <- baseline_iteration + (i-1)*TC_iteration_per_component + k
        D[time_no] <- tmp_D
        TC[time_no] <- tmp_TC
        time <- 1:time_no
      }
    }
    # followup simulation
    for (k in 1:follow_up_iteration) {
      A <- NULL
      P <- NULL
      for (j in 1:length(b)) {
        A[j] <- sum(connectivity*W[,j]*X)
        P[j] <- 1/(1+exp(abs(b[j])-A[j]))
        X[j] <- sample(x=c(1,0),size=1,replace = T,prob=c(P[j],1-P[j]))
      }
      tmp_D <- sum(X[1:num_symptom])/num_symptom
      tmp_TC <- sum(X[(num_symptom+1):(num_symptom+i)])/i
      time_no <-length(D) + 1
      D[time_no] <- tmp_D
      TC[time_no] <- tmp_TC
      time <- 1:time_no
    }

    return(list(D = D, TC = TC, W = W))
  }



  target <- unlist(target)

  total_time <- baseline_iteration + num_TC * TC_iteration_per_component + follow_up_iteration
  D_iteration <- matrix(rep(0,total_time*trial),trial,total_time)
  TC_iteration <- matrix(rep(0,total_time*trial),trial,total_time)


  initial_num_symptoms <- nrow(W_init)
  W_iteration <- matrix(rep(0,((initial_num_symptoms + num_TC)^2)*trial),trial,(initial_num_symptoms + num_TC)^2)


  for (l in 1:trial) {
    # simulation
    result <- simulate_network(W_init,b_init,target,connectivity,edge_between_TC,
                               weight_bias,TB,baseline_iteration,num_TC,
                               TC_iteration_per_component,follow_up_iteration,
                               symptom_name)

    D_iteration[l,] <- result$D
    TC_iteration[l,] <- result$TC
    W_iteration[l,] <- as.vector(result$W)
  }

  W_names <- colnames(result$W)

  D_plot <- colMeans(D_iteration)
  TC_plot <- colMeans(TC_iteration)


  D_sd <- apply(D_iteration, 2, sd)
  TC_sd <- apply(TC_iteration, 2, sd)


  W_iteration_mean <- matrix(colMeans(W_iteration),(initial_num_symptoms + num_TC),(initial_num_symptoms + num_TC))
  colnames(W_iteration_mean) <- W_names
  rownames(W_iteration_mean) <- W_names

  color_node <- c(rep("lightblue", initial_num_symptoms),rep("orange", num_TC))
  Q <- qgraph(W_iteration_mean, DoNotPlot=TRUE)
  v <- (centrality(Q)$OutDegree/max(centrality(Q)$OutDegree))*10

  # plot
  temp_file_path <- file.path(tempdir(), "tmp.tiff")
  tiff(temp_file_path, width = 1200, height = 800, res = 400)
  qgraph(W_iteration_mean, layout="spring", theme="colorblind",vsize = v,color=color_node,posCol = "blue", negCol = "red",negDashed = F)
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
    scale_y_continuous(limits = c(0.0, 1.0),oob = scales::squish) +
    theme_classic()+
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
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

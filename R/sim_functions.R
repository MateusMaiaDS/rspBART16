# Just running the default values so I can go through the function and debug all
#the things
# library(purrr)
# library(tidyverse)
#Auxliar object
# friedman_df <- std.mlbench.friedman1(n = 100,sd = 1) %>% as.data.frame()
# plot_nrow = 2
# break_friedman_df <- std.mlbench.friedman1(n = 250,sd = 1) %>% as.data.frame()


plot_spBART_me <- function(train_me,x_train, plot_nrow){


  # Getting the
  par(mfrow = c(4,3))
  fx1x2 <-colMeans(main_effects_train_list[[11]]) #+ main_effects_train_list[[1]] %>% colMeans()

  interaction_df <- cbind(x_train[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(x_train)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*20,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx1 <- colMeans(main_effects_train_list[[1]])
  plot(x_train[,1],fx1,pch=20, xlab = "x.1", main = "f(x.1)")

  # Plotting x.3
  fx2 <- colMeans(main_effects_train_list[[2]])
  plot(x_train[,2],fx2,pch=20, xlab = "x.2", main = "f(x.2)")

  # Plotting x.3
  fx3 <- colMeans(main_effects_train_list[[3]])
  plot(x_train[,3],fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- colMeans(main_effects_train_list[[4]])
  plot(x_train[,4],fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- colMeans(main_effects_train_list[[5]])
  plot(x_train[,5],fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(x_train[,i],colMeans(main_effects_train_list[[i]]),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }


}

plot_break_friedman <- function(break_friedman_df, plot_nrow){


  # Getting the
  par(mfrow = c(plot_nrow,(NCOL(break_friedman_df)-1)/plot_nrow))
  fx1x2 <- 10* sin(pi * break_friedman_df[, 1] * break_friedman_df[, 2])
  fx1x2[break_friedman_df[,1]<0.5] <- fx1x2[break_friedman_df[,1]<0.5]  - 5
  fx1x2[break_friedman_df[,1]>0.5] <- fx1x2[break_friedman_df[,1]>0.5]  + 5

  interaction_df <- cbind(break_friedman_df[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(break_friedman_df)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*0.3,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx3 <- numeric(nrow(break_friedman_df))
  fx3[break_friedman_df[,3]>0.5] <- fx3[break_friedman_df[,3]>0.5] + 20*(break_friedman_df$x.3[break_friedman_df$x.3>0.5]-0.5)^2 + 5
  fx3[break_friedman_df[,3]<0.5] <- fx3[break_friedman_df[,3]<0.5] - 20*(break_friedman_df$x.3[break_friedman_df$x.3<0.5]-0.5)^2 - 5

  plot(break_friedman_df$x.3,fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- numeric(nrow(break_friedman_df))
  fx4[break_friedman_df[,4]>0.3] <- fx4[break_friedman_df[,4]>0.3] + 10*(break_friedman_df$x.4[break_friedman_df$x.4>0.3]) + 5
  fx4[break_friedman_df[,4]<0.3] <- fx4[break_friedman_df[,4]<0.3] - 15*(break_friedman_df$x.4[break_friedman_df$x.4<0.3]) - 5

  plot(break_friedman_df$x.4,fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- 5*break_friedman_df$x.5
  plot(break_friedman_df$x.5,fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(break_friedman_df[,i],rep(0,nrow(break_friedman_df)),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }

}

# Plotting main effects  friedman1
plot_friedman1 <- function(friedman_df, plot_nrow = 2){

  # Getting the
  par(mfrow = c(plot_nrow,(NCOL(friedman_df)-1)/plot_nrow))
  fx1x2 <- 10* sin(pi * friedman_df[, 1] * friedman_df[, 2])
  interaction_df <- cbind(friedman_df[,1:2],fx1x2) %>% dplyr::arrange(fx1x2) %>% mutate(color = gray.colors(nrow(friedman_df)))

  # Plotting interaction
  plot(interaction_df[,1:2],cex = interaction_df$fx1x2*0.3,col = interaction_df$color, pch=19,main = "f(x.1,x.2)")

  # Plotting x.3
  fx3 <- 20*(friedman_df$x.3-0.5)^2
  plot(friedman_df$x.3,fx3,pch=20, xlab = "x.3", main = "f(x.3)")

  # Plotting x.4
  fx4 <- 10*friedman_df$x.4
  plot(friedman_df$x.4,fx4,pch=20, xlab = "x.4", main = "f(x.4)")

  # Plotting x.5
  fx5 <- 5*friedman_df$x.5
  plot(friedman_df$x.5,fx5,pch=20, xlab = "x.5", main = "f(x.5)")

  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

  for(i in 6:10){
    plot(friedman_df[,i],rep(0,nrow(friedman_df)),pch=20,
         xlab = paste0("x.",i,""), main = paste0("f(x.",i,")"),
         ylab = paste0("f(x.",i,")"))
  }

}


break.mlbench.friedman1 <- function (n, sd = 1)
{
  x <- matrix(runif(10 * n), ncol = 10)
  # Adding a break
  y <- 10 * sin(pi * x[, 1] * x[, 2])
  y[ x[,1] > 0.5] < -  y[ x[,1] > 0.5] + 5
  y[ x[,1] < 0.5] < -  y[ x[,1] < 0.5] - 5

  # Adding a break for x.3
  y[x[,3] > 0.5] <- y[x[,3] > 0.5] + 20* (x[x[,3] > 0.5, 3] - 0.5)^2 + 5
  y[x[,3] < 0.5] <- y[x[,3] < 0.5] - 20* (x[x[,3] < 0.5, 3] - 0.5)^2 - 5

  y[x[,4] > 0.3] <- y[x[,4] > 0.3] + 10* (x[x[,4] > 0.3, 4]) + 5
  y[x[,4] < 0.3] <- y[x[,4] < 0.3] - 15* (x[x[,4] < 0.3, 4]) - 5

  y <- y + 5 * x[, 5]

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

std.mlbench.friedman1 <- function (n, sd = 1)
{
  x <- matrix(runif(10 * n), ncol = 10)
  y <- 10 * sin(pi * x[, 1] * x[, 2])
  y <- y + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}


mlbench.friedman1.nointeraction <- function (n, sd = 1)
{
  x <- matrix(runif(4 * n), ncol = 4)
  y <- 10 * sin(pi * x[, 1])
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.friedman1.interaction.only <- function (n, sd = 1)
{
  x <- matrix(runif(2* n), ncol = 2)
  y <- 10 * sin(pi * x[, 1]*x[,2])
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.friedman1.nointeraction.noise <- function (n, sd = 1)
{
  x <- matrix(runif(8 * n), ncol = 8)
  y <- 10 * sin(pi * x[, 1])
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1 <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1.break <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)
  y[x<0] <- y[x<0] + 5
  y[x>=0] <- y[x>=0] - 5

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

mlbench.d1 <- function(n, sd = 1) {
  x <- matrix(runif(n,min = -pi,max = pi),ncol = 1)
  y <- sin(2*x)

  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

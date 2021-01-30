rm(list = ls())


library(ggplot2)
# library(MASS)
# install.packages("tscount")
library(tscount)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(reshape2)



# get data
setwd("~/covid_19/NB model/")
# original_data = as.data.frame(read.csv("./Korea death estimates.csv", header=TRUE, skip = 1, stringsAsFactors = F))

today <- Sys.Date()
today <- format(today, format="%y%m%d")
data_date <- "20210128"
save_path <- paste0(getwd(), "/", today)
if (!dir.exists(save_path)) {
  dir.create(save_path)
}


color_list <- brewer.pal(n = 8, name = "Dark2")
y_variables <- c("DailyConfirmed", "DailyDeath")

original_data = as.data.frame(read.csv(paste0("~/covid_19/Korea_Prep/Data/COVID_19_Cases_Korea_Merged_Fin_", data_date, ".csv"), stringsAsFactors = F))
original_data$Date <- as.Date(original_data$Date)
colnames(original_data)

data_df <- original_data %>% filter(Country_Region == "Domestic") %>% select(one_of(as.character(c(y_variables, "Date"))))
data_df <- data_df[complete.cases(data_df), ]


# lag_days <- c(7, 14, 21)
lag_days <- c(1:3, 7, 14, 21)
n_pred <- 42


# variables setting
{
  predictor_list <- c()
  for (i in 1:(length(lag_days)-1)) {
    predictor_list <- rbind(predictor_list, t(combn(lag_days, m = i, function(x) as.vector(lag_days %in% x))))
  }
  predictor_list <- rbind(predictor_list, t(combn(lag_days, m = length(lag_days), function(x) as.vector(lag_days %in% x))))
  colnames(predictor_list) <- paste0("L", as.character(lag_days))
}


predictor_list <- as.data.frame(predictor_list)
idx <- intersect(which(predictor_list$L1 | predictor_list$L2 | predictor_list$L3), which(predictor_list$L7 | predictor_list$L14 | predictor_list$L21))

predictor_list <- predictor_list[idx, ]

model_result <- c()

for (each_y in y_variables) {
  # each_y <- y_variables[1]
  for (i in 1:nrow(predictor_list)) {
    # i <- 5
    # predictors <- colnames(predictor_list)[predictor_list[i, ]]
    predictors <- colnames(predictor_list)[as.vector(t(predictor_list[i, ]))]
    
    each_formula <- paste(each_y, paste(predictors, collapse = " + "), sep = " ~ ")
    
    predictor_days <- gsub("L", "", predictors)
    y_date <- predictor_days[which(as.numeric(predictor_days)<7)]
    mean_date <- predictor_days[which(as.numeric(predictor_days)>=7)]
    
    
    
    try({
      # negative binomical model
      ts_model_nb <- tsglm(data_df[, each_y], model = list(past_obs = as.numeric(mean_date), past_mean = as.numeric(y_date)), distr = "nbinom", link = 'log')
      # ts_model_nb <- tsglm(data_df[, each_y], model = list(past_obs = as.numeric(predictor_days)), distr = "nbinom", link = 'log')
      # Cholesky factorization is used to solve system(s) of linear equations where the matrix is symmetric and POSITIVE DEFINITE. 
      # All leading minors of positive definite matrix are positive. 
      # So your error message means that your matrix is not positive definite. 
      # It may be either indefinite (i.e. have both positive and negative eigenvalues) or your matrix may be near singular, i.e. it's smallest eigenvalue is very close to 0 (and so computationally it is 0).

      sink(paste0(save_path, "/ts_NB_model_", each_formula, ".txt"))
      print(ts_model_nb)
      summary(ts_model_nb)
      print(summary(ts_model_nb))
    })
    sink()  # returns output to the console
    
    if (ts_model_nb$sigmasq != Inf) {
      pred <- predict(ts_model_nb, n.ahead = n_pred, level = 0.95)
      
      predicted_mat <- data.frame("Date" = max(data_df$Date)+(1:n_pred), 
                                  "fit" = pred$pred, 
                                  "L" = pred$interval[, "lower"],
                                  "U" = pred$interval[, "upper"])
      
      plot_df <- merge(data_df, predicted_mat, by = "Date", all = T)
      plot_df[1:nrow(data_df), "fit"] <- ts_model_nb$fitted.values
      
      p1 <- ggplot(plot_df) +
        geom_line(aes(x = Date, y = fit)) +
        geom_point(aes(x = Date, y = !!sym(each_y)), color = alpha('blue', 0.5), size = 2) +
        # geom_line(aes(x = Date, y = U), linetype="dashed", color = alpha('red', 0.5)) +
        # geom_line(aes(x = Date, y = L), linetype="dashed", color = alpha('red', 0.5)) +
        # scale_x_date(limits = as.Date(c("2020-01-01", "2021-03-31"))) +
        labs(title = paste0("COVID-19 ", each_y), subtitle = each_formula) +
        theme_bw()
      
      p1
      ggsave(p1, filename = paste0(save_path, "/ts_NB_model_", each_formula, ".png"), width = 12, height = 7, dpi = 600)
      
      model_result <- rbind(model_result, c(each_y, each_formula, AIC(ts_model_nb)))
      print(paste(each_formula, AIC(ts_model_nb)))
    }
  }
  write.csv(model_result, file = paste0(save_path, "/", each_y, "_model_result.csv"))
}

# write model result csv

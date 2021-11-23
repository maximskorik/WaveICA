#' @title WaveICA
#' @description Removing batch effects for metabolomics data.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-matrix metabolomics data.
#' @param wf selecting wavelet functions, the default is "haar".
#' @param batch batch labels.
#' @param group denoting the biological group such as disease vs group.
#' This param is optional. The default is NULL.
#' @param K The maximal component that ICA decomposes.
#' @param t The threshold to consider a component associate with the batch,
#' should be between 0 and 1.
#' @param t2 The threshold to consider a component associate with the group,
#' should be between 0 and 1.
#' @param alpha The trade-off value between the independence of samples and those
#' of variables and should be between 0 and 1.
#' @return A list that contains the clean data.
#' @export
WaveICA <- function(data,
                    wf = "haar",
                    batch = NULL,
                    cutoff = NULL,
                    injection_order = NULL,
                    group = NULL,
                    K = 20,
                    t = 0.05,
                    t2 = 0.05,
                    alpha = 0) {
  ### Wavelet Decomposition
  library(waveslim)
  level <- floor(log(nrow(data), 2))
  coef <- list()
  for (k in 1:(level + 1)) {
    coef[[k]] <- matrix(NA, nrow(data), ncol(data))
  }
  for (j in 1:ncol(data)) {
    # cat(paste("######Decomposition",j,"########\n"))
    data_temp <- data[, j]
    x_modwt <- modwt(data_temp, wf = wf, n.levels = level)
    for (k in 1:(level + 1)) {
      coef[[k]][, j] <- x_modwt[[k]]
    }
  }
  ##### ICA
  index <- level + 1
  data_wave_ICA <- list()
  for (i in (1:index)) {
    # cat(paste("######### ICA",i,"#############\n"))
    data_coef <- coef[[i]]
    data_coef_ICA <- normFact(fact = "stICA", X = t(data_coef), ref = batch, refType = "categorical", k = K, t = t, ref2 = group, refType2 = "categorical", t2 = t2, alpha)
    data_wave_ICA[[i]] <- t(data_coef_ICA$Xn)
  }
  ### Wavelet Reconstruction
  index <- ncol(data)
  index1 <- length(data_wave_ICA)
  data_coef <- matrix(NA, nrow(data_wave_ICA[[1]]), index1)
  data_wave <- matrix(NA, nrow(data_wave_ICA[[1]]), ncol(data_wave_ICA[[1]]))
  for (i in 1:index) {
    # cat(paste("######Reconstruction",i,"########\n"))
    for (j in 1:index1) {
      data_coef[, j] <- data_wave_ICA[[j]][, i]
    }
    data_temp <- data[, i]
    data_coef <- as.data.frame(data_coef)
    colnames(data_coef) <- c(paste("d", 1:(index1 - 1), sep = ""), paste("s", (index1 - 1), sep = ""))
    y <- as.list(data_coef)
    attributes(y)$class <- "modwt"
    attributes(y)$wavelet <- wf
    attributes(y)$boundary <- "periodic"
    data_wave[, i] <- imodwt(y) + mean(data_temp)
  }
  rownames(data_wave) <- rownames(data)
  colnames(data_wave) <- colnames(data)
  return(list(data_wave = data_wave))
}

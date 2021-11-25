wt_decomposition <- function(data, level, wf) {
  ### Wavelet Decomposition
  coef <- list()
  for (k in 1:(level + 1)) {
    coef[[k]] <- matrix(NA, nrow(data), ncol(data))
  }
  cat(paste("Decomposing...\n"))
  for (j in 1:ncol(data)) {
    data_temp <- data[, j]
    x_modwt <- modwt(data_temp, wf = wf, n.levels = level)
    for (k in 1:(level + 1)) {
      coef[[k]][, j] <- x_modwt[[k]]
    }
  }
  return(coef)
}

wt_reconstruction <- function(data, data_wave_ICA, wf) {
  ### Wavelet Reconstruction
  index <- ncol(data)
  index1 <- length(data_wave_ICA)
  data_coef <- matrix(NA, nrow(data_wave_ICA[[1]]), index1)
  data_wave <- matrix(NA, nrow(data_wave_ICA[[1]]), ncol(data_wave_ICA[[1]]))
  cat(paste("Reconstructing...\n"))
  for (i in 1:index) {
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
  return(data_wave)
}

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
#' @importFrom mgcv gam
#' @importFrom parallel mclapply
#' @importFrom waveslim modwt imodwt
#' @export
WaveICA <- function(data,
                    wf = "haar",
                    batch = NULL,
                    Cutoff = NULL,
                    Injection_order = NULL,
                    group = NULL,
                    K = 20,
                    t = 0.05,
                    t2 = 0.05,
                    alpha = 0) {

  level <- floor(log(nrow(data), 2))
  coef <- wt_decomposition(data, level, wf)

  ##### ICA
  index <- level + 1
  data_wave_ICA <- list()
  cat(paste("Performing ICA...\n"))
  for (i in (1:index)) {
    data_coef <- coef[[i]]
    if (is.null(batch)) {
      data_coef_ICA<-unbiased_stICA(X=t(data_coef),k=K,alpha)
      B <- data_coef_ICA$B
      A <- data_coef_ICA$A
      B <- as.data.frame(B)

      ## Gam
      corr <- mclapply(B,function(x){
      corr <- gam(x~s(Injection_Order))
      corr_summary <- summary(corr)
      corr_r <- corr_summary$r.sq
      return(corr_r)
      })
      corr <- unlist(corr)
      label <- which(corr>=Cutoff)
      B_new <- B[,label,drop=F]
      A_new <- A[,label,drop=F]
      Xn = data_coef-t(A_new %*% t(B_new))

      data_wave_ICA[[i]]<-Xn
    } else {
    data_coef_ICA <- normFact(fact = "stICA", X = t(data_coef), ref = batch, refType = "categorical", k = K, t = t, ref2 = group, refType2 = "categorical", t2 = t2, alpha)
    data_wave_ICA[[i]] <- t(data_coef_ICA$Xn)
    }
  }

  data_wave <- wt_reconstruction(data, data_wave_ICA, wf)

  rownames(data_wave) <- rownames(data)
  colnames(data_wave) <- colnames(data)
  return(list(data_wave = data_wave))
}

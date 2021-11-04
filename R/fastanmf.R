#' Fast and Stable NMF
#' 
#' @export
#'
#' @import NMF
#' @import ggplot2
#' @import matrixStats

fastanmf <- function(data, rank, nrun = 200, nmf_seed = 123456, mvg = 1000, cophenetic = 0, ...) {

  # add row and column names if missing
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("sample",as.character(1:ncol(data)))
  }
  if (is.null(rownames(data))) {
    rownames(data) <- paste0("gene",as.character(1:nrow(data)))
  }
  
  # Find ~1000ish most variable genes
  print("Reducing the dataset to the most variable genes...")
  reduced_data <- nmf_p0(data, mvg)
  
  # Initial NMF run on reduced data set
  print("Performing the initial NMF run...")
  nmf_object1 <- NMF::nmf(reduced_data, rank = rank, nrun = nrun, seed = nmf_seed, ...)
  # matrix_analyzer(nmf_object1, 'h')
  
  # find core subset 
  # create initial synthetic seed
  print("Creating the initial synthetic seed...")
  seed_matrix <- nmf_p1(nmf_object1)
  if (cophenetic == 1) {
    return(ncol(seed_matrix))
  }
  
  # run seeded nmf on core
  # expand to reduced data set
  # run seed nmf on reduced data set
  print("Running seeded NMF...")
  nmf_object3 <- nmf_p2(seed_matrix, nmf_object1, reduced_data, ...)
 
  # expand to full data set  
  # run seeded nmf on full data set
  print("Expanding back to the full dataset...")
  nmf_object4 <- nmf_p3(nmf_object3, data, ...) 

  return(nmf_object4)
  
 
}

# Find ~1000ish most variable genes
nmf_p0 <- function(data, mvg) {
  data <- data[order(rowVars(data), decreasing = TRUE),]
  return( t(data[1:min(nrow(data), mvg),]) )
}

# Find the genes in the top ~1000ish that always cluster together
nmf_p1 <- function(nmf_object) {
  rank = ncol(nmf_object@fit@W)
  l <- c()
  w <- matrix(data = NA, nrow = ncol(nmf_object@fit@H), ncol = rank)
  row.names(w) <- colnames(nmf_object@fit@H)
  means <- c(0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50)
  for (i in 1:rank) {
    ind <- which(predict(nmf_object) == i)
    for (j in 1:length(means)) {
      lt <- helper_function(nmf_object, means[j], i)
      lens <- lengths(lapply(lt, unlist))
      if (length(lens) == 0) {
        if (means[j] == 0.50) {
          l <- c(l, ind)
          break
        } else {
          next 
        }
      } else {
        max_lens <- max(lens)
        which_lens <- which(lens == max_lens)
        l <- c(l, ind[unique(unlist(lt[which_lens]))])
        w_i <- w[,i]
        w_i[ind[unique(unlist(lt[which_lens]))]] <- 1
        w[,i] <- w_i
        break
      }
    }
  }
  w[is.na(w)] <- 0.001
  w <- w[which(rowSums(w) > 1),]
  return(t(w))
}

# Run NMF with the genes that always cluster together and then once again after extending it
nmf_p2 <- function(seed_matrix, nmf_object1, original_data, cophenetic, ...) {
  rank = nrow(seed_matrix)
  #reduced_seed_matrix <- seed_matrix[,colSums(seed_matrix[])>(0.001 * rank)]
  cn <- colnames(seed_matrix)
  if(anyDuplicated(cn)){warning("duplicated gene names")}
  reduced_data <- original_data[,cn]
  #w_matrix <- matrix(, nrow = nrow(original_data), ncol = rank)
  seed_model <- NMF::nmfModel(rank = rank, W = nmf_object1@fit@W, H = seed_matrix)
  # now seeded with better than 1s, this is forgetting useful information
  nmf_object <- NMF::nmf(reduced_data, rank = rank, seed = seed_model, nrun = 1, ...)
  # matrix_analyzer(nmf_object, 'h')
  c1 <- c()
  c2 <- c()
  for (i in 1:rank) {
    c1 <- c(c1, nmf_object@fit@H[i,])
    c2 <- c(c2, seed_matrix[i,])
  }
  prdct <- NMF::predict(nmf_object)
  temp_df <- cbind(c1, c2, prdct)
  temp_df <- as.data.frame(temp_df)
  colnames(temp_df) <- c("new", "old", "fac")
  # print(ggplot(temp_df, aes(x = as.numeric(old), y = as.numeric(new), color = as.factor(fac))) + geom_point())
  opposite_reduced_data <- original_data[, !colnames(original_data) %in% cn]
  
  mean.data.bygene <- colMeans(reduced_data) #scale of the real data genes
  mean.nmf.bygene <- colMeans(nmf_object@fit@H) #scale of the gene weights from NMF
  scale.factor <- mean(mean.nmf.bygene/mean.data.bygene) # average ratio of these
  mean.data.bygene <- colMeans(opposite_reduced_data) # scale of the real data genes
  mean.nmf.bygene <- mean.data.bygene * scale.factor # implied scale that the NMF will be on
  h_m <- t(matrix(mean.nmf.bygene, ncol = rank, nrow = ncol(opposite_reduced_data))) # useful seed
  colnames(h_m) <- colnames(opposite_reduced_data)
  
  h_matrix2 <- cbind(nmf_object@fit@H, h_m)
  h_matrix3 <- h_matrix2[ ,colnames(original_data)]
  seed_model2 <- NMF::nmfModel(rank = rank, W = nmf_object@fit@W, H = h_matrix3)
  nmf_object <- NMF::nmf(original_data, rank = rank, seed = seed_model2, nrun = 1, ...)
  # matrix_analyzer(nmf_object, 'h')
  c1 <- c()
  c2 <- c()
  for (i in 1:rank) {
    c1 <- c(c1, nmf_object@fit@H[i,])
    c2 <- c(c2, h_matrix3[i,])
  }
  prdct <- NMF::predict(nmf_object)
  temp_df <- cbind(c1, c2, prdct)
  temp_df <- as.data.frame(temp_df)
  colnames(temp_df) <- c("new", "old", "fac")
  # print(ggplot(temp_df, aes(x = as.numeric(old), y = as.numeric(new), color = as.factor(fac))) + geom_point())
  return(nmf_object)
}

# Extend to the entire data
nmf_p3 <- function(nmf_object, original_data, ...) {
  h <- nmf_object@fit@H
  w <- nmf_object@fit@W
  rank = nrow(h)
  cn <- colnames(h)
  opposite_reduced_data <- original_data[!rownames(original_data) %in% cn,]
  reduced_data <- original_data[cn,]
  
  mean.data.bygene <- rowMeans(reduced_data) #scale of the real data genes
  mean.nmf.bygene <- colMeans(nmf_object@fit@H) #scale of the gene weights from NMF
  scale.factor <- mean(mean.nmf.bygene/mean.data.bygene) # average ratio of these
  mean.data.bygene <- rowMeans(opposite_reduced_data) # scale of the real data genes
  mean.nmf.bygene <- mean.data.bygene * scale.factor # implied scale that the NMF will be on
  h_m <- t(matrix(mean.nmf.bygene, ncol = rank, nrow = nrow(opposite_reduced_data))) # useful seed
  colnames(h_m) <- rownames(opposite_reduced_data)

  h_matrix2 <- cbind(h, h_m)
  h_matrix3 <- h_matrix2[ ,rownames(original_data)]
  seed_model2 <- NMF::nmfModel(rank = rank, W = t(h_matrix3), H = t(w))
  nmf_object <- NMF::nmf(original_data, rank = rank, seed = seed_model2, nrun = 1, ...)
  # matrix_analyzer(nmf_object, 'w')
  c1 <- c()
  c2 <- c()
  c3 <- c()
  for (i in 1:rank) {
    c1 <- c(c1, nmf_object@fit@W[,i])
    c2 <- c(c2, t(h_matrix3)[,i])
    c3 <- c(c3, colnames(h_matrix3) %in% cn)
  }
  temp_df <- cbind(c1, c2, c3)
  temp_df <- as.data.frame(temp_df)
  colnames(temp_df) <- c("new", "old", "seeded")
  # print(ggplot(temp_df, aes(x = as.numeric(old), y = as.numeric(new), color = seeded)) + geom_point())
  return(nmf_object)
}

# Cuts the dendrogram
helper_function <- function(nmf_object, prob, i) {
  lt <- list()
  cons <- NMF::consensus(nmf_object)
  ind <- which(predict(nmf_object) == i)
  subset <- cons[ind,ind]
  x <- hclust(as.dist(1-subset))
  val = 0
  prev = NULL
  while (val < 0.99) {
    val = val + 0.001
    y <- cut(as.dendrogram(x), val)
    z <- sapply(y$lower, nobs)
    if (is.null(prev)) {
      prev <- z
    } else {
      if (length(prev) != length(z)) {
        prev <- z
      } else {
        if (all(prev == z)) {
          next
        } else {
          prev <- z
        } 
      }
    }
    for (j in 1:length(z)) {
      if (z[j] >= 5) {
        a <- y$lower[j]
        unlist_a <- unlist(a)
        test_consensus <- subset[unlist_a, unlist_a]
        if (mean(test_consensus) >= prob) {
          lt[[length(lt) + 1]] <- unlist_a
        }
      } 
    }
  }
  return(lt)
}

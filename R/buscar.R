#' Calculate K for H0
#'
#' `k_h0` calculates the smallest number of relevant documents there would
#' have to be in our urn, for us to have missed a given recall target
#'
#' @export
#' @param r_al The number of relevant documents seen before drawing from the urn began
#' @param r_seen The total number of relevant documents seen.
#' @param recall_target The recall target.
#' @return k The smallest integer compatible with our null hypothesis
k_h0 <- function(r_al, r_seen, recall_target){
  return(floor(r_seen/recall_target-r_al+1))
}

#' Calculate h0
#'
#' Calculates a p-score for our null hypothesis h0, that we have
#' missed our recall target `recall_target`.
#'
#' @export
#' @param df A data.frame that contains the columns `relevant` and `seen`
#' The dataframe should have as many rows as there are documents, and be
#' ordered in the order dictated by the ML prioritisation algorithm.
#' relevant should contain 1s and 0s for relevant and irrelevant documents,
#' and NAs for documents that have not yet been screened. Seen should contain
#' 1s where documents have been screened by a human, and 0s where documents
#' have not yet been screened
#' @param recall_target The recall target (default=0.95). Must be between 0 and 1
#' @param bias a number which represents our estimate of how much more likely we
#' were to select a random relevant document than a random irrelevant
#' document. The higher this is, the better we think the machine learning went.
#' @param seen_docs an integer which overrides the seen column, telling us
#' how many of the first documents have been screened
#' @return p, a p-score for our null hypothesis. We can reject the null
#' hypothesis (and stop screening) if p is below 1 - our confidence level.
#' @examples
#' N <- 60000 # number of documents
#' prevalence <- 0.01 # prevalence of relevant documents
#' r <- N*0.01 # number of relevant documents
#' bias <- 10
#' docs <- rep(0,N)
#' docs[1:r] <- 1
#' weights = rep(1,N)
#' weights[1:r] <- bias
#' set.seed(2023)
#' docs <- sample(
#'   docs, prob=weights, replace=F
#' )
#' df <- data.frame(relevant=docs)
#' df$seen <- 0
#' df$seen[1:1000] <- 1
#' calculate_h0(df)
calculate_h0 <- function(
    df, recall_target=0.95, bias=1, seen_docs=NULL
    ) {
  # Check dataframe in structured properly
  if (!"relevant" %in% colnames(df)) {
    stop("Dataframe must contain a column named 'relevant'")
  }
  if (!"seen" %in% colnames(df)) {
    stop("Dataframe must contain a column named 'seen'")
  }

  N <- nrow(df)
  if (is.null(seen_docs)) {
    docs <- df[df$seen==1,"relevant"]
  } else {
    docs <- df$relevant[1:seen_docs]
  }

  if (length(docs)==N) {
    print(length(docs))
    print(N)
    stop("Dataframe must contain more than zero unseen documents")
  } else if (length(docs)==0) {
    stop("Dataframe must contain more than zero screened documents")
  }

  for (val in unique(docs)) {
    if (!val %in% c(1,0)) {
      stop(paste0("screened docs can only contain values 1 and 0, found ", val))
    }
  }

  r_seen <- sum(docs)
  n_vec <- seq(1:length(docs))
  k_vec <- cumsum(rev(docs))
  r_al_vec <- r_seen - k_vec
  k_hat_vec <- vapply(r_al_vec, k_h0, numeric(1),
                      recall_target=recall_target, r_seen=r_seen)
  red_ball_vec <- N-(length(docs)-n_vec)-k_hat_vec
  if (bias==1) {
    p_vec <- phyper(k_vec, k_hat_vec, red_ball_vec, n_vec)
  } else {
    odds <- rep(bias,length(docs))
    p_vec <- mapply(BiasedUrn::pWNCHypergeo, k_vec, k_hat_vec, red_ball_vec, n_vec, odds)
  }
  n <- which.min(p_vec)
  return(p_vec[n])
}

#' Calculate recall frontier
#'
#' Calculates p scores across different recall targets
#' @export
#' @inheritParams calculate_h0
#' @param plot Boolean describing whether to plot a graph (default=True).
#' @returns A dataframe with a column `p` showing the p score for h0 calculated
#' given each recall target `target`
#' @examples
#' N <- 60000 # number of documents
#' prevalence <- 0.01 # prevalence of relevant documents
#' r <- N*0.01 # number of relevant documents
#' bias <- 10
#' docs <- rep(0,N)
#' docs[1:r] <- 1
#' weights = rep(1,N)
#' weights[1:r] <- bias
#' set.seed(2023)
#' docs <- sample(
#'   docs, prob=weights, replace=F
#' )
#' df <- data.frame(relevant=docs)
#' df$seen <- 0
#' df$seen[1:20000] <- 1
#' recall_df <- recall_frontier(df)
recall_frontier <- function(df, bias=1) {
  recall_targets <- seq(0.01,0.99,0.005)
  p <- sapply(recall_targets, calculate_h0, df=df)
  recall_df <- data.frame(target=recall_targets, p=p)
  plot_recall_df <- recall_df[recall_df$p>0.001,]
  plot(
    plot_recall_df$target,
    plot_recall_df$p,
    type='b',
    main='p score for H0 across recall targets',
    xlab='Recall target',
    ylab='p score',
    pch=21
  )
  return(recall_df)
}

#' Calculate retrospective H0
#'
#' Calculates p scores for our null hypothesis H0, every `batch_size` documents
#' of the documents we have already screened, and plot a graph these p scores,
#' alongside the curve showing the number of relevant documents identified.
#' @export
#' @inheritParams calculate_h0
#' @param batch_size The p score will calculated every `batch_size` documents.
#' Smaller batches will result in greater granularity but larger computation
#' time (default=1000).
#' @param plot Boolean describing whether to plot a graph (default=True).
#' @returns A dataframe with a column `p` showing the p score for h0 calculated
#' at number of screened documents in column `seen`
#' @examples
#' N <- 60000 # number of documents
#' prevalence <- 0.01 # prevalence of relevant documents
#' r <- N*0.01 # number of relevant documents
#' bias <- 10
#' docs <- rep(0,N)
#' docs[1:r] <- 1
#' weights = rep(1,N)
#' weights[1:r] <- bias
#' set.seed(2023)
#' docs <- sample(
#'   docs, prob=weights, replace=F
#' )
#' df <- data.frame(relevant=docs)
#' df$seen <- 0
#' df$seen[1:30000] <- 1
#' h0_df <- retrospective_h0(df)
retrospective_h0 <- function(df, recall_target=0.95, bias=1, batch_size=1000, plot=TRUE) {
  n_seen <- nrow(df[df$seen==1,])
  batches <- seq(batch_size, n_seen, batch_size)
  p <- sapply(
    batches,
    calculate_h0,
    df=df,
    recall_target=recall_target,
    bias=bias
  )
  h0_df <- data.frame(seen=batches, p=p)
  if (plot==TRUE) {
    par(mfrow = c(2, 1))
    par(oma = c(4, 0, 0, 0))
    par(mar = c(2, 2, 1, 1))
    par(
      tck = -.02 # Reduce tick length
    )
    plot(
      cumsum(df[df$seen==1,"relevant"]),
      type='l',
      main="Relevant documents identified",
      xlab="",
      ylab="Relevant  documents seen",
      xlim=c(0,nrow(df))
    )
    plot(
      h0_df$seen,
      h0_df$p,
      main='p score for H0',
      type='b',
      ylab='p H0',
      xlab='Documents seen',
      xlim=c(0,nrow(df)),
    )
    mtext('Documents screened', side = 1, outer = TRUE, line = 2)
  }
  return(h0_df)
}

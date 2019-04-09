################################################################################
#Downloaded from: https://github.com/defleury/adjusted_mutual_information
#and modified under
#GNU General Public License v3.0
#
#Calculate pairwise AMI, NMI and ARI between pairs of cluster sets
#
#AMI		=> Adjusted Mutual Information
#NMI		=> Normalized Mutual Information
#ARI		=> Adjusted Rand Index
#
#All measures are calculated following the formulas given in
#Vinh et al, 2009, Information Theoretic Measures for Clustering Comparison [...]
#Schmidt et al, 2014, Limits to Robustness & Reproducibility [...] (-> Supplementary Text S1)
#
#2016-04-10
#sebastian.schmidt@imls.uzh.ch
################################################################################

# Load Packages
# library("parallel", warn.conflicts=F, quietly=T);
# for now the parallel computations have been disabled

clustering.indices <- function(clustering1,clustering2) {
  #if(is.null(use.cores))
  #  use.cores <- 1 #how many CPUs to use during parallelization
  
  ################################################################################
  #Load data
  # ground.truth data is in format, like so:
  #[otu_idx]	[sequence_indices, comma-separated]
  #=> listing for every cluster all the texts it contains.
  #For example:
  #1	1,10,12,15,20,
  #2	2,7,8,9,11,25,27,
  #...
  # the clustering set being compared is in format:
  #[text_idx]	[cluster_idx]
  # => listing for every text the cluster into which it clustered:
  #For example (using the above example):
  #1	1
  #2	2
  #...
  #7	2
  #8	2
  #9	2
  #10	1
  #...

  
  #Get cluster sizes for set.1
  sizes.1 <- unlist(lapply(clustering1, length));
  #Get total number of clusters for set.1
  clust_count.1 <- length(clustering1);
  
  #Get cluster sizes for set.2
  sizes.2 <- tabulate(clustering2);
  #Get total number of clusters for set.2
  clust_count.2 <- length(sizes.2);
  
  n <- length(clustering2) # get the number of texts

  
  ################################################################################
  #Precalculate sum(n.ij) matrix for expected mutual information (EMI)
  #A lookup is much faster than recalculating EMI values at each step all over again.
  #EMI is calculated for all combinations of cluster sizes from 1 to maximally 1000,
  #providing a (triangular) lookup matrix (or list, in this case).
  #EMI is calculated as described in formula 24a in Vinh et al, 2009.
  
  #Nested listwise computation
  emi.lookup <- lapply(
    seq(1, min(n,1000)),
    function (a.i) {
      unlist(lapply(
        seq(1, min(n,1000)),
        function (b.j) {
          n.ij <- seq(max(a.i + b.j - n, 1), min(a.i, b.j));
          sum((n.ij / n) * log((n.ij * n) / (a.i * b.j)) * exp(lfactorial(a.i) + lfactorial(b.j) + lfactorial(n - a.i) + lfactorial(n - b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(a.i - n.ij) - lfactorial(b.j - n.ij) - lfactorial(n - a.i - b.j + n.ij)))
        }
      ));
    }
  );
  
  ################################################################################
  #Calculate indices
  
  #Compute shared mapping (n.ij) => the sparse "contingency table"
  #The output list "current.shared" contains a sublist for each cluster in set.1
  #(length(current.shared) == clust_count.1)
  #This sublist has two entries, "$n.ij" and "sizes.2". With the example data:
  #
  #> current.shared[[1]]
  # $n.ij
  # [1] 573 118 188  53  52   3  12   1   4   8   1   1   1   1   1   1
  # 
  # $sizes.2
  # [1] 573 118 188  53  52   3  12   1   4   8   1   1   1   1   1   1
  #
  #...this means that the 1018 sequences in cluster 1 of set.1 (check sizes.1[1]) clustered
  #into 16 different clusters in set.2 (try length(current.shared[[1]]$n.ij)). For each of these
  #clusterS, "$n.ij" holds the number of concordant sequences (those in cluster 1 of set.1 and the respective cluster in set.2).
  #At the same time, "$sizes.2" has the size of these clusters in set.2. In the current case, all 16 clusters in set.2
  #contain *only* sequences of cluster 1 in set.1. In other words, cl clustering split the al cluster 1 into 16 smaller clusters
  #which are completely contained in that one cluster from set.1.
  #
  #Note: this is usually the most time-consuming step.
  writeLines(paste(date(), "Calculating current shared mapping..."));
  current.shared <- lapply(clustering1, function (map.1) {
    tmp <- tabulate(clustering2[map.1]);
    out <- list();
    out$n.ij <- as.numeric(tmp[tmp != 0]);
    out$sizes.2 <- sizes.2[which(tmp != 0)];
    out
  });
  
  ################################################################################
  #Compute ARI
  #=> equation (3) in Vinh et al., 2009
  
  writeLines(paste(date(), "Calculating Adjusted Rand Index..."));
  #Calculate size-dependent terms ("sums") for the expected index
  sum.1 <- sum(unlist(lapply(sizes.1, function(X) {choose(X, 2)})));
  sum.2 <- sum(unlist(lapply(sizes.2, function(X) {choose(X, 2)})));
  #Calculate expected Rand Index for these cluster size distributions
  expected.index <- (sum.1 * sum.2) / choose(n, 2);
  #Calculate Rand Index based on shared mapping (contingency table)
  rand.index <- sum(unlist(lapply(current.shared, function (shared) {choose(shared$n.ij, 2)})));
  #Calculate ARI from all of the above
  ari <- (rand.index - expected.index) / (0.5*(sum.1 + sum.2) - expected.index);
  
  ################################################################################
  #Compute Mutual Information & Entropies (for NMI & AMI calculation)
  #=> equations (4) & (5) in Vinh et al., 2009
  
  writeLines(paste(date(), "Calculating Mutual Information and Entropies..."));
  mutual.information <- sum(unlist(lapply(current.shared, function (shared) {sum((shared$n.ij / n) * log((shared$n.ij * n) / (sum(shared$n.ij) * shared$sizes.2)))})));
  entropy.1 <- sum((sizes.1 / n) * log(sizes.1 / n));
  entropy.2 <- sum((sizes.2 / n) * log(sizes.2 / n));
  
  #Compute NMI from mutual information and entropies
  nmi <- (-2 * mutual.information) / (entropy.1 + entropy.2);
  
  ################################################################################
  #Compute Expected Mutual Information (EMI) => for the AMI formula
  
  writeLines(paste(date(), "Calculating Expected Mutual Information..."));
  tab.sizes.1 <- as.numeric(tabulate(sizes.1)); idx.sizes.1 <- which(tab.sizes.1 != 0);
  tab.sizes.2 <- as.numeric(tabulate(sizes.2)); idx.sizes.2 <- which(tab.sizes.2 != 0);
  idx.small.sizes.2 <- idx.sizes.2[idx.sizes.2 <= min(n,1000)];
  idx.large.sizes.2 <- idx.sizes.2[idx.sizes.2 >  min(n,1000)];
  expected.mutual.information <- sum(unlist(lapply(
    idx.sizes.1,
    function (a.i) {
      #Treat a.i <= 1000
      if (a.i <= min(n,1000)) {
        #Sum up n.ij (precalculated) for all b.j <= 1000
        small.sum <- sum(tab.sizes.1[a.i] * tab.sizes.2[idx.small.sizes.2] * emi.lookup[[a.i]][idx.small.sizes.2]);
        #Additionally, calculate n.ij for b.j > 1000
        large.sum <- sum(tab.sizes.1[a.i] * tab.sizes.2[idx.large.sizes.2] * unlist(lapply(
          idx.large.sizes.2,
          function (b.j) {
            n.ij <- as.numeric(seq(max(a.i + b.j - n, 1), min(a.i, b.j)));
            num.a.i <- as.numeric(a.i); num.b.j <- as.numeric(b.j);
            sum((n.ij / n) * log((n.ij * n) / (num.a.i * num.b.j)) * exp(lfactorial(num.a.i) + lfactorial(num.b.j) + lfactorial(n - num.a.i) + lfactorial(n - num.b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(num.a.i - n.ij) - lfactorial(num.b.j - n.ij) - lfactorial(n - num.a.i - num.b.j + n.ij)))
          }
        )))
        #Pass output
        small.sum + large.sum
      }
      else {
        #Sum over n.ij for all a.i >= 1000 and all b.j
        sum(tab.sizes.1[a.i] * tab.sizes.2[idx.sizes.2] * unlist(lapply(
          idx.sizes.2,
          function (b.j) {
            n.ij <- as.numeric(seq(max(a.i + b.j - n, 1), min(a.i, b.j)));
            num.a.i <- as.numeric(a.i); num.b.j <- as.numeric(b.j);
            sum((n.ij / n) * log((n.ij * n) / (num.a.i * num.b.j)) * exp(lfactorial(num.a.i) + lfactorial(num.b.j) + lfactorial(n - num.a.i) + lfactorial(n - num.b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(num.a.i - n.ij) - lfactorial(num.b.j - n.ij) - lfactorial(n - num.a.i - num.b.j + n.ij)))
          }
        )))
      }
    }
  )));
  
  ################################################################################
  #Compute AMI
  
  ami <- (mutual.information - expected.mutual.information) / (sqrt(entropy.1 * entropy.2) - expected.mutual.information);
  
  result_list<-list(ami,ari,nmi)
  names(result_list)<-list("Adjusted Mutual Information","Adjusted Rand Index","Normalised Mutual Information")
  
  return(result_list)
}
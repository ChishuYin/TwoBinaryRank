#' @title PCS0 for Subset Selection
#'
#' @description
#' Compute PCS0 (probability of select the control treatment) for the
#' subset selection procedure based on two binary endpoints.
#'
#' @param c1 Integer. Efficacy cutoff used in the subset rule.
#' @param c2 Integer. Safety cutoff used in the subset rule.
#' @param n Integer. Sample size per treatment arm.
#' @param p Numeric matrix with at least 2 rows. Row 1 contains joint endpoint
#'   probabilities for the control arm; row 2 contains joint endpoint
#'   probabilities for the experimental arm. The ordering must match `joint_prob()`.
#' @param k Integer. Number of experimental treatments (used as an exponent in the PCS0 formula).
#'
#' @return A numeric scalar giving PCS0.
#'
#' @export
#'
#' @examples
#' # Requires joint_prob() to be available
#' # p <- rbind(c(0.1, 0.2, 0.3, 0.4),
#' #            c(0.1, 0.1, 0.3, 0.5))
#' # subset_pcs0(c1 = 3, c2 = 3, n = 10, p = p, k = 2)
subset_pcs0 <- function(c1,c2,n,p,k)
{
  ss<-0
  for(i in 0:min(n-c1,n))
  {
    for(j in 0:min(n-c2,n))
    {
      ss<- ss+(1-  max(1-pbinom(c1+i-1,n,p[2,1]+p[2,2]),1- pbinom(c2+j-1,n,p[2,1]+p[2,3])) )^k * joint_prob(i,j,n,p[1,])
    }
  }
  ss
}

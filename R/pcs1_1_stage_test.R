#' @title Power for One Stage Test
#'
#' @description
#' Compute the probability of reject the null hypothesis for a one-stage test
#' comparing two treatments under joint binary endpoints.
#'
#' @param n Integer. Sample size per treatment arm.
#' @param e Numeric/integer. Efficacy margin/threshold; treatment 2 is considered better on efficacy if (c - a) >= e.
#' @param s Numeric/integer. Safety margin/threshold; treatment 2 is considered better on safety if (d - b) >= s.
#' @param p Numeric matrix with 2 rows. Row 1 contains joint endpoint probabilities for treatment 1;
#'   row 2 contains joint endpoint probabilities for treatment 2. Each row should be a length-4 vector
#'   giving probabilities for the four joint outcomes (ordering consistent with `joint_prob()`).
#' @param k Integer. (Currently unused.) Reserved for compatibility with other procedures.
#'
#' @return A numeric scalar giving the PCS.
#'
#' @export
#'
#' @examples
#' # Example (requires `joint_prob()` to be defined in the package):
#' # p <- rbind(c(0.1, 0.2, 0.3, 0.4),
#' #            c(0.1, 0.1, 0.3, 0.5))
#' # one_stage_test_pcs(n = 10, e = 1, s = 1, p = p, k = 0)
one_stage_test_pcs1 <- function(n,e,s,p,k)
{
  ss<-0
  for(a in 0:n)
  {
    for(b in 0:n)
    {
      for(c in 0:n)
      {
        for(d in 0:n)
        {
          if((c-a>=e)&(d-b>=s))
          {
            ss<-ss+joint_prob(a,b,n,p[1,])*joint_prob(c,d,n,p[2,])
          }
        }
      }
    }
  }
  ss
}

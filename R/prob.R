#' @title Calculate Multinomial Distribution Parameters
#'
#' @description
#' Computes the parameters of a multinomial distribution represented by a 4-element vector \code{p = [p1, p2, p3, p4]}.
#' The function takes the probability of efficacy (\code{pe}), the probability of safety (\code{ps}), and the odds ratio (\code{phi})
#' to determine the joint probabilities under the constraints:
#'
#' \itemize{
#'   \item \code{p1 + p2 = pe}
#'   \item \code{p1 + p3 = ps}
#'   \item \code{p1 + p2 + p3 + p4 = 1}
#' }
#'
#'
#' @param pr Numeric. Probability of efficacy (\code{pe}), where \code{0 <= pe <= 1}.
#' @param pt Numeric. Probability of safety (\code{ps}), where \code{0 <= ps <= 1}.
#' @param fi Numeric. Odds ratio (\code{phi}) representing the relationship between efficacy and safety outcomes, defined as \code{phi = (p1 * p4) / (p2 * p3)}.
#'
#' @returns A numeric vector of length 4 representing the multinomial distribution parameters \code{[p1, p2, p3, p4]}.
#' @export
#'
#' @examples prob(0.7,0.5,1.5)
#'
prob <- function(pr,pt,fi)
{
  if(fi==1)
  {
    p<-c(pr*pt,pr-pr*pt,pt-pr*pt,1-pr-pt+pr*pt)
    p
  }
  else{
    a<- 1+(fi-1)*(pr+pt)
    b <- -4*fi*(fi-1)*pr*pt
    p <- matrix(nrow = 2, ncol = 2)
    p[1,1] <- round(( a - ( a^2+b )^(0.5) ) / (2*(fi-1)),digits = 7)
    p[1,2] <- round(pr-p[1,1],digits = 7)
    p[2,1] <- round(pt-p[1,1],digits = 7)
    p[2,2] <- round(1 - p[1,1] - p[1,2] - p[2,1],digits = 4)
    if(p[2,2]<0)
    {
      p[1,1]<-p[1,1]+1/3*p[2,2]
      p[1,2]<-p[1,2]+1/3*p[2,2]
      p[2,1]<-p[2,1]+1/3*p[2,2]
      p[2,2]<-0
    }
    p <- c(p[1,1],p[1,2],p[2,1],p[2,2])
    p}
}

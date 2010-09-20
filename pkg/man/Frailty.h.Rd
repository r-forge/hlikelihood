\name{Frailty.h}
\Rdversion{1.1}
\alias{Frailty.h}

\docType{package}

\title{Fitting Frailty Models using H-likelihood Approach}

\description{Frailty.h is used to fit frailty models using h-likelihood approach. The package fits Cox's proportional hazards models with random effects (or frailties). Here, for the frailty distribution lognormal or gamma is allowed. The h-likelihood obviates the need for marginalization over the frailty distribution, leading to a statistically efficient procedure for various random-effect models including frailty models.}

\usage{
Frailty.h(formulaMain, censor, DataMain, RandDist = "Normal", mord = 0,
          dord = 1, Maxiter = 200, convergence = 1e-7, contrasts = NULL)
}

\arguments{
  \item{formulaMain}{Specify the formula for the mean structure of the model. e.g. \code{y ~ x + (1|id)}, \code{y}: response, \code{x}: fixed covariate, \code{id}: random effect.}
  \item{censor}{One vector for censoring indicator (1 = uncensored, 0 = censored).}
  \item{DataMain}{Dataframe for formulaMain.}
  \item{RandDist}{Distribution for random effect ("Normal" or "Gamma").}
  \item{mord}{Order for the mean model (0 or 1); default = 0.}
  \item{dord}{Order for the dispersion components (1 or 2); default = 1.}
  \item{Maxiter}{Maximum number of iterations; default = 1.}
  \item{convergence}{Specify the convergence criterion, the default is 1e-7.}
  \item{contrasts}{Caution as it is currently not fully developed.}
}

\details{
\code{Frailty.h} produces estimates of fixed effects and frailty parameters as well as their standard errors. Also, \code{Frailty.h} makes it possible to fit models where the frailty distribution is not necessarily normal and estimate variance components when frailty structure is shared or nested.
}

\references{
Ha, I.D. and Lee, Y. (2003). Estimating frailty models via Poisson Hierarchical generalized linear models. \emph{Journal of Computational and Graphical Statistics}, \bold{12}, 663-681.

Ha, I.D. and Lee, Y. (2005). Comparison of hierarchical likelihood versus orthodox best linear unbiased predictor approaches for frailty models. \emph{Biometrika}, \bold{92}, 717-723.

Ha, I.D., Lee, Y. and Song, J.-K. (2001). Hierarchical likelihood approach for frailty models. \emph{Biometrika}, \bold{88}, 233-243.

Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.
}
\author{
Il Do Ha
}

\examples{
\dontrun{
data(cgd)
data(kidney)

#### Frailty model

#### Analysis of Kidney data
kln11 <- Frailty.h(time ~ sex + age + (1|patient)-1, kidney$censor, kidney,
                   RandDist = "Normal", mord = 1, dord = 1)
kln12 <- Frailty.h(time ~ sex + age + (1|patient) - 1, kidney$censor, kidney,
                   RandDist = "Normal", mord = 1, dord = 2)
kg11 <- Frailty.h(time ~ sex + age + (1|patient) - 1, kidney$censor, kidney, 
                  RandDist = "Gamma", mord = 1, dord = 1)
kg12 <- Frailty.h(time ~ sex + age + (1|patient) - 1, kidney$censor, kidney,
                  RandDist = "Gamma", mord = 1, dord = 2)

#### Analysis of CGD data
cgd11 <- Frailty.h(TIME ~ TRT + (1|HOSPITAL) + (1|PATIENT) - 1, cgd$DEL, cgd,
                   RandDist = "Normal", mord = 1, dord = 1)
cgd12 <- Frailty.h(TIME ~ TRT + (1|HOSPITAL) + (1|PATIENT) - 1, cgd$DEL, cgd,
                   RandDist = "Normal", mord = 1, dord = 2)

}           
}
\keyword{HGLM}
\keyword{h-likelihood}

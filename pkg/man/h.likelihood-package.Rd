\name{h.likelihood-package}
\alias{h.likelihood-package}
\alias{h.likelihood}
\docType{package}
\title{
Statistical Modeling and Inference via Hierarchical Likelihood
}
\description{
The package provides a top interface of hierarchical likelihood (h-likelihood) based models. It covers the estimation of hierarchical generalized linear models (HGLMs), and ... 
}
\details{
\tabular{ll}{
Package: \tab h.likelihood\cr
Type: \tab Package\cr
Version: \tab 2010.8.25\cr
Date: \tab 2010-08-25\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
Depends: \tab hglm, HGLMMM\cr
}
}
\author{
Xia Shen and Marek Molas

Maintainer: Xia Shen <xia.shen@lcb.uu.se>
}
\references{
Lee, Y. and Nelder, J.A. (1996). Hierarchical generalized linear models (with discussion). \emph{Journal of the Royal Statistical Society. Series B (Methological)}
\bold{58}, 619-678.

Lee, Y. and Nelder, J.A. (2001). Hierarchical generalised linear models: A synthesis of generalised linear models, random-effect models and structured dispersions.
\emph{Biometrika} \bold{88}, 987-1006. 

Lee, Y., Nelder, J.A., and Pawitan, Y. (2006). \emph{Generalized Linear Models with Random Effects}. Boca Raton: Chapman & Hall/CRC.

Noh, M. and Lee, Y. (2007). REML estimation for binary data in GLMMs. \emph{Journal of Multivariate Analysis} \bold{98}, 896-915.

Ronnegard, L., Shen, X. and Alam, M. (2010). hglm: A Package for Fitting Hierarchical Generalized Linear Models. \emph{The R Journal. \bold{(to appear)}}

Molas, M. and Lesaffre, E. (2010). Hierarchical Generalized Linear Models: the R Package HGLMMM. \bold{Submitted}.
}
\seealso{
\code{\link{hglm-package}}, \code{\link{HGLMMM-package}}
}
\examples{
\dontrun{

data(semiconductor)

# ----- use 'EQL'

h.gamma.normal <- HGLM(fixed = y ~ x1 + x3 + x5 + x6,
                       random = ~ 1|Device,
                       family = Gamma(link = log),
                       disp = ~ x2 + x3, data = semiconductor)
                       
summary(h.gamma.normal)

plot(h.gamma.normal, cex = .6, pch = 1,
     cex.axis = 1/.6, cex.lab = 1/.6,
     cex.main = 1/.6, mar = c(3, 4.5, 0, 1.5))
     
# ----- use 'HL(0,1)'

RSC <- data.frame(int = rep(1, 16))
h.gamma.normal <- HGLM(fixed = y ~ x1 + x3 + x5 + x6,
                       random = ~ 1|Device,
                       family = Gamma(link = log),
                       disp = ~ x2 + x3, data = semiconductor, 
                       method = 'HL01', disp.start = c(0, 0, 0),
                       disp.random = list(one = ~ 1), data.random = list(RSC))
                       
# ----- use 'HL(1,1)'

RSC <- data.frame(int = rep(1, 16))
h.gamma.normal <- HGLM(fixed = y ~ x1 + x3 + x5 + x6,
                       random = ~ 1|Device,
                       family = Gamma(link = log),
                       disp = ~ x2 + x3, data = semiconductor, 
                       method = 'HL11', disp.start = c(0, 0, 0),
                       disp.random = list(one = ~ 1), data.random = list(RSC))
}
}

\name{HGLM}
\alias{HGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Hierarchical Generalized Linear Models
}
\description{
This function fits hierarchical generalized linear models (HGLMs) using various approximation methods.
}
\usage{
HGLM(y = NULL, X = NULL, Z = NULL, X.disp = NULL, 
     family = gaussian(link = identity), 
     random.family = gaussian(link = identity), method = "EQL", 
     conv = 1e-04, maxit = 20, fixed = NULL, random = NULL, 
     disp = NULL, link.disp = "log", disp.random = NULL, 
     data = NULL, data.random = NULL, fix.disp = NULL, 
     Offset = NULL, Weights = NULL, disp.start = 0, binomial.N = NULL, 
     start.fixed = NULL, start.random = NULL, start.disp = NULL, 
     start.disp.random = NULL, info = TRUE, debug = FALSE, 
     contrasts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the dependent variable, only available when \code{method = 'EQL'}.}
  \item{X}{a design matrix for the fixed effects, only available when \code{method = 'EQL'}.}
  \item{Z}{an design matrix for the random effects, only available when \code{method = 'EQL'}.}
  \item{X.disp}{a design matrix for the fixed effects in the dispersion part of the model, only available when \code{method = 'EQL'}.}
  \item{family}{a description of the error distribution and link function to be used in the mean part of the model. (See \code{\link{family}} for details of family functions.)}
  \item{random.family}{a description of the error distribution and link function to be used in the variance part of the model.}
  \item{method}{estimation method, which can be \code{'EQL'}, \code{'HL01'}, or \code{'HL11'}, where \code{'EQL'} can ONLY be used when ONLY ONE random effect term is specified. \code{'EQL'} is the method of interconnected GLMs presented in Lee et al. (2006), and for \code{'HL01'} and \code{'HL11'}, see Lee and Nelder (2001).}
  \item{conv}{convergence criterion, the default is 1e-4, for models with many random effects could be set less strict.}
  \item{maxit}{maximum number of iterations in the IWLS algorithm, only available when \code{method = 'EQL'}.}
  \item{fixed}{a \code{\link{formula}} specifying the fixed effects part of the model, and the format is \code{Response ~ Fixed.Effect.1 + ... + Fixed.Effect.p}.}
  \item{random}{a one-sided \code{\link{formula}} specifying the random effects part of the model, and the format is \code{~ (Random.Effect.1 | Subject.1) + ... + (Random.Effect.q | Subject.q)}.}
  \item{disp}{a one-sided \code{\link{formula}} specifying the fixed effects in the dispersion part of the model, and the format is \code{~ Effect.1 + ... + Effect.N}.}
  \item{link.disp}{the link function for the dispersion part of the model, only available when \code{method = 'EQL'}.}
  \item{disp.random}{a list of one-sided \code{\link{formula}}e for the dispersion strucutre of each random effects, which has the format of \code{list(one = ~ Effect.1.1 + ..., two = ~ Effect.2.1 + ..., three = ..., ...)}, only available when \code{method = 'HL01'} or \code{'HL11'}.}
  \item{data}{the data frame to be used together with \code{fixed} and \code{random}.}
  \item{data.random}{a list of \code{\link{data.frame}}s for \code{disp.random}, which has the format of \code{list(one = data.Random.1, two = data.Random.2,...)}, only available when \code{method = 'HL01'} or \code{'HL11'}.}
  \item{Weights}{prior weights to be specified in weighted regression, only available when \code{method = 'EQL'}.}
  \item{fix.disp}{a numeric value if the dispersion parameter of the mean model is known for example 1 for binomial and Poisson models.}
  \item{Offset}{an offset for the linear predictor of the mean model.}
  \item{disp.start}{(starting) values for the overdispersion structure - vector of length equal to the number of parameters in the overdispersion structure, only available when \code{fix.disp = NULL} and \code{method = 'HL01'} or \code{'HL11'}.}
  \item{binomial.N}{the number of trials for each observation for binomial models.}
  \item{start.fixed}{optional starting values for fixed effects in the mean structure (one vector of numeric values).}
  \item{start.random}{optional starting values for random effects in the mean structure (one vector of numeric values).}
  \item{start.disp}{optional starting values for parameters of dispersion components of the residuals (one vector of numeric values).}
  \item{start.disp.random}{optional starting values for parameters of dispersion components of random effects (one vector of numeric values).}
  \item{info}{a request to display of iteration information if \code{TRUE}, only available when \code{method = 'HL01'} or \code{'HL11'}.}
  \item{debug}{a request to display of iteration mechanism progress in detail if \code{TRUE}, only available when \code{method = 'HL01'} or \code{'HL11'}.}
%  \item{na.action}{see \code{\link{lm}}, only available when \code{method = 'HL01'} or \code{'HL11'}.}
  \item{contrasts}{see \code{\link{lm}}, caution as it is currently not fully developed, only available when \code{method = 'HL01'} or \code{'HL11'}.}
}
\details{
When \code{method = 'EQL'}, all the model checking functions in the \code{\link{hglm-package}} are available on the object returned; Otherwise, all the model checking functions in the \code{\link{HGLMMM-package}} are available on the object returned.
}
\value{
When \code{method = 'EQL'}, an object of class \code{hglm} is returned, see \code{\link{hglm}}; Otherwise, an object of class \code{\link{HGLM}} is returned, see \code{\link{HGLMfit}}.
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
\author{
Xia Shen and Marek Molas
}
\note{
The function provides a unified interface to the \code{\link{hglm-package}} developed by Moudud Alam, Lars Ronnegard and Xia Shen, and the \code{\link{HGLMMM-package}} developed by Marek Molas.
}
\seealso{
\code{\link{hglm-package}}, \code{\link{HGLMMM-package}}, \code{\link{hglm}}, \code{\link{HGLMfit}}.
}
\examples{
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{HGLM}
\keyword{h-likelihood}

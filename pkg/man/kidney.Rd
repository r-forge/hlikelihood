\name{kidney}
\Rdversion{1.1}
\alias{kidney}

\docType{package}

\title{Kidney Infection Data}

\description{The data presented by McGilchrist and Aisbett (1991) consist of times to the first and second recurrences of infection in 38 kidney patients using a portable dialysis machine.}
\usage{data(kidney)}

\format{
  A data frame with 76 observations on the following 9 variables.
  \describe{
    \item{\code{time}}{Survival time: Time to infection since insertion of the catheter}
    \item{\code{censor}}{Censoring indicator(1 = uncensored, 0 = censored).}
    \item{\code{obs}}{There were 76 observations.}
    \item{\code{patient}}{Patient ID: There were 38 patients with two recurrences.}
    \item{\code{age}}{Age of each patient.}
    \item{\code{sex}}{Sex of each patient(1 = male, 2 = female).}
    \item{\code{dy0}}{GN type of disease(1 = Yes, 0 = No).}
    \item{\code{dy1}}{AN type of disease(1 = Yes, 0 = No).}
    \item{\code{dy2}}{PKD type of disease(1 = Yes, 0 = No).}
  }
}


\references{
McGilchrist and Aisbett(1991). Regression with frailty in survival analysis. \emph{Biometrics}, \bold{47}, 461-466.
}
\examples{
data(kidney)
}


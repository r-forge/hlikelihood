\name{cgd}
\Rdversion{1.1}
\alias{cgd}

\docType{package}

\title{Chronic Granulomatous Disease (CGD) Infection Data}

\description{The CGD data set in Fleming and Harrington (1991) consists of a placebo-controlled randomized trial of gamma interferon in chronic granulomatous disease.}
\usage{data(cgd)}

\format{
  A data frame with 203 observations on the following 16 variables.
  \describe{
    \item{\code{OBS}}{There were 203 observations.}
    \item{\code{HOSPITAL}}{Hospitial ID: There were 13 hospotials where trials were conducted.}
    \item{\code{PATIENT}}{Patient ID: There were 128 patients.}
    \item{\code{TIME}}{The recurrent infection times of each pateient from the different hospitals.}
    \item{\code{DEL}}{Censoring indicator(1 = uncensored, 0 = censored).}
    \item{\code{TRT}}{Treatment Code(1 = gamma interferon, 0 = placebo).}
    \item{\code{INHERIT}}{Pattern of inheritance (0 = autosomal recessive, 1 = X-linked).}
    \item{\code{AGE}}{Age of each patient, years.}
    \item{\code{HEIGHT}}{Height of each patient, in cm.}
    \item{\code{WEIGHT}}{Weight of each patient, in kg.}
    \item{\code{STEROIDS}}{Using corticosteroids at times of study centry(1 = Yes, 0 = No).}
    \item{\code{PROPYLAC}}{Using prophylactic antibiotics at time of study entry(1 = Yes, 0 = No).}    
    \item{\code{SEX}}{Sex of each patient(0 = male, 1 = female).}    
    \item{\code{H}}{Hospital region(0 = U.S.A., 1 = Europe).}    
    \item{\code{LONGI}}{A longitudinal variable representing the accumulated time from the firts infection(in years).}        
    \item{\code{ENUM}}{Sequence number. For each patient, the infection records are in sequnce number.}            
  }
}


\references{Fleming and Harrington (1991). \emph{Counting processes and survival analysis}. Wiley: New York.
}
\examples{
data(cgd)
}


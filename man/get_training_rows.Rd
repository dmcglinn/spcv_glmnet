\name{get_training_rows}
\alias{get_training_rows}
\title{Get the row indices for the training data sets}
\usage{
get_training_rows(coords, dist_thres, longlat = FALSE)
}
\arguments{
  \item{coords}{2 dimensional matrix of spatial
  coordinates}

  \item{dist_thres}{the threshold distance, if longlat =
  TRUE this is in units of kilometers}
}
\description{
For each sample in the dataset a set of indices for the
training data to compare to that sample is provided based
on a threshold distance. If the threshold distance is zero
then all other samples are included in the training set.
}
\examples{
coords = expand.grid(1:4, 1:4)
training_rows = get_training_row(coords, 3)
coords[training_rows[[i]]
}


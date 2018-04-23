\name{plot_peaks}
\alias{plot_peaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the peaks of a MS or NMR dataset.
}
\description{
Function returns a plot where each point represents the intensity of a peak in a sample. Peaks are coloured according to a metadata class. 
}
\usage{
plot_peaks(dataset, column.class, samples = NULL, variable.bounds = NULL,
xlab = NULL, ylab = NULL, legend.place = "topright", cex = 0.8,
reverse.x = F, p.size=0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
    list representing the dataset from a metabolomics experiment.
}
  \item{column.class}{
    string indicating the metadata's variable.
}
  \item{samples}{
    vector with samples' names, if NULL all the samples will be considered.
}
  \item{variable.bounds}{
    numeric vector with two elements indicating the interval of x-values to plot.
}
  \item{xlab}{
    x-axis label.
}
  \item{ylab}{
    y-axis label.
}
  \item{legend.place}{
    string indicating the place that the legend's box will be placed.
}
  \item{cex}{
    numeric value that indicates the amount by which the legend is magnified relative to the default.
}
  \item{reverse.x}{
    boolean value indicating if the x-axis will be shown reversed or not.  
}
  \item{p.size}{
  numeric value indicating the amount by which the plot points are magnified relative to the default.
}
  \item{\dots}{
    additional parameters to matplot.
}
}
\examples{
  data(propolis)
  plot_peaks(propolis, "seasons", variable.bounds = c(0,3), samples=c("XX_au", "XX_sm", "XX_wi"))
}

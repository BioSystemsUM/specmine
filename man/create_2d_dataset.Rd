\name{create_2d_dataset}
\alias{create_2d_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create 2D dataset
}
\description{
  Creates a 2D dataset from existing objects.
}
\usage{
create_2d_dataset(list_2d, type = "undefined", metadata = NULL, description = "", 
                  sample.names = NULL, F1 = NULL, F2 = NULL, label.x = NULL,
                  label.y = NULL, label.values = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list_2d}{
  A list of matrices where each matrix represents a 2D spectra from one sample.
  }
  
  \item{type}{
  Type of data: string that can be one of the following:
  \itemize{
    \item 2d-nmr
  }}
  
  \item{metadata}{
  Data frame with the dataset's metadata: columns represent each metadata variable and rows represent the value of the metadata for the sample.
  }
  
  \item{description}{
  String with a short description of the dataset.
  }
  
  \item{sample.names}{
  Vector with sample names, if NULL then the names of the 2D list or sequential numbers will be used.
  }
  
  \item{F1}{
  Vector of the indirect dimension' ppm values.
  }
  
  \item{F2}{
  Vector of the diret dimension' ppm values.
  }
  
  \item{label.x}{
  Label for the x axis.
  }
  
  \item{label.y}{
  Label for the y axis.
  }
  
  \item{label.values}{
  Label for the variable represented by a pair (x,y).
  }
  
}

\value{
  List representing the 2D dataset:
  \item{data}{
  A list of matrices where each matrix matches one 2D spectra.
  }
  \item{type}{
  The type of the data in the dataset.
  }
  \item{description}{
  A short text describing the dataset.
  }
  \item{metadata}{
  A dataframe with the metadata variables.
  }
  \item{F1_ppm}{
  The ppm values regarding indirect dimension.
  }
  \item{F2_ppm}{
  The ppm values regarding direct dimension.
  }
  \item{labels}{
  A list of labels for the x, y and pairs'(x,y) values.
  }
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ 2d-dataset }


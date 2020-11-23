#convert transmittance to absorbance
transmittance_to_absorbance = function(dataset, percent = TRUE){
  datamat = dataset$data
  if (!percent){
    absorbance.datamat = 2-log10(datamat*100)
  } else {
    absorbance.datamat = 2-log10(datamat)
  }
  dataset$data = absorbance.datamat
  dataset$labels$val = "A"
  dataset
}

#convert absorbance to transmittance
absorbance_to_transmittance = function(dataset){
  datamat = dataset$data
  transmittance.datamat = 100*(10^(-datamat))
  dataset$data = transmittance.datamat
  dataset$labels$val = "%T"
  dataset
}

# Smoothing - specmine (specmine.bin) and hyperspec (spc.loess)
smoothing_interpolation = function(dataset, method = "bin", reducing.factor = 2, x.axis = NULL, p.order = 3, window = 11, deriv = 0, na.rm = TRUE){
  if (method == "bin") {
		dataset = smoothing_spcbin(dataset, reducing.factor, na.rm = na.rm)
	} 
  else if (method == "loess") {
		dataset = smoothing_spcloess(dataset, x.axis)
	} else if (method == "savitzky.golay"){
		dataset = savitzky_golay(dataset, p.order, window, deriv)
	}
	dataset
}

# Specmine binning smoothing interpolation
smoothing_spcbin <- function(dataset, reducing.factor = 2, na.rm = TRUE) {
  res.dataset <- specmine.bin(dataset, reducing.factor, na.rm = na.rm)
  res.dataset$description <- paste(dataset$description, "smoothed with specmine bin", sep="-")
  res.dataset$type <- dataset$type
  res.dataset
}

# Implementation taken from hyperspec spc.bin function
specmine.bin <- function(dataset, reducing.factor = 2, na.rm = TRUE) {
  
  new_dataset <- dataset
  
  n <- ceiling( nrow(new_dataset$data) / reducing.factor)
  
  small <- nrow(new_dataset$data) %% reducing.factor
  
  if (small != 0){
    warning(paste(c("Last data point averages only ", small, " points.")))
  }
  
  bin <- rep(seq_len(n), each = reducing.factor, length.out = nrow(new_dataset$data))
  
  na <- is.na(new_dataset$data)
  
  
  if ((na.rm > 0) && any(na)) {
    if (na.rm == 1) {
      na <- apply(!na, 2, tapply, bin, sum, na.rm = FALSE)
      new_dataset$data <- t (apply(new_dataset$data, 2, tapply, bin, sum, na.rm = TRUE) / na)
      
    } else {
      tmp <- t (apply (new_dataset$data, 2, tapply, bin, sum, na.rm = FALSE))
      tmp <- sweep (tmp, 2, rle(bin)$lengths, "/")
      
      na <- which(is.na(tmp), arr.ind = TRUE)
      bin <- split(wavelength.seq(new_dataset, bin))
      
      for (i in seq_len(ncol(na))) {
        tmp [na [i, 2], na [i, 1]] <- mean(new_dataset$data[bin[[na[i,2]]], na[i, 1]], na.rm = TRUE)
      }
      
      new_dataset$data <- tmp
    }
  } else {
    new_dataset$data <- t (apply(new_dataset$data, 2, tapply, bin, sum, na.rm = FALSE))
    new_dataset$data <- sweep (new_dataset$data, 2, rle(bin)$lengths, "/")
  }
  
  new_dataset$data <- t(new_dataset$data)
  
  new_wl <- wavelengths(dataset, bin, na.rm = na.rm)
  
  rownames(new_dataset$data) <- new_wl
  
  new_dataset$xSet <- NULL
  
  warning("xSet is NULL because bthis process reseted wavelength values.")
  
  new_dataset
}

#Function to provide sequence of wavelength values, taken from hyperspec wl.seq implementation
wavelength.seq <- function(dataset, from = 1, to = nrow(dataset$data), ...){
  if (nrow(dataset$data) == 0) {
    integer(0)
  } else {
    seq (from = from, to = to, ...)
  }
}

#Function to adjust wavelenght values due to binning
wavelengths <- function(dataset, bin, na.rm = TRUE){
  wl <- suppressWarnings(as.numeric(rownames(dataset$data)))
  
  if (all(is.na(wl))) {
    splits <- strsplit(rownames(dataset$data), "/")
    first <- c()
    second <- c()
    for (i in 1:length(splits)){
      first <- c(first, splits[[i]][1])
      second <- c(second, splits[[i]][2])
    }
    new_first <- as.character(as.numeric (tapply(as.numeric(first), bin, mean, na.rm = na.rm > 0)))
    new_second <- as.character(as.numeric (tapply(as.numeric(second), bin, mean, na.rm = na.rm > 0)))
    
    new_wl <- paste(new_first, new_second, sep = "/")
    
  } else {
    new_wl <- as.character(as.numeric (tapply(wl, bin, mean, na.rm = na.rm > 0)))
  }
  
  new_wl
}


# Specmine loess smoothing interpolation 
smoothing_spcloess <- function(dataset, x.axis = NULL){
  # In cases where rownames are strings not convertible to numeric
  wl <- suppressWarnings(as.numeric(rownames(dataset$data)))
  if (all(is.na(wl))) {
    new_wl <- seq(from = 1, to = nrow(dataset$data))
  } else {
    new_wl <- wl
  }
	if (is.null(x.axis)){
		res.dataset <- specmine.loess(dataset, newx = new_wl, na.rm = TRUE)
	} else {
		res.dataset = specmine.loess(dataset, newx = x.axis, na.rm = TRUE)
	}
  res.dataset$description = paste(dataset$description, "smoothed with specmine loess", sep="-")
  res.dataset$type = dataset$type
  res.dataset
}

# Implementation taken from hyperspec spc.loess function
specmine.loess <- function(dataset, newx, enp.target = nrow(dataset$data) / 4, surface = "direct", ...) {
  
  .loess <- function (y, x) {
    if (all (is.na (y))) {
      NA
    } else {
      loess (y ~ x, enp.target = enp.target, surface = surface, ...)
    }
  }
  
  .predict <- function (loess, x) {
    if (! methods::is(loess, "loess") && is.na(loess)) {
      rep (NA_real_, length(x))
    } else {
      predict (loess, x)
    }
  }
  
  wl <- suppressWarnings(as.numeric(rownames(dataset$data)))
  if (all(is.na(wl))) {
    new_wl <- seq(from = 1, to = nrow(dataset$data))
  } else {
    new_wl <- wl
  }
  
  loess <- apply (dataset$data, 2, .loess, new_wl)
  
  dataset$data <- t (sapply(loess, .predict, newx))
  
  dataset$data <- t(dataset$data)
  
  rownames(dataset$data) <- newx
  
  if (any (is.na (dataset$data))) {
    warning ("NAs were generated. Probably newx was outside the spectral range covered by spc.")
  }
  
  dataset
}



savitzky_golay = function(dataset, p.order, window, deriv = 0){
    if (window %%2 != 1 || window < 0) 
        stop("window size (window) must be a positive odd number")
    if (p.order >= window) 
        stop("window size (window) is too small for the polynomial order (p.order)")
    if (p.order < deriv) 
        stop("polynomial order p (p.order) should be geater or equal to differentiation order (deriv)")
    X = t(dataset$data)
    half_window = (window -1)/2
    b = outer(-half_window:half_window, 0:p.order, "^")
    A = MASS::ginv(b)
    result = matrix(data = 0, ncol=ncol(X),nrow=nrow(X))
    for (i in 1:nrow(X)){
        first.values = X[i,1] - abs( X[i,1:(half_window)] - X[i,1] )
        last.values = tail(X[i,], n = 1) + abs(X[i,(ncol(X)-half_window+1):ncol(X)] - tail(X[i,],n=1))
        all = c(first.values, X[i,], last.values)
        result[i,] = factorial(deriv) * convolve(all, A[deriv+1,], type="f")
    }
    colnames(result) = colnames(X)
    dataset$data = t(result)
    rownames(dataset$data) = colnames(X)
    colnames(dataset$data) = rownames(X)
    dataset$description = paste(dataset$description, "smoothed with savitzky-golay filter", sep = "-")
    dataset
} 


# DATA CORRECTION - functions to do spectra correction

"data_correction" = function(dataset, type = "background", method = "modpolyfit", ...){
	if (type == "background"){
		dataset = background_correction(dataset)
	} 
  else if (type == "offset"){
		dataset = offset_correction(dataset)
	} 
  else if (type == "baseline"){
		dataset = baseline_correction(dataset, method, ...)
	} 
	dataset
}

background_correction <- function(dataset) {
  background <- apply(dataset$data, 1, quantile, probs = 0.05)
  dataset$data <- sweep(dataset$data, 1, background, "-")
  dataset$description <- paste(dataset$description, "background correction", sep="; ")
  dataset
}


offset_correction <- function(dataset) {
  offsets <- apply(dataset$data, 2, min)
  dataset$data <- sweep(dataset$data, 2, offsets, "-")
  dataset$description <- paste(dataset$description, "offset correction", sep="; ")
  dataset
}

# ... - extra parameters to baseline function
baseline_correction = function(dataset, method = "modpolyfit", ...){
	rnames = rownames(dataset$data)
	cnames = colnames(dataset$data)
	samples.df = t(dataset$data)
	bl = baseline::baseline(samples.df, method = method, ...)
	samples.df = baseline::getCorrected(bl)
  dataset$data = t(samples.df)
	rownames(dataset$data) = rnames
	colnames(dataset$data) = cnames
  dataset$description = paste(dataset$description, "baseline correction", sep="; ")
	dataset
}

# shifting spectra
# method - "constant" - uses a constant shift that is added to the x.values
# method - "interpolation" - uses interpolation - linear or spline according to "interp.function"
# shift.val - value of the shift (for constant and interpolation methods); can be a single value for all spectra
#			  "auto" - shifts are automatically determined
# or a vector of length = number of samples; can also be the string "auto" for automatic calculation of shifts
"shift_correction" = function(dataset, method = "constant", shift.val = 0, interp.function = "linear",
                              ref.limits = NULL) {
  
  x.vals = get_x_values_as_num(dataset)
  
  if(! length(shift.val) %in% c(1,num_samples(dataset)) ) {
    stop("Shift.val parameter has incorrect size: should be 1 or number of samples in the dataset")
  }
  else if (length(shift.val) == 1) {
	if (shift.val == "auto") {
        if (is.null(ref.limits) | length(ref.limits) != 2) {
			stop("Parameter ref.limits incorrect for automatic determination of shifts")
		}
        else { 
			shift.val = calculate_shifts(dataset, ref.limits)
		}
	}
  }
  if (method == "constant") {
    new.x.values = x.vals + shift.val
    dataset = set_x_values(dataset, new.x.values)
  }
  else if (method == "interpolation") {
    if (interp.function == "spline") {
      interp_fn = function(data, shift, x.values) {
        spline (x.values + shift, data, xout = x.values, method = "natural")$y
      }     
    }
    else if(interp.function == "linear") {
      interp_fn = function(data, shift, x.values) {
        approx(x.values + shift, data, xout = x.values, method = "linear")$y
      }
    }
    else stop("Interpolation function not defined")
    
    if (length(shift.val) == 1)
      newdata = apply(dataset$data, 2, interp_fn, shift = shift.val, x.values = x.vals)
    else {
      newdata = matrix(NA, nrow(dataset$data), ncol(dataset$data))
      for (i in 1:length(shift.val))
        newdata[,i] = interp_fn(dataset$data[,i], shift.val[i], x.vals)
    }
    rownames(newdata) = rownames(dataset$data)
    colnames(newdata) = colnames(dataset$data)
    dataset$data = newdata
  }
  else stop("Method is not defined")
  
  dataset
}

# calculate shifts based on a band of the spectra (see hyperSpec vignette sect. 12.2.1)
"calculate_shifts" = function(dataset, ref.limits = NULL)
{ 
  #x.vals = get.x.values.as.num(dataset)
  dataM = subset_x_values_by_interval (dataset, ref.limits[1], ref.limits[2])
  #xvalsM = x.vals[x.vals >= ref.limits[1] & x.vals <= ref.limits[2]]
  bandpos = apply (t(dataM$data), 1, find_max, get_x_values_as_num(dataM))
  refpos = find_max (colMeans(t(dataM$data)), get_x_values_as_num(dataM))
  refpos - bandpos
}

"find_max" = function (y, x){
  pos = which.max (y) + (-1:1)
  X = x [pos] - x [pos [2]]
  Y = y [pos] - y [pos [2]]
  X = cbind (1, X, X^2)
  coef = qr.solve (X, Y)
  - coef [2] / coef [3] / 2 + x [pos [2]]
}

# multiplicative scatter correction

"msc_correction" = function(dataset) {
  temp = t(dataset$data)
  newdata = pls::msc(temp)
  dataset$data = t(newdata)
  dataset
}

# first derivative

first_derivative = function(dataset) {
  new.data = apply(dataset$data, 2, diff)
  dataset$data = new.data
  dataset
}

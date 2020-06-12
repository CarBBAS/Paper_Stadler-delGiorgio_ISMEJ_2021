# quality control
# we do not trust single observations by bin, correct to 0
# before calculating relative abundances
bin_qual_control <- function(data, x.var, y.var, keep, bins){
  if(!"data.table" %in% (.packages())){
    library(data.table)
  }
  # make data smaller and convert to data.table
  if(any(class(data) == "data.table")){
    data <- data[,c(keep, x.var, y.var), with = F]
  } else if(class(data) == "data.frame"){
    setDT(data)
    data <- data[,c(keep, x.var, y.var), with = F]
  } else {
    stop("'data' needs to be of class 'data.frame' or 'data.table'")
  }
  data$interval <- cut(data[[x.var]], breaks = bins, 
                       include.lowest = TRUE)
  data[, n := .N, by = interval]
  df <- data[, .(n = unique(n)), by = interval]

  setDT(data); setDT(df)
  # extract interval of bins with only one observation
  int <- df[df$n > 1,]$interval
  # run loop if several single observations are per bin
  if(length(int) != 0L){
    for(i in 1:length(int)){
      # extract only the numbers as vector
      n.int <- as.numeric(unlist(regmatches(int[i], gregexpr("-?[[:digit:]]+", int[i]))))
      # extract number of observations
      #bin.n <- df$n[df$interval == int[i]]
  
      # include condition, ( means that it does not include the first interval number for binning
      # set to >
      # include condition, [] means that it does include the first interval number for binning
      # set to >=
      if(substr(int[i],1,1) == "("){
        # number of actual observations
        n.obs <- data[data[[x.var]] > n.int[1] &
                        data[[x.var]] <= n.int[2] &
                        data[[y.var]] != 0,
                      .N]
        if(n.obs == 1){
          # overwrite single observation by 0
          # if read number is smaller than 10
          data[data[[x.var]] > n.int[1] &
                 data[[x.var]] <= n.int[2] &
                 data[[y.var]] < 10,
               eval(y.var) := 0]
        }
      } else {
        # number of actual observations
        n.obs <- data[data[[x.var]] >= n.int[1] &
                        data[[x.var]] <= n.int[2] &
                        data[[y.var]] != 0,
                      .N]
        if(n.obs == 1){
          # overwrite single observation by 0
          # if read number is smaller than 10
          data[data[[x.var]] >= n.int[1] &
                 data[[x.var]] <= n.int[2] &
                 data[[y.var]] < 10,
               eval(y.var) := 0]
        }
      }
    }
  }
  setDF(data); return(data)
}

# Identify positions of local maxima
# this function is very sensitive some quality control has to be done afterwards
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  y
  # return locations of local maxima
}


# Bin values by desired interval
# extracts peaks of binned values
# returns [[1]] raw x y data, [[2]] binned x y data, [[3]] peaks as list
bin_peak <- function(data, x.var, y.var, keep = NULL, bins){
  if(!"data.table" %in% (.packages())){
    library(data.table)
  }
  # make data smaller and convert to data.table
  if(any(class(data) == "data.table")){
    data <- data[,c(keep, x.var, y.var), with = F]
  } else if(class(data) == "data.frame"){
    setDT(data)
    data <- data[,c(keep, x.var, y.var), with = F]
  } else {
    stop("'data' needs to be of class 'data.frame' or 'data.table'")
  }
  # make data smaller
  data$interval <- cut(data[[x.var]], breaks = bins, 
                                              include.lowest = TRUE)
  
  data[, n := .N, by = interval]
  data[get(y.var) > 0 & n >= 3, wg := 9][get(y.var) == 0 & n >= 3, wg := 1]
  data[get(y.var) > 0 & n < 3, wg := 1][get(y.var) == 0 & n < 3, wg := 9]
  
  df <- data[, .(w.mean = weighted.mean(get(y.var), wg),
           mean = mean(get(y.var)),
           min = min(get(y.var)),
           max = max(get(y.var))), by = interval]

  n.int <- regmatches(df$interval, gregexpr("-?[[:digit:]]+", df$interval))
  x <- unlist(lapply(n.int, function(x){
    mean(as.numeric(x))
  }))
  
  y <- df$w.mean # mean of bins
  # find peaks
  i.max <- localMaxima(y)
  
  # return empty data frame if no peaks were found
  if (length(i.max) == 0L | sum(y[i.max]) == 0) {
    list(
      data = data,
      binned = data.frame(
        x = x,
        y = y,
        min = df$min,
        max = df$max
      ),
      pks.df = data.frame(
        no.peaks = 0,
        peaks = NA,
        x = NA,
        y = NA,
        slope.b = NA,
        slope.a = NA
      )
    )
  } else {
    # remove peaks that are 0
    if (any(y[i.max] <= min(y))) {
      i.max <- i.max[which(!y[i.max] <= min(y))]
    }
    
    # Get slopes along x-axis
    slopes <- diff(y) / diff(x)
    # find slopes before and after peak
    if (i.max[1] == 1) {
      before <-
        c(NA, slopes[i.max - 1]) # add NA in case a peak is at the beginning
    } else {
      before <- slopes[i.max - 1]
    }
    after <-
      slopes[i.max + 1] # NA is automatically added if the peak is in the end
    
    list(
      data = data,
      binned = data.frame(
        x = x,
        y = y,
        unweigh.mean = df$mean,
        min = df$min,
        max = df$max
      ),
      pks.df = data.frame(
        no.peaks = rep(length(i.max), times = length(i.max)),
        peaks = c(1:length(i.max)),
        x = x[i.max],
        y = y[i.max],
        slope.b = before,
        slope.a = after
      )
    )
  }
}



## Inegrate function that adjusts to data
adapt_integrate <- function(x, f, abs.tol = .Machine$double.eps^0.25, subdiv = 400){
  
  area <- try(integrate(f, min(x), max(x), stop.on.error = F,
                        abs.tol = .Machine$double.eps^0.25 , subdivisions = subdiv))
  
  if(area$message == "maximum number of subdivisions reached"){
    while(area$message != "OK"){
      subdiv <- subdiv + 50
      area <- try(integrate(f, min(x), max(x), stop.on.error = F, subdivisions = subdiv))
      if(subdiv > 10000 & area$message != "OK"){
        area$message <- "Subdivisions reached 10000. No convergence."
        break
      }
    }
  } else if(area$message == "extremely bad integrand behaviour") {
    while(area$message != "OK"){
      abs.tol <- abs.tol * 10
      area <- try(integrate(f, min(x), max(x), stop.on.error = F,
                            abs.tol = abs.tol, subdivisions = subdiv))
      if(abs.tol > 0.1 & area$message != "OK"){
        area$message <- "Absolute tolerance reached 0.1. No convergence."
        break
      }
    }
  }
  return(area)
}


dist_overlap <- function(data, x.var, y.var, treat.var, .id = NULL){
  # make treat.var as .id if not given
  if(is.null(.id) == T){
    .id <- treat.var
  }
  # extract factor levels of treat.var
  treat <- levels(factor(data[[treat.var]]))
  if(length(treat) > 2L){
    stop("'treat.var' has more than two factors")
  }
  if(length(treat) < 2L){
    warning("'treat.var' has only one factor. No overlap can be calculated.")
    data.frame(Type = .id, 
               between = NA, 
               overlap = NA,
               ratio = NA,
               all = NA,
               mess.bet = NA, mess.ov = NA, mess.all = NA)
  } else {
    # extract x and y variables for two factors
    x1 <- data[[x.var]][data[[treat.var]] == treat[1]]
    x2 <- data[[x.var]][data[[treat.var]] == treat[2]]
    y1 <- data[[y.var]][data[[treat.var]] == treat[1]]
    y2 <- data[[y.var]][data[[treat.var]] == treat[2]]
    
    # handling error if vectors are not same length
    if(length(x1) != length(x2)){
      if(length(x1) > length(x2)){
        y1 <- y1[x1 %in% x2]
        x1 <- x1[x1 %in% x2]
      } else if(length(x2) > length(x1)){
        y2 <- y2[x2 %in% x1]
        x2 <- x2[x2 %in% x1]
      }
    }
    # stop if x axes of both variables do not match
    if(all(x1 %in% x2) == F){
      stop("x-axes of factors do not match")
    }
    
    # area between two curves (non-overlap)
    f1 <- approxfun(x1, y1-y2)
    f2 <- function(z){abs(f1(z))}
    bet <- adapt_integrate(x1, f2)
    
    # area of overlap
    f3 <- approxfun(x1, pmin(y1,y2))
    over <- adapt_integrate(x1, f3)
    
    # all area
    f4 <- approxfun(x1, pmax(y1,y2))
    all <- adapt_integrate(x1, f4)
    
    data.frame(Type = .id, 
                      between = bet$value, 
                      overlap = over$value,
                      ratio = over$value / bet$value,
                      all = all$value,
                      mess.bet = bet$message, mess.ov = over$message, mess.all = all$message)
  }
}



####
expo.fit <- function(data, x.var, y.var) {
  model <- drc::drm(as.formula(paste(y.var, x.var, sep = "~")),
                    fct = aomisc::DRC.expoGrowth(),
                    data = data)
  
  out <- list(
    fitted = data.frame(
      fitted = fitted(model),
      x.var = model$data[[x.var]],
      y.var = model$data[[y.var]],
      stringsAsFactors = F
    ),
    coefficients = data.frame(
      Equation = "a * exp(b * x)",
      a = coef(model)[1],
      b = coef(model)[2],
      stringsAsFactors = F
    )
  )
  colnames(out[["fitted"]])[2:3] <- c(x.var, y.var)
  return(out)
}



select_newest <- function(path, file.pattern, by = NULL) {
  # some error dealing
  if (length(path) > 1) {
    stop("'path' needs to be a single string.")
  } else if(length(path) == 0L){
    stop("'path' is empty.")
  }
  if (length(file.pattern) > 1) {
    stop("'file.pattern' needs to be a single string.")
  } else if(length(file.pattern) == 0L){
    stop("'file.pattern' is empty.")
  }
  
  # read in files in directory
  files <- list.files(path, pattern = file.pattern)
  
  if(length(files) == 0L){
    stop("No file that matches pattern exists.")
  }
  
  if(!is.null(by)){
    # do we have several files per object? -> take newest version
    for (h in 1:length(by)) {
      versions <- files[grepl(by[h], files)]
      versions.sans <-tools::file_path_sans_ext(files)
      if (length(versions) > 1) {
        find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
          col.names <- paste0("V", seq(length = length(x)))
          as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
        }))
        dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
          is.date <- try(as.Date(x), silent = T)
          if(inherits(is.date, "try-error")){
            x <- NULL
          }
          return(x)
        }))))
        
        newest <- max(dates)
        
        files <-
          c(versions[grepl(newest, versions)], files[!grepl(by[h], files)])
      }
    }
  } else {
    versions <- files
    versions.sans <-tools::file_path_sans_ext(files)
    find.date <- dplyr::bind_rows(lapply(stringr::str_split(versions.sans,"_"), function(x){
      col.names <- paste0("V", seq(length = length(x)))
      as.data.frame.list(x, col.names = col.names, stringsAsFactors = F)
    }))
    dates <- as.Date(unlist(Filter(Negate(is.null), apply(find.date, 2, function(x){
      is.date <- try(as.Date(x), silent = T)
      if(inherits(is.date, "try-error")){
        x <- NULL
      }
      return(x)
    }))))
    
    newest <- max(dates)
    files <- versions[grepl(newest, versions)]
  }
  return(files)
}


source("C:/Users/phils/Documents/DSP/R/Common/Common/common.R")

#' combineruns Combines files matching pattern into a single data frame.
#'
#' combineruns
#'
#' @param directory Passed to `list.files` as a directory.
#' @param pattern Passed to `list.files` as `"pattern.(.*).scalars.csv"`
#' @param verbose Will list what it is doing.
#'
#' @return Invisibly returns a data frame with each file's data appended.
#' The form of the data #' frame will depend on the form of the CSV file.
#' @export
#'
#' @examples
#' x = combineruns(  ".", "S1\\.0014\\.") # reads files "S1.0014.*.scalars.csv"
#'
combineruns <- function( directory, pattern, verbose=FALSE){
   pattsplit <- unlist( strsplit( patt, split="\\\\."))
   # cat( "patsplit: ", pattsplit[1], "& ", pattsplit[2], "\n")
   runname <- paste(pattsplit[1], pattsplit[2], sep = ".") # bare run name
   rfdtname <- paste( runname, "rfdt.RDS", sep = ".") # rftd file name
   rfdtpath <- paste( directory, rfdtname, sep="/")
   crname <- paste( runname, "cr.RDS", sep = ".")
   crpath <- paste( directory, crname, sep = "/")
   # cat("fullpath: ", fullpath, "\n")
   # find old data.
   if (file.exists( rfdtpath)) {
      rfdt <- readRDS( rfdtpath) # read the run file data table rfdt
      compare <- TRUE # set as flag to compare file list to rfdt.
      sdf <- readRDS( crpath)
   } else {
      rfdt <- data.table()
      compare <- FALSE
      sdf <- data.table::data.table() # init to empty data.table
   } # end if (file.exists...)

   # get list of all files that match pattern.
   patt <-  paste( pattern,
                   "(.*).scalars.csv", sep=".")
   if (verbose) cat( "Searching Directory: ", directory,
                     " for pattern:", patt, "\n")
   files <- list.files( directory, pattern = patt, full.names = TRUE)
   filetimes <- file.mtime( files)
   stopifnot( !anyNA( filetimes)) # if any return NA, the sort will fail.
   files <- files[ order( filetimes)]
   if (compare) {
      # find files that haven't been read yet.
      # rfdt should have been sorted by mtime, so last record should have
      # latest time
      lt <- last( rfdt)$time
      files <- files[ filetimes > lt]
      filetimes <- filetimes[ filetimes > lt]
   } # end if (compare) ...
   lenf <- length( files)
   stopifnot( lenf > 0) # in ordinary use this should not occur.
   if (verbose) cat( "Files found: ", files, "\n")
   stepoffset <- 0

   stepj <- rep( 0, lenf)
   for (j in 1:lenf) { # walk through all files in group
      data <- data.table::fread( files[ j])
      steps <- nrow( data)
      stepj[ j] <- steps
      #data[ "step"] <- data[ "step"] + stepoffset
      data[, step := step + stepoffset] #dt method of modifying columns
      stepoffset <- stepoffset + steps
      sdf <- data.table::rbindlist( list(sdf, data))
   } # end for (j ...
   newrfdt <- data.table( file = files, time = filetimes, steps = stepj)
   rfdt <- data.table( list( newrfdt, rfdt))
   saveRDS( rfdt, file = rfdtpath) # save for next time.
   saveRDS(  sdf, file = crpath)
   invisible( sdf)
}

#' Call `combineruns` and save RDS file; read the file if it already exists.
#'
#' @param output Filename that will be combined with the directory.
#' @param directory,pattern,verbose Same as `combineruns`.
#'
#' @return
#' @export
#'
#' @examples
crsave <- function( output, directory, pattern, verbose=FALSE){
   stopifnot( endsWith( output, ".RDS"))
   fullpath <- paste( directory, output, sep="/")
   if (file.exists( fullpath)) {
      cat( 'Reading from file ', fullpath)
      sdf <- readRDS( fullpath) }
   else {
      cat( "Combining Files from ", directory, pattern)
      sdf <- combineruns( directory, pattern, verbose)
      saveRDS( sdf, file = fullpath)
   }
   return( sdf)
}

#' plotMD2 Plots a run of data step-by-step.
#'
#' Plots all columns from 4 up.
#'
#' @param data A data frame, usually from `combineruns`.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' x = combineruns(  ".", "S1\\.0014\\.", verbose=TRUE)
#' plotMD2( x)
plotMD2 <- function( data, dt = 1) {
   xlab = if (dt == 1) "Steps" else "Time"
   steps = dt * data[,"step"]
   data$psi <- Mod( complex( real=data$Rpsi, imaginary=data$Impsi))
   nc = ncol(data)
   rows <- ceiling( (nc - 3)/2)
   op <- par( mfrow = c(rows, 2), mar = c(3, 3, 0.5, 0.5), tcl = -0.3,
              mgp = c(1.7, 0.4, 0) )
   #if (nc > 8) paste( "plotMD2 will only plot the columns 4 through 8")
   pl = names( data[ 4:nc])
#   for (v in pl[ c( 1, 4, 2, 3, 5, 6)]) {
   for (v in pl) {
      plot( steps, data[, v], type="l", xlab=xlab, ylab=v, panel.first=grid())
   }
   par( op)
   invisible( NULL)
}

#' redplot Reads an MD2 csv file and plots all data (columns 4 up).
#'
#' @param filename Name of an MD2 produces csv file.
#'
#' @return The data read is returned invisibly.
#' @export
#'
#' @examples
#' setwd("H:/Data/DSP/MD/Data/1980")
#' files = Sys.glob( file.path( "S2.*.scalars.csv"))
#' i <- 40
#' files[i]
#' readplot( files[i])
#'
readplot <- function( filename) {
   data = read.csv( filename)
   plotMD2( data)
   invisible( data)
}

#' colStats Applies a function across all columns adds the function to the name.
#'
#' Executes `lapply( x, f, ...)` to apply `f` to each column of `x`.
#'
#' @param x Typically, `x` is a data frame, but may be anything that `lapply`
#' can process.
#' @param FUN A function that will return one value when handed a vector or list.
#' @param tag A string that will by appended to each column name.
#' @param ... optional arguments for `FUN`.
#'
#' @return A named with each element named.
#' @export
#'
#' @examples
#' colStats(mtcars, mean, 'mean')
#' # mpg.mean   cyl.mean  disp.mean    hp.mean  drat.mean    wt.mean ...
#' # 20.090625   6.187500 230.721875 146.687500   3.596563   3.217250...
#'
colStats <- function( x, FUN, tag, ...) {
  t <- unlist( lapply( x, FUN, ...))
  names(t) <- paste( names( t), tag, sep=".")
  t
}

#' Create block averages and re-sample.
#'
#' @param x A vector to re-sample.
#' @param n An integer number of samples in the block averages.
#'
#' @return A vector of the re-sampled data.
#' @export
#'
#' @examples
#' rsmpavg( rnorm(1000), 10) # results in a length 100 vector.
#'
rsmpavg <- function(x, n){
   # re-sample x by n point averages, returning one sample per average.
   # truncates if length of x is not multiple of n
   outlen <- floor( length(x) / n)
   x <- x[1:(n * outlen)]
   dim(x) <- c( n, outlen)
   apply( x, 2, mean) # essentially, column means
}

#' lmMD2 Computes the regression of the ke vs step
#'
#'
#'
#' @param df a data frame with columns "ke" & "step".
#'
#' @return a named vector with the model's linear coefficients.  Retains the
#' names in `summary.lm`, thus, is of the form
#' `c( Estimate=e, "Std. Error"=s, r.squared=r`.  Note the space in
#' `"Std. Error"`
#' @export
#'
#' @examples
lmMD2 <- function( df) {
   mod <- summary( lm( ke ~ step, df))
   coeff <- unname(mod$coefficients["step", c('Estimate', 'Std. Error')])
   t <- mod$fstatistic
   Fpvalue <-  pf( t['value'], t['numdf'], t['dendf'], lower.tail=FALSE)
   t <- mod$fstatistic
   c( nr = nrow( df), Estimate = coeff[1], Std.Error = coeff[2],
      r.squared = mod$r.squared, Fval =  unname( t['value']),
      Fpvalue = unname( Fpvalue), df = t['dendf'])
}

#' Finds a reasonable length for an LM fit to a full run.
#'
#' `findlmlen` will first look at the whole record and see if it meets the
#' `minp`criteria.  If that fails, then it will look at the last 2/3rds of
#' the record (the tail).  If that works, it will do a regression on the
#' first third
#' to get a direction (rising or falling), and will then find the first
#' point that is above the max of the tail (for falling) or below the min of
#' the tail (for rising), and add all the points up to that point to the
#' regression.  It will return the resulting regression. If it fails to find
#' a valid regression, it will return the the first (whole record)
#' regression (the failure indicated by the `Fpvalue` that is below `minp`).
#'
#' @param df Ordinarily a re-sampled ke data frame with items ke and step.
#' @param minp A threshold to test if the fit is considered valid. `minp = 0.5`
#'    means that the probability that the fit was to purely random (no trend)
#'    data is 50% or greater.
#'
#' @return
#' @export
#'
#' @examples
#'    Rke <- rsmpavg( ke, rsblk)
#'    Rkedf <- data.frame( step = seq_along( Rke), ke = Rke)
#'    kefit <- findlmlen( Rkedf, minp)
#'
findlmlen <- function( df, minp = 0.50) {
   kefit <- lmMD2( df)
   lendf <- nrow( df)
   if (kefit[ 'Fpvalue'] < minp) { # The data appears to be correlated
      third <- floor( lendf / 3) # Errs on the side of a longer 2/3s.
      ttdf <- tail( df, -third) # two thirds df
      kefit2 <- lmMD2( ttdf)
      if (kefit2[ 'Fpvalue'] >= minp) { # The 2nd 2/3s passed
       # Try to expand fit
         ftdf <- head( df, third) #first third df
         approach <- lmMD2( ftdf) # fit the approach to equilibrium
         rftdfke <- rev( ftdf$ke)
         if (sign( approach[ 'Estimate']) > 0) { # ke was growing
            # we are looking for first point (rftdfke is in reverse)
            # that is below the min.
            adder <- which.max( rftdfke < min( ttdf$ke))
         } else {
            # looking for the first point above the max
            adder <- which.max( rftdfke > max( ttdf$ke))
         } # if (sign( approach...))
         newlen <- lendf - third + adder - 1 # exclude the point exceeding limit
         kefit3 <- lmMD2( tail( df, newlen))
         if (kefit3[ 'Fpvalue'] >= minp) return( kefit3) else
            return( kefit2)
         # the above `if` should return either way - no way to get here.
      } # if (kefit2....)
      # if we get here, then neither whole nor 2/3's passed.
   } # if kefit...
   # if we get here, then either the original fit passed the F-test,
   # or neither passed, so the original fit is returned.
   # The F p-value will show if this record is suspect.
   return( kefit)
}

meansdba <- function( x, rsblk = 250) {
   # utility function to create block averages and compute mean and sd.
   z <- rsmpavg( x, rsblk)
   c( mean=mean( z), sd=sd(z))
}

#' procMD2DF Processes the scalar files from Julia MD2 program runs
#'
#' Summarizes the runs by calculating the mean, standard deviation, min, max,
#' of each column in the data frame,
#' and the relative variance of the kinetic energy
#' (i.e., `rtv = ke.sd / ke.mean`), and the estimate of the slope of `ke`, the
#' standard error of that estimate, and R^2 for the fit (`Estimate`,
#' `Std.Error`, and `r.squared` respectively). Also computes the pressure from
#' the average potential and kinetic energies.
#' The complex order parameter
#' magnitude is calculated two ways: step by step (and then averaged), and the
#' magnitude of the means of the real and imaginary components (named
#' `psiabs.mean` and `psiavg`, respectively).
#' @param df A data frame representing the entire series of runs.  Usually the
#' result of `combineruns`
#' @param skip The number of steps to skip in the averages.  `1,000` is an
#' historic default, but `1` is used when running from the command line, as one
#' usually wishes to see the initial effects.
#' @param fftblock Passed to `sdf`, this is the size of the fourier transform
#' block.  The algorithm will be fastest if this is a highly composite number,
#' best if a power of four ($4^N$).
#' @param rsblk The integer block size for ke re-sampling.  If greater than one,
#' the ke data will be averaged in blocks of size `rsblk` and those block
#' averages used to detect correlations.  The block should be long enough to
#' span correlations in the original data, resulting in essentially uncorrelated
#' data samples.
#' @param minp The minimum F-test p-value which will be be accepted as no
#' correlation. The default of 0.5 indicates that half of truly random data
#' will be falsely rejected.
#'
#' @return A record suitable for a data frame which includes the total steps,
#' the steps used in the averages (i.e. `total.steps - skip`).
#' @export
#'
#' @examples
#' d <- procMD2DF(x, skip=1000)
#'
procMD2DF <- function( df, skip = 0, fftblock = 4096, rsblk = 250,
                       minp = 0.50) {
   steps <- nrow( df) # total number of steps in original record.
   ke <- if (skip == 0) df$ke else tail(df$ke, -skip)

   # Spectral Processing
   spec <- sdf( detrend( ke), window="hanning", blocksize=fftblock, overlap=0.5,
                normalize=TRUE, prints=TRUE)
   lensdf <- length( spec)
   len2 <- (lensdf %/% 2)
   indicies <- c( (len2 + 1):lensdf, 1:len2)
   tp <- sum( spec) # total power
   sdfsum <- cumsum( spec[ indicies]) / tp
   width <- which.max( sdfsum > 0.995) - which.max( sdfsum > 0.005)
   width <- 2 * width / fftblock # % of spectrum occupied.

   TE.init = NA # If the entire record includes the first step, then it has
   # the initialization energy.
   if (df[1, "step"] == 1) TE.init = df[1, "TE"]

   # First Regression
   Rke <- rsmpavg( ke, rsblk)
   Rkedf <- data.frame( step = seq_along( Rke), ke = Rke)
   kefit <- findlmlen( Rkedf, minp)

   # skip initial rows & drop step number
   df <- subset( df, (step > skip) & (step < (skip + rsblk * kefit[ 'nr'] + 1)),
                      select = -1)
   # Must drop steps column, as `colstats will be run on all remaining columns.
   df$psiabs <- Mod( complex( real = df$Rpsi, imaginary = df$Impsi))
   meansd <- colStats( df, meansdba, "rs", rsblk)
   #sds   <- colStats( df, function( x) sd( rsmpavg( x, rsblk)), "rs.sd")
   minmax <- colStats( df, function(x) setNames( range(x), c('min', 'max')),
                       "raw")
   #maxs  <- colStats( df, max, "max")
   cv <- unname( (meansd[ "ke.sd.rs"] / meansd[ "ke.mean.rs"])^2)
   c( total.steps = steps, steps.avg = nrow( df), TE.init = TE.init,
      meansd, minmax, kefit, rtv = cv,
      psiavg = Mod( complex( real = meansd[ "Rpsi.mean.rs"],
                             imaginary = meansd[ "Impsi.mean.rs"])),
      specwidth = width)
}

#' createdf Creates the observation data frame for a series of files.
#'
#' `filenames` is a list of files and must come from `directory`.
#' `filenames` assumed to be of the form "series.energy.sequence.scalars.csv",
#' where series, energy, and sequence can be any identifiers (but ".",
#' of course), and "scalars.csv" is fixed. If `df` is supplied, it must have
#' the same structure as currently returned by `procMD2DF`.
#'
#' @param directory A string with a directory, passed unaltered to `list.files`
#' @param filenames A vector of strings with filenames to be added.
#' @param skip Integer with the number of steps to skip before computing
#' statistics.  Passed to `procMD2DF`.
#' @param df Existing data frame to which these new observations are being
#' added.
#' @param clusters The number of clusters to use. This can be a bit time
#' consuming, so file series are processed in parallel.
#'
#' @return Data frame with the summarized data.
#' @export
#'
#' @examples
#' num_cores <- detectCores() - 1
#' files2 = Sys.glob( file.path( "S1.*.scalars.csv"))
#' MD2DF <- createdf(  ".", files, clusters=num_cores)
#'
createdf <- function( directory, filenames, skip = 1000, df=data.frame(),
                      clusters=0, fftblock = 4096, rsblk = 250,
                      minp = 0.50, verbose=FALSE) {
   filenames <- sort( filenames)
   count <- length( filenames)
   # runfiles[[i]] is array of fields
   runfiles <- strsplit( filenames,split="[.]")
   dfrows <- unlist( unique( lapply( runfiles,
                                     function(y) paste(y[1], y[2], sep="\\."))))
   # series and energy will define an init.
   if (verbose) cat( dfrows)
   rowcnt <- length( dfrows)
   #
   # init df column vectors
   #
   series   <- rep( "",  rowcnt)
   energy   <- rep( "",  rowcnt)
   tmp <- procMD2DF( combineruns( directory, dfrows[1])) # init the names
   arr <- matrix( nrow = rowcnt, ncol = length(tmp),
                  dimnames = list(NULL, names(tmp)))
   for (i in 1:rowcnt) {
      serener <- unlist( strsplit( dfrows[i], split="\\\\."))
      series[ i] <- serener[1]
      energy[ i] <- serener[2]
   }

   func <- function( row) { # function to execute for each row
      x <- combineruns( directory, row)
      procMD2DF( x, skip, fftblock = fftblock, rsblk = rsblk, minp = minp)
   } # function( row)...

   if (clusters == 0) {
      parlist <- lapply( dfrows, func)
   }
   else {
      require( parallel) # only load if needed
      cl <- parallel::makeCluster( clusters)
      parallel::clusterExport( cl,
                               c( "procMD2DF", "lmMD2", "colStats", "detrend",
                                  "Vpower","rsmpavg",'meansdba', "findlmlen",
                                  "sdf", "taper", "combineruns"))
      parlist = parallel::parLapply( cl=cl, dfrows, func)
      parallel::stopCluster( cl)
   } # end if (clusters...)
   newdf <- as.data.frame( t( simplify2array(( parlist))))
   newdf <- cbind( data.frame( series = series, energy = energy), newdf)
   newdf <- rbind( df, newdf) # Add newdf to df passed as argument
   attr( newdf, 'call') <- list(directory = directory, filenames = filenames,
                                skip = skip, fftblock = fftblock, rsblk = rsblk,
                                minp = minp)
   return( newdf)
}

#'  Checks if the output file exists, and reads it in, otherwise creates it.
#'
#'  See `createdf` for the rest of the parameters.
#'
#' @param output The filename for the new (or existing) MD2DF file.
#' @param directory The directory to be searched.
#' @param pattern Globbing pattern passed through `glob2rx` before passing to
#' `list.files`.
#' @param skip
#' @param df
#' @param clusters
#' @param fftblock
#' @param rsblk
#' @param minp
#'
#' @return the MD2DF file created or read in.
#' @export
#'
#' @examples
MD2DFfile <- function( output, directory, pattern, skip = 1000, df=data.frame(),
                       clusters = detectCores() - 1, fftblock = 4096,
                       rsblk = 500, minp = 0.50, verbose=FALSE) {
   stopifnot( endsWith( output, ".RDS"))
   if (file.exists( output)) MD2DF <- readRDS( output) else {
      files = list.files( directory, glob2rx( pattern))
      cat( length( files), " files found.\n")
      MD2DF <- createdf( directory, files, skip = skip, df = df,
                             clusters = num_cores, fftblock = fftblock,
                             rsblk = rsblk, minp = minp, verbose=verbose)
      saveRDS( MD2DF, file = output)
   }
   return( MD2DF)
}

#' plotconf Plots a line or points with a shaded confidence band.
#'
#' plots the data with confidence bands. If `Fpvalue` is included, the plot
#' will be colored by points exceeding `minp`.
#' Note that if the data is not sorted
#' by x values, line plots can be confusing.
#'
#' @param x,y,sd,Fpvalue vectors of x, y, the standard
#' deviation of y points (this may be omitted, and no conf interval will be
#' plotted), and the F-test p-value.
#' @param linethres if the number of all points (regardless of `minp` coloring)
#' exceeds this, a line will be drawn
#' rather than points for the y values.
#' @param minp minimum F-test p-value for considering the point to be valid.
#' @param vpch the plot symbol for the y points considered valid.
#' @param ipch the plot symbol for the y points not considered valid.
#' @param vcol plot color for points considered valid.
#' @param icol plot color for points not considered valid.
#' @param ... passed to the initial `plot` call.
#'
#' @return
#' @export
#'
#' @examples
#'  plotconf( MD2DF, "TE.init", "ke.mean","ke.sd")
#'
plotconf <- function( x, y, sd, Fpvalue = rep( 1, length.out = length( x)),
                      linethresh = 100, minp = 0.5, vpch = 1, ipch = vpch + 1,
                      vcol = 'black', icol = 'darkred',
                      xlab="", ylab="", ...) {
   plci <- !(missing( sd) || is.null(sd))
   plot( x, y, xlab=xlab, ylab=ylab, type='n', ...) #set boundary on all points.
   # polygon has to be plotted first or it will cover the points and lines.
   if (plci) {
      so <- sort.list( x)
      x <- x[ so]
      y <- y[ so]
      sd <- sd[ so]
      Fpvalue <- Fpvalue[ so]
      polygon( c( rev(x), x), c( rev( y - sd), y + sd), col = 'grey',
               border = NA)
   }
   vp <- (Fpvalue >= minp) # Vector of valid point indices.
   if (length( x) < linethresh) { # Threshold set on all points.
      points( x[ vp], y[ vp], pch=vpch, col=vcol, ...)
      points( x[ !vp], y[ !vp], pch=ipch, col=icol, ...)
   }
   else {
      lines(  x[ vp], y[ vp], col=vcol, ...)
      lines(  x[ !vp], y[ !vp], col=icol, ...)
   }
   if (plci) {
      lines( x, y + sd, col ='darkblue')
      lines( x, y - sd, col ='darkblue')
   }
   grid()
}

#' plotconfMD2 Calls plotconf with names specific to an MD2 data frame.
#'
#' plotconfMD2 will sort the data according to the `x` variable before calling
#' plotconf.
#'
#' @param df The data frame
#' @param x The x variable.  If omitted, will default to the resampled mean
#' total energy "TE.mean.rs".
#' @param y The main y variable name as a string.  Must be one of the variables
#' that are statistical summaries and thus have a "y.mean" and a "y.sd"
#' column in `df`.
#' @param minp minimum F-test p-value for considering the point to be valid.
#' @param plci If `TRUE`, plot the confidence interval.
#' @param xlab,ylab passed to `plotconf`.
#' @param ... Passed to `plotconf`.
#'
#' @return nothing
#' @export
#'
#' @examples
#' plotconfMD2( MD2DF, "TE.mean", "ke")
#' plotconfMD2( MD2DF, y="ke")
#'
plotconfMD2 <- function( df, x = "TE.mean.rs", y, minp = 0.5, plci = TRUE,
                         linethresh = 100, xlab = x, ylab = y, ...) {
   ord = order( df[, x])
   df <- df[ord, ]
   xdata = df[, x]
   ydata = df[, paste( y, "mean.rs", sep=".")]
   sddat = if (plci) df[, paste( y, "sd.rs", sep=".")] else NULL
   plotconf( x=xdata, y=ydata, sd=sddat, Fpvalue = df$Fpvalue,
             linethresh = linethresh, xlab = xlab, ylab = ylab,
             minp = minp, ...)
}

#' Plot a predetermined collection of Scalars
#'
#' @param MD2DF A dataframe returned from `createdf`
#' @param minp A minimum p-value from the F-test to declare a run invalid
#' @param ... extra parameters passed to `plotconfMD2`
#'
#' @return Returns nothing.
#' @export
#'
#' @examples
#' MD2DF <- MD2DF[ order( MD2DF$TE.mean.rs),] # put df in increasing order
#' op=  par( mfrow = c(4, 2), mar = c(3, 3, 0.5, 0.5), tcl = -0.3,
#'           mgp = c(1.7, 0.4, 0) )
#' plotScalars( subset( MD2DF, TE.mean.rs< 0.2))
#' par(op)
#'
plotScalars <- function( MD2DF, minp = 0.5, ...){
   valid <- MD2DF$Fpvalue > minp
   pch <- rep_len( 1, length.out = nrow( MD2DF))
   pch[ !valid] <- 2
   col <- rep_len( 'black', length.out = nrow( MD2DF))
   col[ !valid] <- 'darkred' # change all non-valid points to red.
   plot( ke.sd.rs ~ TE.mean.rs, data = MD2DF, pch = pch, col = col, ...)
   grid()
   plotconfMD2( y="ke", df = MD2DF, minp=minp, ...)
   plotconfMD2( y="Pressure", df = MD2DF, minp=minp, ...)
   plot( rtv ~ TE.mean.rs, data = MD2DF, pch = pch, col = col, ...)
   grid()
   plot( r.squared ~ TE.mean.rs, data = MD2DF, pch = pch, col = col, ...)
   grid()
   MD2DF$ke.rrange <- (MD2DF$ke.max.raw - MD2DF$ke.min.raw) / MD2DF$ke.mean.rs
   plot( ke.rrange ~ TE.mean.rs, data = MD2DF, pch = pch, col = col, ...)
   plot( psiavg ~ TE.mean.rs, data = MD2DF, pch = pch, col = col, ...)
   grid()
   plot( Std.Error ~ abs(Estimate), data = MD2DF, pch = pch, col = col, ...)
   abline( a = 0, b = 1)
   grid()
}

"%within%" <- function( vector, range)
   (vector > range[ 1]) & (vector < range[ 2])
# Returns a logical true whenever the value of the vector is within the
# range specified by range[1] < vector < range[2].
# Use as infix operator: x[ y %within% C(0,l)] will return x's whose
# y's are between 0 and l.
# stopifnot( length(range)==2)
#' combineruns Combines files matching pattern into a single data frame.
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
#' x = combineruns(  ".", "S1\\.0014\\.") # will read files "S1.0014.*.csv"
#'
combineruns <- function( directory, pattern, verbose=FALSE){
   patt <-  paste( pattern,
                   "(.*).scalars.csv", sep=".")
   if (verbose) cat( "Searching Directory: ", directory,
                     " for pattern:", patt, "\n")
   files <- list.files( directory, pattern = patt)
   # Sys.glob( file.path( directory,
   # paste( pattern, "(.*).scalars.csv", sep=".")))
   files <- sort( files)
   if (verbose) cat( "Files found: ", files, "\n")
   stepoffset <- 0
   sdf <- data.frame() # init to empty dataframe
   for (j in 1:length( files)) { # walk through all files in group
      data <- read.csv( paste( directory, files[ j], sep="\\"))
      steps <- nrow( data)
      # cat( steps, "\n")
      data[ "step"] <- data[ "step"] + stepoffset
      stepoffset <- stepoffset + steps
   sdf <- rbind( sdf, data)
   }
   sdf$psiabs <- sqrt( sdf$Rpsi^2 + sdf$Impsi^2)
   invisible( sdf)
}

#' plotMD2 Plots a run of data step-by-step
#'
#' @param data A data frame, usually from `combineruns`.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' x = combineruns(  ".", "S1\\.0014\\.", verbose=TRUE)
#' plotMD2( x)
plotMD2 <- function( data) {
   op <- par( mfrow = c(3, 2), mar = c(3, 3, 0.5, 0.5), tcl = -0.3,
        mgp = c(1.7, 0.4, 0) )
   steps = data[,"step"]
   nc = ncol(data)
   #if (nc > 8) paste( "plotMD2 will only plot the columns 4 through 8")
   pl = names( data[ 4:nc])
   for (v in pl[ c( 1, 4, 2, 3, 5, 8)]) {
      plot( steps, data[, v], type="l", ylab=v)
      grid()
   }
   par( op)
   return( NULL)
}

readplot <- function( filename) {
   data = read.csv( filename)
   plotMD2( data)
   invisible( data)
}

colStats <- function( x, f, tag) {
  t <- unlist( lapply( x, f))
  names(t) <- paste( names( t), tag, sep=".")
  t
}

#' lmMD2 Computes the regression of the ke vs step/
#'
#' @param df a data frame
#'
#' @return
#' @export
#'
#' @examples
lmMD2 <- function( df) {
   mod <- summary( lm( ke ~ step, df))
   c( coeff <- mod$coefficients["step", 1:2], r.squared = mod$r.squared)
}

#' procMD2DF Processes the scalar files from Julia MD2 program runs
#'
#' Summarizes the runs by calculating the mean, standard deviation, min, max,
#' and the relative variance of the kinetic energy
#' (i.e., `rtv = ke.sd / ke.mean`), and the estimate of the slope of `ke`, the
#' standard error of that estimate, and R^2 for the fit (`Estimate`,
#' `Std.Error`, and `r.squared` respectively).
#'
#' @param df A data frame representing the entire series of runs.  Usually the
#' result of `combineruns`
#' @param skip The number of steps to skip in the averages.  `1,000` is an
#' historic default, but `1` is used as whn running form the command line, one
#' usually wishes to see the initial effects.
#'
#' @return A record suitable for a data frame which includes the total steps,
#' the steps used in the averages (i.e. `total.steps - skip`),
#' @export
#'
#' @examples
#' d <- procMD2DF(x, skip=1000)
#'
procMD2DF <- function( df, skip=1) {
   steps <- nrow(df) # total number of steps in record.
   kefit <- lmMD2( df[skip:steps,]) # must do this first as lm is against step
   TE.init = NA
   if (df[1, "step"] == 1) TE.init = df[1, "TE"]
   # skip initial rows & drop step number
   df <- subset( df, step > skip, xMomentum:Pressure)
   means <- colStats( df, mean, "mean")# vector
   sds <- colStats( df, sd, "sd")
   mins <- colStats( df, min, "min")
   maxs <- colStats( df, max, "max")
   cv <- unname( (sds["ke.sd"] / means[ "ke.mean"])^2)
   c( total.steps = steps, steps.avg = (steps - skip + 1), TE.init = TE.init,
      means, sds, mins, maxs, kefit, rtv = cv)
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
#' statistics.
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
                      clusters=0) {
   filenames = sort( filenames)
   count <- length( filenames)
   runfiles <- strsplit( filenames,split="[.]") # runfiles[[i]] is array of fields
   dfrows <- unlist( unique( lapply( runfiles,
                                     function(y) paste(y[1], y[2], sep="\\."))))
   # series and energy will define an init.
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
   func <- function( row) {
      x <- combineruns( directory, row)
      procMD2DF( x, skip)
   } # function( row)...
   if (clusters == 0) {
      parlist <- lapply( dfrows, func)
   }
   else {
      require( parallel) # only load if needed
      cl <- parallel::makeCluster( clusters)
      # registerDoParallel( cl)
      parallel::clusterExport( cl, c( "procMD2DF", "lmMD2", "colStats",
                               "combineruns"))
      parlist = parallel::parLapply( cl=cl, dfrows, func)
      parallel::stopCluster( cl)
   }
   newdf <- as.data.frame( t( simplify2array(( parlist))))
   newdf <- cbind( data.frame( series = series, energy = energy), newdf)
   rbind( df, newdf)
}

#' plotconf Plots a line or points with a shaded confidence band.
#'
#' @param df A data frame with the variables.
#' @param x,y,sd Strings naming the x and y data columns and the standard
#' deviation data column.
#' @param linethres if the number of points exceeds this, a line will be drawn
#' rather than points for the y values.
#' @param pch The plot symbol for the y points.
#' @param ... Passed to the initial `plot` call.
#'
#' @return
#' @export
#'
#' @examples
#'  plotconf( MD2DF, "TE.init", "ke.mean","ke.sd")
#'
plotconf <- function( x, y, sd, linethresh=100, pch=19, xlab="", ylab="",
                      ...) {
   plot( x, y, xlab=xlab, ylab=ylab, type='n', ...)
   grid()
   polygon( c( rev(x), x), c( rev( y - sd), y + sd),
            col = 'grey', border = NA)
   if (length( x) < 100) {
      points( x, y, pch=pch,  ...)
   }
   else {
      lines( x, y, ...)
   }
   lines( x, y + sd, col ='blue')
   lines( x, y - sd, col ='blue')
}

#' plotconfMD2 Calls plotconf with names specific to an MD2 data frame.
#'
#' plotconfMD2 will sort the data accordding to the `x` variable before calling
#' plotconf.
#'
#' @param df The data frame
#' @param x The x variable.  If omitted, will default to the mean total energy
#'  "TE.mean"
#' @param y The main y variable name.  Must be one of the variables that are
#' statistical summaries and thus have a "y.mean" and a "y.sd" column in `df`.
#' @param ... Passed to `plotconf`.
#'
#' @return nothing
#' @export
#'
#' @examples
#' plotconfMD2( MD2DF, "TE.mean", "ke")
#' plotconfMD2( MD2DF, y="ke")
#'
plotconfMD2 <- function( df, x = "TE.mean", y, ...) {
   ord = order( df[, x])
   xdata = df[ ord, x]
   ydata = df[ ord, paste( y, "mean", sep=".")]
   sddat = df[ ord, paste( y, "sd", sep=".")]
   plotconf( x=xdata, y=ydata, sd=sddat, xlab=x, ylab=y)
}

plotScalars <- function( MD2DF, ...){
   plot( ke.sd ~ TE.mean, data = MD2DF, ...)
   grid()
   plotconfMD2( y="ke", df = MD2DF, ...)
   plotconfMD2( y="Pressure", df = MD2DF, ...)
   plot( rtv ~ TE.mean, data = MD2DF, ...)
   grid()
   plot( r.squared ~ TE.mean, data = MD2DF, ...)
   grid()
   plotconfMD2( y="TE", x="TE.init", df = MD2DF, ...)
}

"%within%" <- function( vector, range)
   (vector > range[ 1]) & (vector < range[ 2])
# Returns a logical true whenever the value of the vector is within the
# range specified by range[1] < vector < range[2].
# Use as infix operator: x[ y %within% C(0,l)] will return x's whose
# y's are between 0 and l.
# stopifnot( length(range)==2)
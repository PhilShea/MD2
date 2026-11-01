dydx <- function( y, x, lag = 1) {
   # find derivative over lag + 1 points by fitting polynomial
   # and using the derivative of the polynomial at the center of
   # points.
   lx <- length( x)
   stopifnot( length( y) == lx)
   data <- data.table( x = x, y = y, keyby = x)
   xm <- (data$x[ 1:(lx - lag)] + data$x[ (lag + 1):lx]) / 2
   ym <- (data$y[ 1:(lx - lag)] + data$y[ (lag + 1):lx]) / 2
   # midpoint of eval sets
   dydxfit <- function( i) { # get dy/dx from points i:(i + lag)
      c <- coef( lm( y ~ stats::poly( x, degree = lag, raw = TRUE), # get poly coef
                     data = data[ i:(i + lag)]))
      dpc <- (1:lag) * tail( c, -1) # derivative of poly coeff
      polyval( rev(dpc), xm[ i])
      # evaluate polynomial at midpoint of data & return
   }
   dd <- vapply( 1:(lx - lag), dydxfit, 1.0)
   list( dd=dd, xm = xm, ym = ym)
}

# test
tx <- c(1, 4, 3, 3.5, 2, 5)
ty <- tx^2
(dd <- dydx(ty, tx, 2))
dd$dd - 2 * dd$xm
rm( tx, ty, dd) # remove these variables

stepdet <- function( x, nsig = 3, nstart = 10, verbose = FALSE) {
   # x is a numeric vector
   # nsig is the number of standard deviations from the running mean x[i] must
   # depart to be considered a step.
   # nstart is the number of points to start the running mean.
   # x must in in relevant order.
   # Finds the index of the point such that
   #     abs( x[nv] - mean( x[1:n])) > nsig * sd(x[1:n])
   #
   # returns a list:
   #     (sx: cumsum(x[1:nv], ss: sum squares, n: points in running mean,
   #      nv: detected point number)
   #
   stopifnot( length( x) > nstart +1)
   if (verbose) print( x)

   cumdet <- function( state, new){
      # state is a list with the following elements:
      #   sx: cumsum of x
      #   ss: the sum of squares
      #   n: the number of points in the running mean
      #   nv: the count when first detected
      # new is a new value to add to the running mean.
      if (state$nv == 0 ) { # only do this if we haven't detected a departure
         meanx <- state$sx / state$n
         sdx <- sqrt( (state$ss - state$n * meanx^2)/ (state$n - 1))
         if (verbose) {
            #print(state)
            print( c( meanx, sdx))
         }
         if( abs( new - meanx) < (nsig * sdx)) { # nothing detected
            state$sx <- state$sx + new
            state$ss <- state$ss + new^2
            state$n <- state$n + 1
            state$nv <- 0
         } else { # departure detected
            state$nv <- state$n + 1
         }
      } # simply return state.
      return( state)
   }
   sxs <- sum( x[ 1:nstart])
   sss <- sum( x[ 1:nstart]^2) #sum squared start.
   s <- list( sx = sxs, ss = sss, n = nstart, nv = 0)
   if (verbose) print( s)
   Reduce( cumdet, x[(nstart + 1):length( x)], init = s)
}

# Test
(testcumdet <- stepdet( tx <- c( rnorm(20), 5 + rnorm(10)))) # n=20, nv=21

linedet <- function( x, y, level = 0.99, nstart = 10) {
   # x & y are a numeric vectors of the same length
   # nsig is the probability that a noise deviation will cause a false detection
   #   over the total run of points.
   # depart to be considered a step.
   # nstart is the number of points to start the running mean.
   #
   # returns a named vector nv is the point number, x is x[nv], y is y[nv],
   # pred is the predicted point, and lwr & upr are the level deviations on
   # pred.
   lx <- length( x)
   stopifnot( lx > nstart +1)
   stopifnot( lx == length( y))
   dt <- data.table( x = x, y = y)
   for ( i in nstart:(lx - 1)) { # step through looking for a departure
      pr <- predict( lm( y ~ x, data = dt[ 1:i,]),
                     newdata = data.frame( x = x[ i + 1]),
                     interval = "prediction", level = level)
      # cat( i, pr, "\n")
      if (!(y[ i + 1] %within% pr[, c("lwr","upr")])) { # if departure detected
         return( c( nv = i + 1, x = x[ i + 1], y = y[ i + 1],
                    pred = pr[1, 1], lwr = pr[1, 2], upr = pr[1, 3]))
      } # if (!(y...)
   } # for ( i...)
   # if we get here, there was no detection (if above was never TRUE)
   return( c( nv = NA, x = NA, y = NA,
              pred = pr[1, 1], lwr = pr[1, 2], upr = pr[1, 3]))
}

# Test
tx <- 1:30
ty <- tx + rnorm(30)
ty[ 21:30] <- ty[ 21:30] + 5
(testcumdet <- linedet( x = tx, y = ty)) # nv=21, x=20

rm( tx, ty, testcumdet) # remove test variables.

detect <- function( data, nsig = 3, nstart = 10, luc, level = 0.99) {
   # data is an MD2DF data.table melted like
   #    id.var = c( "run", "Series", "Temp", "valid")
   #  # or
   #    id.var = c( "run", "Series", "TE", "valid")
   #    SDFd2m3 <- melt( SDFd2ss, id.var = id.var, direction = 'long',
   #                     variable.name = 'variable', variable.factor = TRUE,
   #                     value.name = 'value')
   # nsig is the number of standard deviations for stepdet to detect start
   #     of diffusion
   # nstart is the number of points to include in the start of running means
   #     or regression lines.
   # luc is the look-up column, as a string; luc = "Temp" or "TE"
   # level is the significance threshold to exceed to call a line departure
   #     a detection
   #
   #  Returns a data.table of temperatures (or whatever is passed in `luc`)
   #  at which detections occurred with a string describing each detection.
   #
   setkeyv( data, luc) # sorts according to luc
   nd <- 8L # number of rows to initialize in result dt
   row <- 1L
   dt <- data.table( variable = rep( "", nd), value = rep( as.numeric(NA), nd),
                     ex = rep( "", nd), delx = rep( as.numeric(NA), nd))
   setdt <- function( i, variable="", value, ex = "", delx = 0.0) {
      # this closure captures dt and facilitates setting row i
      set( dt, i, "variable", variable)
      set( dt, i, "value", value)
      set( dt, i, "ex", ex)
      set( dt, i, "delx", delx)
      i + 1L
   }
   subs <- data[ variable == "D" & !is.na( value),]
   # print( nrow( subs))
   # not all records have diffusion estimates.
   if (nrow( subs) > nstart) {
      sd <- stepdet( subs$value, nsig, nstart)
      nv <- sd$nv
      value <- subs[[ luc]][ nv]
      delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      row <- setdt( row, "diff", value, "Start of diffusion (stepdet)", delx)
   }

   subs <- data[ variable == "Cv",]
   if (nrow( subs) > 1) {
      nv <- which.max( subs$value)
      value <- subs[[ luc]][ nv]
      delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      row <- setdt( row, "maxCv", value, "Max Cv via rv", delx)
   }

   subs <- data[ variable == "Log(Cv)",]
   if (nrow( subs) > 1) {
      nv <- which.max( subs$value)
      value <- subs[[ luc]][ nv]
      delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      row <- setdt( row, "maxLCv", value, "Max Cv via 1st diff", delx)
   }

   subs <- data[ variable == "rv",]
   if (nrow( subs) > 1) {
      nv <- which.max( subs$value)
      value <- subs[[ luc]][ nv]
      delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      row <- setdt( row, "Maxrv", value,
                                 "Max relative Temp variance", delx)
   }

   subs <- data[ variable == "LogFpD" & !is.na( value),]
   if (nrow( subs) > 3) {# need min of 3 items
      nv <- which( subs$value < log10(0.01))[1] # first index with lower probability
      value <- subs[[ luc]][ nv]
      if (nv %within% c( 2, length(subs) - 1)) {
         delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      }
      else{
         delx <- subs[[ luc]][ nv]
      }
      row <- setdt( row, "LogFpD", value, "Start of Diff (Fp < 0.01)", delx)
   }

   subs <- data[ variable == "Pressure",]
   if (nrow( subs) > nstart) { # if not long enough, skip all linedet() code
      lv <- linedet( x = subs[[ luc]], y = subs$value, nstart = nstart,
                     level = level)
      nv <- lv[ "nv"]
      value <- lv[ "x"]
      delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
      row <- setdt( row, "Icld1", value, "Isochor Solid1 (linedet)", delx)

      if (nrow( subs) > nv + 1) { # here, we start after Icld1 detection
         lv <- linedet( x = subs[[ luc]], y = subs$value, nstart = nv + 1,
                        level = .99999)
         nv <- lv[ "nv"]
         value <- lv[ "x"]
         delx <- (subs[[ luc]][ nv + 1] - subs[[ luc]][ nv - 1])/4 # 1/2 avg interval
         row <- setdt( row, "Icld2", value, "Isochor Solid2 (linedet 5 nines)",
                       delx)
      } # second linedet() call

      setorderv( subs, luc, -1) # Reverse order
      lv <- linedet( x = subs[[ luc]], y = subs$value, nstart = nstart,
                     level = level)
      nv <- lv[ "nv"]
      value <- lv[ "x"]
      # note that this is reversed, as the sort is reversed.
      delx <- (subs[[ luc]][ nv - 1] - subs[[ luc]][ nv + 1])/4 # 1/2 avg interval
      row <- setdt( row, "Icld3", value, "Isochor Fluid (linedet)", delx)
   } # all linedet calls

   return( dt)
}

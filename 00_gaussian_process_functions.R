norm.solve <- function( y, params )
{
  qnorm( params$q, y[1], y[2] ) - params$x
}

#' Back out Guassian mean/sd from pairs of quantiles
#'
#' @x the values of the quantiles
#' @q the percentiles of the quantiles

CI_to_guassian <- function( x, q )
{
  if( length( x ) != 2 || length( q) != 2 )
    stop( "x's and q's must come in matched pairs")
  
  meaninit <- mean( x )
  sdinit   <- abs( max(x)-min(x))
  init     <- c( meaninit, sdinit)
  
  parms <- list( q=q, x=x)
  mm <- multiroot(  norm.solve, init, parms = parms )
  
  list( 
    mean = mm$root[1],
    sd   = mm$root[2]
  )
}




get_weights <- function( gcomposite )
{
  if( is.null( gcomposite$weights ))
    return( rep(1,length(gcomposite$mean)))
  
  gcomposite$weights
}


##
##  Gaussian mixtures (aka "mixture") -----
##



#' Probability density of Gaussian mixture 
#'
#' @x vector of x's 
#' @gmixture a list with the relevant gmixture parameters

dmixture <- function( x, gmixture )
{
  weights = get_weights( gmixture )
  map_dbl( x, ~ weighted.mean( dnorm( ., gmixture$mean, gmixture$sd ), weights))
}

#' Distribution function of Gaussian mixture 
#'
#' @x vector of x's 
#' @gmixture a list with the relevant gmixture parameters

pmixture <- function( x, gmixture )
{
  weights = get_weights( gmixture )
  map_dbl( x, ~ weighted.mean( pnorm( ., gmixture$mean, gmixture$sd ), weights))
}



#' Helper function for backing out quantiles

x.solve <-  function( x, params )
{
  pmixture( x, params$gmixture )-params$p
}

# ' Quantile of a Guassian mixture for a single q level

qmixture1 <- function( p1, gmixture )
{  
  x.comp <- qnorm( p1, gmixture$mean, gmixture$sd )
  eps = max( gmixture$sd )*1e-6
  params = list( p = p1, gmixture = gmixture )
  uu <- uniroot( x.solve, c(min(x.comp)-eps, max(x.comp)+eps), params )
  uu$root
}

#' Quantile function of Gaussian mixture 
#'
#' @q vector of quantiles
#' @gmixture a list with the relevant gmixture parameters

qmixture <- function( p, gmixture )
{  
  map_dbl( p, ~qmixture1( ., gmixture ))
}


#' Random samples of Gaussian mixture 
#'
#' @n number of samples
#' @gmixture a list with the relevant gmixture parameters

rmixture <- function( n, gmixture)
{
  weights <- get_weights( gmixture )
  cump <- cumsum( weights[-length(weights)] )/sum(weights)
  idx  <- map_int( runif( n ), ~ sum( . > cump )) +1
  rnorm( n, gmixture$mean[idx], gmixture$sd[idx])
}



##
##  Gaussian quantile interpolations (aka "interp") -----
##


#' Quantile function of Gaussian mixture 
#'
#' @q vector of quantiles
#' @gqinterp a list with the relevant gmixture parameters


qinterp <- function( p, gqinterp )
{
  weights = get_weights( gqinterp )
  map_dbl( p, ~ weighted.mean( qnorm( ., gqinterp$mean, gqinterp$sd ), weights))
}


#' Helper function for backing out probabilities in the CDF

p.solve <-  function( p, params )
{
  if( p>1 || p<0)
    return( NaN)
  qinterp( p, params$gqinterp )-params$x1
}

# ' Distribution of a Guassian q-interpolation for a single x

pinterp1 <- function( x1, gqinterp )
{  
  # avoid problems if all distributions are the same
  p.comp <- pnorm( x1, gqinterp$mean, gqinterp$sd )
  
  params = list( x1 = x1, gqinterp = gqinterp )
  uu <- uniroot( p.solve, c(max(0,min(p.comp)-1e-6),min(1,max(p.comp)+1e-6)), params )
  uu$root
}

#' Distribution function of Gaussian quantile interpolations 
#'
#' @x vector of x's 
#' @gqinterp a list with the relevant gqinterp parameters
#' 
#' 

pinterp <- function( x, gqinterp )
{
  map_dbl( x, ~ pinterp1( ., gqinterp ))
}

# ' Density of a Guassian q-interpolation for a single x

dinterp1 <- function( x1, gqinterp )
{
  
  # the derivative of the distribution is 1/the derivative of the quantile function
  # the derivative of the quantile function is 1/the derivative of the distribution
  # so this is the harmonic mean of the component pdfs
  
  q <-pinterp( x1, gqinterp )
  xs <- qnorm( q,  gqinterp$mean, gqinterp$sd)
  pp <- dnorm( xs, gqinterp$mean, gqinterp$sd )
  1/weighted.mean( 1/pp, gqinterp$weights )
}


#' Probability density function of Gaussian quantile interpolations 
#'
#' @x vector of x's 
#' @gqinterp a list with the relevant gqinterp parameters
#' 
#' 

dinterp <- function( x, gqinterp )
{
  if( is.null( gqinterp$weights ))
    gqinterp$weights = rep(1,length(gqinterp$mean))
  
  map_dbl( x, ~ dinterp1( ., gqinterp ))
}


#' Random samples for Gaussian quantile interpolations 
#'
#' @n number of samples
#' @gqinterp a list with the relevant gqinterp parameters
#' 
#' 

rinterp <- function( n, gqinterp)
{
  weights <- get_weights( gqinterp )
  samples <- runif(n)
  qinterp( samples, gqinterp )
}
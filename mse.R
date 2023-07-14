#' Spline using Genetics Algorithm
#' @param x independent variables matrix
#' @param y dependent variable matrix
#' @return spline function using GA
#' @export

##########################################################
#################### MSE #################################
##########################################################
mse = function (y,x,knot,p)
{
  nknot = nrow(knot)
  m = ncol(x)
  n = nrow(x)
  alpa = p + 1 + nknot
  ncolx = alpa * m
  X = matrix (0,nrow = n, ncol = alpa)

  for (k in 1:m)
  {
    X1 = matrix(0, nrow = n, ncol = alpa)
    for (i in 1:alpa)
    {
      for (j in 1:n)
      {
        if (i <= p+1)
        {
          X1[j,i] = x[j,k]^(i-1)
        }
        else
        {
          X1[j,i] = truncate (x[j,k],knot[(i-(p+1)),k],p)
        }
      }
    }
    if (k == 1)
    {
      X = X1
    }
    else
    {
      X=cbind(X,X1)
    }
  }

  I = diag(nrow(x))
  tetha =pinv((t(X)%*%X))%*%t(X)%*%(y)
  H = X%*%pinv((t(X)%*%X))%*%t(X)
  Ytopi = X%*%tetha
  M = y-Ytopi
  HA = 1/n*(euclidnorm(M))
  mse = HA

  return(mse)
}

#' @export
error = function (y,x,knot,p)
{
  nknot = nrow(knot)
  m = ncol(x)
  n = nrow(x)
  alpa = p + 1 + nknot
  ncolx = alpa * m
  X = matrix (0,nrow = n, ncol = alpa)

  for (k in 1:m)
  {
    X1 = matrix(0, nrow = n, ncol = alpa)
    for (i in 1:alpa)
    {
      for (j in 1:n)
      {
        if (i <= p+1)
        {
          X1[j,i] = x[j,k]^(i-1)
        }
        else
        {
          X1[j,i] = truncate (x[j,k],knot[(i-(p+1)),k],p)
        }
      }
    }
    if (k == 1)
    {
      X = X1
    }
    else
    {
      X=cbind(X,X1)
    }
  }

  I = diag(nrow(x))
  tetha =pinv((t(X)%*%X))%*%t(X)%*%(y)
  H = X%*%pinv((t(X)%*%X))%*%t(X)
  Ytopi = X%*%tetha
  M = y-Ytopi
  HA = 1/n*(euclidnorm(M))
  error = M

  return(error)
}

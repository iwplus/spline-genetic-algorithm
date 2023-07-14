#' Splines with genetic algorithms optimization method
#' Spline using Genetics Algorithm
#' @param x independent variables matrix
#' @param y dependent variable matrix
#' @return spline function using GA
#' @export

##########################################################
################## KNOT POINTS ############################
##########################################################

knot = function(x,nknot)
{
	knot = matrix(0, nrow=nknot,ncol=ncol(x))
	for (i in 1:nknot)
	{
		for (j in 1:ncol(x))
		{
			knot[i,j] = sample(max(x[,j]):min(x[,j]),1)
		}
	}
 print(knot)
}


##########################################################
#################### Yhat #################################
##########################################################
#' @export
ytopi = function (y,x,knot,p)
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
  HB = ((1/n)*sum(diag(I-H)))^2
  gcv = HA/HB
  ytopi = Ytopi


  return(ytopi)
}


##########################################################
#################### X matrix #################################
##########################################################
#' @export
ram = function (y,x,knot,p)
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
  HB = ((1/n)*sum(diag(I-H)))^2
  gcv = HA/HB
  ram = X
  return(ram)
}


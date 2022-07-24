#' Splines with genetic algorithms optimization method
#'
#' Fungsi ini dibuat untuk mencari nilai GCV pada Spline Truncated. Dimana dalam
#' fungsi ini terdapat fungsi-fungsi yang lain yang dibutuhkan. data yang akan
#' dimasukkan di fungsi ini harus disimpan terdahulu di dalam data frame, dengan
#' format Y adalah variabel dependen, dan X adalah variabel indevenden.
#'
#' @param x matriks independen
#' @param y matriks dependen
#' @param knot dalam bentuk matriks
#' @param p orde yang diinginkan
#' @return fungsi Spline yang menggunakan metode Algoritma Genetika
#' @export

gcv1 = function (y,x,knot,p)
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
 gcv1 = HA/HB
 mse = HA


 print(gcv1)
 return(gcv1)

}


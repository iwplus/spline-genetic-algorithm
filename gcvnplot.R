#' Fungsi Spline menggunakan metode Algoritma genetika
#'
#' Fungsi ini dibuat untuk mencari nilai GCV pada Spline Truncated. Dimana dalam
#' fungsi ini terdapat fungsi-fungsi yang lain yang dibutuhkan. data yang akan
#' dimasukkan di fungsi ini harus disimpan terdahulu di dalam data frame, dengan
#' format Y adalah variabel dependen, dan X adalah variabel indevenden.
#'
#' @param x matriks independen
#' @param y matriks dependen
#' @return fungsi Spline yang menggunakan metode Algoritma Genetika
#' @export
#############################################################
################ HAL YANG DIPERLUKAN ########################
#############################################################

library(pracma)

truncate = function(t,k,d)
{
  if (t >= k)
   {
   	trun = (t-k)^d
   }
  else
   {
	trun = 0
   }
 return(trun)
}
#' @export
euclidnorm = function(b) #### Norma Euclid
{
 norma = 0
 for (i in 1:length(b))
 {
	norma = norma + b[i]^2
 }
 return(norma)
}
##########################################################
#################### GCV #################################
##########################################################
#' @export
gcv = function (y,x,knot,p)
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


 plot(y, type="o", col="blue", lty=1)
 lines(Ytopi, type="o",col="red",  lty=2)
 legend(1, max(y), legend=c("Data Aktual", "Data Estimasi"), col=c("blue", "red"), lty=1:2, cex=0.8, title="Keterangan", text.font=4, bg='white',box.lty=2, box.lwd=2, box.col="steelblue")
 return(gcv)
}



#' Splines with genetic algorithms optimization method
#'
#' Fungsi ini dibuat untuk mencari nilai GCV pada Spline Truncated. Dimana dalam
#' fungsi ini terdapat fungsi-fungsi yang lain yang dibutuhkan. data yang akan
#' dimasukkan di fungsi ini harus disimpan terdahulu di dalam data frame, dengan
#' format Y adalah variabel dependen, dan X adalah variabel indevenden.
#'
#' @param x matriks independen
#' @param y matriks dependen
#' @param nknot banyaknya titik knot
#' @param p orde yang diinginkan
#' @param a banyak populasi yang diinginkan
#' @param iterasi banyaknya pengulangan
#' @param permu nilai permutasi dalam proses mutasi
#' @return fungsi Spline yang menggunakan metode Algoritma Genetika
#' @export


genetika = function(y,x,a,nknot,iterasi,permu,p)
{

 	b = ncol(x)
 	popA1 = matrix(0,nrow=(nknot*a),ncol=b)
 	for (i in 1:(nknot*a))
 	{
		for (j in 1:b)
		{
			popA1[i,j] = sample(max(x[,j]):min(x[,j]),1)
		}
 	}
 	GCV = matrix (0,nrow=a,ncol=1)
	if (nknot == 1)
	{
		for (i in 1:a)
		{
			knot = matrix(popA1[i,],1,b)
			GCV[i,] = gcv1(y,x,knot,p)
		}
	}else
	{
 		for (i in 1:a)
 		{
			knot = popA1[(((i-1)*nknot)+1),]
			for (j in 1:(nknot-1))
			{
				knot = rbind(knot,popA1[(((i-1)*nknot)+j+1),])
			}
	  	   GCV[i,] = gcv1(y,x,knot,p)
 		}
	}

 	for (z in 1:iterasi)
 	{
		GCV = matrix (0,nrow=((length(popA1)/b)/nknot),ncol=1)
 		if (nknot == 1)
		{
			for (i in 1:((length(popA1)/b)/nknot))
			{
				knot = matrix(popA1[i,],1,b)
				GCV[i,] = gcv1(y,x,knot,p)
			}
		}else
		{
			for (i in 1:((length(popA1)/b)/nknot))
 			{
				knot = popA1[(((i-1)*nknot)+1),]
				for (j in 1:(nknot-1))
				{
					knot = rbind(knot,popA1[(((i-1)*nknot)+j+1),])
				}
	 	  	   GCV[i,] = gcv1(y,x,knot,p)
 			}
		}

 		F = matrix (0,nrow=((length(popA1)/b)/nknot),ncol=1)
 		for (i in 1:((length(popA1)/b)/nknot))
 		{
 			F[i,] = 1/(GCV[i,]+1)
 		}
		
		print(F)
 		Ftotal = sum(F)
		print(Ftotal)
 		F1 = matrix (Ftotal,nrow=((length(popA1)/b)/nknot),ncol=1)
 		pro = F/F1
		print(pro)

 		k = 0
 		kum = matrix (0,nrow=((length(popA1)/b)/nknot),ncol=1)
		for (i in 1:((length(popA1)/b)/nknot))
 		{
 			k = k + (pro[i,])
 			kum[i,] = k
		}
		print(kum)

		if (mod(((length(popA1)/b)/nknot),2) == 1)
		{
			ns = ((length(popA1)/b)/nknot)-1
		}else
		{
			ns = ((length(popA1)/b)/nknot)
		}
 		s = matrix((runif(ns,0,1)),nrow=ns,ncol=1)
		print(s)
 		pop = matrix (0,nrow=(nknot*ns),ncol=b)
 		for (i in 1:ns)
 		{
 			for (j in 1:((length(popA1)/b)/nknot))
 			{
 				if (s[i,]<=kum[j,])
 				{
 					pop[(((i-1)*nknot)+1):(i*nknot),] = popA1[(((j-1)*nknot)+1):(j*nknot),]
 					break()
 				}
 			}
 		}
		
		print(pop)
 		xa = sample(2:b,ns,replace=T)
		print(xa)
 		x1 = matrix (xa,nrow=ns,ncol=1)
		print(x1)
 		pop1 = matrix (0,nrow=nknot*ns,ncol=b)
		acros = ceil(ns/2)
		for (i in 1:acros)
 		{
			ind = acros+i 	# indeks pasangan crossover individu ke i
			for (j in 1:nknot)
			{
				pop1[(((i-1)*nknot)+j),1:x1[i,]] = pop[(((i-1)*nknot)+j),1:x1[i,]]
				pop1[(((ind-1)*nknot)+j),1:x1[i,]] = pop[(((ind-1)*nknot)+j),1:x1[i,]]

				pop1[(((i-1)*nknot)+j),x1[i,]:b] = pop[(((ind-1)*nknot)+j),x1[i,]:b]
				pop1[(((ind-1)*nknot)+j),x1[i,]:b] = pop[(((i-1)*nknot)+j),x1[i,]:b]
			}
 		}

 		xb = runif(ns,0,1)
 		x2 = matrix(xb,nrow=ns,ncol=1)
		print(x2)

 		for (i in 1:ns)
 		{
 			if (x2[i,]<=permu)
 			{
				nmut = sample(1:b,1)
				for (w in 1:nmut)
				{
 					indk = sample(1:b,1)
					indb = sample(1:nknot,1)
					pop1[(((i-1)*nknot)+indb),indk] = sample(max(x[,indk]):min(x[,indk]),1)

				}
 			}
 		}
		print(pop1)
 		popA1 = rbind(popA1,pop1)
 	}
	

 	GCV = matrix (0,nrow=ns,ncol=1)
 	if (nknot == 1)
	{
		for (i in 1:ns)
		{
			knot = matrix(popA1[i,],1,b)
			GCV[i,] = gcv1(y,x,knot,p)
		}
	}else
	{
		for (i in 1:ns)
 		{
			knot = popA1[(((i-1)*nknot)+1),]
			for (j in 1:(nknot-1))
			{
				knot = rbind(knot,popA1[(((i-1)*nknot)+j+1),])
			}
	 	  GCV[i,] = gcv1(y,x,knot,p)
 		}
	}

 	GCVmin = min(GCV)
 	knotA = matrix(0,nrow=nknot,ncol=ncol(x))
	if (nknot == 1)
	{
		for (i in 1:ns)
		{
			if (min(GCV) == GCV[i,])
 			{
				knot = matrix(popA1[i,],1,b)
				knotA = knot
			}
		}
	}else
	{
		for (i in 1:ns)
 		{
			if (min(GCV) == GCV[i,])
 			{
				knot = popA1[(((i-1)*nknot)+1),]
				for (j in 1:(nknot-1))
				{
					knot = rbind(knot,popA1[(((i-1)*nknot)+j+1),])
				}
			    knotA = knot
 			}
 		}
	print(knotA)
	}
 genetika = knotA
 return(genetika)
}





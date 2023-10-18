##################################################################################################################
# Autor: Gabriel Jaume Martin
#
# TFM Analisis de Datos Masivos en Economia y Empresa
# Extension crisp y difusa del problema del agente viajero aplicado al modelado de rutas de una empresa de reparto
#
# Ultima modificacion: 18/10/2023
##################################################################################################################


################################################################################
# LIBRERIAS
################################################################################
library(lpSolve)


################################################################################
# FUNCIONES
################################################################################
# Funcion extraccion coeficientes de matrices
no_diag = function(B){
  
  # Permite expresar los elementos de una matriz, distintos a la diagonal, como vector
  
  #Args:
  # B: matriz de la que obtener los coeficientes
  
  #Returns:
  # elementos: vector con los coeficientes deseados
  
  elementos = c()
  for(i in 1:dim(B)[1]){
    for(j in 1:dim(B)[2]){
      if(i!=j){
        elementos = c(elementos, B[i,j])
      }
    }
  }
  
  return(elementos)
}

# Restriccion (2) de la memoria
rest_2 = function(A, count = 1, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (2).
  # \forall j \in P : \sum_{c=1}^{C}\sum_{q=1}^Q\sum_{k=q}^K\sum_{\begin{subarray}{l} i \in L\\ i\neq j\end{subarray}}prec_{ijkqc} = 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  for(j in 2:n){
    prec=c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[,j] = 1
          prec= c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
    count = count +1
  }
  
  return(list(A, count))
}

# Restriccion (3) de la memoria
rest_3 = function(A, count, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (3).
  # \forall i \in P: \sum_{c=1}^{C}\sum_{q=1}^Q\sum_{k=q}^K\sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}} prec_{ijkqc} = 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  for(i in 2:n){
    prec=c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[i,] = 1
          prec= c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0,(n-1)))
    count = count +1
  }
  
  return(list(A, count))
}

# Restriccion (4) de la memoria
rest_4 = function(A, count, dif, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (4).
  # \forall q, \forall c: \sum_{k=q}^K\sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}}prec_{1jkqc} \leq 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  for(cam in 1:Cam){
    ind_mat=0
    for(s in 1:Q){
      if(dif == FALSE){
        A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
      }else{
        A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
      }
      for(k in s:K){
        ind_mat=ind_mat+1
        M = matrix(0, nrow = n, ncol = n)
        M[1,] = 1
        if(dif == FALSE){
          prec=c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
        }else{
          prec=c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
        }
        A[count,] = A[count,]+prec
      }
      count = count+1
    }
  }
  
  return(list(A, count))
}

# Restriccion (5) de la memoria
rest_5 = function(A, count, dif, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (5).
  # &\sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}}prec_{ijkqc} = \sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}}prec_{jikqc}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  for(cam in 1:Cam){
    for(j in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[,j] = 1
          M[j,] = -1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (6) de la memoria
rest_6 = function(A, count, dif, n, K, Q, Cam, n_mat, C, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (6).
  # cant_{ikqc} \leq DEM_i\cdot prec_{1ikqc} + CAP_c(1-prec_{1ikqc})
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # C: vector de capacidades de los agentes
  # dem: vector de demanda de los nodos
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(j in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[1,j] = C[cam]-dem[j-1]
          cantidades = rep(0,n_mat*(n-1))
          cantidades[(j-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (7) de la memoria
rest_7 = function(A, count, dif, n, K, Q, Cam, n_mat, C, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (7).
  # cant_{ikqc} \geq DEM_i\cdot prec_{1ikqc} - CAP_c(1-prec_{1ikqc})
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # C: vector de capacidades de los agentes
  # dem: vector de demanda de los nodos
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(j in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[1,j] = -(C[cam]+dem[j-1])
          cantidades = rep(0,n_mat*(n-1))
          cantidades[(j-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (8) de la memoria
rest_8 = function(A, count, dif, n, K, Q, Cam, n_mat, C){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (8).
  # cant_{ikqc} \leq CAP_c  \cdot \sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}} prec_{jikqc}\\
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # C: vector de capacidades de los agentes
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[,i] = -C[cam]
          cant = rep(0,n_mat*(n-1))
          cant[(i-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cant, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cant, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (10) de la memoria
rest_10 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (10).
  # acum_{ikqc} \leq (t_{1i}+SAL_k-T_{max})\cdot prec_{1ikqc}+ T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[1,i] = -(D[1,i]+salida[k]-tmax)
          acum = rep(0,n_mat*(n-1))
          acum[(i-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (11) de la memoria
rest_11 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (11).
  # acum_{ikqc} \geq (t_{1i}+SAL_k+T_{max})\cdot prec_{1ikqc}- T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[1,i] = -(D[1,i]+salida[k]+tmax)
          acum = rep(0,n_mat*(n-1))
          acum[(i-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (12) de la memoria
rest_12 = function(A, count, dif, n, K, Q, Cam, n_mat, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (12).
  # acum_{ikqc} \leq T_{max} \cdot \sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}} prec_{jikqc}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[,i] = -tmax
          acum = rep(0,n_mat*(n-1))
          acum[(i-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          count = count+1
          ind_mat = ind_mat+1
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (14) de la memoria
rest_14 = function(A, count, dif, n, K, Q, Cam, n_mat, C, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (14).
  # cant_{jkqc}-cant_{ikqc} \leq (DEM_j-CAP_c)\cdot prec_{ijkqc} +CAP_c
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # C: vector de capacidades de los agentes
  # dem: vector de demanda de los nodos
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      for(j in 2:n){
        if(i!=j){
          ind_mat = 1
          for(k in 1:K){
            for(s in 1:k){
              M = matrix(0, nrow = n, ncol = n)
              M[i,j] = C[cam]-dem[j-1]
              cantidades = rep(0,n_mat*(n-1))
              cantidades[(i-1)+(n-1)*(ind_mat-1)] = -1
              cantidades[(j-1)+(n-1)*(ind_mat-1)] = 1
              if(dif == FALSE){
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)))
              }else{
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
              }
              count = count+1
              ind_mat = ind_mat+1
            }
          }
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (15) de la memoria
rest_15 = function(A, count, dif, n, K, Q, Cam, n_mat, C, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (15).
  # cant_{jkqc}-cant_{ikqc} \geq (DEM_j+CAP_c)\cdot prec_{ijkqc} - CAP_c
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # C: vector de capacidades de los agentes
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      for(j in 2:n){
        if(i!=j){
          ind_mat = 1
          for(s in 1:Q){
            for(k in s:K){
              M = matrix(0, nrow = n, ncol = n)
              M[i,j] = -(C[cam]+dem[j-1])
              cantidades = rep(0,n_mat*(n-1))
              cantidades[(i-1)+(n-1)*(ind_mat-1)] = -1
              cantidades[(j-1)+(n-1)*(ind_mat-1)] = 1
              if(dif == FALSE){
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)))
              }else{
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), cantidades, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
              }
              count = count+1
              ind_mat = ind_mat+1
            }
          }
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (16) de la memoria
rest_16 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (16).
  # acum_{jkqc}-acum_{ikqc} \leq (t_{ij}-T_{max})\cdot prec_{ijkqc}+ T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      for(j in 2:n){
        if(i!=j){
          ind_mat = 1
          for(s in 1:Q){
            for(k in s:K){
              M = matrix(0, nrow = n, ncol = n)
              M[i,j] = -(D[i,j]-tmax)
              acum = rep(0,n_mat*(n-1))
              acum[(i-1)+(n-1)*(ind_mat-1)] = -1
              acum[(j-1)+(n-1)*(ind_mat-1)] = 1
              if(dif == FALSE){
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
              }else{
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
              }
              count = count+1
              ind_mat = ind_mat+1
            }
          }
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (17) de la memoria
rest_17 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (17).
  # acum_{jkqc}-acum_{ikqc} \geq (t_{ij}+T_{max})\cdot prec_{ijkqc}- T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      for(j in 2:n){
        if(i!=j){
          ind_mat = 1
          for(s in 1:Q){
            for(k in s:K){
              M = matrix(0, nrow = n, ncol = n)
              M[i,j] = -(D[i,j]+tmax)
              acum = rep(0,n_mat*(n-1))
              acum[(i-1)+(n-1)*(ind_mat-1)] = -1
              acum[(j-1)+(n-1)*(ind_mat-1)] = 1
              if(dif == FALSE){
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
              }else{
                A[count,] = c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
              }
              count = count+1
              ind_mat = ind_mat+1
            }
          }
        }
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (19) de la memoria
rest_19 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (19).
  # \sum_{s=1}^q \sum_{k=s}^K (acum_{iksc}+t_{i1}\cdot prec_{i1ksc}) \leq \sum_{k=q+1}^K \sum_{j\in P} (SAL_{k}-T_{max}) \cdot prec_{1jk(q+1)c}+T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      for(s in 1:(Q-1)){
        if(dif == FALSE){
          A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
        }else{
          A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
        }
        ind_mat=1
        for(w in 1:s){
          for(k in w:K){
            M = matrix(0, nrow = n, ncol = n)
            M[i,1] = D[i,1]
            acum = rep(0,n_mat*(n-1))
            acum[(i-1)+(n-1)*(ind_mat-1)] = 1
            if(dif == FALSE){
              A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
            }else{
              A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
            }
            ind_mat = ind_mat+1
          }
        }
        
        for(k in (s+1):K){
          M = matrix(0, nrow = n, ncol = n)
          M[1,] = tmax-salida[k]
          if(dif == FALSE){
            A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-(ind_mat))), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-(ind_mat))), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          ind_mat = ind_mat+1
        }
        
        count = count+1
      }
    }
  }
  
  return(list(A, count))
}

# Restriccion (20) de la memoria
rest_20 = function(A, count, dif, n, K, Q, Cam, n_mat, D, salida, tmax){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (20).
  # \forall i\in P, \forall c: \sum_{q=1}^Q \sum_{k=q}^K (acum_{ikqc}+t_{i1}\cdot prec_{i1kqc}) \leq T_{max}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # salida: vector de horas para las posibles salidas
  # tmax: tiempo limite
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(cam in 1:Cam){
    for(i in 2:n){
      if(dif == FALSE){
        A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)))
      }else{
        A[count,] = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
      }
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          M = matrix(0, nrow = n, ncol = n)
          M[i,1] = D[i,1]
          acum = rep(0,n_mat*(n-1))
          acum[(i-1)+(n-1)*(ind_mat-1)] = 1
          if(dif == FALSE){
            A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)))
          }else{
            A[count,] = A[count,] + c(rep(0, (cam-1)*n_mat*n*(n-1)), rep(0,n*(n-1)*(ind_mat-1)), no_diag(M), rep(0,n*(n-1)*(n_mat-ind_mat)), rep(0, (Cam-cam)*n_mat*n*(n-1)), rep(0, Cam*n_mat*(n-1)), rep(0, (cam-1)*n_mat*(n-1)), acum, rep(0, (Cam-cam)*n_mat*(n-1)), rep(0, (n-1)), 0)
          }
          ind_mat = ind_mat+1
        }
      }
      
      count = count+1
    }
  }
  
  return(list(A, count))
}

# Restriccion (21) de la memoria
rest_21 = function(A, count, dif, n, K, Q, Cam, n_mat, h0){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (21).
  # u_i \leq \sum_{q=1}^Q \sum_{k=q}^K \sum_{c=1}^C\frac{acum_{ikqc}}{H_{0i}}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # h0: vector extremos izquierdos intervalo preferido de llegada
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(i in 2:n){
    acum = rep(0, Cam*(n-1)*n_mat)
    for(cam in 1:Cam){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          acum[(cam-1)*(n-1)*n_mat+(i-1)+(n-1)*(ind_mat-1)] = -1/(h0[i-1])
          ind_mat = ind_mat+1
        }
      }
    }
    
    u = rep(0, (n-1))
    u[i-1] = 1
    if(dif == FALSE){
      A[count,] = c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), acum, u)
    }else{
      A[count,] = c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), acum, u, 0)
    }
    count = count+1
  }
  
  return(list(A, count))
}

# Restriccion (22) de la memoria
rest_22 = function(A, count, dif, n, K, Q, Cam, n_mat, tmax, hf){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (22).
  # u_i \leq \sum_{q=1}^Q \sum_{k=q}^K \sum_{c=1}^C\left(\frac{acum_{ikqc}}{H_{fi}-T_{max}}\right) -\frac{T_{max}}{H_{fi}-T_{max}}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # tmax: tiempo limite
  # hf: vector extremos derechos intervalo preferido de llegada
  
  #Returns:
  # list(A, count, indep_der): lista que tiene por elementos la matriz A modificada, el valor actualizado de count y los terminos independientes asociados a la restriccion 22
  
  
  indep_der = c()
  for(i in 2:n){
    acum = rep(0, Cam*(n-1)*n_mat)
    indep_der = c(indep_der, -tmax/(hf[i-1]-tmax))
    for(cam in 1:Cam){
      ind_mat = 1
      for(s in 1:Q){
        for(k in s:K){
          acum[(cam-1)*(n-1)*n_mat+(i-1)+(n-1)*(ind_mat-1)] = -1/(hf[i-1]-tmax)
          ind_mat = ind_mat+1
        }
      }
    }
    
    u = rep(0, (n-1))
    u[i-1] = 1
    if(dif == FALSE){
      A[count,] = c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), acum, u)
    }else{
      A[count,] = c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), acum, u, 0)
    }
    count = count+1
  }
  
  return(list(A, count, indep_der))
}

# Restriccion (23) de la memoria
rest_23 = function(A, count, dif, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (23).
  # \forall i \in P: u_i \leq 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # dif: variable booleana que es cierta si el modelo es difuso
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(i in 2:n){
    u=rep(0, (n-1))
    u[i-1]=1
    if(dif == FALSE){
      A[count,]=c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), u)
    }else{
      A[count,]=c(rep(0, Cam*n*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), rep(0, Cam*(n-1)*n_mat), u, 0)
    }
    count = count+1
  }
  
  return(list(A, count))
}

# Restriccion (31) de la memoria
rest_31 = function(A, count, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (31).
  # \forall j \in P : \sum_{c=1}^{C}\sum_{q=1}^Q\sum_{k=q}^K\sum_{\begin{subarray}{l} i \in L\\ i\neq j\end{subarray}}prec_{ijkqc} \leq 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(j in 2:n){
    prec=c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[,j] = 1
          prec= c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 0)
    count = count +1
  }
  
  return(list(A, count))
}

# Restriccion (32) de la memoria
rest_32 = function(A, count, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (32).
  # \forall i \in P: \sum_{c=1}^{C}\sum_{q=1}^Q\sum_{k=q}^K\sum_{\begin{subarray}{l} j \in L\\ i\neq j\end{subarray}} prec_{ijkqc} \leq 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(i in 2:n){
    prec=c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[i,] = 1
          prec= c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0,(n-1)), 0)
    count = count +1
  }
  
  return(list(A, count))
}

# Restriccion (33) de la memoria
rest_33 = function(A, count, n, K, Q, Cam, n_mat, D, d0){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (33).
  # TT + d_0 \cdot \alpha \leq TT_{max} + d_0
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # d0: margen de tiempo por la derecha
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  A[count,] = c(rep(no_diag(D), Cam*n_mat), rep(0, 2*Cam*(n-1)*n_mat), rep(0, (n-1)), d0)
  count = count+1
  
  return(list(A, count))
}

# Restriccion (34) de la memoria
rest_34 = function(A, count, n, K, Q, Cam, n_mat, D, d1){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (34).
  # -TT + d_1 \cdot \alpha \leq -TT_{min} + d_1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # D: matriz de tiempos
  # d1: margen de tiempo por la izquierda
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  A[count,] = c(rep(no_diag(-D), Cam*n_mat), rep(0, 2*Cam*(n-1)*n_mat), rep(0, (n-1)), d1)
  count = count+1
  
  return(list(A, count))
}

# Restriccion (35) de la memoria
rest_35 = function(A, count, n, K, Q, Cam, n_mat, u_0){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (35).
  # -\sum_{i \in P} u_i +u_0 \cdot \alpha \leq -U_{min} + u_0
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # u_0: margen de utilidad
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  A[count,] = c(rep(0, Cam*n_mat*n*(n-1)), rep(0, 2*Cam*(n-1)*n_mat), rep(-1, (n-1)), u_0)
  count = count+1
  
  return(list(A, count))
}

# Restriccion (36) de la memoria
rest_36 = function(A, count = 1, n, K, Q, Cam, n_mat, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (36).
  #\forall j \in P : \frac{-Dem_j \cdot IN}{\sum_{w \in P} Dem_w} + \alpha \leq \frac{\sum_{w \neq j} Dem_w}{\sum_{w \in P} Dem_w}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # dem: vector de demanda de los nodos
  
  #Returns:
  # list(A, count, indep_dem): lista que tiene por elementos la matriz A modificada, el valor actualizado de count y los terminos independientes asociados a la restriccion 36
  
  
  indep_dem = c()
  
  for(j in 2:n){
    prec = c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[,j] = -dem[(j-1)]/sum(dem)
          prec = c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 1)
    count = count +1
    indep_dem = c(indep_dem, sum(dem[-(j-1)])/sum(dem))
  }
  
  return(list(A, count, indep_dem))
}

# Restriccion (37) de la memoria
rest_37 = function(A, count = 1, n, K, Q, Cam, n_mat, dem){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (37).
  #\forall i \in P : \frac{-Dem_i \cdot OUT}{\sum_{w \in P} Dem_w} + \alpha \leq \frac{\sum_{w \neq i} Dem_w}{\sum_{w \in P} Dem_w}
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # dem: vector de demanda de los nodos
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  for(i in 2:n){
    prec=c()
    for(cam in 1:Cam){
      for(k in 1:K){
        for(s in 1:k){
          M = matrix(0, nrow = n, ncol = n)
          M[i,] = -dem[(i-1)]/sum(dem)
          prec= c(prec, no_diag(M))
        }
      }
    }
    A[count,] = c(prec, rep(0, 2*Cam*n_mat*(n-1)), rep(0,(n-1)), 1)
    count = count +1
  }
  
  return(list(A, count))
}

# Restriccion (39) de la memoria
rest_39 = function(A, count = 1, n, K, Q, Cam, n_mat){
  
  # Permite calcular elementos de la matriz de coeficientes del sistema de ecuaciones del problema, asociados a la restriccion (39).
  # alpha \leq 1
  
  #Args:
  # A: matriz de coeficientes del sistema de ecuaciones del problema
  # count: ordinal de la ecuacion que estamos determinando
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  
  #Returns:
  # list(A, count): lista que tiene por elementos la matriz A modificada y el valor actualizado de count
  
  
  A[count,] = c(rep(0, Cam*n_mat*n*(n-1)), rep(0, 2*Cam*(n-1)*n_mat), rep(0, (n-1)), 1)
  
  return(list(A, count))
}

# Informacion de la solucion
info_sol = function(solution, dif, n, K, Q, Cam, n_mat, n_var){
  
  # Permite obtener toda la informacion deseada de la solucion del problema.
  
  #Args:
  # solution: resultado de la funcion "lp"
  # dif: variable booleana que indica si estamos en la version difusa
  # n: numero de nodos
  # K: numero de salidas posibles
  # Q: cantidad de ciclos permitidos por camion
  # Cam: numero de camiones
  # n_mat: total de ciclos considerados por camion
  # n_var: numero de variables que tiene el problema
  
  # Soluciones
  M_sol = matrix(solution$solution[1:(Cam*(n-1)*n*n_mat)], ncol=n-1, byrow = TRUE)
  cant_sol = matrix(solution$solution[(Cam*(n-1)*n*n_mat+1):(Cam*(n-1)*n*n_mat+Cam*(n-1)*n_mat)], ncol=n-1, byrow = TRUE)
  acum_sol = matrix(solution$solution[(Cam*(n-1)*n*n_mat+Cam*(n-1)*n_mat+1):(Cam*(n-1)*n*n_mat+Cam*2*(n-1)*n_mat)], ncol=(n-1), byrow = TRUE)
  if(dif == TRUE){
    u_sol = solution$solution[(n_var-1-(n-2)):(n_var-1)]
  }else{
    u_sol = solution$solution[(n_var-(n-2)):n_var]
  }
  
  
  # Camino de la solucion
  for(cam in 1:Cam){
    print("CAMION:")
    print(cam)
    ind_mat = 1
    for(s in 1:Q){
      for(k in s:K){
        # Miramos si hay alguna salida
        salidas_origen = which(M_sol[(1+n*(ind_mat-1)+(cam-1)*n*(n_mat)),] == 1) + 1
        if(length(salidas_origen) != 0){
          print("Vuelta:")
          print(s)
          print("salida:")
          print(k)
          for(i in salidas_origen){
            # Nodo al que vamos desde la salida
            nodo_actual = i
            
            # Añadimos el nodo a la ruta
            ruta = c(1, nodo_actual)
            
            # Miramos las cantidades y tiempos acumulados en esa ruta
            cantidades = cant_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
            acumulados = acum_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
            
            # Repetimos el bucle hasta regresar al origen
            while(nodo_actual != 1){
              # Miremos el siguiente nodo de la ruta
              aux = which(M_sol[(nodo_actual+n*(ind_mat-1)+(cam-1)*n*(n_mat)),] != 0)
              
              # Ajsutamos a la notacion utilizada
              if(aux < nodo_actual){
                nodo_actual = aux
                aux_c = cant_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
                aux_a = acum_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
              }else{
                nodo_actual = aux + 1
                aux_c = cant_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
                aux_a = acum_sol[(ind_mat+(cam-1)*(n_mat)), nodo_actual-1]
              }
              
              # Añadimos los datos a la ruta
              ruta = c(ruta, nodo_actual)
              cantidades = c(cantidades, aux_c)
              acumulados = c(acumulados, aux_a)
            }
            
            # Mostramos la ruta calculada
            print(ruta)
            print(cantidades)
            print(acumulados)
          }
        }
        ind_mat = ind_mat + 1
      }
    }
  }
  
  print("Las utilidades conseguidas son:")
  print(u_sol)
  print(sum(u_sol))
  
  if(dif){
    print("El alpha obtenido es:")
    print(solution$objval)
  }
}


################################################################################
# FUNCION SELECCION PARAMETROS
################################################################################
parametros = function(dif, ver){
  # funcion que nos permite los parametros que trabajaremos
  
  #Args:
  # dif: variable booleana que indica si estamos en la version difusa
  # ver: indica si tratamos el caso simple, el caso practico o la comparativa de funciones
  
  
  # SELECCION DE PARAMETROS
  if(ver == 1){
    # Numero de nodos
    n = 4
    
    # Matriz de tiempos
    D = diag(rep(0,n))
    D[lower.tri(D)] = c(2, 4, 5,
                        2, 3,
                        7)
    D = t(D)
    D[lower.tri(D)] = c(2, 4, 5,
                        2, 3,
                        7)
    
    # Vector de demandas
    dem = c(3, 5, 1)
    
    # Numero de salidas posibles
    K = 3
    
    # Numero de ciclos posibles 
    Q = 3
    
    # Cantidad de matrices ij para cada camion
    aux = 0
    for(i in 0:Q){
      aux = aux+K-i
    }
    n_mat = aux
    rm(aux)
    
    # Cantidad de camiones
    Cam = 1
    
    # Capacidad máxima del camion
    C = c(8)
    
    # Posibles horas de salida
    salida = c(0, 10, 20)
    
    # Tiempo maximo
    tmax = 30 
    
    # Horarios preferidos
    h0 = c(12, 1, 12)
    hf = c(16, 5, 16)
    
  }else if(ver == 2){
    # Numero de nodos
    n = 8
    
    # Matriz de tiempos
    D = diag(rep(0,n))
    D[lower.tri(D)] = c(20, 25, 23, 26, 29, 20, 20,
                        18, 31, 36, 44, 32, 10,
                        22, 25, 42, 33, 9,
                        10, 27, 30, 24,
                        25, 34, 27, 
                        16, 42,
                        32)
    D = t(D)
    D[lower.tri(D)] = c(20, 25, 23, 26, 29, 20, 20,
                        18, 31, 36, 44, 32, 10,
                        22, 25, 42, 33, 9,
                        10, 27, 30, 24,
                        25, 34, 27, 
                        16, 42,
                        32)
    
    # Vector de demandas
    dem = c(16,42,76,37,14,12,21)
    
    # Numero de salidas posibles
    K = 4
    
    # Numero de ciclos posibles
    Q = 4
    
    # Cantidad de matrices ij para cada camion
    aux = 0
    for(i in 0:Q){
      aux = aux+K-i
    }
    n_mat = aux
    rm(aux)
    
    # Cantidad de camiones
    Cam = 2
    
    # Capacidad maxima de cada camión
    C = c(150, 70)
    
    # Posibles horas de salida
    # 8.5, 12.5, 15, 18
    salida = c(510, 750, 900, 1080)
    
    # Tiempo maximo
    tmax = 1230 
    
    # Horarios preferidos
    h0 = c(510, 510, 510, 765, 765, 900, 1080)
    hf = c(570, 570, 570, 840, 840, 960, 1140)
    
  }else if(ver == 3){
    # Numero de nodos
    n = 8
    
    # Matriz de tiempos
    D = diag(rep(0,n))
    D[lower.tri(D)] = c(20, 25, 23, 26, 29, 20, 20,
                        18, 31, 36, 44, 32, 10,
                        22, 25, 42, 33, 9,
                        10, 27, 30, 24,
                        25, 34, 27, 
                        16, 42,
                        32)
    D = t(D)
    D[lower.tri(D)] = c(20, 25, 23, 26, 29, 20, 20,
                        18, 31, 36, 44, 32, 10,
                        22, 25, 42, 33, 9,
                        10, 27, 30, 24,
                        25, 34, 27, 
                        16, 42,
                        32)
    
    # Vector de demandas
    dem = c(16,42,76,37,14,12,21)
    
    # Numero de salidas posibles
    K = 2
    
    # Numero de ciclos posibles
    Q = 2
    
    # Cantidad de matrices ij para cada camion
    aux = 0
    for(i in 0:Q){
      aux = aux+K-i
    }
    n_mat = aux
    rm(aux)
    
    # Cantidad de camiones
    Cam = 1
    
    # Capacidad maxima de cada camión
    C = c(150)
    
    # Posibles horas de salida
    # 8.5, 12.5, 15, 18
    salida = c(510, 750, 900, 1080)
    
    # Tiempo maximo
    tmax = 1230 
    
    # Horarios preferidos
    h0 = c(510, 510, 510, 510, 510, 510, 510)
    hf = c(540, 540, 540, 540, 540, 540, 540)
  }
  
  
  if(dif == TRUE){
    if(ver == 1){
      # Tiempo total recorrido maximo
      TT_max = 11
      d0 = 3
      
      # Tiempo total recorrido minimo
      TT_min = 5
      d1 = 2
      
      # Utilidad minima
      U_min = 2.1
      u_0 = 0.5
      
    }else if(ver == 2){
      # Tiempo total recorrido maximo
      TT_max = 240
      d0 = 60
      
      # Tiempo total recorrido minimo
      TT_min = 120
      d1 = 30
      
      # Utilidad minima
      U_min = 7
      u_0 = 0.25
    }
    
  }else{
    TT_max = 0
    d0 = 0
    TT_min = 0
    d1 = 0
    U_min = 0
    u_0 = 0
  }
  
  return(list(n, D, dem, K, Q, n_mat, Cam, C, salida, tmax, h0, hf, TT_max, d0, TT_min, d1, U_min, u_0))
}


################################################################################
# FUNCION PRINCIPAL (MAIN)
################################################################################
main = function(dif, ver, fun_obj){
  # funcion principal que resolvera el problema concreto deseado, de los presentados en la memoria
  
  #Args:
  # dif: variable booleana que indica si estamos en la version difusa
  # ver: indica si tratamos el caso simple, el caso practico o la comparativa de funciones
  # fun_obj: indica cual de las funciones objetivo usaremos para resolver el problema
  
  
  # SELECCION DE PARAMETROS
  param = parametros(dif, ver)
  
  # Numero de nodos
  n = param[[1]]
  
  # Matriz de tiempos
  D = param[[2]]
  
  # Vector de demandas
  dem = param[[3]]
  
  # Numero de salidas posibles
  K = param[[4]]
  
  # Numero de ciclos posibles 
  Q = param[[5]]
  
  # Cantidad de matrices ij para cada camion
  n_mat = param[[6]]

  # Cantidad de camiones
  Cam = param[[7]]
  
  # Capacidad máxima del camion
  C = param[[8]]
  
  # Posibles horas de salida
  salida = param[[9]]
  
  # Tiempo maximo
  tmax = param[[10]] 
  
  # Horarios preferidos
  h0 = param[[11]]
  hf = param[[12]]
  
  # Tiempo total recorrido maximo
  TT_max = param[[13]]
  d0 = param[[14]]
  
  # Tiempo total recorrido minimo
  TT_min = param[[15]]
  d1 = param[[16]]
  
  # Utilidad minima
  U_min = param[[17]]
  u_0 = param[[18]]
    
  
  # SELECCION DE FUNCIONES OBJETIVO 
  if(dif == TRUE){
    # CASO DIFUSO
    if(fun_obj == 1){
      # Asignamos valor 1 solo al coeficiente de alpha
      coef = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0, (n-1)), 1)
      
      # Queremos encontrar un maximo
      tipo_op = "max"
      
    } else if(fun_obj == 2){
      # Asignamos valor 1 solo al coeficiente de alpha y 1e-6 a los coeficientes de la utilidad 
      coef = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(0.000001, (n-1)), 1)
      tipo_op = "max"
    }
    
    # Numero de ecuaciones
    n_ec = Cam*Q+Cam*(n-1)*(Q+n_mat*(7+4*(n-2)))+7*(n-1)+4
    
  }else{
    # CASO CRISP
    if(fun_obj == 1){
      # Asignamos valor 1 solo a los coeficientes de la utilidad
      coef = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), rep(1, (n-1)))
      
      # Queremos encontrar un maximo
      tipo_op = "max"
      
    } else if(fun_obj == 2){
      # Asignamos el valor de la demanda a los coeficientes de la utilidad
      u=c()
      for(i in 2:n){
        u=c(u, dem[i-1])
      }
      coef = c(rep(0, Cam*n_mat*(n)*(n-1)), rep(0, 2*Cam*n_mat*(n-1)), u)
      # Queremos encontrar un maximo
      tipo_op = "max"
    }
    
    # Numero de ecuaciones
    n_ec = Cam*Q+Cam*(n-1)*(Q+n_mat*(7+4*(n-2)))+5*(n-1)
  }
  
  
  # PLANTEAMIENTO DEL PROBLEMA
  if(dif ==FALSE){
    # CONSTRUCCION DE LA MATRIZ DE COEFICIENTES CASO CRISP
    # Matriz de coeficientes
    A = matrix(0, nrow = n_ec, ncol = length(coef))
    
    # Completamos la matriz con los coeficientes obtenidos de las funciones de cada restriccion
    list_aux = rest_2(A, 1, n, K, Q, Cam, n_mat)
    list_aux = rest_3(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat)
    list_aux = rest_5(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    list_aux = rest_4(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    list_aux = rest_6(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_14(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_7(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_15(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_10(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_16(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_11(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_17(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_8(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C)
    list_aux = rest_12(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, tmax)
    list_aux = rest_19(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_20(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_21(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, h0)
    list_aux = rest_22(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, tmax, hf)
    indep_der = list_aux[[3]]
    list_aux = rest_23(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    A = list_aux[[1]]
    count = list_aux[[2]]
    
    # Terminos independientes
    C_indep1=c()
    C_indep2=c()
    for(cam in 1:Cam){
      C_indep1=c(C_indep1, rep(C[cam], (n-1)*n_mat))
      C_indep2=c(C_indep2, rep(-C[cam], (n-1)*n_mat))
    }
    for(cam in 1:Cam){
      C_indep1=c(C_indep1, rep(C[cam], (n-2)*(n-1)*n_mat))
      C_indep2=c(C_indep2, rep(-C[cam], (n-2)*(n-1)*n_mat))
    }
    
    # Vector de terminos independientes
    b = c(rep(1, (n-1)*2), rep(0, Cam*(n-1)*n_mat), rep(1, Q*Cam), C_indep1, C_indep2, rep(tmax, Cam*(n-1)*n_mat), rep(tmax, Cam*(n-2)*(n-1)*n_mat), rep(-tmax, Cam*(n-1)*(n-1)*n_mat), rep(0, 2*Cam*(n-1)*n_mat), rep(tmax, Cam*(n-1)*(Q)),  rep(0, (n-1)), indep_der, rep(1, (n-1)))
    
    # Direccion
    dir = c(rep('=', (n-1)*2), rep("=", Cam*(n-1)*n_mat), rep('<=', Cam*(Q+(n-1)*(n-1)*n_mat)), rep('>=', Cam*(n-1)*(n-1)*n_mat), rep('<=', Cam*(n-1)*(n-1)*n_mat), rep('>=', Cam*(n-1)*(n-1)*n_mat), rep("<=", 2*Cam*(n-1)*n_mat), rep("<=", Cam*(n-1)*Q+3*(n-1)))
    
    # CALCULO DE LA SOLUCION
    set.seed(2023)
    solucion = lp(tipo_op, coef, A, dir, b, binary.vec = 1:(Cam*(n-1)*n*n_mat))
    
    info_sol(solucion, dif, n, K, Q, Cam, n_mat, length(coef))
    
  }else{
    # CONSTRUCCION DE LA MATRIZ DE COEFICIENTES CASO DIFUSO
    # Matriz de coeficientes
    A = matrix(0, nrow = n_ec, ncol = length(coef))
    
    # Completamos la matriz con los coeficientes obtenidos de las funciones de cada restriccion
    list_aux = rest_36(A, 1, n, K, Q, Cam, n_mat, dem)
    indep_dem = list_aux[[3]]
    list_aux = rest_37(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat, dem)
    list_aux = rest_31(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat)
    list_aux = rest_32(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat)
    list_aux = rest_5(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    list_aux = rest_4(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    list_aux = rest_6(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_14(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_7(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_15(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C, dem)
    list_aux = rest_10(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_16(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_11(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_17(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_8(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, C)
    list_aux = rest_12(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, tmax)
    list_aux = rest_19(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_20(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, D, salida, tmax)
    list_aux = rest_21(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, h0)
    list_aux = rest_22(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat, tmax, hf)
    indep_der = list_aux[[3]]
    list_aux = rest_23(list_aux[[1]], list_aux[[2]], dif, n, K, Q, Cam, n_mat)
    list_aux = rest_33(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat, D, d0)
    list_aux = rest_34(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat, D, d1)
    list_aux = rest_35(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat, u_0)
    list_aux = rest_39(list_aux[[1]], list_aux[[2]], n, K, Q, Cam, n_mat)
    A = list_aux[[1]]
    count = list_aux[[2]]
    
    # Terminos independientes
    C_indep1=c()
    C_indep2=c()
    for(cam in 1:Cam){
      C_indep1=c(C_indep1, rep(C[cam], (n-1)*n_mat))
      C_indep2=c(C_indep2, rep(-C[cam], (n-1)*n_mat))
    }
    for(cam in 1:Cam){
      C_indep1=c(C_indep1, rep(C[cam], (n-2)*(n-1)*n_mat))
      C_indep2=c(C_indep2, rep(-C[cam], (n-2)*(n-1)*n_mat))
    }
    
    # Vector de terminos independientes
    b = c(rep(indep_dem, 2), rep(1, 2*(n-1)), rep(0, Cam*(n-1)*n_mat), rep(1, Q*Cam), C_indep1, C_indep2, rep(tmax, Cam*(n-1)*n_mat), rep(tmax, Cam*(n-2)*(n-1)*n_mat), rep(-tmax, Cam*(n-1)*(n-1)*n_mat), rep(0, 2*Cam*(n-1)*n_mat), rep(tmax, Cam*(n-1)*(Q)),  rep(0, (n-1)), indep_der, rep(1, (n-1)), TT_max+d0, -TT_min+d1, -U_min+u_0, 1)
    
    # Direccion
    dir = c(rep('<=', (n-1)*4), rep("=", Cam*(n-1)*n_mat), rep('<=', Cam*(Q+(n-1)*(n-1)*n_mat)), rep('>=', Cam*(n-1)*(n-1)*n_mat), rep('<=', Cam*(n-1)*(n-1)*n_mat), rep('>=', Cam*(n-1)*(n-1)*n_mat), rep("<=", 2*Cam*(n-1)*n_mat), rep("<=", Cam*(n-1)*Q+3*(n-1)+4))
    
    # CALCULO DE LA SOLUCION
    set.seed(2023)
    solucion = lp(tipo_op, coef, A, dir, b, binary.vec = 1:(Cam*(n-1)*n*n_mat))
    
    info_sol(solucion, dif, n, K, Q, Cam, n_mat, length(coef))
  }
}


################################################################################
# PROGRAMA PRINCIPAL (MAIN)
################################################################################
# Para este programa consideramos:

# dif = TRUE: caso difuso
# dif = FALSE: caso crisp

# ver = 1: Caso simple de la memoria
# ver = 2: Caso practico de la memoria
# ver = 3: Comparativa de las funciones objetivo

# dif = FALSE y fun_obj = 1: maximo de la utilidad total (restr. 1 de la memoria)
# dif = FALSE y fun_obj = 2: maximo de la utilidad total considerando la demanda (restr. 25 de la memoria)
# dif = TRUE y fun_obj = 1: maximo de alpha (restr. 38 de la memoria)
# dif = TRUE y fun_obj = 2: maximo de alpha considerando epsilon (restr. 41 de la memoria)

main(dif = FALSE, ver = 2, fun_obj = 1)


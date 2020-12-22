############################################################################
## Código para generar laberintos utilizando el algoritmo de PRIM aleatorio
## Basado en el pseudocódigo encontrado en:
## https://stackoverflow.com/questions/29739751/implementing-a-randomly-generated-maze-using-prims-algorithm
## 19 de noviembre del 2020
## Pablo Barranco, Thomas Bladt, Alejandro De Anda
############################################################################

# Librerias necesarias ----------------------------------------------------
library(animation)
library(plot.matrix)

# Funciones auxiliares ---------------------------------------------------------------

# Crea la base del laberinto de mxn
creaGrid <- function(num_fil,num_col){
  grid_lab <- matrix(0, nrow = num_fil,ncol = num_col)
  return(grid_lab)
}

# Selecciona una cela aleatoriamente
celdaAleatoria <- function(num_fil,num_col){
  fila    <- sample(2:num_fil-1,size = 1, replace = F)
  columna <- sample(2:num_col-1,size = 1, replace = F)
  return(c(fila,columna))
}

# Checa si determinada celda cae dentro del laberinto
esLegal <- function(fila,columna, lab){
  num_fil <- nrow(lab)
  num_col <- ncol(lab)
  return( (fila > 1) &&(fila < num_fil) && (columna > 1)&& (columna < num_col) )
}

# Calcula la frontera de una celda específica
# Frontera: celdas a distancia 2 de tipo "cerradas"
frontera <- function(fila,columna,lab){
  #Inicializamos un data frame
  frontera_df <-  data.frame(x = integer(), y = integer())
  
  #Consideramos los 4 puntos posibles
  if(esLegal(fila+2,columna,lab)){
    if(lab[fila+2,columna] == 0)
      frontera_df <- rbind(frontera_df,data.frame(x= fila+2, y = columna))
  }
  if(esLegal(fila-2,columna,lab)){
    if(lab[fila-2,columna] == 0)
      frontera_df <- rbind(frontera_df,data.frame(x= fila-2, y = columna))
  }
  if(esLegal(fila,columna+2,lab)){
    if(lab[fila,columna+2] == 0)
      frontera_df <-rbind(frontera_df,data.frame(x= fila, y = columna+2))
  }
  if(esLegal(fila,columna-2,lab)){
    if(lab[fila,columna-2] == 0)
      frontera_df <-rbind(frontera_df,data.frame(x= fila , y = columna-2))
  }
  
  return(frontera_df)
}

# Calcula las celdas vecinas a una específica
# Vecina: celdas a distancia 2 de tipo "abierta"
vecinos <- function(fila,columna,lab){
  vecinos_df <-  data.frame(x = integer(), y = integer())
  
  if(esLegal(fila+2,columna,lab)){
    if(lab[fila+2,columna] == 1)
      vecinos_df <- rbind(vecinos_df,data.frame(x= fila +2, y = columna))
  }
  if(esLegal(fila-2,columna,lab)){
    if(lab[fila-2,columna] == 1)
      vecinos_df <-rbind(vecinos_df,data.frame(x= fila-2, y = columna))
  }
  if(esLegal(fila,columna+2,lab)){
    if(lab[fila,columna+2] == 1)
      vecinos_df <-rbind(vecinos_df,data.frame(x= fila, y = columna+2))
  }
  if(esLegal(fila,columna-2,lab)){
    if(lab[fila,columna-2] == 1)
      vecinos_df <-rbind(vecinos_df,data.frame(x= fila , y = columna-2))
  }
  return(vecinos_df)
}

# Función principal -----------------------------------------------------------

genera_laberinto <- function(m,n,opcion){
#opcion 1 = return plotss, 0 = return laberinto
#set.seed(2025)

# Creamos un grid de m x n
num_fil <- m
num_col <- n
laberinto <- creaGrid(num_fil,num_col)

# Seleccionamos la celda inicial (puede ser aleatoria también)
#inicio <- celdaAleatoria(num_fil,num_col)
inicio <- c(num_fil-1,2)

# Definimos celda inicial como abierta
laberinto[inicio[1],inicio[2]] <- 1

# Inicializamos lista de la frontera en un dataframe
listaFrontera <- data.frame(x = integer(), y = integer())
listaFrontera <- rbind(listaFrontera,frontera(inicio[1],inicio[2], laberinto))

# Empezamos el contador de iteraciones
k=0

# Lista para guardar las matrices en cada iteración (para visualizaciones únicamente)
plotss <- list()

#Mientras no esté vacía la lista
while( nrow(listaFrontera) != 0 ){
  k=k+1
  # Inicializamos lista de vecinos vacía 
  vecinos_df <- data.frame(x = integer(), y = integer())
  
  # Seleccionamos aleatoriamente un elemento de la frontera 
  tamañoLista <- nrow(listaFrontera)
  celda <- sample(1:tamañoLista, size = 1, replace = F)
  
  # Calculamos a los vecinos del elemento frontera
  vecinos_df <- rbind(vecinos_df, vecinos(listaFrontera[celda,1],listaFrontera[celda,2],laberinto))
  
  # Seleccionamos un vecino aleatoriamente y lo conectamos con la frontera
  if(nrow(vecinos_df) != 0)
  {
     vecino <- sample(1:nrow(vecinos_df),size = 1)
     laberinto[(vecinos_df[vecino,1]+listaFrontera[celda,1])/2,  (vecinos_df[vecino,2]+listaFrontera[celda,2])/2] = 1
     laberinto[listaFrontera[celda,1],listaFrontera[celda,2]] <- 1
     #laberinto[vecinos_df[vecino,1],vecinos_df[vecino,2]] <- 1
  
  # Calculamos la frontera de la celda frontera actual y agregamos a la lista
  listaFrontera <- rbind(listaFrontera,frontera(listaFrontera[celda,1],listaFrontera[celda,2], laberinto))
  
  # Quitamos la frontera anterior de la lista
  listaFrontera <-  listaFrontera[-celda,]
  
  # Cuidamos que no se repitan elementos en la lista
  listaFrontera <- unique(listaFrontera)
  }
  
  # Guardamos el plot (para animación únicamente)
  plotss[[k]] <- laberinto
}

if(opcion == 1){
  
  # png("plot1_prim.png",width = 1800, height = 1300,res=250)
  # plot( plotss[[1]],col = c("black","white"),key = NULL,breaks = NULL,axis.col=NULL, axis.row=NULL,xlab = "", ylab="", main = "")
  # dev.off()
  # 
  # png("plot2_prim.png",width = 1800, height = 1300,res=250)
  # plot( plotss[[floor(k/2)]],col = c("black","white"),key = NULL,breaks = NULL,axis.col=NULL, axis.row=NULL,xlab = "", ylab="", main = "" )
  # dev.off()
  # 
  # png("plot3_prim.png",width = 1800, height = 1300,res=250)
  # plot( plotss[[k]],col = c("black","white"),key = NULL,breaks = NULL,axis.col=NULL, axis.row=NULL,xlab = "", ylab="", main = "")
  # dev.off()
  # En caso de querer guardar un gif de la creación
  saveGIF({
    for (i in 1:k) plot( plotss[[i]],col = c("black","white"),key = NULL,breaks = NULL,axis.col=NULL, axis.row=NULL,xlab = "", ylab="")
  },interval = 0.01)
  return(plotss)
}
else
  return(laberinto)
}




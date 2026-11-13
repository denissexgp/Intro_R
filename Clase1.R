
# suma

a <- 2
b <- 4

a + b


R.version.string

Sys.info()[["sysname"]]

"ggplot2" %in% rownames(installed.packages())

"dplyr" %in% rownames(installed.packages())

#library(ggplot2)
#library(dplyr)

#packageVersion("ggplot2")
#packageVersion("dplyr")



identificador <- "gene_101"
class(identificador)

longitud <- 450L 
class(450L)

proporcion <- 0.53
class(proporcion)


# Para conocer los datos
class(1.5)
typeof(1.5)
str(identificador)




genes<- c("gA", "gB", "gC", "gE")

length(genes)
genes[2]


longitudes <- c(150, 400, 220, 150)

genes[longitudes >= 300]

sum(longitudes >= 300)


#Ejercicio1 Vectores y selecciónPregunta

ids <- c("g1", "g2", "g3", "g4", "g5")
tam <- c(90, 450, 310, 120, 600)

ids[3]

ids[tam >= 300]

sum(tam >= 300)
#----------

valores <- c("12", "18", "sin_dato")
valores <- as.numeric(valores)

#valores faltantes
anyNA(valores)
is.na(valores)

# Matrices 

conteos <- matrix(
  c(12, 8, 4, 6, 10, 2),
  nrow = 3,
  dimnames = list(c("g1", "g2", "g3"), c("M1", "M2"))
)
conteos

conteos[2, 1]
conteos[, "M2"]
dim(conteos)


# arreglos 

arreglo <- array(1:8, dim = c(2, 2, 2))
arreglo

#Dataframe ( puedes tener columnas de diferente tipo)

anotaciones <- data.frame(
  gen = c("g1", "g2", "g3"),
  longitud = c(450, 120, 800),
  producto = c("transportador", "hipotetica", "enzima")
)

anotaciones

anotaciones[2, ]

# lista

muestra <- list(
  id = "M01",
  secuencias = c("ATGCGT", "CCATGA"),
  anotaciones = anotaciones
)

muestra$id
muestra[["secuencias"]] #Extraer la lista, es un vector 
muestra["secuencias"] # Conserva la estructura de una lista

# factores  (se le asiga un orden a las variables categoricas )

condicion <- factor(c("control", "tratamiento", "control"),
                    levels = c("control", "tratamiento"))
levels(condicion)

# Ejercicio 2: Factores

altura <- factor(c("A", "M", "A", "B", "M", "M", "B", "A", "B", "B"),
            levels = c("A", "M", "B" ))
levels(altura)

# condicionales 

secuencia <- "ATGCGT"
longitud <- nchar(secuencia)

#if 
if (longitud <= 10) {
  print("La secuencia tiene al menos 10 nucleótidos")
}

#else 
if (longitud >= 10) {
  print("La secuencia tiene al menos 10 nucleótidos")
} else {
  print("La secuencia tiene menos de 10 nucleótidos")
}

#else if
if (longitud >= 20) {
  print("Secuencia larga")
} else if (longitud >= 10) {
  print("Secuencia mediana")
} else {
  print("Secuencia corta")
}


# ejercicio 
conteo <- NA

cif (is.na(conteo)) {
  print("Dato faltante")
} else if (conteo < 0) {
  print("Conteo inválido")
} else if (conteo == 0) {
  print("Sin lecturas")
} else {
  print("Con lecturas")
}


# 

conteos <- c(0, 12, 4, 0, 25)

estado <- ifelse(conteos > 0, "Detectado", "No detectado")
estado

any(conteos == 0)  # Existe al menos un cero?

which(conteos == 0) #  que posiciones hay ceros?


# For 

for (elemento in coleccion) {
  # Instrucciones que se ejecutan para cada elemento
}


for (i in 1:5) {
  print(i)
}


secuencias <- c("ATCGC", "ATGCGT", "CCATGAA")

for (secuencia in secuencias) {
  print(nchar(secuencia))
}

longitudes <- numeric(length(secuencias))

for (i in seq_along(secuencias)) {
  longitudes[i] <- nchar(secuencias[i])
}

longitudes


#  Ejercicio 3. Longitudes de secuencias

secuencias <- c("ATG", "ATGCGT", "AATTCCGG","CG") 
longuitudes <- c()
clasificacion <- c()

for (i in seq_along(secuencias)) {
  longitudes[i] <- nchar(secuencias[i])
  
  if (longitudes[i] < 5) {
    clasificacion[i] <- "Corta"
  } else {
    clasificacion[i] <- "Larga"
  }
}

longitudes
clasificacion


#while

limite <- 50
numero <- 0

while (numero^2 <= limite) {
  numero <- numero + 1
}

numero

# ejercicio while 

secuencia
paste0(secuencia, "TA")


secuencia <- ""

while (nchar(secuencia) < 8) {
  secuencia <- paste0(secuencia, "A")
}

secuencia


 ##############################################################################
# clase 2 
# funciones 

longitud_adn <- function(secuencia) {
  nchar(secuencia)
}

longitud_adn("ATGCGTAA")

#_______________________________________________________________________________
contiene_motivo <- function(secuencia, motivo) {
  grepl(motivo, secuencia, fixed = TRUE) # fixted hace que busque el texto tal cual.
}

contiene_motivo("ATGCGTAA", "CGT")

#_______________________________________________________________________________

empieza_con <- function(secuencia, prefijo){
  startsWith(secuencia, prefijo)
}

empieza_con("ATGCGT", "ATG")


#_______________________________________________________________________________
preparar_secuencia <- function(secuencia, tipo = "ADN", mayusculas = TRUE) {
  if (mayusculas) {
    secuencia <- toupper(secuencia)
  }
  
  if (tipo == "ARN") {
    secuencia <- chartr("T", "U", secuencia)
  }
  
  secuencia
}

preparar_secuencia("atgcgt")


#_______________________________________________________________________________
resumir_adn <- function(id, secuencia) {
  secuencia <- toupper(secuencia)
  bases <- strsplit(secuencia, split = "", fixed = TRUE)[[1]]

  list(
    id = id,
    longitud = nchar(secuencia),
    conteo_gc = sum(bases %in% c("G", "C"))
  )
}

resultado <- resumir_adn("GEN_01", "ATGCGTAA")
resultado

#_______________________________________________________________________________

#ejercicio funcion 2 


resumir_adn <- function(id, secuencia) {
  secuencia <- toupper(secuencia)
  bases <- strsplit(secuencia, "", fixed = TRUE)[[1]]
  longitud <-  nchar(secuencia)
  conteo_gc <- sum(bases %in% c("G", "C"))
  
  list(
    id = id,
    longitud = longitud,
    conteo_gc = conteo_gc,
    proporcion_gc <- conteo_gc / longitud
  )
}

resumir_adn("GEN_01", "ATGCGTAA")

#_______________________________________________________________________________

validar_adn <- function(secuencia) {
  if (!is.character(secuencia) ||
      length(secuencia) != 1 ||
      is.na(secuencia) ||
      !grepl("^[ACGT]+$", toupper(secuencia))) {
    stop("Se requiere una secuencia no vacía con bases A, C, G o T.")
  }
  
  TRUE
}

validar_adn("ATGCGT")

#______________________________________________________________________________-
longitud_validada <- function(secuencia) {
  stopifnot(
      is.character(secuencia),
      length(secuencia) == 1,
      !is.na(secuencia)
  )
  
  nchar(secuencia)
}


longitud_validada("ATGCGT")

#______________________________________________________________________________



clasificacion_secuencias <- function(secuencias) {
  resultado <- logical(length(secuencias))
  
  for (i in seq_along(secuencias)) {
    secuencia <- secuencias[i]
    
    if (is.na(secuencia)) {
      resultado[i] <- FALSE
    } else {
      resultado[i] <- grepl("^[ACGT]+$", toupper(secuencia))
    }
  } 
  
  # Calculamos el conteo al finalizar el bucle
  validas <- sum(resultado)
  invalidas <- length(secuencias) - validas
  
  # Devolvemos tanto el vector booleano como el resumen en una lista
  list(
    vector_resultado = resultado,
    validas = validas,
    invalidas = invalidas
  )
}


secuencias <- c("ATGCGT", "ACNTGA", "ttacga", NA)
clasificacion_secuencias(secuencias)




#-------------

 


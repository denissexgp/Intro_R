





getwd()  # Directorio de trabajo actual

blast <- read.delim("data_raw/datos_sesion_02_blast/resultados_blast.tsv")
metadata <- read.csv("data_raw/datos_sesion_02_blast/metadata_secuencias.csv",
                     colClasses = c(codigo_muestra = "character"))

class(blast)


library(readr)

blast <- readr::read_tsv("data_raw/datos_sesion_02_blast/resultados_blast.tsv")


readLines("data_raw/datos_sesion_02_blast/blast_comentarios.tsv", n = 4)

blast_con_comentarios <- readr::read_tsv(
  "data_raw/datos_sesion_02_blast/blast_comentarios.tsv",
  comment = "#",
  show_col_types = FALSE 
)

head(blast_con_comentarios)
#_______________________________________________________________________________

blast_na <- readr::read_tsv(
  "data_raw/datos_sesion_02_blast/blast_valores_faltantes.tsv",
  na = c("", "NA"),
  show_col_types = FALSE
)
anyNA(blast_na)

#---
readr::spec(blast)



#____________


# Simulamos el contenido de un archivo TSV


blast_error <- read_tsv(
  I("secuencia\tidentidad
seq1\t98.5
seq2\t95.2
seq3\terror
seq4\t99.1"),
  col_types = cols(
    secuencia = col_character(),
    identidad = col_double()
  )
)


# Revisamos los datos importados
blast_error
# Consultamos los problemas de importación
problems(blast_error)




#_______________________________________________________________________________

#¿Cuántos registros y columnas tiene la tabla?
nrow(blast)      # Número de filas
dim(blast)       # Filas y columnas

  #¿Cuántos query_id diferentes aparecen?
length(blast$query_id)


  #¿Cuántos resultados tiene cada consulta?
table(blast$query_id)

  #¿Qué filas tienen pident mayor o igual a 95?
blast[blast$pident >= 95, ]

  #¿Hay valores faltantes?
anyNA(blast)

################################################################################

# Problemas comunes 

# Separador incorrecto
ruta <- "data_raw/datos_sesion_02_blast/blast_separador_incorrecto.csv"
mal <- readr::read_csv(ruta, show_col_types = FALSE)
names(mal)

bien <- readr::read_tsv(ruta, show_col_types = FALSE)
names(bien)

# ausencias y registros repetidos 


blast_na <- readr::read_tsv(
  "data_raw/datos_sesion_02_blast/blast_valores_faltantes.tsv",
  show_col_types = FALSE
)
colSums(is.na(blast_na))

# hay duplicados?
anyDuplicated(blast_na)

anyDuplicated(blast_na$query_id)


#Fasta 

lineas <- readLines("data_raw/datos_sesion_02_blast/secuencias.fasta")
head(lineas, 6)


encabezados <- grep("^>", lineas, value = TRUE)
encabezados

length(encabezados)




# Ejercicio 6 De los archivos a una descripción de los datos






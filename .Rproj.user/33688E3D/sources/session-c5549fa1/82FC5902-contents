# Ejercicio 6: De los archivos a una descripción de los datos
#Crea un script llamado 02_exploracion.R dentro de tu proyecto. El script debe:

  #1. Comprobar que existen resultados_blast.tsv, metadata_secuencias.csv y secuencias.fasta.
cat("Resultados_blast.tsv existe?", file.exists("data_raw/datos_sesion_02_blast/resultados_blast.tsv"), "\n")

cat("Resultados_blast.tsv existe?", file.exists("data_raw/datos_sesion_02_blast/metadata_secuencias.csv"), "\n")

cat("Resultados_blast.tsv existe?", file.exists("data_raw/datos_sesion_02_blast/secuencias.fasta"), "\n")


  #2. Importar las dos tablas con readr y conservar codigo_muestra como texto.
library(readr)

blast <- readr::read_tsv("data_raw/datos_sesion_02_blast/resultados_blast.tsv", 
                         col_types = readr::cols(codigo_muestra = readr::col_character()),
                         show_col_types = FALSE)

metadata <- readr::read_csv("data_raw/datos_sesion_02_blast/metadata_secuencias.csv", 
                            col_types = readr::cols(codigo_muestra = readr::col_character()),
                            show_col_types = FALSE)


  #3. Mostrar dimensiones, nombres y tipos de columnas.

dim(blast)  
dim(metadata)  

str(blast)
str(metadata)

  #4. Identificar los query_id distintos y contar los resultados por consulta.
unique(blast$query_id)
unique(metadata$query_id)

table(blast$query_id)
table(metadata$query_id)

  #5. Mostrar las filas con pident >= 90 y explicar qué representa este filtro.
blast[blast$pident >= 90, ]

# dice que  filas tienen pident mayor o igual a 90


#6. Leer el FASTA e indicar cuántos encabezados contiene.
lineas <- readLines("data_raw/datos_sesion_02_blast/secuencias.fasta")
encabezados <- grep("^>", lineas, value = TRUE)

cat("el numero de encabezados es:", length(encabezados), "\n")

# 7. Mostrar si alguna de las tablas tiene valores faltantes.

colSums(is.na(blast))

cat("Blast tiene valores faltantes?", anyNA(blast), "\n")

cat("Metadatatiene valores faltantes?", anyNA(metadata), "\n")







    






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


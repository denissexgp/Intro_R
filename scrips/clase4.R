


# Integración de datos y automatización de análisis en R -------------------


library(dplyr)
library(readr)
library(purrr)
getwd()

carpeta_datos <- "data_raw/datos_sesion_04_integracion"
carpeta_blast <- file.path(carpeta_datos, "blast_muestras")

if (!dir.exists("results")) {
  dir.create("results")
}

blast_A <- read_tsv(
  file.path(carpeta_blast, "muestra_A.tsv"),
  show_col_types = FALSE
)

metadata <- read_csv(
  file.path(carpeta_datos, "metadata_secuencias.csv"),
  col_types = cols(codigo_muestra = col_character(), .default = col_guess())
)

head(blast_A)
head(metadata)


blast_A_con_metadata <- blast_A %>%
  left_join(metadata, by = "query_id")

blast_A_con_metadata %>%
  select(query_id, subject_id, pident, organismo, origen) %>%
  head()


blast_C <- read_tsv(
  file.path(carpeta_blast, "muestra_C.tsv"),
  show_col_types = FALSE
)

nrow(blast_C)



blast_C %>%
  inner_join(metadata, by = "query_id") %>%
  nrow()

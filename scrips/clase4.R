
rm(list = ls())

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


blast_C %>%
  anti_join(metadata, by = "query_id")


blast_C %>%
  full_join(metadata, by = "query_id") %>%
  select(query_id, subject_id, organismo) %>%
  arrange(query_id)



blast_B <- read_tsv(
  file.path(carpeta_blast, "muestra_B.tsv"),
  show_col_types = FALSE
)

blast_total <- bind_rows(
  muestra_A = blast_A,
  muestra_B = blast_B,
  muestra_C = blast_C,
  .id = "muestra"
)

blast_total %>%
  count(muestra)


#_______________________________________________________________________________

library(dplyr)
library(readr)
library(tidyr)



carpeta_rizosfera <- "data_raw/datos_sesion_04_integracion/datos_rizosfera"

conteos_riz <- read_tsv(
  file.path(carpeta_rizosfera, "conteos_taxonomia_rizosfera.tsv"),
  show_col_types = FALSE
)

metadata_riz <- read_csv(
  file.path(carpeta_rizosfera, "metadata_rizosfera.csv"),
  show_col_types = FALSE
)

dim(conteos_riz)

cuantos <- metadata_riz %>%
  count(host_name, sort = TRUE)


niveles <- c(
  "kingdom", "phylum", "class", "order",
  "family", "genus", "species"
)

conteos_limpios <- conteos_riz %>%
  mutate(
    across(
      all_of(niveles),
      ~ na_if(sub("^[a-z]__", "", .x), "")
    )
  )

conteos_limpios %>%
  select(taxon_id, all_of(niveles)) %>%
  head()


ids_muestras <- metadata_riz$sample_id

stopifnot(
  !anyDuplicated(ids_muestras),
  all(ids_muestras %in% names(conteos_limpios))
)

conteos_largos <- conteos_limpios %>%
  pivot_longer(
    cols = all_of(ids_muestras),
    names_to = "sample_id",
    values_to = "conteo"
  )

datos_riz <- conteos_largos %>%
  left_join(
    metadata_riz %>% select(sample_id, host_name, country, year),
    by = "sample_id",
    relationship = "many-to-one"
  )

glimpse(datos_riz)



#-------------------------------------------------------------------------

#sumar toso lo que hay de un mismo sample_id 



totales_muestra <- datos_riz %>%
  group_by(sample_id) %>%
  summarise(total_lecturas = sum(conteo), .groups = "drop")

stopifnot(all(totales_muestra$total_lecturas > 0))
head(totales_muestra)


write_tsv(totales_muestra, "results/totales_muestra.tsv")


abundancias_genero <- datos_riz %>%
  mutate(genus = if_else(is.na(genus), "Sin clasificar", genus)) %>%
  group_by(sample_id, host_name, genus) %>%
  summarise(conteo = sum(conteo), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(abundancia_relativa = conteo / sum(conteo)) %>%
  ungroup()

abundancias_genero %>%
  select(sample_id, host_name, genus, conteo, abundancia_relativa) %>%
  head(abundancias_genero, n=11)



generos_5 <- abundancias_genero %>%
  mutate(
    taxon_final = if_else(
      genus != "Sin clasificar" & abundancia_relativa < 0.05,
      "Otros",
      genus
    )
  ) %>%
  group_by(sample_id, host_name, taxon_final) %>%
  summarise(
    abundancia_relativa = sum(abundancia_relativa),
    .groups = "drop"
  )



# Cada muestra debe sumar 1 (salvo pequeñas diferencias numéricas).
comprobacion <- generos_5 %>%
  group_by(sample_id) %>%
  summarise(total = sum(abundancia_relativa), .groups = "drop")

stopifnot(all(abs(comprobacion$total - 1) < 1e-8))
head(comprobacion) 





preparar_abundancias <- function(tabla, nivel = "genus", umbral = 0.05) {
  niveles_validos <- c(
    "kingdom", "phylum", "class", "order",
    "family", "genus", "species"
  )
  
  stopifnot( 
    nivel %in% niveles_validos,
    is.numeric(umbral), length(umbral) == 1,
    !is.na(umbral), umbral >= 0, umbral <= 1
  )
  
  etiqueta_otros <- if (nivel == "phylum") {
    "Filos < 5%"
  } else if (nivel == "genus") {
    "Géneros < 5%"
  } else {
    "Taxones < 5%"
  }
  
  tabla %>%
    mutate(
      taxon = .data[[nivel]],
      taxon = if_else(is.na(taxon) | taxon == "",
                      "Sin clasificar", taxon)
    ) %>%
    group_by(sample_id, host_name, taxon) %>%
    summarise(conteo = sum(conteo), .groups = "drop") %>%
    group_by(sample_id) %>%
    mutate(abundancia_relativa = conteo / sum(conteo)) %>%
    ungroup() %>%
    mutate(
      taxon_final = if_else(
        abundancia_relativa < umbral,
        etiqueta_otros,
        taxon
      )
    ) %>%
    group_by(sample_id, host_name, taxon_final) %>%
    summarise(
      abundancia_relativa = sum(abundancia_relativa),
      .groups = "drop"
    )
}



generos_5 <- preparar_abundancias(datos_riz, "genus", 0.05)
filos_5 <- preparar_abundancias(datos_riz, "phylum", 0.05)
familias_1 <- preparar_abundancias(datos_riz, "family", 0.01)

filos_5 %>% head()



write_tsv(generos_5, "results/generos_5.tsv")
write_tsv(filos_5, "results/filos_5.tsv")
write_tsv(familias_1, "results/familias_1.tsv")

#Actividad crear comentarios 

preparar_abundancias <- function(tabla, nivel = "genus", umbral = 0.05) {
  niveles_validos <- c(
    "kingdom", "phylum", "class", "order",
    "family", "genus", "species"
  )
  
  # Validar los argumentos de entrada antes de procesar
  stopifnot( 
    nivel %in% niveles_validos,
    is.numeric(umbral), length(umbral) == 1,
    !is.na(umbral), umbral >= 0, umbral <= 1
  )
  
  # Asignar la etiqueta para los taxones que queden por debajo del 5%
  etiqueta_otros <- if (nivel == "phylum") { 
    "Filos < 5%"
  } else if (nivel == "genus") {
    "Géneros < 5%"
  } else {
    "Taxones < 5%"
  }
  
  # Seleccionar dinámicamente la columna del nivel elegido y etiquetar valores faltantes
  tabla %>%
    mutate(
      taxon = .data[[nivel]],
      taxon = if_else(is.na(taxon) | taxon == "",
                      "Sin clasificar", taxon)
    ) %>% #Sumar los conteos por muestra, hospedero y taxón
    group_by(sample_id, host_name, taxon) %>%
    summarise(conteo = sum(conteo), .groups = "drop") %>%
    group_by(sample_id) %>% 
    #Calcular la abundancia relativa respecto al total de lecturas de cada muestra
    mutate(abundancia_relativa = conteo / sum(conteo)) %>%
    ungroup() %>%
    #Reemplazar los taxones con abundancia menor al umbral por la etiqueta agrupadora
  
    mutate(
      taxon_final = if_else(
        abundancia_relativa < umbral,
        etiqueta_otros,
        taxon
      )
    ) %>%
    # Sumar las abundancias relativas de los taxones agrupados en la categoría residual
    group_by(sample_id, host_name, taxon_final) %>%
    summarise(
      abundancia_relativa = sum(abundancia_relativa),
      .groups = "drop"
    )
}


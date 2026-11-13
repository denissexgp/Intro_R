################################################################################
#Clase 3 Manipulación y transformación de datos bioinformáticos con dplyr y tidyr
rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(tibble)


ruta_blast <- "data_raw/datos_sesion_02_blast/resultados_blast.tsv"
ruta_conteos <- "data_raw/datos_sesion_03_dplyr/conteos_genes.csv"


stopifnot(file.exists(ruta_blast), file.exists(ruta_conteos))

blast <- read_tsv(ruta_blast, show_col_types = FALSE)
conteos <- read_csv(ruta_conteos, show_col_types = FALSE)

glimpse(blast)
glimpse(conteos)

class(blast)
class(as.data.frame(blast))


# El operador %>% es igual que una pipeline

# seleccionar solo algunas columnas 
blast %>%
  select(query_id, subject_id, pident)

# seleccionar todas menos una en especifico 
blast %>%
  select(-evalue)  # Todas las columnas excepto evalue

blast %>% select(starts_with("query"))



# filter(): elegir filas


blast %>% 
  filter(pident >= 90)

blast %>% 
  filter(pident >= 90, length >= 150)

blast %>% 
  filter(query_id == "seq01" | query_id == "seq02")

blast %>% 
  filter(query_id %in% c("seq01", "seq03", "seq05"))



blast %>% 
  arrange(bitscore)         # Menor a mayor

blast %>% 
  arrange(desc(bitscore))   # Mayor a menor

#primero ordena respecto a una y despues a la otra 
blast %>% 
  arrange(query_id, desc(bitscore))

# renombrar 

blast %>%
  rename(
    identidad_porcentaje = pident,
    longitud_alineamiento = length
  ) %>%
  select(query_id, identidad_porcentaje, longitud_alineamiento)

# distinct(): identificar combinaciones únicas --------------------------------

blast %>% distinct(query_id)

blast %>% distinct(query_id, subject_id) # se ve que es distinto de las dos columnas y te da los valores unicos 





#Crear y modificar columnas con mutate() --------------------------------------

blast %>%
  mutate(
    identidad_proporcion = pident / 100,
    longitud_kb = length / 1000
  ) %>%
  select(query_id, pident, identidad_proporcion, length, longitud_kb)


blast %>%
  mutate(
    supera_90 = pident >= 90,
    etiqueta = if_else(pident >= 90, "90 o más", "Menos de 90")
  ) %>%
  select(query_id, pident, supera_90, etiqueta)

# Más de dos categorías con case_when() -----------------------------------
blast %>%
  mutate(
    grupo_didactico = case_when(
      pident >= 95 ~ "Grupo A",
      pident >= 80 ~ "Grupo B",
      TRUE ~ "Grupo C"
    )
  ) %>%
  count(grupo_didactico)


blast %>% count(query_id, sort = TRUE)
  

# group_by() y summarise() ------------------------------------------------

resumen <- blast %>%
  group_by(query_id) %>%
  summarise(
    numero_alineamientos = n(),
    referencias_distintas = n_distinct(subject_id),
    bitscore_maximo = max(bitscore),
    identidad_promedio = mean(pident),
    .groups = "drop"
  )
resumen


blast %>%
  group_by(query_id) %>%
  mutate(alineamientos_por_consulta = n()) %>% # n() cuenta cuantos hay 
  ungroup() %>%
  select(query_id, subject_id, alineamientos_por_consulta)


# pivot_longer(): de ancho a largo ----------------------------------------

conteos_largos <- conteos %>%
  pivot_longer(
    cols = starts_with("muestra_"),
    names_to = "muestra",
    values_to = "conteo"
  )
conteos_largos


# pivot_wider(): de largo a ancho -------------------------------------------
conteos_recuperados <- conteos_largos %>%
  pivot_wider(
    names_from = muestra,
    values_from = conteo
  )
conteos_recuperados

# drop --------------------------------------------------------------------


library(dplyr)

datos <- tibble(
  especie = c("A", "A", "A", "A", "B", "B", "B", "B"),
  tratamiento = c("Control", "Control", "Tratado", "Tratado",
                  "Control", "Control", "Tratado", "Tratado"),
  valor = c(10, 20, 30, 40, 50, 60, 70, 80)
)

datos

resumen1 <- datos %>%
  group_by(especie, tratamiento) %>%
  summarise(
    promedio = mean(valor)
  )
resumen1


group_vars(resumen1) # dice en donde esta haciendo los grupos 

resumen1 %>%
  summarise(
    media = mean(promedio)
  )

# ahora con drop  eso desase los grupos, tiene el mismo efecto de ungrup()
resumen2 <- datos %>%
  group_by(especie, tratamiento) %>%
  summarise(
    promedio = mean(valor),
    .groups = "drop"
  )

resumen2

resumen2 %>%
  summarise(
    media = mean(promedio)
  )


# Guardar resultados sin sobrescribir los datos originales ----------------

# si no existe el directorio results, crea el directorio
if (!dir.exists("results")) dir.create("results")

# guardar  como archivo algo de una variable 
readr::write_csv(resumen, "results/03_resumen_blast.csv") 
readr::write_tsv(conteos_largos, "results/03_conteos_largos.tsv")


file.exists("results/03_resumen_blast.csv")



# Ejercicios  -------------------------------------------------------------

library(readr)
library(dplyr)

datos <- read_csv(
  "data_raw/expresion_genica_estres.csv",
  show_col_types = FALSE
)

glimpse(datos)


# Ejercicio1 --------------------------------------------------------------
#Queremos explorar qué genes podrían estar relacionados con la respuesta bacteriana al estrés.

#Utilizando dplyr:
  
#Calcula el conteo promedio de cada gen en cada condición experimental.
resumen11 <- datos %>%
  group_by(gene_id, condicion) %>%
  summarise(
    promedio = mean(conteo)
  )

resumen11

#Organiza los resultados de manera que puedas comparar los promedios de Control y Estrés para cada gen.

comparar <- resumen11 %>%
  group_by(gene_id) %>% 
  arrange(desc(promedio))

comparar


#Calcula la razón entre el promedio en Estrés y el promedio en Control.

razon_promedios <- datos %>%
  group_by(gene_id, condicion, funcion) %>%
  summarise(promedio = mean(conteo), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condicion, values_from = promedio) %>%
  mutate(razon_estres_control = Estres / Control)

print(razon_promedios)


#Identifica los cinco genes con mayor aumento relativo y muestra su anotación funcional.

top5_genes <- razon_promedios %>%
  slice_max(order_by = razon_estres_control, n = 5) 


print(top5_genes)


#Pregunta de interpretación: ¿Un gen con un aumento relativo grande necesariamente presenta un conteo alto de lecturas? Justifica tu respuesta.

















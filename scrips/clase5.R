

library(ggplot2)

# Resultados procesados en la sesión anterior
#generos_5 
#filos_5 
#totales_muestra 
#metadata_riz 


profundidad <- totales_muestra %>%
  left_join(
    metadata_riz %>%
      select(sample_id, host_name, country, crop_system),
    by = "sample_id",
    relationship = "one-to-one"
  )

stopifnot(!anyNA(profundidad$host_name))

glimpse(profundidad)

p_profundidad <- ggplot(profundidad, aes(x = sample_id, y = total_lecturas)) +
  geom_col(fill = "#367A90") +
  labs(title = "Conteos totales por muestra", x = "Muestra",
       y = "Lecturas contabilizadas") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
p_profundidad
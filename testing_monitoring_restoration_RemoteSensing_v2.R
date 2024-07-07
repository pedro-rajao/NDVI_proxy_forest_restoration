################# 
#### Analyse, Dataframes, Figures and Maps 
#### Author: Pedro Rajão
#### Data: 22/03/2024
#################
################# 
#### packages
#################

library(ggplot2)
library(sf)
library(dplyr)
library(terra)
library(stringr)
library(raster)
library(rgdal)
library(cowplot)
library(lattice)
library(tidyr)
library(sp)
library(stars)
library(leaps)
library(lme4)
library(ggplot2)
library(cluster)
library(vegan)


#########################
################# 
#### Mapas
##

shapefile_pontos <- st_read(file.choose()) #shape_pointsSIRGAS2000
shapefile_buffer <- st_buffer(shapefile_pontos, 40)

## approach ONE : multiple raster, mixed-linear model
import_rasters <- function(directory) {
  files <- list.files(directory, full.names = TRUE)
  raster_list <- list()
  for (file in files) {
    raster <- rast(file)
    raster_name <- tools::file_path_sans_ext(basename(file))
    raster_list[[raster_name]] <- raster
  }
  
  return(raster_list)
}

# Exemplo de uso da função
directory <- "C:/Users/pedro.rajao/Desktop/AGEVAP/13. NDVI proxy sucesso da restauração/tiff/2019"
rasters <- import_rasters(directory)
str(rasters)

# Lista para armazenar os valores extraídos
valores <- list()
for (nome_raster in names(rasters)) {
  raster <- rasters[[nome_raster]]
  valores[[nome_raster]] <- terra::extract(raster, shapefile_buffer)
}


dados_extraidos_A <- bind_rows(valores, .id = "raster_nome")

dados_extraidos_A$EVI <- 2.5 * (dados_extraidos_A$B8 - dados_extraidos_A$B4) / (dados_extraidos_A$B8 + 6 * dados_extraidos_A$B4 - 7.5 * dados_extraidos_A$B1 + 1)
dados_extraidos_A$NDWI <- (dados_extraidos_A$B3 - dados_extraidos_A$B8) / (dados_extraidos_A$B3 + dados_extraidos_A$B8)
dados_extraidos_A$SAVI <- (1.5 * (dados_extraidos_A$B8 - dados_extraidos_A$B4)) / (dados_extraidos_A$B8 + dados_extraidos_A$B4 + 0.5)
dados_extraidos_A$MSAVI <- (2 * dados_extraidos_A$B8 + 1 - sqrt((2 * dados_extraidos_A$B8 + 1) ^ 2 - 8 * (dados_extraidos_A$B8 - dados_extraidos_A$B4))) / 2
dados_extraidos_A$NDMI <- (dados_extraidos_A$B3 - dados_extraidos_A$B8) / (dados_extraidos_A$B3 + dados_extraidos_A$B8)
dados_extraidos_A$NDBI <- (dados_extraidos_A$B4 - dados_extraidos_A$B8) / (dados_extraidos_A$B4 + dados_extraidos_A$B8)
dados_extraidos_A$LAI <- 0.5 * log((1 - dados_extraidos_A$NDVI) / (0.69 * dados_extraidos_A$NDVI))
dados_extraidos_A$TCARI <- ((dados_extraidos_A$B7 - dados_extraidos_A$B6) - 0.2 * (dados_extraidos_A$B7 - dados_extraidos_A$B5) * (dados_extraidos_A$B7 / dados_extraidos_A$B6)) / ((dados_extraidos_A$B7 - dados_extraidos_A$B6) + 0.2 * (dados_extraidos_A$B7 - dados_extraidos_A$B5) * (dados_extraidos_A$B7 / dados_extraidos_A$B6))
dados_extraidos_A$OSAVI <- (1 + 0.16) * (dados_extraidos_A$B8 - dados_extraidos_A$B4) / (dados_extraidos_A$B8 + dados_extraidos_A$B4 + 0.16)
dados_extraidos_A$MTVI <- (1.5 * (2.5 * (dados_extraidos_A$B8 - dados_extraidos_A$B4) - 1.3 * (dados_extraidos_A$B8 - dados_extraidos_A$B3))) / sqrt(((2 * dados_extraidos_A$B8 + 1)^2 - (6 * dados_extraidos_A$B8 - 5 * sqrt(dados_extraidos_A$B4))))
NIR <- dados_extraidos_A$B8
RED <- dados_extraidos_A$B4
dados_extraidos_A$GEMI <- (2 * ((NIR^2 - RED^2) + 1.5 * NIR + 0.5 * RED)) / (NIR + RED + 0.5)


media_sd_por_ID <- dados_extraidos_A %>%
  group_by(raster_nome, ID) %>%
  summarise(across(B1:GEMI, list(media = mean, sd = sd), na.rm = TRUE))

dados_extraidos_com_atributos <- cbind(st_drop_geometry(shapefile_buffer), media_sd_por_ID, row.names = NULL)
dados_extraidos_df_mixed <- as.data.frame(dados_extraidos_com_atributos)


####################
### approach TWO : ONE raster, linear model
## Loading TIFFs
s2_2020 <- rast(file.choose()) #NDVI_median_2019sirgas2000_v1
shapefile_pontos <- st_read(file.choose()) #shape_pointsSIRGAS2000

dados_extraidos <- terra::extract(s2_2020, shapefile_pontos)
atributos_shape <- st_drop_geometry(shapefile_pontos)
dados_extraidos_com_atributos <- cbind(atributos_shape, dados_extraidos)
str(dados_extraidos_com_atributos)


shapefile_buffer <- st_buffer(shapefile_pontos, 10)

# Extrair os valores dos pixels do raster dentro dos buffers
dados_extraidos_buffer <- terra::extract(s2_2020, shapefile_buffer)
media_por_ID <- dados_extraidos_buffer %>%
  group_by(ID) %>%
  summarise_all(mean, na.rm = TRUE)
dados_extraidos_com_buffer_atributos <- cbind(atributos_shape, media_por_ID)

###################
## ### approach THREE : ONE shape with raster values extract in QGIS, linear model
shapefile_pontos_df_AA <- as.data.frame(readOGR(file.choose())) #shapefile_pontos_NDVI_2019
colnames(shapefile_pontos_df_AA)[7:21] <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12", "NDVI", "NDVI_1")


##################
##########
### trocar: 1. dados_extraidos_com_atributos 2. dados_extraidos_com_buffer_atributos  3. shapefile_pontos_df_AA  4. dados_extraidos_df_mixed

shapefile_pontos_df <- dados_extraidos_df_mixed
shapefile_pontos_df$copa <- as.numeric(gsub(",", ".", shapefile_pontos_df$copa))
str(shapefile_pontos_df)



## Gráficos exploratórios de regressão MIXED-EFFECTS
cores <- rainbow(length(unique(shapefile_pontos_df$raster_nome)))
cores_aleatorizadas <- sample(cores)
simbolos <- c(sample(1:30, 20, replace = TRUE))  
grafico <- ggplot(shapefile_pontos_df, aes(x = NDWI_media, y = copa, color = raster_nome, shape = raster_nome)) +
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +  
  scale_color_manual(values = cores_aleatorizadas) +  # Usar vetor de cores
  scale_shape_manual(values = simbolos) +  # Usar vetor de símbolos
  labs(x = "NDWI", y = "copa", color = "Raster", shape = "Raster") +  # Definir rótulos dos eixos e legendas
  theme_classic() +
  guides() + # Exibir legendas
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = TRUE, size = 0.7, fill = "lightgrey") +
  theme(legend.text = element_text(size = 6),  # Ajusta o tamanho da fonte da legenda
        legend.position = "bottom",  # Posição da legenda
        legend.spacing.y = unit(0.1, "cm"),  # Espaçamento vertical da legenda
        legend.title = element_blank()) +  # Remover título da legenda
  guides(color = guide_legend(ncol = 3)) #+ xlim(0.45, 0.8) 
  
print(grafico)
summary(lm(shapefile_pontos_df$copa~shapefile_pontos_df$NDWI_media, data=shapefile_pontos_df))
ggsave(filename='C:/Users/pedro.rajao/Desktop/AGEVAP/16. NDVI proxy Indicadores da restauração/figure/NDWI_media_copa_mix-effects.png', plot=grafico, dpi=300)


# média dos INDEX da vegetação entre as imagens. criando uma "imagem média" e definindo os indicadores da imagem média.
media_por_cod_pol_UA <- shapefile_pontos_df %>%
  group_by(cod_pol, cod_UA) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))



############
######
##

fito <- read.table(file.choose(), sep='\t', dec=',', header=TRUE)
df_mixed <- merge(fito, shapefile_pontos_df, by = 'cod_UA') # por raster_nome. para rodar relações de efeito-misto.
df <- merge(fito, media_por_cod_pol_UA, by = 'cod_UA') #  imagem média.
df <- df %>% filter(cod_UA != 61)

# agrupando tudo para polígono.
df_poligono <- df %>%
  group_by(cod_pol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


# Crie a matriz de espécies
matriz_especies <- table(df$cod_UA, df$especie)
#matriz_especies <- table(df_poligono$cod_UA, df_poligono$especie)


dist_euclidiana <- vegdist(matriz_especies, method = "euclidean")
nmds_resultado_euclidiana <- metaMDS(dist_euclidiana)
summary(nmds_resultado_euclidiana)
coordenadas_nmds <- scores(nmds_resultado_euclidiana, choices = c(1, 2))
df_nmds <- data.frame(cod_UA = rownames(coordenadas_nmds), coordenadas_nmds)

# Calculando a raiz quadrada das contagens e normalizando
sqrt_matriz <- sqrt(matriz_especies)
matriz_hellinger <- sqrt_matriz / sqrt(rowSums(sqrt_matriz^2))
dist_bray <- vegdist(matriz_hellinger, method = "bray")
nmds_resultado_bray <- metaMDS(dist_bray)
summary(nmds_resultado_bray)
coordenadas_nmds <- scores(nmds_resultado_bray, choices = c(1, 2))
df_nmds <- data.frame(cod_UA = rownames(coordenadas_nmds), coordenadas_nmds)
plot(df_nmds[,1], df_nmds[,2])
summary(matriz_hellinger)
num_grupos <- 3
kmeans_resultado <- kmeans(coordenadas_nmds, centers = num_grupos)
df_nmds$grupo <- factor(kmeans_resultado$cluster)
points(kmeans_resultado$centers[,1], kmeans_resultado$centers[,2], col = 1:num_grupos, pch = 8, cex = 2)
plot(df_nmds[,1], df_nmds[,2], col = df_nmds$grupo, pch = 19, main = "NMDS plot com grupos")
legend("bottomright", legend = levels(df_nmds$grupo), col = 1:num_grupos, pch = 19, title = "Grupos")

df_composicao <- merge(df_mixed, df_nmds, by = 'cod_UA')

###########
### POR UA

# no-mixed
novo_dataframe_media <- df_composicao %>%
  group_by(cod_UA, grupo) %>% #, raster_nome
  summarise(across(c(graminea:NMDS2), mean, na.rm = TRUE),
            altura_media = mean(altura, na.rm = TRUE),
            riqueza_especies = n_distinct(especie),
            diversidade_shannon = -sum((table(especie) / n()) * log(table(especie) / n())),
            densidade_individuos = n() / 100*1000,
            equidade_pielou = exp(diversidade_shannon / log(riqueza_especies)),
            percentual_Bio = mean(sd == "Bio") * 100)
novo_dataframe_media <- novo_dataframe_media %>%
  dplyr::select(-raster_nome)

# mixed
novo_dataframe <- df_composicao %>%
  group_by(cod_UA, raster_nome, grupo) %>% #
  summarise(across(c(graminea:NMDS2), mean, na.rm = TRUE),
            altura_media = mean(altura, na.rm = TRUE),
            riqueza_especies = n_distinct(especie),
            diversidade_shannon = -sum((table(especie) / n()) * log(table(especie) / n())),
            densidade_individuos = n() / 100*1000,
            equidade_pielou = exp(diversidade_shannon / log(riqueza_especies)),
            percentual_Bio = mean(sd == "Bio") * 100)

str(novo_dataframe)
str(novo_dataframe_media)

#no-mixed
novo_dataframe_modelo <- novo_dataframe_media[, c("B8_media", "B4_media", "TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media", "TCARI_sd", "OSAVI_sd", "MTVI_sd", "GEMI_sd", "NDVI_sd","EVI_sd", "NDWI_sd", "SAVI_sd", "MSAVI_sd", "NDMI_sd","NDBI_sd")]
#novo_dataframe_modelo <- novo_dataframe_media[, c("TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media" )]

# mixed
novo_dataframe_modelo2 <- novo_dataframe[, c("B1_media", "B2_media", "B3_media", "B4_media", "B5_media", "B6_media", "B7_media", "B8_media", "B8A_media", "B9_media", "B10_media", "B11_media", "B12_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media", "B1_sd", "B2_sd", "B3_sd", "B4_sd", "B5_sd", "B6_sd", "B7_sd","B8_sd", "B8A_sd", "B9_sd", "B10_sd", "B11_sd", "B12_sd", "NDVI_sd","EVI_sd", "NDWI_sd", "SAVI_sd", "MSAVI_sd", "NDMI_sd","NDBI_sd")]

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$altura_media ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B8_media + B4_media + TCARI_media + OSAVI_media + MTVI_media + GEMI_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media ), data = novo_dataframe_modelo)
initial_model2 <- lm(novo_dataframe$altura_media ~ 1, data = novo_dataframe_modelo2)
final_model2 <- step(initial_model2, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd), data = novo_dataframe_modelo2)
# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$riqueza_especies ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B8_media + B4_media + TCARI_media + OSAVI_media + MTVI_media + GEMI_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media), data = novo_dataframe_modelo)
initial_model2 <- lm(novo_dataframe$riqueza_especies ~ 1, data = novo_dataframe_modelo2)
final_model2 <- step(initial_model2, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd), data = novo_dataframe_modelo2)
# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$diversidade_shannon ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo)
initial_model2 <- lm(novo_dataframe_media$diversidade_shannon ~ 1, data = novo_dataframe_modelo2)
final_model2 <- step(initial_model2, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd), data = novo_dataframe_modelo2)

# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$densidade_individuos ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo)
# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$equidade_pielou ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo)
# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$copa ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media + B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo)
initial_model2 <- lm(novo_dataframe_media$copa ~ 1, data = novo_dataframe_modelo2)
final_model2 <- step(initial_model2, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo2)
# Sumário do modelo final
summary(final_model)

# Ajustar o modelo inicial
initial_model <- lm(novo_dataframe_media$NMDS1 ~ 1, data = novo_dataframe_modelo)
final_model <- step(initial_model, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd ), data = novo_dataframe_modelo)
initial_model2 <- lm(novo_dataframe_media$copa ~ 1, data = novo_dataframe_modelo2)
final_model2 <- step(initial_model2, direction = "forward", scope = formula(~ B1_media + B2_media + B3_media + B4_media + B5_media + B6_media + B7_media + B8_media + B8A_media + B9_media + B10_media + B11_media + B12_media + NDVI_media + EVI_media + NDWI_media + SAVI_media + MSAVI_media + NDMI_media + NDBI_media +  B1_sd + B2_sd + B3_sd + B4_sd + B5_sd + B6_sd + B7_sd + B8_sd + B8A_sd + B9_sd + B10_sd + B11_sd + B12_sd + NDVI_sd + EVI_sd + NDWI_sd + SAVI_sd + MSAVI_sd + NDMI_sd + NDBI_sd), data = novo_dataframe_modelo2)
# Sumário do modelo final
summary(final_model)


######
### gráficos

cores <- rainbow(length(unique(novo_dataframe$raster_nome)))
cores_aleatorizadas <- sample(cores, length(unique(novo_dataframe$raster_nome)))
simbolos <- c(sample(1:30, 20, replace = TRUE))  # Adicione mais símbolos conforme necessário
gerar_grafico <- function(data, y_var) {
    ggplot(data, aes(x = NDVI_media, y = !!sym(y_var), color = raster_nome, shape = raster_nome)) +
    geom_point(size = 2.5) + 
    geom_smooth(method = "lm", se = FALSE, size = 1.2) +  
    scale_color_manual(values = cores_aleatorizadas) +  # Usar vetor de cores
    scale_shape_manual(values = simbolos) +  # Usar vetor de símbolos
    labs(x = "NDVI", y = y_var, color = "Raster", shape = "Raster") +  # Definir rótulos dos eixos e legendas
    theme_classic() +
    theme(legend.text = element_text(size = 6),  # Ajusta o tamanho da fonte da legenda
          legend.position = "bottom",  # Posição da legenda
          legend.spacing.y = unit(0.1, "cm"),  # Espaçamento vertical da legenda
          legend.title = element_blank()) +  # Remover título da legenda
    guides(color = guide_legend(ncol = 3)) +
    geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = TRUE, fill = 'grey', size = 0.7) 
}
variaveis_dependentes <- c("grupo", "NMDS1", "NMDS2", "copa", "graminea", "altura_media", 
                           "riqueza_especies", "diversidade_shannon", "densidade_individuos", 
                           "equidade_pielou", "percentual_Bio")
for (var in variaveis_dependentes) {
  print(gerar_grafico(novo_dataframe, var))
}

###########
### POR POLIGONO

novo_dataframe2 <- df_poligono %>%
  group_by(cod_pol, raster_nome) %>%
  summarise(across(c(graminea :NMDS2), mean, na.rm = TRUE),
            altura_media = mean(altura, na.rm = TRUE),
            riqueza_especies = n_distinct(especie),
            densidade_individuos = n()/100*10000,
            diversidade_shannon = -sum((table(especie) / n()) * log(table(especie) / n())),
            equidade_pielou = exp(diversidade_shannon / log(riqueza_especies)),
            percentual_Bio = mean(sd == "Bio") * 100)

novo_dataframe2_modelo <- novo_dataframe_media[, c("B1_media", "B2_media", "B3_media", "B4_media", "B5_media", "B6_media", "B7_media", 
                                            "B8_media", "B8A_media", "B9_media", "B10_media", "B11_media", "B12_media", "EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media", "NDBI_media", "NDVI_media", "NDVI_sd", "LAI_media")]


regressao_altura <- lm(novo_dataframe2$altura_media ~ . , data = novo_dataframe2_modelo)
regressao_riqueza <- lm(novo_dataframe2$riqueza_especies ~ . , data = novo_dataframe2_modelo)
regressao_diversidade <- lm(novo_dataframe2$diversidade_shannon ~ . , data = novo_dataframe2_modelo)
regressao_equidade <- lm(novo_dataframe2$equidade_pielou ~ . , data = novo_dataframe2_modelo)
regressao_copa <- lm(novo_dataframe2$copa ~ . , data = novo_dataframe2_modelo)

summary(regressao_altura)
summary(regressao_riqueza)
summary(regressao_diversidade)
summary(regressao_equidade)
summary(regressao_copa)

# Função para gerar modelos lineares e gráficos de dispersão com linha de regressão
# Criar vetor de cores baseado nos atributos únicos da coluna raster_nome
cores <- rainbow(length(unique(novo_dataframe2$raster_nome)))
cores_aleatorizadas <- sample(cores, length(unique(novo_dataframe2$raster_nome)))
simbolos <- c(sample(1:30, 20, replace = TRUE))  # Adicione mais símbolos conforme necessário
gerar_grafico <- function(data, y_var) {
  # Ajuste do modelo linear
  modelo <- lm(paste(y_var, "~ NDVI_media"), data = data)
  # Gerando o gráfico
  ggplot(data, aes(x = NDVI_media, y = !!sym(y_var), color = raster_nome, shape = raster_nome)) +
    geom_point(size = 2.5) + 
    geom_smooth(method = "lm", se = FALSE, size = 1.2) +  
    scale_color_manual(values = cores_aleatorizadas) +  # Usar vetor de cores
    scale_shape_manual(values = simbolos) +  # Usar vetor de símbolos
    labs(x = "NDVI", y = y_var, color = "Raster", shape = "Raster") +  # Definir rótulos dos eixos e legendas
    theme_classic() +
    guides() + # Exibir legendas
    theme(legend.text = element_text(size = 6),  # Ajusta o tamanho da fonte da legenda
          legend.position = "bottom",  # Posição da legenda
          legend.spacing.y = unit(0.1, "cm"),  # Espaçamento vertical da legenda
          legend.title = element_blank()) +  # Remover título da legenda
    guides(color = guide_legend(ncol = 3))  # Definir número de colunas na legenda de cores
}
variaveis_dependentes <- c("copa", "graminea", "altura_media" ,"riqueza_especies","diversidade_shannon","densidade_individuos", "equidade_pielou", "percentual_Bio")
for (var in variaveis_dependentes) {
  print(gerar_grafico(novo_dataframe2, var))
}
#####
## salvar graficos


# Lista para armazenar os resultados das regressões
resultados <- list()

# Função para extrair informações da regressão
extrair_informacoes_regressao <- function(regressao) {
  print(summary(regressao))  # Verificar o resumo da regressão
  print(summary(regressao)$adj.r.squared)  # Verificar R²
  print(summary(regressao)$coefficients[2, 4])  # Verificar valor p
  list(
    equacao = as.character(summary(regressao)$call),
    r_squared = summary(regressao)$adj.r.squared,
    p_valor = summary(regressao)$coefficients[2, 4]
  )
}

# Regressões lineares
regressoes <- list(
  regressao_altura = regressao_altura,
  regressao_riqueza = regressao_riqueza,
  regressao_diversidade = regressao_diversidade,
  regressao_equidade = regressao_equidade,
  regressao_percentual = regressao_percentual,
  regressao_copa = regressao_copa
)

# Extrair informações das regressões
for (nome in names(regressoes)) {
  print(nome)  # Verificar o nome da regressão
  resultados[[nome]] <- extrair_informacoes_regressao(regressoes[[nome]])
}

# Converter a lista de resultados em um dataframe
tabela_resultados <- as.data.frame(do.call(rbind, resultados))

summary(lm(novo_dataframe_media$densidade_individuos~novo_dataframe_media$NDVI_media))

classificar_INEA <- function(copa, densidade, ind_zoocoricos, cobertura_copa, equidade, riqueza, altura_media, infestacao_gramineas) {
  # Inicializando a variável para armazenar a soma das notas
  soma_notas <- 0
  
  # Critério 1: Densidade (n° ind./ha)
  if (densidade < 1111) {
    soma_notas <- soma_notas + 0
  } else if (densidade >= 1111 & densidade < 1250) {
    soma_notas <- soma_notas + 0.65
  } else {
    soma_notas <- soma_notas + 1
  }
  
  # Critério 2: Ind. Zoocóricos (%)
  if (ind_zoocoricos < 40) {
    soma_notas <- soma_notas + 0
  } else if (ind_zoocoricos >= 40 & ind_zoocoricos < 60) {
    soma_notas <- soma_notas + 0.65
  } else {
    soma_notas <- soma_notas + 1
  }
  
  # Critério 3: Cobertura de copa (%)
  if (!is.na(cobertura_copa)) {
    if (cobertura_copa < 50) {
      soma_notas <- soma_notas + 0
    } else if (cobertura_copa >= 50 & cobertura_copa < 70) {
      soma_notas <- soma_notas + 0.65
    } else {
      soma_notas <- soma_notas + 1
    }
  }
  
  # Critério 4: Equidade J'
  if (!is.na(equidade)) {
    if (equidade < 0.6) {
      soma_notas <- soma_notas + 0
    } else if (equidade >= 0.6 & equidade < 0.8) {
      soma_notas <- soma_notas + 0.65
    } else {
      soma_notas <- soma_notas + 1
    }
  }
  
  # Critério 5: Riqueza S’
  if (riqueza < 15) {
    soma_notas <- soma_notas + 0
  } else if (riqueza >= 15 & riqueza < 25) {
    soma_notas <- soma_notas + 0.65
  } else {
    soma_notas <- soma_notas + 1
  }
  
  # Critério 6: Altura média (m)
  if (altura_media < 2) {
    soma_notas <- soma_notas + 0
  } else if (altura_media >= 2 & altura_media < 3) {
    soma_notas <- soma_notas + 0.65
  } else {
    soma_notas <- soma_notas + 1
  }
  
  # Critério 7: Infestação de gramíneas (%)
  if (!is.na(infestacao_gramineas)) {
    if (infestacao_gramineas > 30) {
      soma_notas <- soma_notas + 0
    } else if (infestacao_gramineas <= 30 & infestacao_gramineas > 20) {
      soma_notas <- soma_notas + 0.65
    } else {
      soma_notas <- soma_notas + 1
    }
  } 
  
  # Retornando a classificação final
  return(soma_notas)
}



novo_dataframe_media$classificacao <- apply(novo_dataframe_media[, c("copa", "densidade_individuos", "percentual_Bio", "altura_media", "equidade_pielou", "riqueza_especies", "graminea")], 
                                     1, 
                                     function(x) classificar_INEA(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]))

# Visualizando o novo dataframe com a coluna de classificação
str(novo_dataframe_media)
cores <- rainbow(length(unique(novo_dataframe2$raster_nome)))
cores_aleatorizadas <- sample(cores, length(unique(novo_dataframe2$raster_nome)))
simbolos <- c(sample(1:10, 10, replace = TRUE))  # Adicione mais símbolos conforme necessário
ggplot(novo_dataframe2, aes(x = NDVI_sd, y = classificacao, color = raster_nome, shape = raster_nome)) + #, shape = raster_nome
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +  
  scale_color_manual(values = cores_aleatorizadas) +  # Usar vetor de cores
  scale_shape_manual(values = simbolos) +  # Usar vetor de símbolos
  labs(x = "NDVIsd", y = "inea", color = "Raster", shape = "Raster") +  # Definir rótulos dos eixos e legendas
  theme_classic() +
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = FALSE, size = 0.7) +  
  guides(color = guide_legend(title = "Raster"), shape = guide_legend(title = "Raster"))  # Exibir legendas

summary(lm(novo_dataframe2$classificacao ~ ., novo_dataframe2_modelo))

regressao_class <- glm(classificacao ~ NDVI_media , data = novo_dataframe2, family = "binomial")
summary(regressao_class)
regressao_class2 <- lm(novo_dataframe2$classificacao ~ B1_media + B3_media + B5_media + B9_media + NDVI_sd, data = novo_dataframe2_modelo)
summary(regressao_class2)
residuos<-resid(regressao_class2)
plot(residuos)
shapiro.test(residuos)






grupo1 <- c("diversidade_shannon", "copa", "densidade_individuos", "percentual_Bio", "altura_media", "equidade_pielou", "riqueza_especies", "graminea")
grupo2 <- novo_dataframe_media[,c("B8_media", "B4_media", "TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media", "EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media", "NDBI_media") ]
correlacoes <- cor(na.omit(grupo2))
# Plotar o mapa de calor da matriz de correlação
heatmap(correlacoes, 
        col = colorRampPalette(c("blue", "white", "red"))(100),
        symm = TRUE,
        main = "Mapa de Calor da Matriz de Correlação")

# Matriz para armazenar as correlações
correlacoes <- matrix(NA, nrow = length(grupo2), ncol = length(grupo1), dimnames = list(grupo2, grupo1))

# Calcular correlações de Pearson
for (var_grupo2 in grupo2) {     for (var_grupo1 in grupo1) {
  correlacao <- cor.test(novo_dataframe_media[[var_grupo1]], novo_dataframe_media[[var_grupo2]])
  if (correlacao$p.value < 0.05) {
    correlacoes[var_grupo2, var_grupo1] <- correlacao$estimate
  } else {
    correlacoes[var_grupo2, var_grupo1] <- NA
  }
}
}

# Visualizar as correlações
print(correlacoes)





library(randomForest)

# Definir as colunas de atributos de entrada e o atributo de destino
atributos_entrada <- novo_dataframe_media[, c("B8_media", "B4_media", "TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media")]  # Atributo de entrada
atributo_destino <- novo_dataframe_media$classificacao    # Atributo de destino

# Dividir os dados em treino e teste
set.seed(123)  # Define uma semente para reprodução dos resultados
indice_treino <- sample(1:nrow(novo_dataframe_media), 0.7 * nrow(novo_dataframe_media))  # 70% dos dados para treino
dados_treino <- novo_dataframe_media[indice_treino, ]
dados_teste <- novo_dataframe_media[-indice_treino, ]

# Selecionar apenas as colunas relevantes para o modelo
dados_treino <- dados_treino[, c("B8_media", "B4_media", "TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media","classificacao","classificacao", "copa", "riqueza_especies", "altura_media", "densidade_individuos")]
dados_teste <- dados_teste[, c("B8_media", "B4_media", "TCARI_media", "OSAVI_media", "MTVI_media", "GEMI_media", "NDVI_media","EVI_media", "NDWI_media", "SAVI_media", "MSAVI_media", "NDMI_media","NDBI_media","classificacao", "copa", "riqueza_especies", "altura_media", "densidade_individuos")]


# Treinar o modelo de RandomForest
modelo_rf <- randomForest(copa ~ B4_media + B8_media +  NDVI_media + EVI_media + SAVI_media, data = dados_treino)
previsoes <- predict(modelo_rf, newdata = dados_teste)
mae <- mean(abs(previsoes - dados_teste$copa))
print(paste("Erro médio absoluto:", mae))

# Para calcular o erro médio quadrático:
mse <- mean((previsoes - dados_teste$copa)^2)
print(paste("Erro médio quadrático:", mse))

# Para calcular o coeficiente de determinação (R²):
r_squared <- cor(previsoes, dados_teste$copa)^2
print(paste("R²:", r_squared))

# Criar dataframe com previsões e valores observados
df_plot <- data.frame(Previsto = previsoes, Observado = dados_teste$copa)
# Calcular o intervalo de confiança para os pontos
df_plot <- df_plot %>%
mutate(IC_sup = Previsto + 1.1 * sd(Previsto), IC_inf = Previsto - 1.1 * sd(Previsto), Dentro_IC = ifelse(Observado >= IC_inf & Observado <= IC_sup, "Dentro", "Fora"))
ggplot(df_plot, aes(x = Previsto, y = Observado, color = Dentro_IC)) +
  geom_point(aes(color = Dentro_IC)) +
  geom_smooth(method = NULL, color = "lightgrey", fill = "grey", alpha = 0.3) +
  geom_smooth(method = 'lm', se = FALSE, color = "black", linetype = 'dashed', alpha = 0.3) +
  scale_color_manual(values = c("blue", "red")) +  # Definir cores para pontos dentro e fora do IC
  labs(x = "Classe Prevista", y = "Classe Observada", title = "Modelo SVM: Modelado x Observado") +
  theme_classic()



# Treinar o modelo SVM
library(e1071)
modelo_svm <- svm(copa ~B4_media+B8_media+NDVI_media+EVI_media+SAVI_media, data = dados_treino)
previsoes_svm <- predict(modelo_svm, newdata = dados_teste)
mae <- mean(abs(previsoes_svm - dados_teste$copa))
print(paste("Erro médio absoluto:", mae))
mse <- mean((previsoes_svm - dados_teste$copa)^2)
print(paste("Erro médio quadrático:", mse))
r_squared <- cor(previsoes_svm, dados_teste$copa)^2
print(paste("R²:", r_squared))

df_plot <- data.frame(Previsto = previsoes_svm, Observado = dados_teste$copa)

ggplot(df_plot, aes(x = Previsto, y = Observado)) +
  geom_point() +
  geom_smooth(method = NULL, color = "lightgrey", fill = "grey", alpha = 0.3) +  # Excluir linha de tendência
  geom_smooth(method = 'lm', se = FALSE, color = "black", linetype = 'dashed', alpha = 0.3) +  # Desligar o 'se'
  labs(x = "Classe Prevista", y = "Classe Observada", title = "Modelo SVM: Modelado x Observado") +
  theme_classic()


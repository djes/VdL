#******************************************************************************
#
# Expérimentation de Classification des arbres à partir des données LiDAR
#
# Export du modèle de hauteur de canopée (CHM) et délimitation des houppiers 
# avec métriques
#
# (c) Jésahel BENOIST
#     Service des Risques Urbains et Sanitaires 
#     Ville de Lille - 2022
#
# Basé sur le package R "lidR" par Jean-Romain Roussel 
# https://github.com/r-lidar/lidR
#
# Testé avec les données LiDAR de la Métropole Européenne de Lille
#
# Merci à J.R. Roussel, Jean-Roc Morréale et à Serge Hombert
#
# Installer R 
#
# Eventuellement installer RTools pour compiler les dernières versions des packages
#  à partir des sources (voir le site de RTools)
#
# Installer et charger les packages lidR, RCSF, sf, mapview, remotes, raster, 
# rasterVis, spatialEco, RSQLite, EBImage, progress, terra 
# -> install.packages("lidR", "RSCF", "sf", "mapview", "remotes", ...)
#
# Eventuellement installer lidRplugins, EBImage (https://github.com/aoles/EBImage)
#  -> 'remotes::install_github("Jean-Romain/lidRplugins")
#
# Quelques fonctions utiles...
# las_check(las)
# summary(las)
# rm(obj)
# str(obj)
# las@data : pour connaître les attributs
# las@data$treeID
# mesure du temps de traitement avec system.time()
#
#******************************************************************************

#******************************************************************************

# Chargement des packages

#******************************************************************************

library(lidR)
library(RCSF)
library(sf)
library(mapview)
library(remotes)
library(lidRplugins)
library(raster)
library(spatialEco)
library(rasterVis)
library(progress)
library(RSQLite)
library(terra)
#library(EBImage)

paths <- ""

#******************************************************************************

# Optionnel : traitement multi-processeurs

#******************************************************************************

library(future)

# 2 coeurs
plan(multisession, workers = 2L)

#set_lidr_threads(2L)

#******************************************************************************

# Chargement des tuiles et retuilage

#******************************************************************************

#--- Dossier où sont stockées les tuiles d'origine

ctg <- readLAScatalog("tuiles/")

#--- En cas de plantage, permet de redémarrer à la nème tuile (chunk)

#opt_restart(ctg) <- 1

#--- Taille des futures tuiles (côté en m.)

opt_chunk_size(ctg) <- 200

#--- Pas de bordure

opt_chunk_buffer(ctg) <- 0

#--- Active la compression du nuage résultant

opt_laz_compression(ctg) <- TRUE

#--- Supprime les bordures (buffer) existantes dans le nuage original

opt_filter(ctg) <- "-drop_withheld"

#--- Dossier de sortie et syntaxe de nommage des fichiers

opt_output_files(ctg) <-  paste0("1_retile200/", "/tuile_{XLEFT}_{YBOTTOM}")
ctg_retile<- catalog_retile(ctg)
lidR:::catalog_laxindex(ctg_retile)

#******************************************************************************

# Premier traitement, classification du sol
# Algorithme CSF (voir note)

#******************************************************************************

#--- Recharge le catalogue si nécessaire

ctg_retile <- readLAScatalog("1_retile200/")

opt_laz_compression(ctg_retile) <- TRUE
#opt_restart(ctg_retile) <- 9
opt_output_files(ctg_retile) <- paste0("2_classified/", "{*}_classified")

#--- pré-filtrage des points en dessous ou au dessus des données utiles

# pour éliminer le plus gros des points hors norme. A adapter en fonction
# du contexte local

opt_filter(ctg_retile) <- "-drop_z_below 15 -drop_z_above 70"

ctg_classified <- classify_ground(ctg_retile, algorithm = csf())
lidR:::catalog_laxindex(ctg_classified)

#******************************************************************************

# Exemple de filtrage d'objets vectoriels divers avec Z réel(non normalisé)

#******************************************************************************

ofilter <- function(chunk)
{
	#--- Bâtiments avec une marge Z de 2 m.
	
	nlas <- merge_spatial(chunk, bdtopo_batiment, "altitude_maximale_toit")
	nlas <- merge_spatial(nlas, bdtopo_batiment, "altitude_minimale_sol")
	nlas <- merge_spatial(nlas, bdtopo_batiment, "hauteur")
	nlas <- filter_poi(nlas, is.na(altitude_minimale_sol) | (Z>(altitude_maximale_toit+2))|(Z>(altitude_minimale_sol+hauteur+2)))

	#--- Pylones électriques
	
	nlas <- merge_spatial(nlas, pylones, "in_pylone")
	nlas <- filter_poi(nlas, in_pylone == FALSE)

	return(nlas)
}

#--- Chargements de couches dans des fichiers geopackage

pylones <- st_read("pylones-tampon8m.gpkg", layer="pylones-tampon8m")
bdtopo_batiment <- st_read("bdtopo_batiment_2D_tampon2m.gpkg", layer="bdtopo_batiment_2D_tampon2m")

ctg_filterb <- readLAScatalog("2_classified/")
#opt_restart(ctg_filterb) <- 96
opt_output_files(ctg_filterb) <-  paste0("3_filtero/", "/{*}_filtero")
opt_laz_compression(ctg_filterb) <- TRUE
options <- list(need_output_file = TRUE, automerge = TRUE, overwrite = TRUE)
ctg_filtero <- catalog_map(ctg_filterb, ofilter, .options = options)
lidR:::catalog_laxindex(ctg_filtero)

#******************************************************************************

# Différents filtrages

#******************************************************************************

myfirstfilter <- function(chunk)
{
	#--- Quitte si le nuage est vide
	
	if (is.empty(chunk)) return(NULL)

	#--- Filtre les données solides et les duplicats
	
	# Pour un maximum de végétation
	nlas <- filter_poi(chunk, NumberOfReturns > 1)
	
	#nlas <- filter_poi(chunk, NumberOfReturns > 2)
	
	if (is.empty(nlas)) return(NULL)
	nlas <- filter_duplicates(nlas)
	if (is.empty(nlas)) return(NULL)

	#--- Enlève les retours intermédiaires (non utilisé ici)
	
	nlas <- filter_firstlast(nlas)
	if (is.empty(nlas)) return(NULL)
	#check if all pulses have only two returns i.e. first and last
	#each gps time should have only two (first and last) returns 
	no <- 0.0
	nlas <- add_lasattribute(nlas, no, "no", "no")
	nlas@data[, no := .N , by = c("gpstime")]
	#keep only those points which are the first and last returns 
	#nlas <- filter_poi(nlas, no==2)
	#if (is.empty(nlas)) return(NULL)
	#retain all the intermediate returns for the 'good' gpstimes 
	#nlas <- filter_poi(nlas, gpstime %in% nlas@data$gpstime)
	#if (is.empty(nlas)) return(NULL)

	nlas <- retrieve_pulses(nlas)

	#--- Calcule le NDVI
	
	ndvi <- 0.0
	nlas <- add_lasattribute(nlas, ndvi, "NDVI", "NDVI")
	nlas@data[, NDVI := (NIR-R)/(NIR+R), by = pulseID]

	#--- Filtre les données en fonction du NDVI (garde le végétal)
	
	#nlas <- filter_poi(nlas, NDVI > 0.25)
	#if (is.empty(nlas)) return(NULL)

	#--- Calcule l'EVI (non utilisé ici)
	
	#evi <- 0.0
	#C1 <- 6
	#C2 <- 7.5
	#L <- 1
	#nlas <- add_lasattribute(nlas, evi, "EVI", "EVI")
	#nlas@data[, EVI := G*((NIR-R)/(NIR+C1*R-C2*B+L)), by = pulseID]

	#--- Calcule le NEVI
	
	nevi <- 0.0
	nlas <- add_lasattribute(nlas, nevi, "NEVI", "NEVI")
	nlas@data[, NEVI := G*((NIR-R)/(NIR+B)), by = pulseID]

	#--- Filtre les données en fonction de cet indice (garde le végétal)
	
	nlas <- filter_poi(nlas, NEVI > 5000)
	if (is.empty(nlas)) return(NULL)

	#--- Filtre certains objets géométriques
	
	#nlas <- segment_shapes(nlas, shp_hplane(k = 100), "CoplanarH")
	#nlas <- segment_shapes(nlas, shp_plane(k = 100), "Coplanar")
	#nlas <- segment_shapes(nlas, shp_hline(k = 100), "ColinearH")
	#nlas <- segment_shapes(nlas, shp_vline(k = 100), "ColinearV")
	nlas <- segment_shapes(nlas, shp_line(th1 = 5, k = 100), "Colinear")
	#nlas <- segment_shapes(nlas, shp_line(th1 = 1, k = 100), "Colinear")
	nlas <- filter_poi(nlas, Colinear == FALSE)
	if (is.empty(nlas)) return(NULL)

	#--- Lignes électriques
	
	st_crs(nlas) <- 2154
	nlas <- merge_spatial(nlas, ligne_elec, "in_elec")
	nlas <- filter_poi(nlas, (in_elec == FALSE) | (Z < 10))
	if (is.empty(nlas)) return(NULL)

	#--- Filtre les points et objets isolés
	
	nlas <- classify_noise(nlas , sor(50, 10))
	nlas <- filter_poi(nlas, Classification != LASNOISE)
	if (is.empty(nlas)) return(NULL)
	nlas <- classify_noise(nlas , sor(10, 5))
	nlas <- filter_poi(nlas, Classification != LASNOISE)
	if (is.empty(nlas)) return(NULL)

	#--- Filtrages plus forts
	
	#nlas <- classify_noise(nlas , sor(50, 3))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)

	# 2000 points dans les 27 voxels (5 m de large)
	#nlas <- classify_noise(nlas , ivf(5, 2000))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	# 500 points in the 27 voxels (3 m wide)
	#nlas <- classify_noise(nlas , ivf(3, 500))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	# 50 points in the 27 voxels (1.5 m wide)
	#nlas <- classify_noise(nlas , ivf(1.5, 50))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	# 20 points in the 27 voxels (1 m wide)
	#nlas <- classify_noise(nlas , ivf(1, 20))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)

	nlas <- filter_poi(nlas, Classification != LASUNCLASSIFIED)

	#--- Si le nuage est vide
	
	if (is.empty(nlas)) return(NULL)

	return(nlas)
}

ligne_elec <- st_read("lignes-elec-tampon-8m.gpkg", layer="lignes-elec-tampon-8m2")
ctg_classified <- readLAScatalog("3_filtero/")
#opt_restart(ctg_classified) <- 20
opt_output_files(ctg_classified) <-  paste0("4_filter/", "/{*}_f1")
opt_chunk_buffer(ctg_classified) <- 0
opt_laz_compression(ctg_classified) <- TRUE
options <- list(need_output_file = TRUE, overwrite = TRUE)
ctg_filter <- catalog_map(ctg_classified, myfirstfilter, .options = options)
lidR:::catalog_laxindex(ctg_filter)

#******************************************************************************

# Normalisation du nuage
# Algorithme kniddw (voir note)

#******************************************************************************

ctg_filter<- readLAScatalog("4_filter/")
#opt_restart(ctg_filter1) <- 174
opt_output_files(ctg_filter) <-  paste0("5_norm/", "/{*}_norm")
opt_laz_compression(ctg_filter) <- TRUE
ctg_norm <- normalize_height(ctg_filter, knnidw())
lidR:::catalog_laxindex(ctg_norm)

#******************************************************************************

# Filtrage après la normalisation

#******************************************************************************

bfilter <- function(chunk)
{
	#--- Enlève les points au sol
	
	nlas <- filter_poi(chunk, Classification == 0L)
	if (is.empty(nlas)) return(NULL)

	#--- Filtrage des bâtiments avec le MNTB disponible
	# (tout ce qui est supérieur à 2.50 est supprimé)
	
	nlas <- merge_spatial(nlas, bati , "zbati")
	nlas <- filter_poi(nlas, zbati < 2.50)
	if (is.empty(nlas)) return(NULL)

	#--- Petit filtrage du bruit complémentaire
	
	nlas <- classify_noise(nlas , sor(50, 10))
	nlas <- filter_poi(nlas, Classification != LASNOISE)
	nlas <- classify_noise(nlas , sor(10, 5))
	nlas <- filter_poi(nlas, Classification != LASNOISE)

	return(nlas)
}
bati <- rast("LHL_2018_lidar_mel_mntb_norm.tif")
ctg_filter2 <- readLAScatalog("5_norm/")
#opt_restart(ctg_filter2) <- 96

#--- Après un premier test, aucun arbre n'a été trouvé plus haut que 42 m.
# -> supprime donc les points largement au-dessus 
# supprime également les points en dessous de 2 m. (essentiellement pour enlever véhicules et piétons)
# Cette méthode est plus rapide qu'un filtre

opt_filter(ctg_filter2) = "-drop_z_below 2.0 -drop_z_above 42"

opt_output_files(ctg_filter2) <-  paste0("6_filterb/2.0-42", "/{*}_f2")
opt_laz_compression(ctg_filter2) <- TRUE
opt_chunk_buffer(ctg_filter2) <- 0
options <- list(need_output_file = TRUE, overwrite = TRUE)
ctg_filterb <- catalog_map(ctg_filter2, bfilter, .options = options)
lidR:::catalog_laxindex(ctg_filterb)

#******************************************************************************

# Calcul du modèle de hauteur de canopée (CHM)
# Algorithme pitfree
# A noter que le fichier rasterize_canopy.vrt résultant sur Windows est buggé, 
# -> il faut l'éditer et enlever les "/" en début de chemin

#******************************************************************************

ctg_filtero <- readLAScatalog("6_filterb/2.0-42")
#opt_restart(ctg_filtero) <- 9
opt_output_files(ctg_filtero) <-  paste0("7_chm/2.0-42/0.25/", "/{*}_chm")
opt_laz_compression(ctg_filtero) <- TRUE
opt_chunk_buffer(ctg_filtero) <- 20
options <- list(need_output_file = TRUE, automerge = TRUE)

#--- option non documentée pour écraser les fichiers

ctg_filtero@output_options$drivers$Raster$param$overwrite <- TRUE

#--- Résolution 0.25 m.

chm <- rasterize_canopy(ctg_filtero, 0.25, pitfree(c(0, 2, 5, 10, 15), c(3,1), subcircle = 0.15))

#--- Autres exemples

# Résolution 0.50 m.
#chm <- grid_canopy(ctg_filtero, 0.50, pitfree(c(0, 2, 5, 10, 15), c(3, 1), subcircle = 0.2))
#chm <- rasterize_canopy(ctg_filtero, res = 0.50, algorithm = p2r(subcircle = 0.15), pkg = "terra")

#******************************************************************************

# Lissage du CHM (optionnel mais recommandé)

#******************************************************************************

ker = matrix(1, 5, 5)
r <- raster("7_chm/2.0-42/0.25/rasterize_canopy.vrt")
#nchm <- focal(r, w = ker, fun = median, na.rm = TRUE, filename = "7_chm/chm_focal.gid", overwrite = TRUE)
nchm <- focal(r, w = ker, fun = median, na.rm = TRUE, filename = "7_chm/chm_focal.tif", overwrite = TRUE)

#--- Variante
# Attention, celle-ci ne prend pas en charge les buffers, il peut en 
# résulter des artefacts au niveau des bordures des tuiles

#f = list.files("7_chm/2.0-42/0.5/", pattern = "*.tif")
f = list.files("7_chm/2.0-42/0.25/", pattern = "*.tif")

print(paste0("Nb: ", length(f)))

for( i in 1:length(f)){
	ker = matrix(1, 5, 5)
  	r <- rast(paste0("7_chm/2.0-42/0.25/", f[i]))
	nchm <- focal(r, w = ker, fun = median, na.rm = TRUE)
	writeRaster(nchm, paste0("7_chm/2.0-42/0.25/focal5x5/", f[i]))
	print(i)
	flush.console()
}

#******************************************************************************

# Rééchantillonnage du CHM (si besoin)

#******************************************************************************

f = list.files("7_chm/0.25/focal5x5/", pattern = "*.tif")

print(paste0("Nb: ", length(f)))

for( i in 1:length(f)){
  	r <- rast(paste0("7_chm/0.25/focal5x5/", f[i]))
	r_resampled <- r
	res(r_resampled) <- res(r)*2
	r_resampled <- resample(r, 
                                r_resampled, 
                                method = "cubic")
	writeRaster(r_resampled, paste0("7_chm/0.5/focal5x5/", f[i]))
	print(i)
	flush.console() 
}

#******************************************************************************

# Si le calcul du CHM a planté, il faut recréer et convertir un VRT avec un outil quelconque 
# (QGIS par exemple) puis recréer un GéoTIFF

#******************************************************************************

# chm <- raster("7_chm/rasterize_canopy.vrt")
# writeRaster(nchm, "7_chm/chm_focal.tif")

#******************************************************************************

# Localise les arbres en se basant sur le nuage de points (optionnel)
# Algorithme lmf (voir note)

#******************************************************************************

f <- function(x) {
  y <- 5 * sin(x*pi/30) + 3
  y[x < 3.5] <- 3
  y[x > 30] <- 3
  return(y)
}

heights <- seq(-5, 30, 0.5)
ws <- f(heights)

ctg_retile <- readLAScatalog("6_filterb/2.0-42")
#ctg_retile <- readLAScatalog("6_filterb/3.5-42")
#opt_restart(ctg_retile) <- 14
opt_output_files(ctg_retile) <-  paste0("8_ttops/", "/{*}_ttops")
opt_chunk_buffer(ctg_retile) <- 20
opt_laz_compression(ctg_retile) <- TRUE

#--- Utilise le format geopackage plutot que shapefile

ctg_retile@output_options$drivers$sf$extension <- '.gpkg'

ctg_ttops <- locate_trees(ctg_retile , lmf(ws = f, hmin = 3.5, shape = "circular"), uniqueness = "bitmerge")

#--- Autre exemple

# ctg_ttops <- locate_trees(ctg_retile , lmf(4), uniqueness = "bitmerge")

#--- Exemple d'affichage

# plot(ctg_retile) |> add_treetops3d(ctg_ttops)

#******************************************************************************

# Segmentation du nuage de points par arbre
# Algo watershed (voir note)
# D'autres algorithmes basés sur le nuage de points sont disponibles 

#******************************************************************************

ctg_retile <- readLAScatalog("6_filterb/2.0-42")
#ctg_retile <- readLAScatalog("6_filterb/3.5-42")
#opt_restart(ctg_retile) <- 29
#chm <- raster("7_chm/2.0-42/0.25/focal5x5/rasterize_canopy.vrt")
#chm <- raster("7_chm/2.0-42/0.5/rasterize_canopy.vrt")
chm <- raster("7_chm/chm_focal.tif")
opt_laz_compression(ctg_retile) <- TRUE
opt_filter(ctg_retile) <- "-drop_withheld"
opt_chunk_buffer(ctg_retile) <- 30
opt_output_files(ctg_retile) <- paste0("9_segmented/2.0-42", "/{*}_segmented")
algo <- lidR::watershed(chm, th_tree = 3)
ctg_segmented<- segment_trees(ctg_retile, algo, uniqueness = "bitmerge")

#******************************************************************************

# Autres nettoyage post segmentation
# avant vectorisation des couronnes (optionnel)

#******************************************************************************

treefilter <- function(chunk)
{
	#--- Enlève les points non associés à un arbre
	
	nlas <- filter_poi(chunk, !is.na(treeID))
	if (is.empty(nlas)) return(NULL)
	nlas <- filter_poi(nlas, !is.null(treeID))
	if (is.empty(nlas)) return(NULL)

	#nlas <- classify_noise(nlas , sor(50, 3))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#nlas <- classify_noise(nlas , sor(50, 3))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#nlas <- classify_noise(nlas , sor(50, 3))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)

	#--- 2000 points in the 27 voxels (5 m wide)
	
	#nlas <- classify_noise(nlas , ivf(5, 2000))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	
	#--- 500 points in the 27 voxels (3 m wide)
	
	#nlas <- classify_noise(nlas , ivf(3, 500))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	
	#--- 50 points in the 27 voxels (1.5 m wide)
	
	#nlas <- classify_noise(nlas , ivf(1.5, 50))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	
	#--- 20 points in the 27 voxels (1 m wide)
	
	#nlas <- classify_noise(nlas , ivf(1, 20))
	#nlas <- filter_poi(nlas, Classification != LASNOISE)
	#if (is.empty(nlas)) return(NULL)
	
	nlas <- filter_poi(nlas, Classification != LASUNCLASSIFIED)
	if (is.empty(nlas)) return(NULL)

	#--- Nécessaire car l'aggregation ne fonctionne pas sur des données vides
	
	#if (is.empty(nlas)) return(NULL)

	#--- Exemple de suppression d'arbres en fonction de paramètres
	
	# Enlève les points qui appartiennent à des arbres de moins de 35 points
	
	M2 <- aggregate(X ~ treeID, data = nlas@data, FUN = length)
	M3 <- M2[which(M2$X>35),]
	nlas <- filter_poi(nlas, treeID %in% M3$treeID)

	if (is.empty(nlas)) return(NULL)

	return(nlas)
}

ctg_segmented <- readLAScatalog("9_segmented/")
#opt_restart(ctg_segmented) <- 14
opt_output_files(ctg_segmented) <-  paste0("9_segmented/filtered/", "/{*}_filtered")
opt_laz_compression(ctg_segmented) <- TRUE
opt_filter(ctg_segmented) <- "-drop_withheld"
opt_chunk_buffer(ctg_segmented) <- 30
options <- list(need_output_file = TRUE, automerge = TRUE, overwrite = TRUE)
ctg_trees <- catalog_map(ctg_segmented, treefilter, .options = options)
lidR:::catalog_laxindex(ctg_trees)

#******************************************************************************

# Crée les couronnes vectorielles avec quelques 
# métriques

#******************************************************************************

#ctg_trees <- readLAScatalog("9_segmented/filtered")
ctg_trees <- readLAScatalog("9_segmented/2.0-42")
#opt_restart(ctg_trees) <- 14
opt_output_files(ctg_trees) <-  paste0("A_final/2.0-42/", "/{*}_final")
opt_laz_compression(ctg_trees) <- TRUE
opt_chunk_buffer(ctg_trees) <- 30
opt_filter(ctg_trees) <- "-drop_withheld"
options <- list(need_output_file = TRUE, automerge = TRUE, autoread = FALSE)
ctg_trees@output_options$drivers$sf$extension <- '.gpkg'

#--- Différents exemples de vectorisation des couronnes (polygones convexes, concaves, boites englobantes...)

pol = crown_metrics(ctg_trees, ~list(0, z_max = max(Z), z_mean = mean(Z), r_mean = mean(R), g_mean = mean(G), b_mean = mean(B), NIR_mean = mean(NIR), NEVI_mean = mean(NEVI), z_sd = sd(Z), i_mean = mean(Intensity)), geom = "convex")
#pol = crown_metrics(ctg_trees, ~list(0, z_max = max(Z), z_mean = mean(Z), r_mean = mean(R), g_mean = mean(G), b_mean = mean(B), NIR_mean = mean(NIR), NEVI_mean = mean(NEVI), z_sd = sd(Z), i_mean = mean(Intensity)), geom = "concave", concaveman = c(3, 5))
#pol = crown_metrics(ctg_trees, ~list(z_max = max(Z), z_mean = mean(Z), r_mean = mean(R), g_mean = mean(G), b_mean = mean(B), NIR_mean = mean(NIR), NEVI_mean = mean(NEVI), z_sd = sd(Z), i_mean = mean(Intensity)), geom = "bbox")
#pol = delineate_crowns(nlas, func = ~list(z_max = max(Z), z_mean = mean(Z), r_mean = mean(R), g_mean = mean(G), b_mean = mean(B), NIR_mean = mean(NIR), NEVI_mean = mean(NEVI), z_sd = sd(Z), i_mean = mean(Intensity)), type = "concave", concavity = 3, length_threshold = 5)

#******************************************************************************

# Autres petites fonctions utiles

#******************************************************************************

bbox <- st_bbox(ctg_trees)
st_write(bbox, "bbox.gpkg")
st_write(ctg_trees$data, "bboxes.gpkg")

#--- Possibilité de sauver directement le résultat dans un géopackage
st_write(pol , "crowns.gpkg")


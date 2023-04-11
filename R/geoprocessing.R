library(sf)
library(raster)
library(leaflet)

# Web mapping
pts <- sf::st_as_sf(sites_meta, coords = c("x_lv03", "y_lv03"), crs = 21781)
pts <- sf::st_transform(pts, 4326)

leaflet(pts) %>% addTiles() %>% addCircleMarkers(sf::st_coordinates(pts)[,1], sf::st_coordinates(pts)[,2], )

sf::write_sf(pts, "data/sampling_points.shp")


# Subcatchments
pts <- sf::st_as_sf(sites_meta, coords = c("x_lv03", "y_lv03"), crs = 21781)
catchments <- read_sf("/home/ecoadmin/Documents/postdoc/Swiss_catchments/EZG_Gewaesser.gpkg", layer = "Teileinzugsgebiet") %>%
  st_transform(21781)

catch_pts <- pts %>%
  st_intersects(catchments) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(up_catchments = map(col.id, function(x) {
    filter(catchments,
           H1 >= catchments[x, "H1", drop = TRUE],
           H2 <= catchments[x, "H2", drop = TRUE]) %>%
      st_union()
  }))

catch_pts <- catch_pts %>%
  unnest(up_catchments) %>%
  st_as_sf()

catch_pts <- pts %>%
  st_buffer(dist = 10000) %>%
  rename(buffer = geometry) %>%
  mutate(points = pts$geometry, up_catchments = catch_pts$up_catchments, .before = buffer) %>%
  st_set_geometry("points") %>%
  group_by(site_code) %>%
  mutate(upstr_area = st_intersection(up_catchments, buffer)) %>%
  ungroup()

catch_pts$upstr_area %>% write_sf("/home/ecoadmin/Documents/temp/upstr_area.shp")

catch_pts <- dplyr::select(catch_pts, -up_catchments, -buffer)

# Land use

var_names <- read_csv("data/landuse_variables.csv")
lc <- readRDS("/home/ecoadmin/Documents/temp/Swiss_landcover.rds")

lc_vals <- names(lc)[4:12] %>%
  map(function(x) {
    rasterFromXYZ(lc[,c("E", "N", x)], crs = "+init=epsg:2056") %>%
      projectRaster(crs = "+init=epsg:21781") %>%
      raster::extract(catch_pts %>% st_set_geometry("upstr_area"))
  })

lc_all <- map(lc_vals, function(x) {
  x %>%
    set_names(catch_pts$site_code) %>%
    map(as_tibble) %>%
    bind_rows(.id = "site_code") %>%
    mutate(value = as.integer(round(value))) %>%
    group_by(site_code) %>%
    count(value) %>%
    ungroup()
}) %>%
  set_names(names(lc)[4:12]) %>%
  bind_rows(.id = "layer") %>%
  rename(category = value) %>%
  left_join(var_names)

write_csv(lc_all, "data/landuse.csv")






lu_val %>%
  left_join(var_names) %>%
  filter(layer == "AS18_4") %>%
  filter(description != "Surfaces improductives") %>%
  group_by(site_code) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  spread_cdm(site_code, description, prop) %>%
  rda()

lu_val %>%
  left_join(var_names) %>%
  filter(layer == "AS18_4") %>%
  filter(description != "Surfaces improductives") %>%
  group_by(site_code) %>%
  mutate(lu_prop = n/sum(n)) %>%
  ungroup() %>%
  mutate(description = case_when(description == "Surfaces d’habitat et d’infrastructure" ~ "Urban",
                                 description == "Surfaces agricoles" ~ "Agriculture",
                                 description == "Surfaces boisées" ~ "Forest")) %>%
  dplyr::select(site_code, description, lu_prop) %>%
  write_csv("data/landuse_3cat.csv")


landuse <- read_csv("data/landuse_3cat.csv")
landuse %>%
  #left_join(sites_meta) %>%
  left_join(samples_meta) %>%
  filter(!is.na(project)) %>%
  pivot_wider(names_from = description, values_from = lu_prop) %>%
  ggtern::ggtern(aes(x = Urban, y = Agriculture, z = Forest, color = project)) +
  geom_point() +
  ggtern::geom_crosshair_tern(data = tibble(Urban = 1/3, Agriculture = 1/3, Forest = 1/3, project = NA))



landuse_nmds <- landuse_all %>%
  filter(layer == "LU18_10") %>%
  filter(!is.na(description)) %>%
  group_by(layer, site_code, description) %>%
  summarise(n = sum(n)) %>%
  group_by(layer, site_code) %>%
  mutate(lu_prop = n/sum(n)) %>%
  ungroup() %>%
  select(site_code, description, lu_prop) %>%
  spread_cdm(site_code, description, lu_prop, fill.missing = 0) %>%
  metaMDS()

landuse_nmds %>%
  scores() %>%
  as_tibble(rownames = "site_code") %>%
  ggplot() +
  geom_point(aes(NMDS1, NMDS2))


landuse_rda <- landuse_all %>%
  filter(layer == "LU18_10") %>%
  filter(!is.na(description)) %>%
  group_by(layer, site_code, description) %>%
  summarise(n = sum(n)) %>%
  group_by(layer, site_code) %>%
  mutate(lu_prop = n/sum(n)) %>%
  ungroup() %>%
  select(site_code, description, lu_prop) %>%
  spread_cdm(site_code, description, lu_prop, fill.missing = 0) %>%
  rda()

ggvegan:::autoplot.rda(landuse_rda)


# Get Main basins

pts <- sf::st_as_sf(sites_meta, coords = c("x_lv03", "y_lv03"), crs = 21781)
catchments <- read_sf("/home/ecoadmin/Documents/postdoc/Swiss_catchments/EZG_Gewaesser.gpkg", layer = "Teileinzugsgebiet") %>%
  st_transform(21781)

pts_match <- pts %>%
  st_intersects(catchments) %>%
  as.data.frame() %>%
  as_tibble()

pts_match %>%
  mutate(pts[row.id, "site_code"],
         catchments[col.id, "FLUSSGB"])



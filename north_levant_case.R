####North Levant case study

##Import data
#Vector
north_sites <- st_read("data/levant_sites_north_3395.shp")
north_roads <- st_read("data/levant_roads_north.shp")

north_roads <- st_zm(north_roads, drop=TRUE, what = "ZM")

#Raster
#Conductance surface for leastcostpath (0)
north_cs90 <- terra::rast("data/north_conductance_90.tif")

#Conductance surface for gdistance (NA)
north_cs_na <- raster::raster("data/north_conductance_90_na.tif")
target_crs <- CRS("EPSG:3395")
north_cs_na_3395 <- projectRaster(north_cs_na, crs = target_crs, method = "bilinear")

#DEM
north_dem90 <- terra::rast("data/north_90.tif")
north_dem_rast <- raster::raster("data/north_90.tif")

###Create network, get edge list
start_points <- st_sfc(lapply(st_geometry(north_roads), get_endpoints, "start"), crs = st_crs(north_roads))
end_points   <- st_sfc(lapply(st_geometry(north_roads), get_endpoints, "end"), crs = st_crs(north_roads))

start_sf <- st_sf(geometry = start_points)
end_sf   <- st_sf(geometry = end_points)

start_nearest <- st_nearest_feature(start_sf, north_sites)
end_nearest   <- st_nearest_feature(end_sf, north_sites)

north_roads$from_id <- north_sites$idAll[start_nearest]
north_roads$to_id   <- north_sites$idAll[end_nearest]

#Convert idAll, from_id and to_id fields to character. They are currently numeric type which sfnetworks does not accept. Also rename old 'id' field to 'PleiadesId'
north_roads$from_id <- as.character(north_roads$from_id)
north_roads$to_id   <- as.character(north_roads$to_id)
north_sites$idAll   <- as.character(north_sites$idAll)
colnames(north_sites)[5] <- "PleiadesID"

#Reorder columns in roads
north_roads_ordered <- north_roads[, c("from_id", "to_id", "lengthGeo", "Type", "typeWeight", "Avg_Slope", "pace", "timeWeight", "geometry")]

#Build the network
north_road_network <- as_sfnetwork(x = north_sites, edges = north_roads_ordered, node_key = "idAll", from = "from_id", to = "to_id", directed = FALSE, edges_as_lines = TRUE, length_as_weight = FALSE)

#As igraph
north_road_network_ig <- as.igraph(north_road_network)

##Get edge list
#Get vertex attributes (especially idAll)
node_ids <- north_road_network %>%
  activate("nodes") %>%
  as_tibble() %>%
  pull(idAll)

#Extract the edge list using vertex indices
edge_df <- igraph::as_data_frame(north_road_network_ig, what = "edges")

#Add from_id and to_id by indexing into node_ids
edge_df$from_id <- node_ids[as.numeric(edge_df$from)]
edge_df$to_id   <- node_ids[as.numeric(edge_df$to)]

#Join with original edge attributes
edge_attributes <- north_road_network %>%
  activate("edges") %>%
  as_tibble()

#Combine from_id and to_id
north_edge_list <- bind_cols(
  edge_df[, c("from_id", "to_id")],
  edge_attributes
)

##Cost distance matrix
#Convert sites to sp object (for gdistance)
north_sites_sp <- as(north_sites, "Spatial")
north_sites_sp <- spTransform(north_sites_sp, crs(north_cs_na_3395))

# Create transition layer using mean(1/x) function since our raster layer values represent conductance, raster must have NA values where not passable (water)
tr_n <- transition(north_cs_na_3395, transitionFunction = mean, directions = 16)
tr_n <- geoCorrection(tr_n, type = "c")

#Calculate cost distance matrix
cd_n <- gdistance::costDistance(tr_n, north_sites_sp)
cd_n_matrix <- as.matrix(cd_n)

###Create Gabriel Graph
gg_north <- cccd::gg(cd_n_matrix, r = 1)

#Extract edge list
gg_north_edges <- igraph::as_data_frame(gg_north, what = "edges")

#Map indices to idAll
id_map_n <- north_sites$idAll

gg_north_edges$from_id <- id_map_n[as.integer(gg_north_edges$from)]
gg_north_edges$to_id   <- id_map_n[as.integer(gg_north_edges$to)]

#Get GG edge list where from_id, to _id corresponds to idAll
gg_north_edge_list <- gg_north_edges %>%
  dplyr::select(from_id, to_id)

###Compare edges (common edges, only in road network, only in modeled network)
north_roads_edges_key <- sort_edges(north_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

gg_north_edges_key <- sort_edges(gg_north_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

#Edge lists
common_roads_gg_north <- intersect(north_roads_edges_key, gg_north_edges_key)
only_roads_gg_north <- setdiff(north_roads_edges_key, gg_north_edges_key)
only_gg_roads_north <- setdiff(gg_north_edges_key, north_roads_edges_key)

#To data frames
common_roads_gg_north_df <- do.call(rbind, strsplit(common_roads_gg_north, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_roads_gg_north_df <- do.call(rbind, strsplit(only_roads_gg_north, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_gg_roads_north_df <- do.call(rbind, strsplit(only_gg_roads_north, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

##Add type categories
common_roads_gg_north_df$type <- "RoadsGG"
only_roads_gg_north_df$type <- "RoadsNoGG"
only_gg_roads_north_df$type <- "GGNoRoads"

#Add to south_edge_list as boolean
north_edge_list <- north_edge_list %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id),
    key = paste(edge_min, edge_max, sep = "-"),
    commonGG  = key %in% gg_north_edges_key) %>%
  dplyr::select(-edge_min, -edge_max, -key)

##Export
north_roads_edges <- st_as_sf(north_edge_list)
st_write(north_roads_edges, "output/north_roads_edges.shp")

###Export constructed graphs as straight line sf objects
#Extract point coordinates into a data frame
sites_coord_north <- north_sites %>%
  mutate(
    idAll = as.character(idAll),
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

#Extract edges and add coordinates GG
gg_north_coord <- gg_north_edge_list %>%
  left_join(sites_coord, by = c("from_id" = "idAll")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sites_coord, by = c("to_id" = "idAll")) %>%
  rename(x_to = x, y_to = y)

gg_north_sf <- gg_north_coord %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2, byrow = FALSE)),
      crs = st_crs(sites)
    )
  ) %>%
  ungroup() %>%
  st_as_sf()

st_write(gg_north_sf, "output/gg_north.shp")

###Gabriel graphs LCP
#Create conductivity surface with leastcostpath
north_cs <- leastcostpath::create_cs(north_cs90, neighbours = 16)

#Drop some obviously wrongly calculated edges
drop_gg_n_pairs <- data.frame(
  from_id = c(4, 5, 289, 293, 293, 316, 340, 346, 348, 355, 366, 670),
  to_id = c(352, 391, 346, 330, 640, 348, 346, 412, 670, 412, 376, 775)
)

drop_gg_n_pairs$from_id <- as.character(drop_gg_n_pairs$from_id)
drop_gg_n_pairs$to_id   <- as.character(drop_gg_n_pairs$to_id)

gg_north_edge_list_cl <- gg_north_edge_list %>%
  anti_join(drop_gg_n_pairs, by = c("from_id", "to_id"))

##LCP network for GG
gg_north_lcps_list <- list()

calculate_lcp_north <- function(row, x, sites) {
  # Get origin and destination points by matching from_id and to_id to idAll
  origin <- north_sites %>% filter(idAll == row$from_id)
  destination <- north_sites %>% filter(idAll == row$to_id)
  #Calculate LCP
  lcp <- create_lcp(x, origin, destination)
  # Add from_id to_id to the result
  lcp$from_id <- row$from_id
  lcp$to_id <- row$to_id
  
  return(lcp)
}

batch_size <- 1
indices <- seq_len(nrow(gg_north_edge_list_cl))
batches <- split(indices, ceiling(indices / batch_size))

for (i in seq_along(batches)) {
  cat("Processing batch", i, "of", length(batches), "\n")
  
  batch_results <- map(batches[[i]], function(j) {
    tryCatch({
      calculate_lcp_north(gg_north_edge_list_cl[j, ], x = north_cs, sites = north_sites)
    }, error = function(e) {
      message("Error in row ", j, ": ", e$message)
      return(NULL)
    })
  })
  
  # Store results and clean memory
  gg_north_lcps_list <- c(gg_north_lcps_list, batch_results)
  gc(verbose = FALSE)
}

gg_north_lcps <- do.call(rbind, gg_north_lcps_list)
st_write(gg_north_lcps, "output/gg_north_lcps.shp")

####Sites/settlements only <- for L-networks? Filtering stations 291 and 717 - uncertain identification.
north_sites_sel <- north_sites %>%
  filter(featureTyp %in% c("city", "settlement", "station")) %>%
  filter(idAll != 291, idAll !=717)

#Convert sites to sp object (for gdistance)
north_sites_sel_sp <- as(north_sites_sel, "Spatial")

#Calculate cost distance matrix
cd_n_a <- gdistance::costDistance(tr_n, north_sites_sel_sp)
cd_n_a_matrix <- as.matrix(cd_n_a)

#Calculate time distance matrix
north_cs_tobler <- leastcostpath::create_slope_cs(north_dem90, cost_function = "tobler", neighbours = 16)

#Create an empty transitivity layer that will be populated with values from south_cs_tobler
north_tobler_tr <- gdistance::transition(north_dem_rast, transitionFunction = function(x) 1, directions = 16)
north_tobler_tr@transitionMatrix <- north_cs_tobler$conductanceMatrix
north_tobler_tr@matrixValues <- "transitionMatrix"

cd_n_time <- gdistance::costDistance(north_tobler_tr, north_sites_sel_sp)
cd_n_time_matrix <- as.matrix(cd_n_time)

#Time threshold at 57600 s (16 hours), to connect Germanicea (remote site)
north_L_network <- function(cd_time_matrix, cd_matrix, north_sites_sel, time_cutoff = 57600, L_range = 4) {
  result_list <- list()
  
  for (L in L_range) {
    message("======================")
    message("Processing for L = ", L)
    message("======================")
    
    network_list <- list()
    start_time <- Sys.time()
    
    for (i in seq_len(nrow(cd_time_matrix))) {
      message("Processing site ", i, "/", nrow(cd_time_matrix), " ...")
      
      reachable_idx <- which(cd_time_matrix[i, ] <= time_cutoff & !is.infinite(cd_time_matrix[i, ]))
      message("  - Found ", length(reachable_idx), " reachable sites")
      if (length(reachable_idx) < 2) {
        message("  - Skipped: not enough reachable sites")
        next
      }
      
      subset_sites <- north_sites_sel[reachable_idx, ]
      subset_ids <- reachable_idx
      subset_sites$node_id <- north_sites_sel$idAll[subset_ids]
      
      subset_cd <- as.matrix(cd_matrix)[subset_ids, subset_ids]
      
      edge_list <- t(combn(seq_len(nrow(subset_cd)), 2))
      edges_df <- data.frame(
        from = edge_list[, 1],
        to = edge_list[, 2],
        weight = subset_cd[edge_list]
      )
      edges_df_rev <- edges_df %>% rename(from = to, to = from)
      edges_all <- bind_rows(edges_df, edges_df_rev)
      
      net <- sfnetwork(
        nodes = subset_sites,
        edges = edges_all,
        directed = FALSE
      )
      
      message("  - Pruning edges...")
      nodes_df <- net %>% activate("nodes") %>% as_tibble()
      
      net <- net %>%
        activate("edges") %>%
        mutate(keep = map_lgl(row_number(), function(e_idx) {
          edge <- .E()[e_idx, ]
          X <- edge$from
          Z <- edge$to
          XZ <- edge$weight
          
          node_ids <- seq_len(nrow(subset_cd))
          intermediates <- setdiff(node_ids, c(X, Z))
          
          # Keep edge if no triangle detour is shorter than XZ * (1 + 1/L)
          keep_edge <- !any(sapply(intermediates, function(Y) {
            XY <- subset_cd[X, Y]
            YZ <- subset_cd[Y, Z]
            if (is.infinite(XY) || is.infinite(YZ)) return(FALSE)
            return(XY + YZ < XZ * (1 + 1/L))
          }))
          
          return(keep_edge)
        })) %>%
        filter(keep) %>%
        mutate(
          from_id = nodes_df$node_id[from],
          to_id   = nodes_df$node_id[to]
        )
      
      network_list[[i]] <- net
      
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      avg_time <- as.numeric(elapsed) / i
      remaining <- round(avg_time * (nrow(cd_time_matrix) - i))
      message("  - Final edge count: ", nrow(as_tibble(net, "edges")))
      message("  - Done. Est. time remaining: ", round(remaining / 60, 1), " min\n")
    }
    
    network_list <- network_list[!sapply(network_list, is.null)]
    
    all_nodes <- map_dfr(network_list, ~ 
                           .x %>% activate("nodes") %>% as_tibble() %>% select(node_id, geometry))
    
    unique_nodes <- all_nodes %>% distinct(node_id, .keep_all = TRUE) %>% st_as_sf()
    
    all_edges <- map_dfr(network_list, ~ 
                           .x %>% activate("edges") %>% as_tibble() %>% select(from_id, to_id, weight))
    
    all_edges <- all_edges %>%
      rowwise() %>%
      mutate(
        from_clean = min(from_id, to_id),
        to_clean   = max(from_id, to_id)
      ) %>%
      ungroup() %>%
      select(from_id = from_clean, to_id = to_clean, weight) %>%
      distinct(from_id, to_id, weight)
    
    unique_nodes <- unique_nodes %>% arrange(node_id) %>% mutate(new_index = row_number())
    
    edges_indexed <- all_edges %>%
      left_join(unique_nodes %>% select(node_id, from_idx = new_index), by = c("from_id" = "node_id")) %>%
      left_join(unique_nodes %>% select(node_id, to_idx = new_index), by = c("to_id" = "node_id")) %>%
      select(from = from_idx, to = to_idx, weight)
    
    network_final <- sfnetwork(
      nodes = unique_nodes,
      edges = edges_indexed,
      directed = FALSE
    )
    
    result_list[[paste0("L", L)]] <- network_final
  }
  
  return(result_list)
}

north_network_L3 <- north_L_network(cd_n_time_matrix, cd_n_matrix, north_sites_sel)

north_L_network2 <- function(cd_time_matrix, cd_matrix, north_sites_sel, time_cutoff = 45000, L_range = 4) {
  result_list <- list()
  
  for (L in L_range) {
    message("======================")
    message("Processing for L = ", L)
    message("======================")
    
    network_list <- list()
    start_time <- Sys.time()
    
    for (i in seq_len(nrow(cd_time_matrix))) {
      message("Processing site ", i, "/", nrow(cd_time_matrix), " ...")
      
      reachable_idx <- which(cd_time_matrix[i, ] <= time_cutoff & !is.infinite(cd_time_matrix[i, ]))
      message("  - Found ", length(reachable_idx), " reachable sites")
      if (length(reachable_idx) < 2) {
        message("  - Skipped: not enough reachable sites")
        next
      }
      
      subset_sites <- north_sites_sel[reachable_idx, ]
      subset_ids <- reachable_idx
      subset_sites$node_id <- north_sites_sel$idAll[subset_ids]
      
      subset_cd <- as.matrix(cd_matrix)[subset_ids, subset_ids]
      
      edge_list <- t(combn(seq_len(nrow(subset_cd)), 2))
      edges_df <- data.frame(
        from = edge_list[, 1],
        to = edge_list[, 2],
        weight = subset_cd[edge_list]
      )
      edges_df_rev <- edges_df %>% rename(from = to, to = from)
      edges_all <- bind_rows(edges_df, edges_df_rev)
      
      net <- sfnetwork(
        nodes = subset_sites,
        edges = edges_all,
        directed = FALSE
      )
      
      message("  - Pruning edges...")
      nodes_df <- net %>% activate("nodes") %>% as_tibble()
      
      net <- net %>%
        activate("edges") %>%
        mutate(keep = map_lgl(row_number(), function(e_idx) {
          edge <- .E()[e_idx, ]
          X <- edge$from
          Z <- edge$to
          XZ <- edge$weight
          
          node_ids <- seq_len(nrow(subset_cd))
          intermediates <- setdiff(node_ids, c(X, Z))
          
          # Keep edge if no triangle detour is shorter than XZ * (1 + 1/L)
          keep_edge <- !any(sapply(intermediates, function(Y) {
            XY <- subset_cd[X, Y]
            YZ <- subset_cd[Y, Z]
            if (is.infinite(XY) || is.infinite(YZ)) return(FALSE)
            return(XY + YZ < XZ * (1 + 1/L))
          }))
          
          return(keep_edge)
        })) %>%
        filter(keep) %>%
        mutate(
          from_id = nodes_df$node_id[from],
          to_id   = nodes_df$node_id[to]
        )
      
      network_list[[i]] <- net
      
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      avg_time <- as.numeric(elapsed) / i
      remaining <- round(avg_time * (nrow(cd_time_matrix) - i))
      message("  - Final edge count: ", nrow(as_tibble(net, "edges")))
      message("  - Done. Est. time remaining: ", round(remaining / 60, 1), " min\n")
    }
    
    network_list <- network_list[!sapply(network_list, is.null)]
    
    all_nodes <- map_dfr(network_list, ~ 
                           .x %>% activate("nodes") %>% as_tibble() %>% select(node_id, geometry))
    
    unique_nodes <- all_nodes %>% distinct(node_id, .keep_all = TRUE) %>% st_as_sf()
    
    all_edges <- map_dfr(network_list, ~ 
                           .x %>% activate("edges") %>% as_tibble() %>% select(from_id, to_id, weight))
    
    all_edges <- all_edges %>%
      rowwise() %>%
      mutate(
        from_clean = min(from_id, to_id),
        to_clean   = max(from_id, to_id)
      ) %>%
      ungroup() %>%
      select(from_id = from_clean, to_id = to_clean, weight) %>%
      distinct(from_id, to_id, weight)
    
    unique_nodes <- unique_nodes %>% arrange(node_id) %>% mutate(new_index = row_number())
    
    edges_indexed <- all_edges %>%
      left_join(unique_nodes %>% select(node_id, from_idx = new_index), by = c("from_id" = "node_id")) %>%
      left_join(unique_nodes %>% select(node_id, to_idx = new_index), by = c("to_id" = "node_id")) %>%
      select(from = from_idx, to = to_idx, weight)
    
    network_final <- sfnetwork(
      nodes = unique_nodes,
      edges = edges_indexed,
      directed = FALSE
    )
    
    result_list[[paste0("L", L)]] <- network_final
  }
  
  return(result_list)
}

north_networks_L <- north_L_network2(cd_n_time_matrix, cd_n_matrix, north_sites_sel)

north_network_L4 <- north_networks_L[["L4"]]

#Export network edges
net <- north_network_L4

#Extract edges and nodes
edges_df <- net %>% activate("edges") %>% as_tibble()
nodes_tbl <- net %>% activate("nodes") %>% as_tibble() %>%
  mutate(index = row_number())                # network's internal node index

#Make nodes an sf object
nodes_sf <- st_sf(nodes_tbl, crs = st_crs(net))

# Ensure same CRS between nodes and source sites
if (!identical(st_crs(nodes_sf), st_crs(north_sites_sel))) {
  north_sites_sel <- st_transform(north_sites_sel, st_crs(nodes_sf))
}

#Map node -> idAll by geometry (exact), fallback to nearest
eq_list <- st_equals(nodes_sf, north_sites_sel)  # list of integer vectors
id_map <- map_chr(eq_list, function(idx) {
  if (length(idx) == 1) {
    as.character(north_sites_sel$idAll[idx])
  } else {
    NA_character_
  }
})

#Fallback for unmatched via nearest feature
na_idx <- which(is.na(id_map))
if (length(na_idx) > 0) {
  nearest <- st_nearest_feature(nodes_sf[na_idx, ], north_sites_sel)
  id_map[na_idx] <- as.character(north_sites_sel$idAll[nearest])
}

#Optional sanity check: warn if duplicates (shouldn’t happen unless duplicate coords)
dup_ids <- id_map[duplicated(id_map)]
if (length(dup_ids) > 0) {
  warning("Duplicated idAll mapped to multiple network nodes: ",
          paste(dup_ids, collapse = ", "))
}

#Build node -> (idAll, geometry) map
node_map <- tibble(
  index     = nodes_sf$index,
  idAll     = id_map,
  geometry  = nodes_sf$geometry
)

#Join node ids/geoms onto edges
edges_enriched <- edges_df %>%
  left_join(node_map %>% select(index, from_id = idAll, geom_from = geometry),
            by = c("from" = "index")) %>%
  left_join(node_map %>% select(index, to_id = idAll, geom_to = geometry),
            by = c("to" = "index"))

#Build edge geometries
edge_geom <- st_sfc(
  map2(edges_enriched$geom_from, edges_enriched$geom_to, ~ {
    st_linestring(rbind(st_coordinates(.x)[, 1:2], st_coordinates(.y)[, 1:2]))
  }),
  crs = st_crs(net)
)

edge_lines <- st_sf(
  edges_enriched %>% select(from_id, to_id, weight),
  geometry = edge_geom
)

st_write(edge_lines, "output/edges_L4_north.shp", delete_layer = TRUE)

#Get edge keys
north_edge_keys_L4 <- edges_df %>%
  mutate(
    from_id = id_map[from],
    to_id   = id_map[to]
  )

#Normalize edges
north_edge_keys_L4 <- sort_edges(edges_with_ids) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

compare_edges_n <- function(network_edges_key, label_prefix, other_key, other_label = "Roads") {
  common <- intersect(other_key, network_edges_key)
  only_other <- setdiff(other_key, network_edges_key)
  only_network <- setdiff(network_edges_key, other_key)
  
  edge_df <- function(keys, type) {
    if (length(keys) == 0) return(data.frame(from_id=character(), to_id=character(), type=character()))
    do.call(rbind, strsplit(keys, "-")) %>%
      as.data.frame() %>%
      setNames(c("from_id", "to_id")) %>%
      mutate(type = type)
  }
  
  list(
    common     = edge_df(common, paste0(other_label, label_prefix)),
    only_other = edge_df(only_other, paste0(other_label, "No", label_prefix)),
    only_net   = edge_df(only_network, paste0(label_prefix, "No", other_label))
  )
}

#Edges common between south_roads and L networks, only contained in roads, only contained in L networks
L4_north_compare <- compare_edges_n(north_edge_keys_L4, "L4", north_roads_edges_key)

L4_roads_north_comparison_table <- data.frame(
  common     = nrow(L4_north_compare$common),
  only_roads = nrow(L4_north_compare$only_other),
  only_L4    = nrow(L4_north_compare$only_net)
)

write.csv(L4_roads_north_comparison_table, "output/L4_north_roads_comparison.csv")

#Edges common between GG and L networks, only contained in GG, only contained in L networks
L4_GG_north <- compare_edges_n(north_edge_keys_L4, "L4", gg_north_edges_key)

L4_GG_north_comparison_table <- data.frame(
  common   = nrow(L4_GG_north$common),
  only_GG  = nrow(L4_GG_north$only_other),
  only_L4  = nrow(L4_GG_north$only_net)
)

write.csv(L4_GG_north_comparison_table, "output/L4_GG_north_comparison.csv")

#Edges in north_roads contained in L4 or GG (Boolean)
L4_GG_roads_north <- edges_union(north_edge_keys_L4, "L4", north_roads_edges_key, gg_north_edges_key)

L4_GG_roads_north_comparison_table <- data.frame(
  common = nrow(L4_GG_roads_north$common)
)

write.csv(L4_GG_roads_north_comparison_table, "output/L4_GG_roads_north_comparison.csv")

####Model L=4 LCPS
##Since L=4 network shares edges with GG graph, we need to model only the unique edges.Then combine L=4 LCPs with GG LCPs to get complete network.

#Filter GG LCPs for LCPs contained in L=4
gg_north_lcps <- gg_north_lcps %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id),
    edge_key = paste(edge_min, edge_max, sep = "-")
  )

L4_GG_north_lcps_common <- gg_north_lcps %>%
  filter(edge_key %in% north_edge_keys_L4)

#Create edge key unique to L=4 (i.e., excluding edges already contained in GG).
L4_GG_north_common <- L4_GG_north[["common"]] %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id),
    edge_key = paste(edge_min, edge_max, sep = "-")
  ) %>%
  pull(edge_key)

L4_north_unique_edges <- setdiff(north_edge_keys_L4, L4_GG_north_common)

L4_north_unique_edges_df <- data.frame(edge_key = L4_north_unique_edges) %>%
  tidyr::separate(edge_key, into = c("from_id", "to_id"), sep = "-", convert = TRUE)

#Calculate L=4 unique LCPs
L4_north_unique_lcps_list <- list()

batch_size <- 1
indices <- seq_len(nrow(L4_north_unique_edges_df))
batches <- split(indices, ceiling(indices / batch_size))

for (i in seq_along(batches)) {
  cat("Processing batch", i, "of", length(batches), "\n")
  
  batch_results <- map(batches[[i]], function(j) {
    tryCatch({
      calculate_lcp_north(L4_north_unique_edges_df[j, ], x = north_cs, sites = north_sites_sel)
    }, error = function(e) {
      message("Error in row ", j, ": ", e$message)
      return(NULL)
    })
  })
  
  #Store results and clean memory
  L4_north_unique_lcps_list <- c(L4_north_unique_lcps_list, batch_results)
  gc(verbose = FALSE)
}

L4_north_unique_lcps <- do.call(rbind, L4_north_unique_lcps_list)

L4_GG_north_lcps_common <- L4_GG_north_lcps_common %>%
  select(-edge_min, -edge_max, -edge_key)

L4_north_lcps <- rbind(L4_north_unique_lcps, L4_GG_north_lcps_common)

#Add GG edges not in L=4 for connectivity
GG_north_unique_edges <- setdiff(gg_north_edges_key, L4_GG_north_common)

GG_north_unique_lcps <- gg_north_lcps %>%
  filter(edge_key %in% GG_north_unique_edges) %>%
  select(-edge_key, -edge_min, -edge_max)

L4_GG_north_lcps <- rbind(L4_north_lcps, GG_north_unique_lcps)

st_write(L4_GG_north_lcps, "output/L4_GG_north_lcps.shp")

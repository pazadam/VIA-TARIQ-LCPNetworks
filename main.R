###Data preparation, network creation
##Load libraries
library(sf)
library(sfnetworks)
library(tidygraph)
library(dplyr)
library(igraph)
library(ggplot2)

##Load data, sites=nodes, roads=edges
sites <- sf::st_read("data/levant_sites.shp")
roads <- sf::st_read("data/levant_roads.shp")

##Drop Z dimension from the road dataset
roads <- st_zm(roads, drop=TRUE, what = "ZM")

##In order to link the nodes to the edges, we need a node key associated with the edges. We get start and end points of the edges and associate them with the nodes IDs, then add 2 new fields <from_ID and to_ID> to the edge dataset.
get_endpoints <- function(line, position = "start") {
  coords <- st_coordinates(line)
  if (position == "start") {
    st_point(coords[1, 1:2])
  } else {
    st_point(coords[nrow(coords), 1:2])
  }
}

start_points <- st_sfc(lapply(st_geometry(roads), get_endpoints, "start"), crs = st_crs(roads))
end_points   <- st_sfc(lapply(st_geometry(roads), get_endpoints, "end"), crs = st_crs(roads))

start_sf <- st_sf(geometry = start_points)
end_sf   <- st_sf(geometry = end_points)

start_nearest <- st_nearest_feature(start_sf, sites)
end_nearest   <- st_nearest_feature(end_sf, sites)

roads$from_id <- sites$idAll[start_nearest]
roads$to_id   <- sites$idAll[end_nearest]

#Convert idAll, from_id and to_id fields to character. They are currently numeric type which sfnetworks does not accept. Also rename old 'id' field to 'PleiadesId'
roads$from_id <- as.character(roads$from_id)
roads$to_id   <- as.character(roads$to_id)
sites$idAll   <- as.character(sites$idAll)
colnames(sites)[4] <- "PleiadesID"

#Reorder columns in roads
roads_ordered <- roads[, c("from_id", "to_id", "lengthGeo", "Type", "typeWeight", "Avg_Slope", "pace", "timeWeight", "geometry")]

#Build the network
road_network <- as_sfnetwork(x = sites, edges = roads_ordered, node_key = "idAll", from = "from_id", to = "to_id", directed = FALSE, edges_as_lines = TRUE, length_as_weight = FALSE)

###Analysis
##Global measures (size of the largest component,  density, global clustering coefficient, average local clustering coefficient, gamma index, detour)
#Size of the largest component
size_largest_comp <- components(road_network)
size_largest_comp <- data.frame(
  component = which.max(size_largest_comp$csize),
  size = max(size_largest_comp$csize)
)

#Separate subgraph of connected roads and sites
road_network_1 <- igraph::subgraph_from_edges(road_network, 1:964)

#Density
density <- data.frame(
    density_global = edge_density(road_network_1)
)

#Clustering coefficient global
clustering_coefficient_global <- data.frame(
  clustering_coefficient_global = transitivity(road_network_1, type = "global")
)

#Clustering coefficient local
clustering_coefficient_local <- data.frame(
  clustering_coefficient_local = transitivity(road_network_1, type = "average")
)

#Gamma index
gamma_index <- data.frame(
  gamma_index = ecount(road_network_1)/(3*(vcount(road_network_1)-2))
)

##Detour
#First, we need the subgraph of all connected sites and roads to be a sfnetwork object
road_network_2_sf <- road_network %>%
  activate(nodes) %>%
  mutate(component = group_components(type = "weak")) %>%
  filter(component == which.max(table(component)))

#Calculate Euclidean distance. Since lengthGeo is geodesic distance, we need to project nodes to WGS 84 to get geodesic distance in R.
edge_list <- as_data_frame(road_network_2_sf, what = "edges")
nodes_coord <- st_geometry(road_network_2_sf, what = "nodes")
nodes_coord <- st_transform(nodes_coord, 4326)
from_node <- nodes_coord[edge_list$from]
to_node <- nodes_coord[edge_list$to]
distance_euc <- st_distance(from_node, to_node, by_element = TRUE)
distance_euc <- units::drop_units(distance_euc)

#Get actual distances in the network
distance_network <- (edge_list$lengthGeo)

#Calculate Detour index for edges
detour_index <- data.frame (
  detour_index = distance_euc/distance_network)

##Calculate detour centrality
#Get node distance matrix
distance_euc_matrix <- st_distance(nodes_coord, nodes_coord)

#Assign geoLength as weight in the network 
distance_network_matrix <- distances(
  graph = road_network_2_sf %>% as.igraph(),
  weights = E(road_network_2_sf)$lengthGeo,
  algorithm = "automatic"
  )
distance_network_matrix <- units::set_units(distance_network_matrix, "m")

#Calculate detour matrix
detour_matrix <- distance_network_matrix/distance_euc_matrix

#Detour centrality
detour_centrality_normalised <- apply(detour_matrix, 1, 
                             function(row) {sum(row[is.finite(row)], na.rm = TRUE)/sum(is.finite(row))})

median_detour_centrality <- data.frame(
  median_detour_centrality = median(detour_centrality_normalised, na.rm = TRUE)
)

##Centrality measures

#Degree
road_network_degree <- degree(road_network_2_sf %>% as.igraph())

#Add degree as a node attribute attribute
road_network_degree_df <- data.frame(
  node_id = seq_along(road_network_degree),
  degree = as.numeric(road_network_degree)
)

road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(node_id = row_number()) %>%
  left_join(road_network_degree_df, by = "node_id")%>%
  select(-node_id)

#Mean and max degree
mean_degree <- data.frame(
  mean_degree = mean(road_network_degree)
)

max_degree <- data.frame(
  max_degree = max(road_network_degree)
)

#Degree distribution
ggplot(road_network_degree_df, aes(x = factor(degree))) +
  geom_bar(fill = "darkgrey", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  scale_x_discrete(drop = FALSE) +  # show all factor levels
  labs(
    title = "Node Degree Distribution",
    x = "Degree",
    y = "Frequency"
  ) +
  theme_minimal()

#By categories and periods
nodes_df <- road_network_2_sf %>%
  activate("nodes") %>%
  as.data.frame() %>%
  select(-degree.y) %>%
  rename(degree = degree.x)
  mutate(degree = as.integer(degree))

ggplot(nodes_df, aes(x = factor(degree))) +
  geom_bar(fill = "darkgrey", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  facet_wrap(~featureTyp, scales = "free_y") +
  labs(
    title = "Node Degree Distribution by Feature Type",
    x = "Degree",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

nodes_df_hel <- nodes_df %>%
  filter(hel == "yes")

nodes_df_rom <- nodes_df %>%
  filter(rom == "yes")

nodes_df_byz <- nodes_df %>%
  filter(byz == "yes")

ggplot(nodes_df_hel, aes(x = factor(degree))) +
  geom_bar(fill = "darkgrey", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  facet_wrap(~featureTyp, scales = "free_y") +
  labs(
    title = "Node Degree Distribution by Feature Type, Hellenistic period",
    x = "Degree",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggplot(nodes_df_rom, aes(x = factor(degree))) +
  geom_bar(fill = "darkgrey", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  facet_wrap(~featureTyp, scales = "free_y") +
  labs(
    title = "Node Degree Distribution by Feature Type, Roman period",
    x = "Degree",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggplot(nodes_df_byz, aes(x = factor(degree))) +
  geom_bar(fill = "darkgrey", color = "white") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  facet_wrap(~featureTyp, scales = "free_y") +
  labs(
    title = "Node Degree Distribution by Feature Type, Byzantine period",
    x = "Degree",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

#Mean degree by type and period
mean_degree_all <- nodes_df %>%
  group_by(featureTyp) %>%
  summarise(mean_deg = mean(degree, na.rm = TRUE)) %>%
  ungroup()

mean_degree_hel <- nodes_df %>%
  filter(hel == "yes") %>%
  group_by(featureTyp) %>%
  summarise(mean_deg = mean(degree, na.rm = TRUE)) %>%
  ungroup()

mean_degree_rom <- nodes_df %>%
  filter(rom == "yes") %>%
  group_by(featureTyp) %>%
  summarise(mean_deg = mean(degree, na.rm = TRUE)) %>%
  ungroup()

mean_degree_byz <- nodes_df %>%
  filter(byz == "yes") %>%
  group_by(featureTyp) %>%
  summarise(mean_deg = mean(degree, na.rm = TRUE)) %>%
  ungroup()

mean_degree_comparison <- mean_degree_byz %>%
  full_join(mean_degree_rom, by = "featureTyp") %>%
  full_join(mean_degree_hel, by = "featureTyp")%>%
  rename(hel = mean_deg) %>%
  rename(rom = mean_deg.y) %>%
  rename(byz = mean_deg.x) %>%
  select(featureTyp, hel, rom, byz)

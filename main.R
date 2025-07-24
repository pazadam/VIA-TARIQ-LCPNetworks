####Data preparation, network creation
##Load libraries
library(sf)
library(sfnetworks)
library(tidygraph)
library(dplyr)
library(igraph)
library(ggplot2)
library(spatgraphs)
library(spatstat)
library(deldir)
library(raster)
library(sp)
library(cccd)
library(gdistance)
library(leastcostpath)
library(purrr)

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

#### Analysis
### Global measures (size of the largest component,  density, global clustering coefficient, average local clustering coefficient, gamma index, alpha index, detour)
##Size of the largest component
size_largest_comp <- components(road_network)
size_largest_comp <- data.frame(
  component = which.max(size_largest_comp$csize),
  size = max(size_largest_comp$csize)
)

#Separate subgraph of connected roads and sites
road_network_1 <- igraph::subgraph_from_edges(road_network, 1:964)

##Diameter
network_diameter <- diameter(road_network_1, weights = NA)
network_distance_diameter <- diameter(road_network_1, weights = E(road_network_1)$lengthGeo)

diameter <- data.frame(
  diameter = network_diameter
)

distance_diameter <- data.frame(
  distance_diameter = network_distance_diameter
)

##Density
density <- data.frame(
    density_global = edge_density(road_network_1)
)

##Average shortest path
mean_shortest <- mean_distance(road_network_1, directed = FALSE, weights = NA)
average_shortest_path <- data.frame(
  average_shortest_path = mean_shortest)

##Clustering coefficient global
clustering_coefficient_global <- data.frame(
  clustering_coefficient_global = transitivity(road_network_1, type = "global")
)

##Clustering coefficient local
clustering_coefficient_local <- data.frame(
  clustering_coefficient_local = transitivity(road_network_1, type = "average")
)

##Gamma index
gamma_index <- data.frame(
  gamma_index = ecount(road_network_1)/(3*(vcount(road_network_1)-2))
)

##Alpha Index
alpha_index <- (ecount(road_network_1)-vcount(road_network_1)+1)/(2*vcount(road_network_1)-5)
alpha <- data.frame(
  alpha_in = alpha_index
)

##Detour index (edge based)
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
detour_index_df <- data.frame (
  detour_index = distance_euc/distance_network)

detour_index <- distance_euc/distance_network

#Add detour index as an edge attribute
road_network_2_sf <- road_network_2_sf%>%
  activate("edges")%>%
  mutate(detour_in = detour_index)

##Detour centrality (node based)
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

#Calculate detour centrality
detour_centrality_normalised <- apply(detour_matrix, 1, 
                             function(row) {sum(row[is.finite(row)], na.rm = TRUE)/sum(is.finite(row))})

median_detour_centrality <- data.frame(
  median_detour_centrality = median(detour_centrality_normalised, na.rm = TRUE)
)

###Centrality measures

##Degree
road_network_degree <- degree(road_network_2_sf %>% as.igraph())

#Add degree as a node attribute
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

##Closeness
#Distance-weighted closeness
closeness_distance <- closeness(road_network_1, weights = E(road_network_1)$lengthGeo)
closeness_distance_norm <- closeness(road_network_1, weights = E(road_network_1)$lengthGeo, normalized = TRUE)

#Add as a node attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(close = closeness_distance) %>%
  mutate(closeN = closeness_distance_norm)

#Time-weighted closeness
closeness_time <- closeness(road_network_1, weights = E(road_network_1)$timeWeight)
closeness_time_norm <- closeness(road_network_1, weights = E(road_network_1)$timeWeight, normalized = TRUE)

#Add as a node attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(closeTime = closeness_time) %>%
  mutate(closeTimeN = closeness_time_norm)

##Betweenness - nodes
#Distance weighted betweenness
betweenness_distance <- betweenness(road_network_1, v = V(road_network_1), weights = E(road_network_1)$lengthGeo)
betweenness_distance_n <- betweenness(road_network_1, v = V(road_network_1), weights = E(road_network_1)$lengthGeo, normalized = TRUE)

#Add as a node attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(betDist = betweenness_distance) %>%
  mutate(betDistN = betweenness_distance_n)

#Time weighted betweenness
betweenness_time <- betweenness(road_network_1, v = V(road_network_1), weights = E(road_network_1)$timeWeight)
betweenness_time_n <- betweenness(road_network_1, v = V(road_network_1), weights = E(road_network_1)$timeWeight, normalized = TRUE)

#Add as a node attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(betTime = betweenness_time) %>%
  mutate(betTimeN = betweenness_time_n)

##Betweenness - edges
#Distance weighted betweenness
betweennessE_distance <- edge_betweenness(road_network_1, e = E(road_network_1), weights = E(road_network_1)$lengthGeo)

#Add as an edge attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("edges") %>%
  mutate(betDist = betweennessE_distance)

#Time weighted betweenness
betweennessE_time <- edge_betweenness(road_network_1, e = E(road_network_1), weights = E(road_network_1)$timeWeight)

#Add as an edge attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("edges") %>%
  mutate(betTime = betweennessE_time)

##Eigenvector
#Unweighted
eigenvector <- eigen_centrality(road_network_1, weights = NA)

#Add as a node attribute
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(eigen = eigenvector$vector)

#Road type weighted
eigen_type <- eigen_centrality(road_network_1, weights = E(road_network_1)$typeWeight)

#Add as a node attribute
eigen_vec <- eigen_type$vector

road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(eigenType = eigen_vec)

##Export network properties table
network_properties <- bind_cols(density, diameter, mean_degree, average_shortest_path, clustering_coefficient_global, clustering_coefficient_local, gamma_index, alpha, median_detour_centrality)
write.csv(network_properties, file = "output/network_properties.csv")

##Export roads and sites with edge and node properties
roads_edges <- road_network_2_sf %>%
  activate("edges") %>%
  st_as_sf()

write_sf(roads_edges, "output/roads_edges.shp")

sites_nodes <- road_network_2_sf %>%
  activate("nodes") %>%
  st_as_sf()

sites_nodes <- sites_nodes %>%
  select(-component) %>%
  select(-degree.y) %>%
  rename(degree = degree.x)

write_sf(sites_nodes, "output/sites_nodes.shp")

###Network comparison
###We compare the network with theoretical planar networks (Gabriel Graph, Relative Neighbourhood Graph, Minimum Spanning Tree, Delaunay Triangulation, K=4 Nearest Neighbours)

#Data preparation
#Sites dataset contains also placeholder nodes at intersections that could not be tied to any archaeological site. For network construction we will exclude these nodes.
#Prepare observation window for spatstat to create point pattern (ppp)
b_box <- st_read("data/b_box.shp")
b_box <- st_zm(b_box, drop=TRUE, what = "ZM")
b_box_geo <- st_geometry(b_box)
b_box_owin <- as.owin(b_box_geo)

#Extract point coordinates from meaningful nodes in the connected network to create ppp
nodes <- road_network_2_sf %>%
  activate("nodes") %>%
  filter(featureTyp %in% c("city", "settlement", "fort", "station", "site", "bridge")) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  setNames(c("x","y")) %>%
  na.omit() %>%
  unique()

nodes_ppp <- ppp(nodes$x, nodes$y, window = b_box_owin)

##Gabriel graph
gg_c <- spatgraph(nodes_ppp, type = "gabriel")
gg_c_igraph <- graph_from_adj_list(gg_c$edges, mode = "all", duplicate = TRUE)

##Relative Neighbourhood Graph
rng_c <- spatgraph(nodes_ppp, type = "RNG")
rng_c_igraph <- graph_from_adj_list(rng_c$edges, mode = "all", duplicate = TRUE
                                    )
##Minimum Spanning Tree
mst_c <- spatgraph(nodes_ppp, type = "MST")
mst_c_igraph <- graph_from_adj_list(mst_c$edges, mode = "all", duplicate = TRUE)

##Delaunay Triangulation
dt_c <- deldir(nodes_ppp)
dt_c_edges <- dt_c$delsgs[, c("ind1", "ind2")]
dt_c_elist <- as.matrix(dt_c_edges)
dt_c_igraph <- graph_from_edgelist(dt_c_elist, directed = FALSE)

#Remove connections across the sea and the desert. This is done manually in the GIS, therefore the edges are exported to .shp, then the graph is recreated again with the removed links. 
dt_c_edges_sf <- dt_c$delsgs %>%
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(matrix(c(x1, x2, y1, y2), ncol = 2, byrow = FALSE)), crs = 3395)) %>%
  st_as_sf()

dt_cleaned <- st_read("data/dt_cleaned.shp")

dt_cleaned_edgelist <- dt_cleaned %>%
  select(ind1, ind2)

dt_igraph_cleaned <- graph_from_data_frame(dt_cleaned_edgelist, directed = FALSE)

##K=4 Nearest Neighbours
k4_c <- spatgraph(nodes_ppp, type = "knn", par = 4)
k4_c_igraph <- graph_from_adj_list(k4_c$edges, mode = "all", duplicate = FALSE)

##Add all to list
constructed_networks <- list(gg = gg_c_igraph, rng = rng_c_igraph, mst = mst_c_igraph, dt = dt_igraph_cleaned, k4 = k4_c_igraph)

##Calculate properties of the constructed networks
#Number of edges
number_of_edges_df <- constructed_networks %>% 
  lapply(ecount) %>%
  stack() %>%
  dplyr::select(ind, values) %>%
  setNames(c("network","number_of_edges"))

#Degree
degree_df <- constructed_networks %>% 
  lapply(degree) %>%
  lapply(mean) %>%
  stack() %>%
  dplyr::select(ind, values) %>%
  setNames(c("network","avg_degree"))

#Density
density_df <- constructed_networks %>%
  lapply(edge_density) %>%
  stack() %>%
  dplyr::select(ind, values) %>%
  setNames(c("network","density"))

#Gamma index
gamma_df <- constructed_networks %>%
  lapply(function(g) {
    ecount(g)/(3*(vcount(g)-2))
    }) %>%
  stack() %>%
  dplyr::select(ind, values) %>%
  setNames(c("network","gamma_index"))

#Alpha index
#Since alpha index for MSP is always 0, it is excluded
alpha_values <- sapply(names(constructed_networks), function(nm) {
  g <- constructed_networks[[nm]]
  if (grepl("mst", nm, ignore.case = TRUE)) {
    NA
  } else {
    (ecount(g)-vcount(g)+1)/(2*vcount(g)-5)
    }
  }
)

alpha_df <- data.frame(
  network = names(alpha_values),
  alpha_index = alpha_values
)

#Clustering coefficient global
ccg_values <- sapply(names(constructed_networks), function(nm) {
    if (grepl("mst|rng", nm, ignore.case = TRUE)) {
    NA_real_
  } else {
    transitivity(constructed_networks[[nm]], type = "global")
  }
})

ccg_df <- data.frame(
  network = names(ccg_values),
  ccg = ccg_values
)

#Clustering coefficient local
ccl_values <- sapply(names(constructed_networks), function(nm) {
  if (grepl("mst|rng", nm, ignore.case = TRUE)) {
    NA_real_
  } else {
    transitivity(constructed_networks[[nm]], type = "average")
  }
})

ccl_df <- data.frame(
  network = names(ccl_values),
  ccl = ccl_values
)

#Combine tables
constructed_networks_table <- number_of_edges_df %>%
  left_join(degree_df, by = "network") %>%
  left_join(density_df, by = "network") %>%
  left_join(gamma_df, by = "network") %>%
  left_join(alpha_df, by = "network") %>%
  left_join(ccg_df, by = "network") %>%
  left_join(ccl_df, by = "network")

write.csv(constructed_networks_table, "output/constructed_networks_table.csv")

###Robustness
###In order to evaluate robustness of the results of the centrality metrics (degree, time-weighted edge betweenness), a simple method adapted from 'Network Science in Archaeology' using Spearman's rho is used here.
sim_missing_edges <- function(net,
                              nsim = 1000,
                              props = c(0.9, 0.8, 0.7, 0.6, 0.5,
                                        0.4, 0.3, 0.2, 0.1),
                              met = NA,
                              missing_probs = NA) {
  # Initialize required library
  require(reshape2)
  
  props <- as.vector(props)
  
  if (FALSE %in% (is.numeric(props) & (props > 0) & (props <= 1))) {
    stop("Variable props must be numeric and be between 0 and 1",
         call. = F)
  }
  
  # Select measure of interest based on variable met and calculate
  if (!(met %in% c("degree", "betweenness"))) {
    stop(
      "Argument met must be either degree, betweenness, or eigenvector.
      Check function call.",
      call. = F
    )
  }
  else {
    if (met == "degree") {
      met_orig <- igraph::degree(net)
    }
    else  {
      if (met == "betweenness") {
        met_orig <- igraph::edge_betweenness(net, directed = FALSE, weights = E(net)$timeWeight)  # <<< UPDATED LINE
      }
    }
  }
  
  # Create data frame for out put and name columns
  output <- matrix(NA, nsim, length(props))
  colnames(output) <- as.character(props)
  
  # Iterate over each value of props and then each value from 1 to nsim
  for (j in seq_len(length(props))) {
    for (i in 1:nsim) {
      # Run code in brackets if missing_probs is NA
      if (is.na(missing_probs)[1]) {
        sub_samp <- sample(seq(1, ecount(net)),
                           size = round(ecount(net) * props[j], 0))
        sub_net <- igraph::delete_edges(net, which(!(seq(1, ecount(net))
                                                     %in% sub_samp)))
      }
      # Run code in brackets if missing_probs contains values
      else {
        sub_samp <- sample(seq(1, ecount(net)), prob = missing_probs,
                           size = round(ecount(net) * props[j], 0))
        sub_net <- igraph::delete_edges(net, which(!(seq(1, ecount(net))
                                                     %in% sub_samp)))
      }
      
      # Select measure of interest based on met and calculate
      if (met == "degree") {
        temp_stats <- igraph::degree(sub_net)
        output[i, j] <- suppressWarnings(cor(temp_stats,
                                             met_orig,
                                             method = "spearman"))
      }
      else   {
        if (met == "betweenness") {
          temp_stats <- igraph::edge_betweenness(sub_net, directed = FALSE, weights = E(sub_net)$timeWeight)  # <<< UPDATED LINE
          kept_edges <- which(seq_len(ecount(net)) %in% sub_samp)  # <<< ADDED LINE
          
          if (length(temp_stats) == length(kept_edges) &&
              length(unique(temp_stats)) >= 2 &&
              length(unique(met_orig[kept_edges])) >= 2) {
            output[i, j] <- suppressWarnings(
              cor(temp_stats, met_orig[kept_edges], method = "spearman")  # <<< UPDATED LINE
            )
          } else {
            output[i, j] <- NA  # <<< ADDED LINE
          }
        }
      }
    }
  }
  
  # Return output as data.frame
  df_output <- suppressWarnings(melt(as.data.frame(output)))
  return(df_output)
}

##Node degree robustness
set.seed(99)
node_degree_rob <- sim_missing_edges(net = road_network_1, met = "degree")

##Plot the results
ggplot(data = node_degree_rob) +
  geom_boxplot(aes(x = variable, y = value)) +
  xlab("Sub-Sample Size as Proportion of Original") +
  ylab(expression("Spearman's" ~ rho)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = rel(1)),
    axis.text.y = element_text(size = rel(1)),
    axis.title.x = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1)),
    legend.text = element_text(size = rel(1))
  )

ggsave(filename = "node_degree_spearmansrho.tiff", path = "output/", device = "tiff", dpi = 300)

##Time-weighted edge betweenness robustness
set.seed(99)
edge_btw_rob <- sim_missing_edges(net = road_network_1, met = "betweenness")

##Plot the results
ggplot(data = edge_btw_rob) +
  geom_boxplot(aes(x = variable, y = value)) +
  xlab("Sub-Sample Size as Proportion of Original") +
  ylab(expression("Spearman's" ~ rho)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = rel(1)),
    axis.text.y = element_text(size = rel(1)),
    axis.title.x = element_text(size = rel(1)),
    axis.title.y = element_text(size = rel(1)),
    legend.text = element_text(size = rel(1))
  )

ggsave(filename = "edge_betweeness_spearmansrho.tiff", path = "output/", device = "tiff", dpi = 300)

###Community detection
##Louvain algorithm
community_louvain <- cluster_louvain(road_network_1, weights = NA, resolution = 1)

#Add as a node attribute
community_louvain_vec <- membership(community_louvain)
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(clusLou = community_louvain_vec)

##Edge betweenness method (Girvan-Newman algorithm)
community_edbtw <- cluster_edge_betweenness(road_network_1, weights = road_network_1$timeWeight)

#Add as a node attribute
community_edgbtw_vec <- membership(community_edbtw)
road_network_2_sf <- road_network_2_sf %>%
  activate("nodes") %>%
  mutate(clusEdg = community_edgbtw_vec)

#######################################################################
####Least cost path networks modelling and comparison (southern Levant)
#######################################################################

###Create network from the Roman road data for the southern case study
#Add data
south_sites <- st_read("data/levant_sites_south.shp")
south_roads <- st_read("data/levant_roads_south.shp")
south_roads <- st_zm(south_roads, drop=TRUE, what = "ZM")

start_points <- st_sfc(lapply(st_geometry(south_roads), get_endpoints, "start"), crs = st_crs(south_roads))
end_points   <- st_sfc(lapply(st_geometry(south_roads), get_endpoints, "end"), crs = st_crs(south_roads))

start_sf <- st_sf(geometry = start_points)
end_sf   <- st_sf(geometry = end_points)

start_nearest <- st_nearest_feature(start_sf, south_sites)
end_nearest   <- st_nearest_feature(end_sf, south_sites)

south_roads$from_id <- south_sites$idAll[start_nearest]
south_roads$to_id   <- south_sites$idAll[end_nearest]

#Convert idAll, from_id and to_id fields to character. They are currently numeric type which sfnetworks does not accept. Also rename old 'id' field to 'PleiadesId'
south_roads$from_id <- as.character(south_roads$from_id)
south_roads$to_id   <- as.character(south_roads$to_id)
south_sites$idAll   <- as.character(south_sites$idAll)
colnames(south_sites)[5] <- "PleiadesID"

#Reorder columns in roads
south_roads_ordered <- south_roads[, c("from_id", "to_id", "lengthGeo", "Type", "typeWeight", "Avg_Slope", "pace", "timeWeight", "geometry")]

#Build the network
south_road_network <- as_sfnetwork(x = south_sites, edges = south_roads_ordered, node_key = "idAll", from = "from_id", to = "to_id", directed = FALSE, edges_as_lines = TRUE, length_as_weight = FALSE)

#As igraph
south_road_network_ig <- as.igraph(south_road_network)

##Get edge list
#Get vertex attributes (especially idAll)
node_ids <- south_road_network %>%
  activate("nodes") %>%
  as_tibble() %>%
  pull(idAll)

#Extract the edge list using vertex indices
edge_df <- as_data_frame(south_road_network_ig, what = "edges")

#Add from_id and to_id by indexing into node_ids
edge_df$from_id <- node_ids[as.numeric(edge_df$from)]
edge_df$to_id   <- node_ids[as.numeric(edge_df$to)]

#Join with original edge attributes
edge_attributes <- south_road_network %>%
  activate("edges") %>%
  as_tibble()

#Combine from_id and to_id
south_edge_list <- bind_cols(
  edge_df[, c("from_id", "to_id")],
  edge_attributes
)

###Create cost distance matrix
#Convert sites to sp object (for gdistance)
south_sites_sp <- as(south_sites, "Spatial")

#Add conductivity surface (for gdistance, NA values for not passable)
cs_na <- raster::raster("data/south_na.tif")

# Create transition layer using mean(1/x) function since our raster layer values represent conductance, raster must have NA values where not passable (water)
tr <- transition(cs_na, transitionFunction = mean, directions = 16)
tr <- geoCorrection(tr, type = "c")

#Calculate cost distance matrix
cd <- gdistance::costDistance(tr, south_sites_sp)
cd_matrix <- as.matrix(cd)

###Create Gabriel Graph
gg_south <- cccd::gg(cd_matrix, r = 1)

#Extract edge list
gg_south_edges <- as_data_frame(gg_south, what = "edges")

#Map indices to idAll
id_map <- south_sites$idAll

gg_south_edges$from_id <- id_map[as.integer(gg_south_edges$from)]
gg_south_edges$to_id   <- id_map[as.integer(gg_south_edges$to)]

#Get GG edge list where from_id, to _id corresponds to idAll
gg_south_edge_list <- gg_south_edges %>%
  dplyr::select(from_id, to_id)

###Create Relative Neighbourhood Graph
rng_south <- cccd::rng(cd_matrix, r = 1)

#Extract edge list
rng_south_edges <- as_data_frame(rng_south, what = "edges")

#Get RNG edge list
rng_south_edges$from_id <- id_map[as.integer(rng_south_edges$from)]
rng_south_edges$to_id   <- id_map[as.integer(rng_south_edges$to)]

rng_south_edge_list <- rng_south_edges %>%
  dplyr::select(from_id, to_id)

###Create K=4 Nearest Neighbour Graph
k4_south <- cccd::nng(cd_matrix, k = 4, mutual = FALSE)

#Extract edge list
k4_south_edges <- as_data_frame(k4_south, what = "edges")

#Get K4 edge list
k4_south_edges$from_id <- id_map[as.integer(k4_south_edges$from)]
k4_south_edges$to_id   <- id_map[as.integer(k4_south_edges$to)]

k4_south_edge_list <- k4_south_edges %>%
  dplyr::select(from_id, to_id)

###Create K=5 Nearest Neighbour Graph
k5_south <- cccd::nng(cd_matrix, k = 5, mutual = FALSE)

#Extract edge list
k5_south_edges <- as_data_frame(k5_south, what = "edges")

#Get K4 edge list
k5_south_edges$from_id <- id_map[as.integer(k5_south_edges$from)]
k5_south_edges$to_id   <- id_map[as.integer(k5_south_edges$to)]

k5_south_edge_list <- k5_south_edges %>%
  dplyr::select(from_id, to_id)

###Compare edge lists (common edges, only in road network, only in modeled networks)
##Create character keys from edge lists
#Normalize edge ordering
sort_edges <- function(df) {
  df %>%
    mutate(
      edge_min = pmin(from_id, to_id),
      edge_max = pmax(from_id, to_id)
    ) %>%
    dplyr::select(edge_min, edge_max)
}

south_roads_edges_key <- sort_edges(south_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

gg_edges_key <- sort_edges(gg_south_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

rng_edges_key <- sort_edges(rng_south_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

k4_edges_key <- sort_edges(k4_south_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

k5_edges_key <- sort_edges(k5_south_edge_list) %>%
  mutate(key = paste(edge_min, edge_max, sep = "-")) %>%
  pull(key)

##GG lists
common_roads_gg <- intersect(south_roads_edges_key, gg_edges_key)
only_roads_gg <- setdiff(south_roads_edges_key, gg_edges_key)
only_gg_roads <- setdiff(gg_edges_key, south_roads_edges_key)

#To data frames
common_roads_gg_df <- do.call(rbind, strsplit(common_roads_gg, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_roads_gg_df <- do.call(rbind, strsplit(only_roads_gg, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_gg_roads_df <- do.call(rbind, strsplit(only_gg_roads, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

##RNG lists
common_roads_rng <- intersect(south_roads_edges_key, rng_edges_key)
only_roads_rng <- setdiff(south_roads_edges_key, rng_edges_key)
only_rng_roads <- setdiff(rng_edges_key, south_roads_edges_key)

#To data frames
common_roads_rng_df <- do.call(rbind, strsplit(common_roads_rng, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_roads_rng_df <- do.call(rbind, strsplit(only_roads_rng, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_rng_roads_df <- do.call(rbind, strsplit(only_rng_roads, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

##K=4 NN lists
common_roads_k4 <- intersect(south_roads_edges_key, k4_edges_key)
only_roads_k4 <- setdiff(south_roads_edges_key, k4_edges_key)
only_k4_roads <- setdiff(k4_edges_key, south_roads_edges_key)

#To data frames
common_roads_k4_df <- do.call(rbind, strsplit(common_roads_k4, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_roads_k4_df <- do.call(rbind, strsplit(only_roads_k4, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_k4_roads_df <- do.call(rbind, strsplit(only_k4_roads, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

##K=5 NN lists
common_roads_k5 <- intersect(south_roads_edges_key, k5_edges_key)
only_roads_k5 <- setdiff(south_roads_edges_key, k5_edges_key)
only_k5_roads <- setdiff(k5_edges_key, south_roads_edges_key)

#To data frames
common_roads_k5_df <- do.call(rbind, strsplit(common_roads_k5, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_roads_k5_df <- do.call(rbind, strsplit(only_roads_k5, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

only_k5_roads_df <- do.call(rbind, strsplit(only_k5_roads, "-")) %>%
  as.data.frame() %>%
  setNames(c("from_id", "to_id"))

##Add type categories
common_roads_gg_df$type <- "RoadsGG"
only_roads_gg_df$type <- "RoadsNoGG"
only_gg_roads_df$type <- "GGNoRoads"

common_roads_rng_df$type <- "RoadsRNG"
only_roads_rng_df$type <- "RoadsNoRNG"
only_rng_roads_df$type <- "RNGNoRoads"

common_roads_k4_df$type <- "RoadsK4"
only_roads_k4_df$type <- "RoadsNoK4"
only_k4_roads_df$type <- "K4NoRoads"

common_roads_k5_df$type <- "RoadsK5"
only_roads_k5_df$type <- "RoadsNoK5"
only_k5_roads_df$type <- "K5NoRoads"

#Add to south_edge_list as boolean
south_edge_list <- south_edge_list %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id),
    key = paste(edge_min, edge_max, sep = "-"),
    commonGG  = key %in% gg_edges_key,
    commonRNG = key %in% rng_edges_key,
    commonK4  = key %in% k4_edges_key,
    commonK5 = key %in% k5_edges_key
  ) %>%
  dplyr::select(-edge_min, -edge_max, -key)

##Export
south_roads_edges <- st_as_sf(south_edge_list)
st_write(south_roads_edges, "output/south_roads_edges.shp")

##Export constructed graphs as straight line sf objects
#Extract point coordinates into a data frame
sites_coord <- sites %>%
  mutate(
    idAll = as.character(idAll),
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

#Extract edges and add coordinates GG
gg_south_coord <- gg_south_edge_list %>%
  left_join(sites_coord, by = c("from_id" = "idAll")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sites_coord, by = c("to_id" = "idAll")) %>%
  rename(x_to = x, y_to = y)

gg_south_sf <- gg_south_coord %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2, byrow = FALSE)),
      crs = st_crs(sites)
    )
  ) %>%
  ungroup() %>%
  st_as_sf()

st_write(gg_south_sf, "output/gg_south.shp")

#Extract edges and add coordinates RNG
rng_south_coord <- rng_south_edge_list %>%
  left_join(sites_coord, by = c("from_id" = "idAll")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sites_coord, by = c("to_id" = "idAll")) %>%
  rename(x_to = x, y_to = y)

rng_south_sf <- rng_south_coord %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2, byrow = FALSE)),
      crs = st_crs(sites)
    )
  ) %>%
  ungroup() %>%
  st_as_sf()

st_write(rng_south_sf, "output/rng_south.shp")

#Extract edges and add coordinates K4
k4_south_coord <- k4_south_edge_list %>%
  left_join(sites_coord, by = c("from_id" = "idAll")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sites_coord, by = c("to_id" = "idAll")) %>%
  rename(x_to = x, y_to = y)

k4_south_sf <- k4_south_coord %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2, byrow = FALSE)),
      crs = st_crs(sites)
    )
  ) %>%
  ungroup() %>%
  st_as_sf()

st_write(k4_south_sf, "output/k4_south.shp")

#Extract edges and add coordinates K4
k5_south_coord <- k5_south_edge_list %>%
  left_join(sites_coord, by = c("from_id" = "idAll")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(sites_coord, by = c("to_id" = "idAll")) %>%
  rename(x_to = x, y_to = y)

k5_south_sf <- k5_south_coord %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2, byrow = FALSE)),
      crs = st_crs(sites)
    )
  ) %>%
  ungroup() %>%
  st_as_sf()

st_write(k5_south_sf, "output/k5_south.shp")

###Create LCP networks

#Add conductivity surface (for lestcostpath package)
south_cs90 <- terra::rast("data/levant_conductance_90_south.tif")

#Create conductivity surface with leastcostpath
south_cs <- leastcostpath::create_cs(south_cs90, neighbours = 16)

#Create a function for calculating LCP by filtering the south_sites row by row with origin and destinations in edge lists
calculate_lcp <- function(row, x, sites) {
  # Get origin and destination points by matching from_id and to_id to idAll
  origin <- south_sites %>% filter(idAll == row$from_id)
  destination <- south_sites %>% filter(idAll == row$to_id)
  #Calculate LCP
  lcp <- create_lcp(x, origin, destination)
  # Add from_id to_id to the result
  lcp$from_id <- row$from_id
  lcp$to_id <- row$to_id
  
  return(lcp)
}

##LCP network for GG
gg_lcps_list <- list()

batch_size <- 1
indices <- seq_len(nrow(gg_south_edge_list))
batches <- split(indices, ceiling(indices / batch_size))

for (i in seq_along(batches)) {
  cat("Processing batch", i, "of", length(batches), "\n")
  
  batch_results <- map(batches[[i]], function(j) {
    tryCatch({
      calculate_lcp(gg_south_edge_list[j, ], x = south_cs, sites = south_sites)
    }, error = function(e) {
      message("Error in row ", j, ": ", e$message)
      return(NULL)
    })
  })
  
  # Store results and clean memory
  gg_lcps_list <- c(gg_lcps_list, batch_results)
  gc(verbose = FALSE)
}

gg_lcps <- do.call(rbind, gg_lcps_list)
st_write(gg_lcps, "output/gg_lcps.shp")

##LCP network for RNG
#Filter GG LCPs (RNG is a subgraph of GG)
gg_lcps$edge_key <- gg_edges_key

rng_lcps <- gg_lcps %>%
  filter(edge_key %in% rng_edges_key)

st_write(rng_lcps, "output/rng_lcps.shp")

##LCP network for K=4
#Filter GG LCPs for LCPs contained in K=4
k4_gg_lcps_common <- gg_lcps %>%
  filter(edge_key %in% k4_edges_key)

#Create edge key unique to K=4 NN (i.e., excluding edges already contained in GG). K=4 NN contains duplicate edges, so the number of unique edges is lower than would be expected just from looking at the numbers.
k4_gg_common_key <- k4_gg_lcps_common %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id),
    edge_key = paste(edge_min, edge_max, sep = "-")
  ) %>%
  pull(edge_key)

k4_unique_edges <- setdiff(k4_edges_key, k4_gg_common_key)

k4_unique_edges_df <- data.frame(edge_key = k4_unique_edges) %>%
  tidyr::separate(edge_key, into = c("from_id", "to_id"), sep = "-", convert = TRUE)

#Calculate K=4 NN unique LCPs
k4_unique_lcps_list <- list()

batch_size <- 1
indices <- seq_len(nrow(k4_unique_edges_df))
batches <- split(indices, ceiling(indices / batch_size))

for (i in seq_along(batches)) {
  cat("Processing batch", i, "of", length(batches), "\n")
  
  batch_results <- map(batches[[i]], function(j) {
    tryCatch({
      calculate_lcp(k4_unique_edges_df[j, ], x = south_cs, sites = south_sites)
    }, error = function(e) {
      message("Error in row ", j, ": ", e$message)
      return(NULL)
    })
  })
  
  # Store results and clean memory
  k4_unique_lcps_list <- c(k4_unique_lcps_list, batch_results)
  gc(verbose = FALSE)
}

k4_unique_lcps <- do.call(rbind, k4_unique_lcps_list)

k4_gg_lcps_common <- k4_gg_lcps_common %>%
  select(-edge_key)

k4_lcps <- rbind(k4_unique_lcps, k4_gg_lcps_common)

st_write(k4_lcps, "output/k4_lcps.shp")

###Road network and LCP networks comparison
##NPDI validation

#RNG NPDI
#Normalize edge directions (from_id and to_id are often switched between the two compared sf objects)
common_roads_rng_df <- common_roads_rng_df %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )

rng_lcps <- rng_lcps %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )

south_roads <- south_roads %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )

#Initialize RNG list to store the results
rng_PDIs <- list()

#Loop through the sf objects using normalized order of the from_id and to_id columns while retaining the original column order and enforcing consistent naming of columns
desired_cols <- c("area", "pdi", "max_distance", "normalised_pdi", "geometry", "from_id", "to_id")

standardize_cols <- function(df, desired_cols) {
  #Add missing columns as NA
  missing <- setdiff(desired_cols, names(df))
  for (col in missing) {
    if (col == "geometry") {
      #Add empty geometry (compatible with sf)
      df[[col]] <- st_sfc(lapply(1:nrow(df), function(x) st_geometrycollection()), crs = st_crs(df))
    } else {
      df[[col]] <- NA
    }
  }
  
  #Drop extra columns not in desired_cols
  df <- df[, intersect(desired_cols, names(df)), drop = FALSE]
  
  #Reorder to desired_cols order (some columns may still be missing so use intersect again)
  df <- df[, intersect(desired_cols, names(df))]
  
  #If geometry column exists, ensure it's class "sfc"
  if ("geometry" %in% names(df)) {
    if (!inherits(df$geometry, "sfc")) {
      df$geometry <- st_sfc(df$geometry, crs = st_crs(df))
    }
  }
  
  return(df)
}

for (i in 1:nrow(common_roads_rng_df)) {
  from_id <- common_roads_rng_df$from_id[i]
  to_id <- common_roads_rng_df$to_id[i]
  origin_ <- common_roads_rng_df$edge_min[i]
  destination <- common_roads_rng_df$edge_max[i]
  
  subset_rng <- rng_lcps %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  subset_south_roads <- south_roads %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  valid_input <- (
    nrow(subset_rng) > 0 &&
      nrow(subset_south_roads) > 0 &&
      all(st_is_valid(subset_rng)) &&
      all(st_is_valid(subset_south_roads)) &&
      length(st_geometry(subset_rng)) > 0 &&
      length(st_geometry(subset_south_roads)) > 0
  )
  
  if (valid_input) {
    tryCatch({
      rng_pdi_results <- leastcostpath::PDI_validation(
        lcp = subset_rng,
        comparison = subset_south_roads
      )
      
      #Add original IDs to the result
      rng_pdi_results$from_id <- from_id
      rng_pdi_results$to_id <- to_id
      
      #Standardize columns and reorder
      rng_pdi_results <- standardize_cols(rng_pdi_results, desired_cols)
      
      rng_PDIs[[paste(from_id, to_id, sep = "_")]] <- rng_pdi_results
    }, error = function(e) {
      message("Error in PDI_validation for ", from_id, "-", to_id, ": ", e$message)
    })
  } else {
    message("Skipping ", from_id, "-", to_id, ": missing or invalid input")
  }
}

rng_PDI_validation <- do.call(rbind, rng_PDIs)
sf::st_write(rng_PDI_validation, "output/rng_PDI_validation.shp")

#GG NPDI validation
common_roads_gg_df <- common_roads_gg_df %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )

gg_lcps <- gg_lcps %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )


gg_PDIs <- list()

for (i in 1:nrow(common_roads_gg_df)) {
  from_id <- common_roads_gg_df$from_id[i]
  to_id <- common_roads_gg_df$to_id[i]
  origin_ <- common_roads_gg_df$edge_min[i]
  destination <- common_roads_gg_df$edge_max[i]
  
  subset_gg <- gg_lcps %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  subset_south_roads <- south_roads %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  valid_input <- (
    nrow(subset_gg) > 0 &&
      nrow(subset_south_roads) > 0 &&
      all(st_is_valid(subset_gg)) &&
      all(st_is_valid(subset_south_roads)) &&
      length(st_geometry(subset_gg)) > 0 &&
      length(st_geometry(subset_south_roads)) > 0
  )
  
  if (valid_input) {
    tryCatch({
      gg_pdi_results <- leastcostpath::PDI_validation(
        lcp = subset_gg,
        comparison = subset_south_roads
      )
      
      #Add original IDs to the result
      gg_pdi_results$from_id <- from_id
      gg_pdi_results$to_id <- to_id
      
      #Standardize columns and reorder
      gg_pdi_results <- standardize_cols(gg_pdi_results, desired_cols)
      
      gg_PDIs[[paste(from_id, to_id, sep = "_")]] <- gg_pdi_results
    }, error = function(e) {
      message("Error in PDI_validation for ", from_id, "-", to_id, ": ", e$message)
    })
  } else {
    message("Skipping ", from_id, "-", to_id, ": missing or invalid input")
  }
}

gg_PDI_validation <- do.call(rbind, gg_PDIs)
sf::st_write(gg_PDI_validation, "output/gg_PDI_validation.shp")

#K=4 NN NPDI
common_roads_k4_df <- common_roads_k4_df %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )

k4_lcps <- k4_lcps %>%
  mutate(
    edge_min = pmin(from_id, to_id),
    edge_max = pmax(from_id, to_id)
  )


k4_PDIs <- list()

for (i in 1:nrow(common_roads_k4_df)) {
  from_id <- common_roads_k4_df$from_id[i]
  to_id <- common_roads_k4_df$to_id[i]
  origin_ <- common_roads_k4_df$edge_min[i]
  destination <- common_roads_k4_df$edge_max[i]
  
  subset_k4 <- k4_lcps %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  subset_south_roads <- south_roads %>%
    filter(edge_min == origin_ & edge_max == destination)
  
  valid_input <- (
    nrow(subset_k4) > 0 &&
      nrow(subset_south_roads) > 0 &&
      all(st_is_valid(subset_k4)) &&
      all(st_is_valid(subset_south_roads)) &&
      length(st_geometry(subset_k4)) > 0 &&
      length(st_geometry(subset_south_roads)) > 0
  )
  
  if (valid_input) {
    tryCatch({
      k4_pdi_results <- leastcostpath::PDI_validation(
        lcp = subset_k4,
        comparison = subset_south_roads
      )
      
      #Add original IDs to the result
      k4_pdi_results$from_id <- from_id
      k4_pdi_results$to_id <- to_id
      
      #Standardize columns and reorder
      k4_pdi_results <- standardize_cols(k4_pdi_results, desired_cols)
      
      k4_PDIs[[paste(from_id, to_id, sep = "_")]] <- k4_pdi_results
    }, error = function(e) {
      message("Error in PDI_validation for ", from_id, "-", to_id, ": ", e$message)
    })
  } else {
    message("Skipping ", from_id, "-", to_id, ": missing or invalid input")
  }
}

k4_PDI_validation <- do.call(rbind, k4_PDIs)
sf::st_write(k4_PDI_validation, "output/k4_PDI_validation.shp")

##Table comparing mean, median, min, max NPDI
#Distinct npdi column and transform to data frame
rng_npdi <- rng_PDI_validation %>%
  st_drop_geometry() %>%
  rename(npdi = normalised_pdi) %>%
  dplyr::select(npdi) %>%
  distinct()

gg_npdi <- gg_PDI_validation %>%
  st_drop_geometry() %>%
  rename(npdi = normalised_pdi) %>%
  dplyr::select(npdi) %>%
  distinct()

k4_npdi <- k4_PDI_validation %>%
  st_drop_geometry() %>%
  rename(npdi = normalised_pdi) %>%
  dplyr::select(npdi) %>%
  distinct()

#Calculate values
rng_npdi_stats <- rng_npdi %>%
  summarise( 
            mean = mean(npdi, na.rm = TRUE),
            median = median(npdi, na.rm = TRUE),
            min = min(npdi, na.rm = TRUE),
            max = max(npdi, na.rm = TRUE))

gg_npdi_stats <- gg_npdi %>%
  summarise( 
    mean = mean(npdi, na.rm = TRUE),
    median = median(npdi, na.rm = TRUE),
    min = min(npdi, na.rm = TRUE),
    max = max(npdi, na.rm = TRUE))

k4_npdi_stats <- k4_npdi %>%
  summarise( 
    mean = mean(npdi, na.rm = TRUE),
    median = median(npdi, na.rm = TRUE),
    min = min(npdi, na.rm = TRUE),
    max = max(npdi, na.rm = TRUE))

#To table
npdi_table <- bind_rows(rng_npdi_stats, gg_npdi_stats, k4_npdi_stats)
rownames(npdi_table) <- c("RNG", "GG", "K4")

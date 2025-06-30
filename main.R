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
###We compare the network with theoretical planar networks (Gabriel Graph, Relative Neighbourhood Graph, Minimum Spanning Tree, Delaunay Triangulation)

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

##Add all to list
constructed_networks <- list(gg = gg_c_igraph, rng = rng_c_igraph, mst = mst_c_igraph, dt = dt_igraph_cleaned)

##Calculate properties of the constructed networks
#Number of edges
number_of_edges_df <- constructed_networks %>% 
  lapply(ecount) %>%
  stack() %>%
  dplyr::select(ind, values) %>%
  setNames(c("network","number_of_edges"))

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

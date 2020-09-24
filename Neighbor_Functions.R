require(dplyr)
#create functions to generate neighbors
#### Step 1: creating a function for generating nearest neighbors ####
makeMooreNeighbors <- function(coords) {
  
  moore_neighbors <- data.frame()
  ## First do left
  left_neighbors <- data.frame(id = coords$id, 
                               x = coords$x-1,
                               y = coords$y)
  
  moore_neighbors <- rbind(moore_neighbors, left_neighbors)
  
  ## Then do right
  right_neighbors <- data.frame(id = coords$id, 
                                x = coords$x+1,
                                y = coords$y) 
  
  moore_neighbors <- rbind(moore_neighbors, right_neighbors)
  
  ## Above
  above_neighbors <- data.frame(id = coords$id, 
                                x = coords$x,
                                y = coords$y+1) 
  
  moore_neighbors <- rbind(moore_neighbors, above_neighbors)
  
  ## Below
  below_neighbors <- data.frame(id = coords$id, 
                                x = coords$x,
                                y = coords$y-1) 
  
  moore_neighbors <- rbind(moore_neighbors, below_neighbors)
  
  ## Above-left
  above_left_neighbors <- data.frame(id = coords$id, 
                                     x = coords$x-1,
                                     y = coords$y+1) 
  
  moore_neighbors <- rbind(moore_neighbors, above_left_neighbors)
  
  
  ## Above-right
  above_right_neighbors <- data.frame(id = coords$id, 
                                      x = coords$x+1,
                                      y = coords$y+1) 
  
  moore_neighbors <- rbind(moore_neighbors, above_right_neighbors)
  
  ## Below-left
  below_left_neighbors <- data.frame(id = coords$id, 
                                     x = coords$x-1,
                                     y = coords$y-1) 
  
  moore_neighbors <- rbind(moore_neighbors, below_left_neighbors)
  
  ## Below-right
  below_right_neighbors <- data.frame(id = coords$id, 
                                      x = coords$x+1,
                                      y = coords$y-1) 
  
  moore_neighbors <- rbind(moore_neighbors, below_right_neighbors)
  
  moore_neighbors <- moore_neighbors %>%
    dplyr::filter(x >= min(coords$x), 
                  y >= min(coords$y), 
                  y <= max(coords$y), 
                  x <= max(coords$x)) %>%
    dplyr::arrange(id)
  
  
  return(moore_neighbors)
  
}

#### Step 2: Define function that creates the # of neighbors (by population) by modifying ####
# makeMooreNeighbors function. Will use this to generate non-density dependent rho (later)

num_neighbors_adjacent <- function(coords) {
  
  num_neighbors_adj <- data.frame()
  ## First do left
  left_neighbors_adj <- data.frame(id = coords$id, 
                                   x = coords$x-1,
                                   y = coords$y,
                                   pop = coords$population)
  
  num_neighbors_adj <- rbind(num_neighbors_adj, left_neighbors_adj)
  
  ## Then do right
  right_neighbors_adj <- data.frame(id = coords$id, 
                                    x = coords$x+1,
                                    y = coords$y,
                                    pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, right_neighbors_adj)
  
  ## Above
  above_neighbors_adj <- data.frame(id = coords$id, 
                                    x = coords$x,
                                    y = coords$y+1,
                                    pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, above_neighbors_adj)
  
  ## Below
  below_neighbors_adj <- data.frame(id = coords$id, 
                                    x = coords$x,
                                    y = coords$y-1,
                                    pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, below_neighbors_adj)
  
  ## Above-left
  above_left_neighbors_adj <- data.frame(id = coords$id, 
                                         x = coords$x-1,
                                         y = coords$y+1,
                                         pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, above_left_neighbors_adj)
  
  
  ## Above-right
  above_right_neighbors_adj <- data.frame(id = coords$id, 
                                          x = coords$x+1,
                                          y = coords$y+1,
                                          pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, above_right_neighbors_adj)
  
  ## Below-left
  below_left_neighbors_adj <- data.frame(id = coords$id, 
                                         x = coords$x-1,
                                         y = coords$y-1,
                                         pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, below_left_neighbors_adj)
  
  ## Below-right
  below_right_neighbors_adj <- data.frame(id = coords$id, 
                                          x = coords$x+1,
                                          y = coords$y-1,
                                          pop = coords$population) 
  
  num_neighbors_adj <- rbind(num_neighbors_adj, below_right_neighbors_adj)
  
  num_neighbors_adj <- num_neighbors_adj %>%
    dplyr::filter(x >= min(coords$x), 
                  y >= min(coords$y), 
                  y <= max(coords$y), 
                  x <= max(coords$x)) %>%
    dplyr::arrange(id)
  
  
  return(num_neighbors_adj)
  
}
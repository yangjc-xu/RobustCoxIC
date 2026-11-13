## ---------------------------------------------------------------------------
## Helper function shared by the simulation and application scripts.
## ---------------------------------------------------------------------------

## Unique finite interval endpoints, sorted in decreasing order. The second and
## third columns of simpdat are the left and right interval endpoints.
find_time = function(simpdat){
  Time_unique = unique(c(simpdat[which(simpdat[,2] != 0),2], simpdat[which(simpdat[,3] != Inf),3])) %>%
    sort(decreasing = TRUE)
}

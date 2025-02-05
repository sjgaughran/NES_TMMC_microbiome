library(ggplot2) # Use ggplot2 to add layer for visualization
library(ggspatial) # Special ggplot2 package for interacting with spatial data 
library(tigris) # Package that has shape files from US Census Bureau

# Pulling data for all California counties and the six counties in our study
ca_counties <- counties(state = "CA", cb = TRUE, class = "sf")

counties_of_interest = ca_counties[ca_counties$NAME %in% c("Marin", "Monterey", 
                                                           "San Francisco",
                                                           "San Mateo", 
                                                           "Santa Cruz", 
                                                           "San Luis Obispo"),]

# Creating map using ggplot  
ggplot() +
  geom_sf(data=ca_counties, fill = "white") +
  geom_sf(data=counties_of_interest) + # default is filling with gray 
  coord_sf(xlim=c(-124, -119.5), ylim=c(34.75, 38.3)) + # setting lat/long 
  theme_bw() + 
  theme(axis.text=element_text(size=8)) + 
  annotation_scale() + # adding scale bar and north arrow 
  annotation_north_arrow(location = "bl", style = north_arrow_fancy_orienteering, 
                         height = unit(1, "cm"), width = unit(1, "cm"),
                         pad_x = unit(1,"cm"), pad_y = unit(0.75,"cm"))

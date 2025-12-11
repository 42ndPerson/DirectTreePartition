### Information on Dual Graph Provided

All these graphs are connected and planar: 

- t_wy.json: The census tracts of Wyoming, 160 vertices total
- t_ia.json: The census tracts of Iowa, 896 vertices total
- t_tx.json: The census tracts of Texas, 6896 vertices total
- bg_oh.json: The census block groups of Ohio, 9472 vertices total

Each node in each graph has several attributes, but the most relevant are: 

- "boundary_node": whether the census unit lies on the boundary of the state
- "INTPTLAT20": the latitude of a point in the interior of the census unit
- "INTPTLON20": the longitude of a point in the interior of the census unit
- "P0010001": total population of the census unit

(It appears there's also attributes for INTPTLAT, INTPTLON, and POP100 that match the attributes above, not sure why or if they always have the same values)
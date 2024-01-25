sizes <- calculateNodeSizes(edges)
# Create a data frame for nodes with sizes
nodes$size <- sizes * sizes

# Create the forceNetwork plot
ColourScale <- 'd3.scaleOrdinal().domain(["origin", "citation"]).range(["#FF6900", "#694489"]);'
forceNetwork(
         Links = edges,
         Nodes = nodes,
         Source = "source",
         Target = "target",
         NodeID = "name",
         Group = "group",
         Value = "width",
         Nodesize = "size",
         radiusCalculation = "Math.sqrt(d.nodesize)+4",
         zoom = TRUE,
         legend = TRUE,
         bounded = TRUE,
         opacity = 0.9,
         colourScale = JS(ColourScale)
)

## Graph without size if needed
# forceNetwork(Links = edges, Nodes = nodes, 
#          Source = "source",
#          Target = "target",
#          NodeID ="name",
#          Group = "group",
#          Value = "width",
#          opacity = 0.9,
#          zoom = TRUE,
#          colourScale = JS(ColourScale))



library(DiagrammeR)

daily.gr = 
grViz(" digraph {
      graph [layout = dot, rankdir = BT]
      node [shape = oval]
      
      1 [label = 'Phototrophs', style = bold]
      2 [label = 'Mixotrophs']
      3 [label = 'Heterotrophs']
      4 [label = 'N,P', style = bold]
      5 [label = 'PAR']
      6 [label = 'Bacteria']
      7 [label = 'Browning']
      
      4 -> 1 [style = bold]
      4 -> 2
      4 -> 6
      6 -> {2 3}
      7 -> 5 [style = dashed]
      7 -> 6
      5 -> 2 [style = dashed]
      5 -> 1 [style = dotted]
      
}", height = 200)


extreme.gr = 
  grViz(" digraph {
      graph [layout = dot, rankdir = BT]
      node [shape = oval]
      
      1 [label = 'Phototrophs', style = dashed]
      2 [label = 'Mixotrophs', style = bold]
      3 [label = 'Heterotrophs', style = bold]
      4 [label = 'N,P']
      5 [label = 'PAR']
      6 [label = 'Bacteria', style = bold]
      7 [label = 'Browning']
      
      4 -> {1 2}
      4 -> 6
      6 -> {2 3} [style = bold]
      7 -> 5 [style = dashed]
      7 -> 6
      5 -> {1 2} [style = dashed]
      
}", height = 200)

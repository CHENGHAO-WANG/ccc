
library(readxl)
library(httr)
library(dplyr)

url1 <- "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms8866/MediaObjects/41467_2015_BFncomms8866_MOESM611_ESM.xlsx"

GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
df <- read_excel(tf, 2L)
df <- as.data.frame(df)
ramilowski <- df %>% 
  filter(Pair.Evidence=="literature supported") %>% 
  select(Ligand.ApprovedSymbol, Receptor.ApprovedSymbol) %>% 
  rename(ligand = Ligand.ApprovedSymbol, receptor = Receptor.ApprovedSymbol)

usethis::use_data(ramilowski, overwrite = TRUE)

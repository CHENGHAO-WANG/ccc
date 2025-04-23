## code to prepare `omnipathr` dataset goes here

intercell_network <- OmnipathR::intercell_network(ligand_receptor = TRUE, high_confidence = TRUE, simplify = TRUE)

omnipathr <- data.frame(ligand=intercell_network$source_genesymbol,receptor=intercell_network$target_genesymbol)

usethis::use_data(omnipathr, overwrite = TRUE)

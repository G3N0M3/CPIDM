######
# Script that identify the overlapping or unique genes, Immune
######


### need directory "location"

# Input Data --------------------------------------------------------------

# intersection of all 3
inter_all <- intersect(DEG_ImPC, DEG_ImPctr)
inter_all <- intersect(inter_all, DEG_ImCctr)

# intersection of two categories
Pctr_PC <- intersect(DEG_ImPctr, DEG_ImPC)
Pctr_PC <- setdiff(Pctr_PC, inter_all)

Pctr_Cctr <- intersect(DEG_ImPctr, DEG_ImCctr)
Pctr_Cctr <- setdiff(Pctr_Cctr, inter_all)

PC_Cctr <- intersect(DEG_ImPC, DEG_ImCctr)
PC_Cctr <- setdiff(PC_Cctr, inter_all)

# sole elements
flower <- c(inter_all, Pctr_PC, Pctr_Cctr, PC_Cctr)
sole_Pctr <- setdiff(DEG_ImPctr, flower)
sole_PC <- setdiff(DEG_ImPC, flower)
sole_Cctr <- setdiff(DEG_ImCctr, flower)


# Print Results -----------------------------------------------------------

# intersection of all 3
inter_all

# intersection of two categories
"Gzma" %in% Pctr_PC
"Cxcl10" %in% Pctr_Cctr
"Cxcl10" %in% PC_Cctr

# sole elements
"Ifng" %in% sole_Pctr
sole_PC
sole_Cctr

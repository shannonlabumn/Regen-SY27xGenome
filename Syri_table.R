setwd("~/Desktop/Harpreet/")
library(readr)
library(dplyr)
library(tidyr)


# ---- read & quick sanity ----
raw <- read_tsv("XJDY_homo34_svtype_len.tsv", col_names = c("alt","len"), show_col_types = FALSE)
raw %>% count(alt, sort = TRUE)  # see which SV types are present

# ---- bin setup (≥100 bp like your example) ----
bins   <- c(100, 1e3, 1e4, 1e5, 1e6, Inf)
labels <- c("100 to <1 kb","1 to <10 kb","10 to <100 kb","100 to <1 Mb",">1 Mb")

# ---- map display names (merge INVTR into Translocations; INVDP -> Inverted duplicates) ----
map_type <- c(
  "INV"   = "Inversions",
  "TRANS" = "Translocations",
  "INVTR" = "Translocations",
  "INS"   = "Insertions",
  "DEL"   = "Deletions",
  "DUP"   = "Duplicated regions",
  "INVDP" = "Inverted duplicates",
  "TDM"   = "Tandem repeat",
  "CPG"   = "Copy gain in hapB",
  "CPL"   = "Copy lost in hapB"
)

row_order <- c("Inversions","Translocations","Insertions","Deletions",
               "Duplicated regions","Inverted duplicates","Tandem repeat",
               "Copy gain in hapB","Copy lost in hapB")

# ---- build table ----
df <- raw %>%
  mutate(
    disp = recode(alt, !!!map_type),
    bin  = cut(len, breaks = bins, labels = labels, right = FALSE)
  ) %>%
  filter(!is.na(disp), !is.na(bin)) %>%
  mutate(disp = factor(disp, levels = row_order))

wide <- df %>%
  count(disp, bin, name = "count") %>%
  complete(disp = factor(row_order, levels = row_order), bin = labels, fill = list(count = 0)) %>%
  pivot_wider(names_from = bin, values_from = count, values_fill = 0) %>%
  arrange(disp) %>%
  mutate(`Total count` = rowSums(across(all_of(labels)))) %>%
  select(Type = disp, `Total count`, all_of(labels))

# ---- add "All" row ----
all_row <- wide %>%
  summarise(
    Type = "All",
    `Total count` = sum(`Total count`),
    across(all_of(labels), sum)
  )

final_tbl <- bind_rows(wide, all_row)
final_tbl
write.csv(final_tbl, "SV_size_bins_XJDY_hap34.csv", row.names = FALSE, quote=FALSE)

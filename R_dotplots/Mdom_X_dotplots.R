# Mdom_X_dotplots.R 
# makes pairwise dotplots of X (or prior X) scaffolds for primaries of all strains

# source FUNCTIONS file, includes tidyverse, RColorBrewer, plotly packages
source("alignment_dotplot-FUNCTIONS.R")
library(cowplot)

MIN_MATCH = 2000

# strain/assembly labels
MdomMYp_SHORTLAB = "MY pri v8"
MdomM1p_SHORTLAB = "M1 pri v3"
MdomM2p_SHORTLAB = "M2 pri v5b"
MdomM3p_SHORTLAB = "M3 pri v3"
MdomM5p_SHORTLAB = "M5 pri v2"

# scaffolds containing X (Y in MY) sequence
l.Mdom.XY_seqs = list(
  MY = paste0("scaffold_", 6),
  M1 = paste0("scaffold_", 6),
  M2 = paste0("scaffold_", 6),
  M3 = paste0("scaffold_", 6),
  M5 = paste0("scaffold_", 6:8)
)

# AGP files for contig boundaries
l.Mdom.AGP_files = list(
  MY = "v8_aabys_primary.FINAL.agp",
  M1 = "MdomM1_pctg_wAG.review.FINAL.agp",
  M2 = "MdomM2-pctg_keepMdmd.review.FINAL.agp",
  M3 = "MdomM3-pctg_v3.review.FINAL.agp",
  M5 = "M5_v2_primary.FINAL.agp"
)

# alignment data for all pairwise primaries
# MY primary has Y chromosome as scaffold 6, others have "X" (more accurately post-X)
l.c.MY_M1 = make.l.coords("v8_aabys_p-vs-MdomM1-pctg_v3-c2000.coords",
                          "v8_aabys_p.fa.fai", "MdomM1-pctg_v3.fa.fai")
l.c.MY_M2 = make.l.coords("v8_aabys_p-vs-MdomM2-pctg_v5b-c2000.coords",
                          "v8_aabys_p.fa.fai", "MdomM2-pctg_v5b.fa.fai")
l.c.MY_M3 = make.l.coords("v8_aabys_p-vs-MdomM3-pctg_v3-c2000.coords",
                          "v8_aabys_p.fa.fai", "MdomM3-pctg_v3.fa.fai")
l.c.MY_M5 = make.l.coords("v8_aabys_p-vs-M5_v2_primary-c2000.coords",
                          "v8_aabys_p.fa.fai", "M5_v2_primary.fa.fai")
l.c.M1_M2 = make.l.coords("MdomM1-pctg_v3-vs-MdomM2-pctg_v5b-c2000.coords",
                          "MdomM1-pctg_v3.fa.fai", "MdomM2-pctg_v5b.fa.fai")
l.c.M1_M3 = make.l.coords("MdomM1-pctg_v3-vs-MdomM3-pctg_v3-c2000.coords",
                          "MdomM1-pctg_v3.fa.fai", "MdomM3-pctg_v3.fa.fai")
l.c.M1_M5 = make.l.coords("MdomM1-pctg_v3-vs-M5_v2_primary-c2000.coords",
                          "MdomM1-pctg_v3.fa.fai", "M5_v2_primary.fa.fai")
l.c.M2_M3 = make.l.coords("MdomM2-pctg_v5b-vs-MdomM3-pctg_v3-c2000.coords",
                          "MdomM2-pctg_v5b.fa.fai", "MdomM3-pctg_v3.fa.fai")
l.c.M2_M5 = make.l.coords("MdomM2-pctg_v5b-vs-M5_v2_primary-c2000.coords",
                          "MdomM2-pctg_v5b.fa.fai", "M5_v2_primary.fa.fai")
l.c.M3_M5 = make.l.coords("MdomM3-pctg_v3-vs-M5_v2_primary-c2000.coords",
                          "MdomM3-pctg_v3.fa.fai", "M5_v2_primary.fa.fai")

# coordinates for scaffold gaps
l.df.Mdom.XY_gap_coords = lapply(
  seq_along(l.Mdom.AGP_files),
  function(ASMi){
    read.table(l.Mdom.AGP_files[[ASMi]]) %>%
      filter(V1 %in% l.Mdom.XY_seqs[[ASMi]],
             V5 == "U") %>% 
      select(c("V1", "V2"))
  })
names(l.df.Mdom.XY_gap_coords) = names(l.Mdom.AGP_files)
# M5 has multiple scaffolds, need to adjust coordinate for stacked dotplot 
l.df.Mdom.XY_gap_coords$M5 = l.df.Mdom.XY_gap_coords$M5 %>%
  # no gaps in scaffold_8 (thus not in df.XY_gaps) so this way is easiest
  mutate(V2 = if_else(V1 == "scaffold_6", V2, 
                      V2 + l.c.M1_M5$breaksQry$SeqLength[6]))

# X alignment plots, all pairwise and in both directions
# standard axes
p.X.MY_M1 = plot.df.coords(l.c.MY_M1, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M1,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM1p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.MY_M2 = plot.df.coords(l.c.MY_M2, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M2,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM2p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.MY_M3 = plot.df.coords(l.c.MY_M3, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.MY_M5 = plot.df.coords(l.c.MY_M5, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M1_M2 = plot.df.coords(l.c.M1_M2, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M2,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM2p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M1_M3 = plot.df.coords(l.c.M1_M3, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M1_M5 = plot.df.coords(l.c.M1_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M2_M3 = plot.df.coords(l.c.M2_M3, 
                           REFERENCE = l.Mdom.XY_seqs$M2, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomM2p_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M2_M5 = plot.df.coords(l.c.M2_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M2, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM2p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH)
p.X.M3_M5 = plot.df.coords(l.c.M3_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M3, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM3p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 90, MIN_MATCH_LENGTH = MIN_MATCH) + 
              theme(axis.text.x = element_text(angle = 0))
# flipped axes
p.X.M1_MY = plot.df.coords(l.c.MY_M1, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M1,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM1p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M2_MY = plot.df.coords(l.c.MY_M2, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M2,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM2p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M3_MY = plot.df.coords(l.c.MY_M3, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M5_MY = plot.df.coords(l.c.MY_M5, 
                           REFERENCE = l.Mdom.XY_seqs$MY, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomMYp_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M2_M1 = plot.df.coords(l.c.M1_M2, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M2,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM2p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M3_M1 = plot.df.coords(l.c.M1_M3, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M5_M1 = plot.df.coords(l.c.M1_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M1, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM1p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M3_M2 = plot.df.coords(l.c.M2_M3, 
                           REFERENCE = l.Mdom.XY_seqs$M2, QUERY = l.Mdom.XY_seqs$M3,
                           LAB_REF = MdomM2p_SHORTLAB, LAB_QRY = MdomM3p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M5_M2 = plot.df.coords(l.c.M2_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M2, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM2p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)
p.X.M5_M3 = plot.df.coords(l.c.M3_M5, 
                           REFERENCE = l.Mdom.XY_seqs$M3, QUERY = l.Mdom.XY_seqs$M5,
                           LAB_REF = MdomM3p_SHORTLAB, LAB_QRY = MdomM5p_SHORTLAB,
                           MAR_T = 5, MAR_R = 5, MAR_B = 5, MAR_L = 5, 
                           POINTSIZE = 0.1, ALPHA = 1,
                           SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH, FLIP_AXES = T)

# all plots in 5X5 grid (NULL diagonal, remove redundant labels)
p.X.all = plot_grid(
  # ROW 1
  NULL, 
  p.X.M1_MY + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M1$V2), 
  p.X.M2_MY + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2), 
  p.X.M3_MY + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2), 
  p.X.M5_MY + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  # ROW 2
  p.X.MY_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M1$V2),
  NULL,
  p.X.M2_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  p.X.M3_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M5_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  # ROW 3
  p.X.MY_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  p.X.M1_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  NULL,
  p.X.M3_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M5_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  # ROW 4
  p.X.MY_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M1_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M2_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  NULL,
  p.X.M5_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M3$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  # ROW 5
  p.X.MY_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M1_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M2_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M3_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M3$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  NULL,
  nrow = 5, ncol = 5
)

# add labels in empty diagonal and coordinates along left edge
# positioning is finicky, looks ok with pdf at 5x6
p.out = ggdraw() +
  draw_plot(p.X.all, x = 0.05+0.03, y = 0.05, width = 0.9, height = 0.9) +
  draw_label("MY ChrY", x=0.5-0.4*0.9+0.03, y=0.5+0.4*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.MY_M1$breaksRef[l.Mdom.XY_seqs$MY,2]/1e6, 2), " Mb"),
             x=0.5-0.285*0.9+0.03, y=0.95, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.5-0.285*0.9+0.03, y=0.95-0.09*2, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M1 ChrX", x=0.5-0.2*0.9+0.03, y=0.5+0.2*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M1_M2$breaksRef[l.Mdom.XY_seqs$M1,2]/1e6, 2), " Mb"),
             x=0.063+0.03, y=0.95-0.09*2, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.063+0.03, y=0.95-0.09*4, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M2 ChrX", x=0.5+0.03, y=0.5, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M2_M3$breaksRef[l.Mdom.XY_seqs$M2,2]/1e6, 2), " Mb"),
             x=0.063+0.03, y=0.95-0.09*4, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.063+0.03, y=0.95-0.09*6, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M3 ChrX", x=0.5+0.2*0.9+0.03, y=0.5-0.2*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M3_M5$breaksRef[l.Mdom.XY_seqs$M3,2]/1e6, 2), " Mb"),
             x=0.063+0.03, y=0.95-0.09*6, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.063+0.03, y=0.95-0.09*8, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M5 ChrX", x=0.5+0.4*0.9+0.03, y=0.5-0.4*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(sum(l.c.M1_M5$breaksQry[l.Mdom.XY_seqs$M5,2])/1e6, 2), " Mb"),
             x=0.063+0.03, y=0.95-0.09*8, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.063+0.03, y=0.95-0.09*10, size = 8, hjust = 1, vjust = -0.65)

pdf(file = "Mdom_X_dotplots-grid_gaps.pdf", width = 6, height = 5)
p.out
dev.off()

# each pairwise with contig boundaries
pdf(file = "Mdom_X_dotplots-gaps.pdf")
p.X.MY_M1 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M1$V2)
p.X.MY_M2 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M2$V2)
p.X.MY_M3 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M3$V2)
p.X.MY_M5 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M5$V2)
p.X.M1_M2 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M2$V2)
p.X.M1_M3 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M3$V2)
p.X.M1_M5 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M5$V2)
p.X.M2_M3 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M3$V2)
p.X.M2_M5 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M5$V2)
p.X.M3_M5 + 
  geom_vline(xintercept = l.df.Mdom.XY_gap_coords$M3$V2) +
  geom_hline(yintercept = l.df.Mdom.XY_gap_coords$M5$V2)
dev.off()


# all plots in 5X5 grid, opposite diagonal, only lower right side
p.X.lower_diag = plot_grid(
  # ROW 1
  NULL, NULL, NULL, NULL, NULL,
  # ROW 2
  NULL, NULL, NULL, NULL, 
  p.X.MY_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M1$V2),
  # ROW 3
  NULL, NULL, NULL, 
  p.X.M1_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  p.X.MY_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  # ROW 4
  NULL, NULL,
  p.X.M2_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M1_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.MY_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  # ROW 5
  NULL,
  p.X.M3_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M3$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M2_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M1_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.MY_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  
  nrow = 5, ncol = 5
)

# add labels in empty diagonal and coordinates along left edge
# positioning is finicky, looks ok with pdf at 5x6
p.out.lower_diag = ggdraw() +
  draw_plot(p.X.lower_diag, x = 0.05+0.03, y = 0.05, width = 0.9, height = 0.9) +
  draw_label("MY ChrY", x=0.5+0.4*0.9+0.03, y=0.5+0.4*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.MY_M1$breaksRef[l.Mdom.XY_seqs$MY,2]/1e6, 2), " Mb"),
             x=0.97, y=0.95-0.09*2-0.01, size = 8, hjust = 0, vjust = 0, angle = 90) +
  draw_label("0 Mb",
             x=0.5+0.315*0.9+0.04, y=0.95-0.09*2-0.01, size = 8, hjust = 0, vjust = 0, angle = 90) +
  draw_label("M1 ChrX", x=0.5+0.2*0.9+0.03, y=0.5+0.2*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M1_M2$breaksRef[l.Mdom.XY_seqs$M1,2]/1e6, 2), " Mb"),
             x=0.5+0.315*0.9+0.025, y=0.95-0.09*2, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.5+0.315*0.9+0.025, y=0.95-0.09*4, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M2 ChrX", x=0.5+0.03, y=0.5, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M2_M3$breaksRef[l.Mdom.XY_seqs$M2,2]/1e6, 2), " Mb"),
             x=0.5+0.115*0.9+0.025, y=0.95-0.09*4, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.5+0.115*0.9+0.025, y=0.95-0.09*6, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M3 ChrX", x=0.5-0.2*0.9+0.03, y=0.5-0.2*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(l.c.M3_M5$breaksRef[l.Mdom.XY_seqs$M3,2]/1e6, 2), " Mb"),
             x=0.5-0.085*0.9+0.025, y=0.95-0.09*6, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.5-0.085*0.9+0.025, y=0.95-0.09*8, size = 8, hjust = 1, vjust = -0.65) +
  draw_label("M5 ChrX", x=0.5-0.4*0.9+0.03, y=0.5-0.4*0.9, hjust = 0.5, vjust = 0.5) +
  draw_label(paste0(round(sum(l.c.M1_M5$breaksQry[l.Mdom.XY_seqs$M5,2])/1e6, 2), " Mb"),
             x=0.5-0.285*0.9+0.025, y=0.95-0.09*8, size = 8, hjust = 1, vjust = 1.4) +
  draw_label("0 Mb",
             x=0.5-0.285*0.9+0.025, y=0.95-0.09*10, size = 8, hjust = 1, vjust = -0.65)

pdf(file = "Mdom_X_dotplots-lower_diag.pdf", width = 6, height = 5)
p.out.lower_diag
dev.off()


# version with grid scaled to sequence lengths
SEQ_LENGTHS = c(l.c.MY_M1$breaksRef[l.Mdom.XY_seqs$MY,2],
                l.c.M1_M2$breaksRef[l.Mdom.XY_seqs$M1,2],
                l.c.M2_M3$breaksRef[l.Mdom.XY_seqs$M2,2],
                l.c.M3_M5$breaksRef[l.Mdom.XY_seqs$M3,2],
                sum(l.c.M1_M5$breaksQry[l.Mdom.XY_seqs$M5,2]))
# annotating endpoints and labels will be a pain, so placing all into the grid, parameters:
MAIN_LABEL_SIZE = 4
SUB_LABEL_SIZE = 3
LABEL_X = 0
# offset for X position, position is LABEL_X + X_OFFSET*SEQ_LENGTHS[i]/sum(SEQ_LENGTHS)
X_OFFSET = 0.2
LABEL_Y = 0.4
# offsets for second line, need to scale with panel dimensions
# second line is LABEL_Y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[i]/sum(SEQ_LENGTHS)
Y_OFFSET_B = 0.3
Y_OFFSET_M = 0.7

p.X.scaled.lower_diag = plot_grid(
  # ROW 1
  NULL, NULL, NULL, NULL, 
  ggplot(data = data.frame(x = LABEL_X, y = LABEL_Y)) + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[1]/sum(SEQ_LENGTHS),
                  y), label = "MY ChrY", hjust = 0, vjust = 0, 
              size = MAIN_LABEL_SIZE, fontface = "bold") + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[1]/sum(SEQ_LENGTHS), 
                  y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[1]/sum(SEQ_LENGTHS)
              ), 
              label = paste0("(", round(SEQ_LENGTHS[1]/1e6, 2), " Mb)"), 
              hjust = 0, vjust = 0, size = SUB_LABEL_SIZE, fontface = "plain") + 
    xlim(0,1) + ylim(0,1) + theme_void(),
  # ROW 2
  NULL, NULL, NULL, 
  ggplot(data = data.frame(x = LABEL_X, y = LABEL_Y)) + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[2]/sum(SEQ_LENGTHS),
                  y), label = "M1 ChrX", hjust = 0, vjust = 0, 
              size = MAIN_LABEL_SIZE, fontface = "bold") + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[2]/sum(SEQ_LENGTHS), 
                  y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[2]/sum(SEQ_LENGTHS)
              ), 
              label = paste0("(", round(SEQ_LENGTHS[2]/1e6, 2), " Mb)"), 
              hjust = 0, vjust = 0, size = SUB_LABEL_SIZE, fontface = "plain") + 
    xlim(0,1) + ylim(0,1) + theme_void(), 
  p.X.MY_M1 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M1$V2),
  # ROW 3
  NULL, NULL, 
  ggplot(data = data.frame(x = LABEL_X, y = LABEL_Y)) + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[3]/sum(SEQ_LENGTHS),
                  y), label = "M2 ChrX", hjust = 0, vjust = 0, 
              size = MAIN_LABEL_SIZE, fontface = "bold") + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[3]/sum(SEQ_LENGTHS), 
                  y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[3]/sum(SEQ_LENGTHS)
              ), 
              label = paste0("(", round(SEQ_LENGTHS[3]/1e6, 2), " Mb)"), 
              hjust = 0, vjust = 0, size = SUB_LABEL_SIZE, fontface = "plain") + 
    xlim(0,1) + ylim(0,1) + theme_void(), 
  p.X.M1_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  p.X.MY_M2 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M2$V2),
  # ROW 4
  NULL, 
  ggplot(data = data.frame(x = LABEL_X, y = LABEL_Y)) + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[4]/sum(SEQ_LENGTHS),
                  y), label = "M3 ChrX", hjust = 0, vjust = 0, 
              size = MAIN_LABEL_SIZE, fontface = "bold") + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[4]/sum(SEQ_LENGTHS),  
                  y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[4]/sum(SEQ_LENGTHS)
              ), 
              label = paste0("(", round(SEQ_LENGTHS[4]/1e6, 2), " Mb)"), 
              hjust = 0, vjust = 0, size = SUB_LABEL_SIZE, fontface = "plain") + 
    xlim(0,1) + ylim(0,1) + theme_void(),
  p.X.M2_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.M1_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  p.X.MY_M3 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M3$V2),
  # ROW 5
  ggplot(data = data.frame(x = LABEL_X, y = LABEL_Y)) + 
    # shift X for last row (rightmost column)
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[5]/sum(SEQ_LENGTHS) + 0.3,
                  y), label = "M5 ChrX", hjust = 0, vjust = 0, 
              size = MAIN_LABEL_SIZE, fontface = "bold") + 
    geom_text(aes(x + X_OFFSET*SEQ_LENGTHS[5]/sum(SEQ_LENGTHS) + 0.3,  
                  y - Y_OFFSET_B + Y_OFFSET_M*SEQ_LENGTHS[5]/sum(SEQ_LENGTHS)
              ), 
              label = paste0("(", round(SEQ_LENGTHS[5]/1e6, 2), " Mb)"), 
              hjust = 0, vjust = 0, size = SUB_LABEL_SIZE, fontface = "plain") + 
    xlim(0,1) + ylim(0,1) + theme_void(),
  p.X.M3_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M3$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M2_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M2$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.M1_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$M1$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  p.X.MY_M5 + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), 
                    axis.title.x = element_blank(), axis.title.y.right = element_blank()) + 
    geom_vline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", xintercept = l.df.Mdom.XY_gap_coords$MY$V2) +
    geom_hline(alpha = 0.7, linewidth = 0.1, linetype = "dashed", yintercept = l.df.Mdom.XY_gap_coords$M5$V2),
  
  nrow = 5, ncol = 5, 
  rel_widths = rev(SEQ_LENGTHS), rel_heights = SEQ_LENGTHS
)

pdf(file = "Mdom_X_dotplots-scaled-lower_diag.pdf", width = 6, height = 5)
p.X.scaled.lower_diag
dev.off()





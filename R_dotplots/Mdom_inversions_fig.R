# Mdom_inversions_fig.R plots inverted regions between Mdom strains, with coord annotation
# inversion between M2 hap1 and hap2 on Chr2
# inversion between MY hap1 and hap2 on Chr5 (MY aka abyss)
# inversions on Chr1 of M5 pri vs other strains

library(scales) # for comma() formatting on coordinate labels

# source FUNCTIONS file, includes tidyverse, RColorBrewer, plotly packages
source("alignment_dotplot-FUNCTIONS.R")

MIN_MATCH = 2000

# alignment data for M2 haps
# without resolved Chr1 inversion (both copies in h1)
l.c.M2.h1v2_h2v3 = make.l.coords("MdomM2h1v2-vs-MdomM2h2v3-c2000.coords",
                              "MdomM2h1v2.fa.fai", "MdomM2h2v3.fa.fai")
MdomM2h1v2_LAB = "M2 strain, haplotype 1 v2"
MdomM2h2v3_LAB = "M2 strain, haplotype 2 v3"
# Chr1 inversion phasing resolved
l.c.M2.h1v3_h2v4 = make.l.coords("MdomM2-hap1_v3-vs-MdomM2-hap2_v4-c2000.coords",
                                 "MdomM2-hap1_v3.fa.fai", "MdomM2-hap2_v4.fa.fai")
MdomM2h1v3_LAB = "M2 strain, haplotype 1 v3"
MdomM2h2v4_LAB = "M2 strain, haplotype 2 v4"


# alignment data for aabys (MY) haps
l.c.MYh1_MYh2 = make.l.coords("v7_aabys_h1-vs-v7_aabys_h2-c2000.coords",
                              "v7_aabys_h1.fa.fai", "v7_aabys_h2.fa.fai")
MdomMYh1_LAB = "MY (aabys) strain, haplotype 1"
MdomMYh2_LAB = "MY (aabys) strain, haplotype 2"

# alignment data for M5 vs other strains (primaries plus M2 haps)
l.c.M1_M5 = make.l.coords("MdomM1-pctg_v3-vs-M5_v2_primary-c2000.coords",
                          "MdomM1-pctg_v3.fa.fai", "M5_v2_primary.fa.fai")
l.c.M2v4_M5 = make.l.coords("MdomM2-pctg_v4-vs-M5_v2_primary-c2000.coords",
                          "MdomM2-pctg_v4.fa.fai", "M5_v2_primary.fa.fai")
l.c.M2v5b_M5 = make.l.coords("MdomM2-pctg_v5b-vs-M5_v2_primary-c2000.coords",
                            "MdomM2-pctg_v5b.fa.fai", "M5_v2_primary.fa.fai")
l.c.M2h1v3_M5 = make.l.coords("MdomM2-hap1_v3-vs-M5_v2_primary-c2000.coords",
                                 "MdomM2-hap1_v3.fa.fai", "M5_v2_primary.fa.fai")
l.c.M2h2v4_M5 = make.l.coords("MdomM2-hap2_v4-vs-M5_v2_primary-c2000.coords",
                              "MdomM2-hap2_v4.fa.fai", "M5_v2_primary.fa.fai")
l.c.M3_M5 = make.l.coords("MdomM3-pctg_v3-vs-M5_v2_primary-c2000.coords",
                          "MdomM3-pctg_v3.fa.fai", "M5_v2_primary.fa.fai")
l.c.MY_M5 = make.l.coords("v8_aabys_p-vs-M5_v2_primary-c2000.coords",
                          "v8_aabys_p.fa.fai", "M5_v2_primary.fa.fai")
MdomM1p_LAB = "M1 strain, primary v3"
MdomM2pV4_LAB = "M2 strain, primary v4"
MdomM2pV5_LAB = "M2 strain, primary v5b"
MdomM3p_LAB = "M3 strain, primary v3"
MdomM5p_LAB = "M5 strain, primary v2"
MdomMYp_LAB = "MY (aabys) strain, primary v8"


# breakspoints for inversion, edge of "forward orientation" blocks pulled from df.coords
breaks.M2_Chr2_h1 = c(139250797, 174968441)
breaks.M2_Chr2_h2 = c(109898361, 145711716)
breaks.M2_Chr5_h1 = c(10246246, 48619955)
breaks.M2_Chr5_h2 = c(10357632, 48282746)
breaks.MY_Chr5_h1 = c(10801295, 48740278)
breaks.MY_Chr5_h2 = c(11793087, 50283522)

# coordinates for Mdmd in M2 Chr1 hap1
Mdmd_loc.M2_Chr2_h1 = c(111313000, 111317000)

# breakpoint coords for primaries: Chr1 M5 vs others
# inversion A, small low divergence inversion surrounded by deletion, coordinates the edges of inverted block
breaks.M5_Chr1A = c(59657191, 62680988) # inclusive of alignments to all strains
breaks.M1_Chr1A = c(64545922, 67644658)
breaks.M2_Chr1A = c(57990272, 60936889)
breaks.M3_Chr1A = c(63869591, 66984816)
breaks.MY_Chr1A = c(63699187, 66831144)
# inversion B, large high divergence inversion, coordinates edge of flanking syntenic sequence
breaks.M5_Chr1B = c(71347106, 84608382) # inclusive of all alignments (narrowest region)
breaks.M1_Chr1B = c(89924657, 105233525)
breaks.M2_Chr1B = c(81490361, 95098335)
breaks.M3_Chr1B = c(92338906, 104955418)
breaks.MY_Chr1B = c(90313137, 102865104)

# breakpoint coords for M2 haps Chr1B inversion 
breaks.M2h1v3_Chr1B = c(71492521,84936273)
breaks.M2h2v4_Chr1B = c(54092957,66352199)

# main plot objects via plot.df.coords()
# old haps with unresolved inversion misphased (both in h1)
p.M2.h1v2_h2v3.Chr1 = plot.df.coords(l.c.M2.h1v2_h2v3, REFERENCE = "Chromosome1", QUERY = "Chromosome1",
                                     SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                                     ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                                     LAB_REF = MdomM2h1v2_LAB, LAB_QRY = MdomM2h2v3_LAB)
# new haps with one copy of inversion moved to h2
p.M2.h1v3_h2v4.Chr1 = plot.df.coords(l.c.M2.h1v3_h2v4, REFERENCE = "Chromosome1", QUERY = "Chromosome1",
                                     SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                                     ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                                     LAB_REF = MdomM2h1v3_LAB, LAB_QRY = MdomM2h2v4_LAB)
p.M2.h1v3_h2v4.Chr2 = plot.df.coords(l.c.M2.h1v3_h2v4, REFERENCE = "Chromosome2", QUERY = "Chromosome2", 
                                     SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                                     ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                                     LAB_REF = MdomM2h1v3_LAB, LAB_QRY = MdomM2h2v4_LAB)

p.M2.h1v3_h2v4.Chr5 = plot.df.coords(l.c.M2.h1v3_h2v4, REFERENCE = "Chromosome5", QUERY = "Chromosome5", 
                        ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                        SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                        LAB_REF = MdomM2h1v3_LAB, LAB_QRY = MdomM2h2v4_LAB)

p.MYChr5 = plot.df.coords(l.c.MYh1_MYh2, REFERENCE = "scaffold_5", QUERY = "scaffold_5",
                          ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                          SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                          LAB_REF = MdomMYh1_LAB, LAB_QRY = MdomMYh2_LAB)

p.M1M5Chr1 = plot.df.coords(l.c.M1_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                            ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                            SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                            LAB_REF = MdomM1p_LAB, LAB_QRY = MdomM5p_LAB)

p.M2v4M5Chr1 = plot.df.coords(l.c.M2v4_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                            ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                            SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                            LAB_REF = MdomM2pV4_LAB, LAB_QRY = MdomM5p_LAB)
p.M2v5bM5Chr1 = plot.df.coords(l.c.M2v5b_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                              ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                              SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                              LAB_REF = MdomM2pV5_LAB, LAB_QRY = MdomM5p_LAB)
p.M2h1v3M5Chr1 = plot.df.coords(l.c.M2h1v3_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                                ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                                SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                                LAB_REF = MdomM2h1v3_LAB, LAB_QRY = MdomM5p_LAB)
p.M2h2v4M5Chr1 = plot.df.coords(l.c.M2h2v4_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                                ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                                SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                                LAB_REF = MdomM2h2v4_LAB, LAB_QRY = MdomM5p_LAB)

p.M3M5Chr1 = plot.df.coords(l.c.M3_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
                            ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                            SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                            LAB_REF = MdomM3p_LAB, LAB_QRY = MdomM5p_LAB)

p.MYM5Chr1 = plot.df.coords(l.c.MY_M5, REFERENCE = "scaffold_1", QUERY = "scaffold_1",
                            ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
                            SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
                            LAB_REF = MdomMYp_LAB, LAB_QRY = MdomM5p_LAB)

# zooms for Jae, plotting his guess at breakpoints based on microsynteny in annotation
PlotDFCoords(l.c.M2.h1v3_h2v4, REFERENCE = "Chromosome1", QUERY = "Chromosome1",
             SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
             ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
             LAB_REF = MdomM2h1v3_LAB, LAB_QRY = MdomM2h2v4_LAB, 
             TIC_LABELS = T, TIC_COUNT = 25, COORD_OFFSET = .05, 
             QRY_LIM = c(53000000 , 69000000), REF_LIM = c(69000000, 87000000)) + 
  geom_hline(yintercept = c(54.9e6, 65.4e6), linetype = 2) + 
  geom_vline(xintercept = c(72.3e6, 84.6e6), linetype = 2)
PlotDFCoords(l.c.M2h2v4_M5, REFERENCE = "Chromosome1", QUERY = "scaffold_1",
             ALPHA = 1, LINEWIDTH = 1, POINTSIZE = 0.5,
             SEQLABANGLE = 0, MIN_MATCH_LENGTH = MIN_MATCH,
             LAB_REF = MdomM2h2v4_LAB, LAB_QRY = MdomM5p_LAB, 
             TIC_LABELS = T, TIC_COUNT = 25, COORD_OFFSET = 0.05, 
             REF_LIM = c(53000000 , 69000000), QRY_LIM = c(70000000, 87000000)) + 
  geom_vline(xintercept = c(54.9e6, 65.4e6), linetype = 2) + 
  geom_hline(yintercept = c(72.2e6, 84.51e6), linetype = 2)


pdf(file = "Fig-M2_inversions-NoOrientation.pdf", width = 5.5, height = 4.5)

p.M2Chr2 + 
  geom_vline(xintercept = breaks.M2_Chr2_h1, linetype = "longdash") + 
  geom_text(data = data.frame(
      lab = sapply(breaks.M2_Chr2_h1, comma), 
      pos = breaks.M2_Chr2_h1
    ), aes(label = lab, x = pos), 
    y = 0, angle = 90,
    size = 3, hjust = -0.1, vjust = -0.5) + 
  geom_hline(yintercept = breaks.M2_Chr2_h2, linetype = "longdash") + 
  geom_text(data = data.frame(
      lab = sapply(breaks.M2_Chr2_h2, comma),
      pos = breaks.M2_Chr2_h2
    ), aes(label = lab, y = pos),
    x = 0, angle = 0,
    size = 3, hjust = -0.1, vjust = -0.5)

p.M2Chr5 +
  geom_vline(xintercept = breaks.M2_Chr5_h1, linetype = "longdash") +
  geom_text(data = data.frame(
      lab = sapply(breaks.M2_Chr5_h1, comma),
      pos = breaks.M2_Chr5_h1
    ), aes(label = lab, x = pos),
    y = l.c.M2h1_M2h2$breaksQry["Chromosome5", 2], angle = 90,
    size = 3, hjust = 1.1, vjust = -0.5) +
  geom_hline(yintercept = breaks.M2_Chr5_h2, linetype = "longdash") +
  geom_text(data = data.frame(
      lab = sapply(breaks.M2_Chr5_h2, comma),
      pos = breaks.M2_Chr5_h2
    ), aes(label = lab, y = pos),
    x = l.c.M2h1_M2h2$breaksRef["Chromosome5", 2], angle = 0,
    size = 3, hjust = 1.1, vjust = -0.5)

dev.off()


pdf(file = "Fig-MY_Chr5_inversion.pdf", width = 5.5, height = 4.5)

p.MYChr5 +
  geom_vline(xintercept = breaks.MY_Chr5_h1, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.MY_Chr5_h1, comma),
    pos = breaks.MY_Chr5_h1
  ), aes(label = lab, x = pos),
  y = l.c.MYh1_MYh2$breaksQry["scaffold_5", 2], angle = 90,
  size = 3, hjust = 1.1, vjust = -0.5) +
  geom_hline(yintercept = breaks.MY_Chr5_h2, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.MY_Chr5_h2, comma),
    pos = breaks.MY_Chr5_h2
  ), aes(label = lab, y = pos),
  x = l.c.MYh1_MYh2$breaksRef["scaffold_5", 2], angle = 0,
  size = 3, hjust = 1.1, vjust = -0.5)

dev.off()



pdf(file = "Fig-M5_Chr1_inversion-v2.pdf", width = 5, height = 4.5)

p.M1M5Chr1 + 
  geom_vline(xintercept = breaks.M1_Chr1A, linetype = "dotted") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M1_Chr1A, comma), 
    pos = breaks.M1_Chr1A
  ), aes(label = lab, x = pos), 
  y = 0, angle = 90, size = 2, 
  hjust = -0.1, vjust = c(-0.5, 1.3)) +
  geom_hline(yintercept = breaks.M5_Chr1A, linetype = "dotted") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1A, comma),
    pos = breaks.M5_Chr1A
  ), aes(label = lab, y = pos),
  x = l.c.M1_M5$breaksRef["Chromosome1",2], angle = 0,
  size = 2, hjust = 1.1, vjust = c(1.3, -0.5)) +
  geom_vline(xintercept = breaks.M1_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M1_Chr1B, comma), 
    pos = breaks.M1_Chr1B
    ), aes(label = lab, x = pos), 
    y = l.c.M1_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
    hjust = 1.1, vjust = c(-0.5, -0.5)) +
  geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.M2v5bM5Chr1 + 
  geom_vline(xintercept = breaks.M2_Chr1A, linetype = "dotted") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M2_Chr1A, comma), 
    pos = breaks.M2_Chr1A
  ), aes(label = lab, x = pos), 
  y = 0, angle = 90, size = 2, 
  hjust = -0.1, vjust = c(-0.5, 1.3)) +
  geom_hline(yintercept = breaks.M5_Chr1A, linetype = "dotted") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1A, comma),
    pos = breaks.M5_Chr1A
  ), aes(label = lab, y = pos),
  x = l.c.M2v5b_M5$breaksRef["Chromosome1",2], angle = 0,
  size = 2, hjust = 1.1, vjust = c(1.3, -0.5)) +
  geom_vline(xintercept = breaks.M2_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M2_Chr1B, comma), 
    pos = breaks.M2_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.M2v5b_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.M3M5Chr1 + 
  geom_vline(xintercept = breaks.M3_Chr1A, linetype = "dotted") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M3_Chr1A, comma), 
    pos = breaks.M3_Chr1A
  ), aes(label = lab, x = pos), 
  y = 0, angle = 90, size = 2, 
  hjust = -0.1, vjust = c(-0.5, 1.3)) +
  geom_hline(yintercept = breaks.M5_Chr1A, linetype = "dotted") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1A, comma),
    pos = breaks.M5_Chr1A
  ), aes(label = lab, y = pos),
  x = l.c.M3_M5$breaksRef["Chromosome1",2], angle = 0,
  size = 2, hjust = 1.1, vjust = c(1.3, -0.5)) +
  geom_vline(xintercept = breaks.M3_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M3_Chr1B, comma), 
    pos = breaks.M3_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.M3_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.MYM5Chr1 + 
  geom_vline(xintercept = breaks.MY_Chr1A, linetype = "dotted") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.MY_Chr1A, comma), 
    pos = breaks.MY_Chr1A
  ), aes(label = lab, x = pos), 
  y = 0, angle = 90, size = 2, 
  hjust = -0.1, vjust = c(-0.5, 1.3)) +
  geom_hline(yintercept = breaks.M5_Chr1A, linetype = "dotted") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1A, comma),
    pos = breaks.M5_Chr1A
  ), aes(label = lab, y = pos),
  x = l.c.MY_M5$breaksRef["scaffold_1",2], angle = 0,
  size = 2, hjust = 1.1, vjust = c(1.3, -0.5)) +
  geom_vline(xintercept = breaks.MY_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.MY_Chr1B, comma), 
    pos = breaks.MY_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.MY_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
  geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.M2h1v3M5Chr1 + 
  geom_vline(xintercept = breaks.M2h1v3_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M2h1v3_Chr1B, comma), 
    pos = breaks.M2h1v3_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.M2v5_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
  geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.M2h2v4M5Chr1 + 
  geom_vline(xintercept = breaks.M2h2v4_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M2h2v4_Chr1B, comma), 
    pos = breaks.M2h2v4_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.M2v5_M5$breaksQry["scaffold_1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
  geom_hline(yintercept = breaks.M5_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M5_Chr1B, comma),
    pos = breaks.M5_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(-0.5, -0.5))

p.M2.h1v3_h2v4.Chr1 + 
  geom_vline(xintercept = breaks.M2h1v3_Chr1B, linetype = "longdash") + 
  geom_text(data = data.frame(
    lab = sapply(breaks.M2h1v3_Chr1B, comma), 
    pos = breaks.M2h1v3_Chr1B
  ), aes(label = lab, x = pos), 
  y = l.c.M2.h1v3_h2v4$breaksQry["Chromosome1",2], angle = 90, size = 2, 
  hjust = 1.1, vjust = c(-0.5, -0.5)) +
  geom_hline(yintercept = breaks.M2h2v4_Chr1B, linetype = "longdash") +
  geom_text(data = data.frame(
    lab = sapply(breaks.M2h2v4_Chr1B, comma),
    pos = breaks.M2h2v4_Chr1B
  ), aes(label = lab, y = pos),
  x = 0, angle = 0,
  size = 2, hjust = -0.1, vjust = c(1.3, 1.3))

dev.off()

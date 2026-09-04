library(HiCExperiment)
library(HiContacts)

### 5M
hic_M5@interactions@regions@seqnames
hic_M5_add <- import("5M_50kb_norm_corrected.cool", format='cool')
plotMatrix(hic_M5_add[c("scaffold1","scaffold2","scaffold3","scaffold4","scaffold5","scaffold6","scaffold7","scaffold8")], caption=F)

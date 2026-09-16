# alignment_dotplot-FUNCTIONS.R 
# takes nucmer alignment coordinates file (from mummer show-coords),
#  and ref/qry fasta indexes (from samtoolds faidx)
# produces dotplot of alignment regions 
# FUNCTIONS file to source by other scripts

library(tidyverse)

# FUNCTIONS

# new version using mummer's coord output format, along with fai index files for seq lengths
#   FAI inputs need sequence ID and length columns, 
#   can use/reuse any version so long as all aligned sequences represented
# function to list dfs of match coordinates and sequence ref/qry lengths ("breaks")
# copies format of older version l.coord header names etc
make.l.coords <- function(COORDSFILE, REF_FAI, QRY_FAI){
  # read coords format, skip first 5 header lines
  df.coords <- read.table(COORDSFILE, stringsAsFactors = F, skip = 5)
  # drop formatting separating columns
  df.coords <- df.coords[,-c(3,6,9,11,14)]
  # rename to work with downstream subsetting/plot functions
  colnames(df.coords) <- c("R1", "R2", "Q1", "Q2", 
                           "ref_length", "qry_length", 
                           "pctID",
                           "ref_SeqName", "qry_SeqName")
  # assign orientation, in coords format reverse matches shown in query start/stop reversal
  df.coords$orientation <- apply(df.coords, 1, FUN = function(x){
    if( x[3] <= x[4] ){return("forward")}
    if( x[3] >= x[4] ){return("reverse")}
  })
  
  # "breaks" dfs for tracking sequence lengths
  df.breaksRef <- read.table(REF_FAI)
  df.breaksRef <- df.breaksRef[,c(1,2)] 
  colnames(df.breaksRef) <- c("SeqName", "SeqLength")
  rownames(df.breaksRef) <- df.breaksRef$SeqName
  
  df.breaksQry <- read.table(QRY_FAI)
  df.breaksQry <- df.breaksQry[,c(1,2)]
  colnames(df.breaksQry) <- c("SeqName", "SeqLength")
  rownames(df.breaksQry) <- df.breaksQry$SeqName
  
  return( list(coords = df.coords, breaksRef = df.breaksRef, breaksQry = df.breaksQry, filename = COORDSFILE) )
}


# Function to build a df with subset of coordinates from given sequence names in order given, by default returns df.coords
make.scaf.df.coords <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, 
                                DROP_EMPTY_SCAFFOLDS=T,
                                MIN_MATCH_LENGTH=0){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }
  
  # get matches for qry and reference
  DF.COORDS_SUBSET <- COORDS[ (COORDS$ref_SeqName %in% REFERENCE) & (COORDS$qry_SeqName %in% QUERY) ,]
  
  # filter based on minimum match length
  DF.COORDS_SUBSET <- DF.COORDS_SUBSET[ (abs(DF.COORDS_SUBSET$ref_length) > MIN_MATCH_LENGTH) &
                                          (abs(DF.COORDS_SUBSET$qry_length) > MIN_MATCH_LENGTH),]
  
  if( DROP_EMPTY_SCAFFOLDS ){
    REFERENCE <- REFERENCE[ REFERENCE %in% DF.COORDS_SUBSET$ref_SeqName ]
    QUERY <- QUERY[ QUERY %in% DF.COORDS_SUBSET$qry_SeqName ]
  }
  
  # lengths of remaining qry/scaffs
  LENGTHS_REF <- sapply(REFERENCE, function(x) { as.numeric(BREAKS_REF[x,]$SeqLength) }) 
  LENGTHS_QRY <- sapply(QUERY, function(x) { as.numeric(BREAKS_QRY[x,]$SeqLength) })
  
  STARTPOS_REF <- c(1, 1+cumsum(LENGTHS_REF)[-length(LENGTHS_REF)] )
  names(STARTPOS_REF) <- REFERENCE
  STARTPOS_QRY <- c(1, 1+cumsum(LENGTHS_QRY)[-length(LENGTHS_QRY)] )
  names(STARTPOS_QRY) <- QUERY
  
  DF.COORDS_SUBSET$ref_StartPos <- STARTPOS_REF[DF.COORDS_SUBSET$ref_SeqName]
  DF.COORDS_SUBSET$qry_StartPos <- STARTPOS_QRY[DF.COORDS_SUBSET$qry_SeqName]
  
  DF.COORDS_SUBSET$ref_SeqLength <- LENGTHS_REF[DF.COORDS_SUBSET$ref_SeqName]
  DF.COORDS_SUBSET$qry_SeqLength <- LENGTHS_QRY[DF.COORDS_SUBSET$qry_SeqName]
  
  DF.COORDS_SUBSET[,c('R1','R2')] <- DF.COORDS_SUBSET[,c('R1','R2')] + #- 
    DF.COORDS_SUBSET$ref_StartPos 
  
  DF.COORDS_SUBSET[,c('Q1','Q2')] <- DF.COORDS_SUBSET[,c('Q1','Q2')] + #- 
    DF.COORDS_SUBSET$qry_StartPos
  
  return(DF.COORDS_SUBSET)
  
}


plot.df.coords <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL,
                           MIN_MATCH_LENGTH = 0, 
                           REF_LIM = NULL, QRY_LIM = NULL, # akin to xlim, ylim - takes c(START, STOP) coordinates
                           DROP_EMPTY_SCAFFOLDS = T, FLIP_AXES = F, COORD_OFFSET = 0.02,
                           LINEWIDTH = 1, ALPHA = 0.25, POINTSIZE = 0.8,
                           SEQLABANGLE = 45,
                           MAR_T = 25, MAR_R = 5, MAR_B = 5, MAR_L = 10,
                           LAB_REF = "", LAB_QRY = ""){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }
  
  REFERENCE = as.character(REFERENCE)
  QUERY = as.character(QUERY)
  COORDS <- COORDS[((COORDS$ref_SeqName %in% REFERENCE) & (COORDS$qry_SeqName %in% QUERY)),]
  
  DF_COORDS <- make.scaf.df.coords(L.COORDS = L.COORDS, REFERENCE = REFERENCE, QUERY = QUERY, 
                                   DROP_EMPTY_SCAFFOLDS = DROP_EMPTY_SCAFFOLDS,
                                   MIN_MATCH_LENGTH = MIN_MATCH_LENGTH)
  
  if( nrow(DF_COORDS) > 0 ){  
    
    # By default only plot scaffolds with matches
    if( DROP_EMPTY_SCAFFOLDS ){  
      
      REF_STARTPOS <- unique(DF_COORDS$ref_StartPos)
      REF_SEQNAMES <- unique(DF_COORDS$ref_SeqName)
      QRY_STARTPOS <- unique(DF_COORDS$qry_StartPos)
      QRY_SEQNAMES <- unique(DF_COORDS$qry_SeqName)
      REF_MAX = max(DF_COORDS$ref_StartPos + DF_COORDS$ref_SeqLength)
      QRY_MAX = max(DF_COORDS$qry_StartPos + DF_COORDS$qry_SeqLength)
      
    } else {
      # otherwise need to replicate some of GetScaffoldCoordsDF() but without filtering for matches
      DF_BREAKSREF <- L.COORDS$breaksRef[REFERENCE,]
      DF_BREAKSQRY <- L.COORDS$breaksQry[QUERY,]
      
      REF_STARTPOS <- c(1, 1+cumsum(DF_BREAKSREF$SeqLength)[-length(REFERENCE)])
      REF_SEQNAMES <- REFERENCE
      QRY_STARTPOS <- c(1, 1+cumsum(DF_BREAKSQRY$SeqLength)[-length(QUERY)])
      QRY_SEQNAMES <- QUERY
      REF_MAX <- REF_STARTPOS[length(REF_STARTPOS)] + DF_BREAKSREF$SeqLength[nrow(DF_BREAKSREF)]
      QRY_MAX <- QRY_STARTPOS[length(QRY_STARTPOS)] + DF_BREAKSQRY$SeqLength[nrow(DF_BREAKSQRY)]
    }
    
    # adjust boundaries for REF_LIM and QRY_LIM
    if( !is.null(REF_LIM) ){
      REF_SEQNAMES = REF_SEQNAMES[intersect(which(c(REF_STARTPOS,REF_MAX) > REF_LIM[1]) - 1, 
                                            which(REF_STARTPOS < REF_LIM[2]))]
      REF_STARTPOS = c(REF_LIM[1], 
                       REF_STARTPOS[(REF_STARTPOS > REF_LIM[1]) & (REF_STARTPOS < REF_LIM[2])])
      REF_MIN = REF_LIM[1]
      REF_MAX = REF_LIM[2]
    } else {
      REF_MIN = 0
    }
    if( !is.null(QRY_LIM) ){
      QRY_SEQNAMES = QRY_SEQNAMES[intersect(which(c(QRY_STARTPOS,QRY_MAX) > QRY_LIM[1]) - 1, 
                                            which(QRY_STARTPOS < QRY_LIM[2]))]
      QRY_STARTPOS = c(QRY_LIM[1], 
                       QRY_STARTPOS[(QRY_STARTPOS > QRY_LIM[1]) & (QRY_STARTPOS < QRY_LIM[2])])
      QRY_MIN = QRY_LIM[1]
      QRY_MAX = QRY_LIM[2]
    } else {
      QRY_MIN = 0
    }
    
    if( FLIP_AXES ){
      XPOS = "top"
      YPOS = "left"
    } else {
      XPOS = "bottom"
      YPOS = "right"
    }
    
    p1 <- ggplot(DF_COORDS, aes(x=R1, y=Q1)) +
      geom_point(alpha=ALPHA, size=POINTSIZE, color = "grey30") +
      geom_point(aes(x=R2, y=Q2), color = "grey30", alpha=ALPHA, size=POINTSIZE) +
      geom_segment(aes(xend=R2, yend=Q2), color = "grey30", linewidth = LINEWIDTH) +
      labs(x=LAB_REF, y=LAB_QRY) +
      geom_segment(data = data.frame(x = c(REF_STARTPOS, REF_MAX), xend=c(REF_STARTPOS, REF_MAX), y=QRY_MIN, yend=QRY_MAX), 
                   aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25, inherit.aes = F) + 
      geom_segment(data = data.frame(x=REF_MIN, xend=REF_MAX, y = c(QRY_STARTPOS, QRY_MAX), yend = c(QRY_STARTPOS, QRY_MAX)), 
                   aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25, inherit.aes = F) +
      scale_x_continuous(breaks = c(REF_STARTPOS), 
                         labels = c(REF_SEQNAMES), position = XPOS,
                         minor_breaks = NULL) +
      scale_y_continuous(breaks = c(QRY_STARTPOS),
                         labels = c(QRY_SEQNAMES), position = YPOS,
                         minor_breaks = NULL ) + 
      guides(color = guide_legend(override.aes = list(linewidth = LINEWIDTH))) + 
      theme_minimal()
    if( FLIP_AXES ){
      p2 <- p1 +
        coord_flip(xlim = c( (REF_MIN + COORD_OFFSET*(REF_MAX-REF_MIN)), 
                             (REF_MAX - COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                   ylim = c( (QRY_MIN + COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                             (QRY_MAX - COORD_OFFSET*(QRY_MAX-QRY_MIN)) )) 
    } else {
      p2 <- p1 +
        coord_cartesian(xlim = c( (REF_MIN + COORD_OFFSET*(REF_MAX-REF_MIN)), 
                                  (REF_MAX - COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                        ylim = c( (QRY_MIN + COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                                  (QRY_MAX - COORD_OFFSET*(QRY_MAX-QRY_MIN)) )) 
      
    }
    return(p2 +
             theme(plot.margin = margin(MAR_T, MAR_R, MAR_B, MAR_L), 
                   axis.title.y.right = element_text(angle = 90),
                   axis.text.x = element_text(angle = -SEQLABANGLE, 
                                              hjust = 0, vjust = 0, 
                                              color = "slateblue", 
                                              face = "bold"), 
                   axis.text.y = element_text(angle = 90-SEQLABANGLE, 
                                              hjust = 0, vjust = 0, 
                                              color = "slateblue", 
                                              face = "bold") )
           
    )
  } else {
    # no coordinates returned
    return(NULL)
  }
}

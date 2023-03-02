dgrp2_892_indel_2L <- read.table(file="dgrp2_892_indel_2L.txt")
dgrp2_437_indel_2L <- read.table(file="dgrp2_437_indel_2L.txt")
dgrp2_437_indel_2L$V3 <- gsub("2L_","", dgrp2_437_indel_2L$V3)
dgrp2_892_indel_2L$V3 <- gsub("2L_","", dgrp2_892_indel_2L$V3)
dgrp2_437_indel_2L$V3 <- gsub(".*_","", dgrp2_437_indel_2L$V3)
dgrp2_892_indel_2L$V3 <- gsub(".*_","", dgrp2_892_indel_2L$V3)
dgrp2_437_indel_2L <- dgrp2_437_indel_2L[dgrp2_437_indel_2L$V10 %in% "1/1",]
dgrp2_892_indel_2L <- dgrp2_892_indel_2L[dgrp2_892_indel_2L$V10 %in% "1/1",]
dgrp2_437_indel_2L$size <- abs(nchar(dgrp2_437_indel_2L$V4) - nchar(dgrp2_437_indel_2L$V5))
dgrp2_892_indel_2L$size <- abs(nchar(dgrp2_892_indel_2L$V4) - nchar(dgrp2_892_indel_2L$V5))                    
                                
dgrp2_437_indel_2La <- dgrp2_437_indel_2L[dgrp2_437_indel_2L$V2 < 3717622,]
dgrp2_892_indel_2La <- dgrp2_892_indel_2L[dgrp2_892_indel_2L$V2 < 3717622,]

3715788 - sum(dgrp2_437_indel_2La[dgrp2_437_indel_2La$V3 %in% "INS",]$size ) + sum(dgrp2_437_indel_2La[dgrp2_437_indel_2La$V3 %in% "DEL",]$size )
3712824 - sum(dgrp2_892_indel_2La[dgrp2_892_indel_2La$V3 %in% "INS",]$size ) + sum(dgrp2_892_indel_2La[dgrp2_892_indel_2La$V3 %in% "DEL",]$size )

dgrp2_437_indel_2Lb <- dgrp2_437_indel_2L[dgrp2_437_indel_2L$V2 < 3717533,]
dgrp2_892_indel_2Lb <- dgrp2_892_indel_2L[dgrp2_892_indel_2L$V2 < 3717533,]

3715718 - sum(dgrp2_437_indel_2Lb[dgrp2_437_indel_2Lb$V3 %in% "INS",]$size ) + sum(dgrp2_437_indel_2Lb[dgrp2_437_indel_2Lb$V3 %in% "DEL",]$size )
3712735 - sum(dgrp2_892_indel_2Lb[dgrp2_892_indel_2Lb$V3 %in% "INS",]$size ) + sum(dgrp2_892_indel_2Lb[dgrp2_892_indel_2Lb$V3 %in% "DEL",]$size )


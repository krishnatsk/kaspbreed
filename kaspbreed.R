

#setwd("/home/irri-hub/Downloads/rp")
#library(seqinr)
#library(stringr)
#library(dplyr)
#all <- read.fasta("all_pr.fas")
#final_primer_filter2 <- read.csv("preeti_snpeff.csv",header = T)
#gene <- final_primer_filter2$gene_id
#gene <- unique(gene)
#final_prim<- list()
#coln<- colnames(final_primer_filter2)

kaspprimer <- function(ref.fasta, snp.table, stream){
  all <- read.fasta(ref.fasta)
  final_primer_filter2 <- read.csv(snp.table,header = T)
  gene <- snp.table$gene_id
  gene <- unique(gene)

for (seq in gene){
  fasta <- all[seq]
  fasta <- unlist(fasta)
  fasta <- paste0(fasta,sep="",collapse = "")
  asd<- subset(final_primer_filter2,final_primer_filter2$gene_id %in% seq)

id<- asd$snp_pos
prime <- matrix(nrow = 0,ncol = 23)
prime <- data.frame(prime)
colnames(prime) <- coln
updown <- stream

   for (pos in id) {
      dsa<- subset(asd,asd$snp_pos %in% pos)
      dsa$strand <- as.character(dsa$strand)
      if (dsa$strand == "+"){

      if (str_count(dsa$REF) == 1 &  str_count(dsa$ALT) == 1){
      snp_pos <- pos
      upstream <- substr(fasta, snp_pos - updown, snp_pos - 1)
        downstream <- substr(fasta, snp_pos+1, snp_pos + updown)
         names(upstream) <- "upstream"
         names(downstream) <- "downstream"
       primer <- c(upstream,downstream)
      primer <- t(data.frame(primer))
     #prim<- rbind(prim,primer)
     #prim[,"gene_id"] <- seq
     pri <- cbind(dsa,primer)
     prime <- rbind(prime,pri)
     }
     else {
       ref <- str_count(dsa$REF)
       alt <- str_count(dsa$ALT)
       snp_pos <- pos
       if (ref > 1){
       upstream <- substr(fasta, snp_pos - updown, snp_pos - 1)
       downstream <- substr(fasta, snp_pos+ref-1, snp_pos + updown)
       names(upstream) <- "upstream"
       names(downstream) <- "downstream"
       primer <- c(upstream,downstream)
       primer <- t(data.frame(primer))
       #prim<- rbind(prim,primer)
       #prim[,"gene_id"] <- seq
       pri <- cbind(dsa,primer)
       prime <- rbind(prime,pri)

       }

       else if (alt >1){
         upstream <- substr(fasta, snp_pos - updown, snp_pos - 1)
         downstream <- substr(fasta, snp_pos+alt-1, snp_pos + updown)
         names(upstream) <- "upstream"
         names(downstream) <- "downstream"
         primer <- c(upstream,downstream)
         primer <- t(data.frame(primer))
         #prim<- rbind(prim,primer)
         #prim[,"gene_id"] <- seq
         pri <- cbind(dsa,primer)
         prime <- rbind(prime,pri)

         }
      else {
        print("good")
      }

        }

  }
  else {
    #revfas<-c2s(rev(comp(s2c(fasta))))
    revfas <-fasta
    if (str_count(dsa$REF) == 1 &  str_count(dsa$ALT) == 1){
      snp_pos <- pos
      upstream <- substr(revfas, snp_pos - updown, snp_pos - 1)
      downstream <- substr(revfas, snp_pos+1, snp_pos + updown)
      names(upstream) <- "upstream"
      names(downstream) <- "downstream"
      primer <- c(upstream,downstream)
      primer <- t(data.frame(primer))
      #prim<- rbind(prim,primer)
      #prim[,"gene_id"] <- seq
      pri <- cbind(dsa,primer)
      prime <- rbind(prime,pri)

    }
    else {
      ref <- str_count(dsa$REF)
      alt <- str_count(dsa$ALT)
      snp_pos <- pos
      if (ref > 1){
        upstream <- substr(revfas, snp_pos - updown, snp_pos - 1)
        downstream <- substr(revfas, snp_pos+ref-1, snp_pos + updown)
        names(upstream) <- "upstream"
        names(downstream) <- "downstream"
        primer <- c(upstream,downstream)
        primer <- t(data.frame(primer))
        #prim<- rbind(prim,primer)
        #prim[,"gene_id"] <- seq
        pri <- cbind(dsa,primer)
        prime <- rbind(prime,pri)

      }

      else if (alt >1){
        upstream <- substr(revfas, snp_pos - updown, snp_pos - 1)
        downstream <- substr(revfas, snp_pos+alt-1, snp_pos + updown)
        names(upstream) <- "upstream"
        names(downstream) <- "downstream"
        primer <- c(upstream,downstream)
        primer <- t(data.frame(primer))
        #prim<- rbind(prim,primer)
        #prim[,"gene_id"] <- seq
        pri <- cbind(dsa,primer)
        prime <- rbind(prime,pri)

      }
      else {
        print("good")
      }

       }
  }

     final_prim[[seq]] <- prime
   }
}

}
#filter_primer <- bind_rows(final_prim)
#write.csv(filter_primer,file = "preeti_kasp_res.csv",row.names = F)



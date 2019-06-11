
# ><>< ========================================================================= ><>< #
# ><><                         Sadiq, Hassan - SDQHAS001                         ><>< #
# ><><                         PhD Statistics - STA6001W                         ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><        ASSESSING THE ROBUSTNESS OF PHYLOGENETIC AND PHYLOGEOGRAPHIC       ><>< #
# ><><             MODELS TO CHANGES IN RESIDUE FITNESSES OVER TIME:             ><>< #
# ><><                            A Simulation Study.                            ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                                  R Code:                                  ><>< #
# ><><           Illustrative Simulation for Literature Review Chapter           ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
# ><><                              9th June, 2019.                              ><>< #
# ><>< ========================================================================= ><>< #

# ><>< ========================================================================= ><>< #
# ><><                           COMPILATION RESOURCES                           ><>< #
# ><><                           ~~~~~~~~~~~~~~~~~~~~~                           ><>< #
#      Computer System (MacBook) :-                                                   #
#        Intel Core i7, CPU 2.90GHz, 16GB RAM                                         #
#        MacBook Pro OS Mojave (Version 10.14.5)                                      #
#      R :-                                                                           #
#        Version 3.4.2 -- "Short Summer"                                              #
# ><>< ========================================================================= ><>< #

# ><><                                 Preamble.                                 ><>< #
rm(list=ls())
library(Matrix)

# ><><             Goldman-Yang Codon Substitution Matrix Definition             ><>< #
t <- 1/61;   k <- 1e-06;   w <- 1e-02
GY94 <- matrix(0, 61, 61); GY94[1,2] = t*w;
      GY94[1,3] = k*t;     GY94[1,4] = t*w;     GY94[1,5] = t*w;   GY94[1,9] = k*t*w
     GY94[1,13] = t*w;    GY94[1,17] = t*w;  GY94[1,33] = k*t*w;     GY94[2,1] = t*w
      GY94[2,3] = t*w;     GY94[2,4] = k*t;     GY94[2,6] = t*w;  GY94[2,10] = k*t*w
     GY94[2,14] = t*w;    GY94[2,18] = t*w;  GY94[2,34] = k*t*w;    GY94[2,49] = t*w
      GY94[3,1] = k*t;     GY94[3,2] = t*w;     GY94[3,4] = t*w;     GY94[3,7] = t*w
   GY94[3,11] = k*t*w;    GY94[3,15] = t*w;    GY94[3,19] = t*w;  GY94[3,35] = k*t*w
      GY94[4,1] = t*w;     GY94[4,2] = k*t;     GY94[4,3] = t*w;     GY94[4,8] = t*w
   GY94[4,12] = k*t*w;    GY94[4,16] = t*w;    GY94[4,20] = t*w;  GY94[4,36] = k*t*w
     GY94[4,50] = t*w;     GY94[5,1] = t*w;       GY94[5,6] = t;     GY94[5,7] = k*t
        GY94[5,8] = t;     GY94[5,9] = t*w;  GY94[5,13] = k*t*w;    GY94[5,21] = t*w
   GY94[5,37] = k*t*w;    GY94[5,51] = t*w;     GY94[6,2] = t*w;       GY94[6,5] = t
        GY94[6,7] = t;     GY94[6,8] = k*t;    GY94[6,10] = t*w;  GY94[6,14] = k*t*w
     GY94[6,22] = t*w;  GY94[6,38] = k*t*w;    GY94[6,52] = t*w;     GY94[7,3] = t*w
      GY94[7,5] = k*t;       GY94[7,6] = t;       GY94[7,8] = t;    GY94[7,11] = t*w
   GY94[7,15] = k*t*w;    GY94[7,23] = t*w;  GY94[7,39] = k*t*w;    GY94[7,53] = t*w
      GY94[8,4] = t*w;       GY94[8,5] = t;     GY94[8,6] = k*t;       GY94[8,7] = t
     GY94[8,12] = t*w;  GY94[8,16] = k*t*w;    GY94[8,24] = t*w;  GY94[8,40] = k*t*w
     GY94[8,54] = t*w;   GY94[9,1] = k*t*w;     GY94[9,5] = t*w;    GY94[9,10] = t*w
     GY94[9,11] = k*t;    GY94[9,12] = t*w;    GY94[9,13] = t*w;      GY94[9,25] = t
   GY94[9,41] = k*t*w;  GY94[10,2] = k*t*w;    GY94[10,6] = t*w;    GY94[10,9] = t*w
    GY94[10,11] = t*w;   GY94[10,12] = k*t;   GY94[10,14] = t*w;   GY94[10,26] = t*w
  GY94[10,42] = k*t*w;   GY94[10,55] = t*w;  GY94[11,3] = k*t*w;    GY94[11,7] = t*w
     GY94[11,9] = k*t;   GY94[11,10] = t*w;   GY94[11,12] = t*w;   GY94[11,15] = t*w
      GY94[11,27] = t; GY94[11,43] = k*t*w;   GY94[11,56] = t*w;  GY94[12,4] = k*t*w
     GY94[12,8] = t*w;    GY94[12,9] = t*w;   GY94[12,10] = k*t;   GY94[12,11] = t*w
    GY94[12,16] = t*w;   GY94[12,28] = t*w; GY94[12,44] = k*t*w;   GY94[12,57] = t*w
     GY94[13,1] = t*w;  GY94[13,5] = k*t*w;    GY94[13,9] = t*w;     GY94[13,14] = t
  GY94[13,15] = k*t*w;     GY94[13,16] = t;   GY94[13,29] = t*w; GY94[13,45] = k*t*w
    GY94[13,58] = t*w;    GY94[14,2] = t*w;  GY94[14,6] = k*t*w;   GY94[14,10] = t*w
      GY94[14,13] = t;   GY94[14,15] = t*w;   GY94[14,16] = k*t;   GY94[14,30] = t*w
  GY94[14,46] = k*t*w;   GY94[14,59] = t*w;    GY94[15,3] = t*w;  GY94[15,7] = k*t*w
    GY94[15,11] = t*w; GY94[15,13] = k*t*w;   GY94[15,14] = t*w;   GY94[15,16] = t*w
    GY94[15,31] = t*w; GY94[15,47] = k*t*w;   GY94[15,60] = t*w;    GY94[16,4] = t*w
   GY94[16,8] = k*t*w;   GY94[16,12] = t*w;     GY94[16,13] = t;   GY94[16,14] = k*t
    GY94[16,15] = t*w;   GY94[16,32] = t*w; GY94[16,48] = k*t*w;   GY94[16,61] = t*w
     GY94[17,1] = t*w;   GY94[17,18] = t*w;   GY94[17,19] = k*t;   GY94[17,20] = t*w
    GY94[17,21] = t*w; GY94[17,25] = k*t*w;   GY94[17,29] = t*w;   GY94[17,33] = t*w
     GY94[18,2] = t*w;   GY94[18,17] = t*w;   GY94[18,19] = t*w;   GY94[18,20] = k*t
    GY94[18,22] = t*w; GY94[18,26] = k*t*w;   GY94[18,30] = t*w;   GY94[18,34] = t*w
  GY94[18,49] = k*t*w;    GY94[19,3] = t*w;   GY94[19,17] = k*t;   GY94[19,18] = t*w
    GY94[19,20] = t*w;   GY94[19,23] = t*w; GY94[19,27] = k*t*w;   GY94[19,31] = t*w
    GY94[19,35] = t*w;    GY94[20,4] = t*w;   GY94[20,17] = t*w;   GY94[20,18] = k*t
    GY94[20,19] = t*w;   GY94[20,24] = t*w; GY94[20,28] = k*t*w;   GY94[20,32] = t*w
    GY94[20,36] = t*w; GY94[20,50] = k*t*w;    GY94[21,5] = t*w;   GY94[21,17] = t*w
      GY94[21,22] = t;   GY94[21,23] = k*t;     GY94[21,24] = t;   GY94[21,25] = t*w
  GY94[21,29] = k*t*w;   GY94[21,37] = t*w; GY94[21,51] = k*t*w;    GY94[22,6] = t*w
    GY94[22,18] = t*w;     GY94[22,21] = t;     GY94[22,23] = t;   GY94[22,24] = k*t
    GY94[22,26] = t*w; GY94[22,30] = k*t*w;   GY94[22,38] = t*w; GY94[22,52] = k*t*w
     GY94[23,7] = t*w;   GY94[23,19] = t*w;   GY94[23,21] = k*t;     GY94[23,22] = t
      GY94[23,24] = t;   GY94[23,27] = t*w; GY94[23,31] = k*t*w;   GY94[23,39] = t*w
  GY94[23,53] = k*t*w;    GY94[24,8] = t*w;   GY94[24,20] = t*w;     GY94[24,21] = t
    GY94[24,22] = k*t;     GY94[24,23] = t;   GY94[24,28] = t*w; GY94[24,32] = k*t*w
    GY94[24,40] = t*w; GY94[24,54] = k*t*w;      GY94[25,9] = t; GY94[25,17] = k*t*w
    GY94[25,21] = t*w;     GY94[25,26] = t;   GY94[25,27] = k*t;     GY94[25,28] = t
    GY94[25,29] = t*w;   GY94[25,41] = t*w;   GY94[26,10] = t*w; GY94[26,18] = k*t*w
    GY94[26,22] = t*w;     GY94[26,25] = t;     GY94[26,27] = t;   GY94[26,28] = k*t
    GY94[26,30] = t*w;   GY94[26,42] = t*w; GY94[26,55] = k*t*w;     GY94[27,11] = t
  GY94[27,19] = k*t*w;   GY94[27,23] = t*w;   GY94[27,25] = k*t;     GY94[27,26] = t
      GY94[27,28] = t;   GY94[27,31] = t*w;   GY94[27,43] = t*w; GY94[27,56] = k*t*w
    GY94[28,12] = t*w; GY94[28,20] = k*t*w;   GY94[28,24] = t*w;     GY94[28,25] = t
    GY94[28,26] = k*t;     GY94[28,27] = t;   GY94[28,32] = t*w;   GY94[28,44] = t*w
  GY94[28,57] = k*t*w;   GY94[29,13] = t*w;   GY94[29,17] = t*w; GY94[29,21] = k*t*w
    GY94[29,25] = t*w;     GY94[29,30] = t;   GY94[29,31] = k*t;     GY94[29,32] = t
    GY94[29,45] = t*w;   GY94[29,58] = k*t;   GY94[30,14] = t*w;   GY94[30,18] = t*w
  GY94[30,22] = k*t*w;   GY94[30,26] = t*w;     GY94[30,29] = t;     GY94[30,31] = t
    GY94[30,32] = k*t;   GY94[30,46] = t*w; GY94[30,59] = k*t*w;   GY94[31,15] = t*w
    GY94[31,19] = t*w; GY94[31,23] = k*t*w;   GY94[31,27] = t*w;   GY94[31,29] = k*t
      GY94[31,30] = t;     GY94[31,32] = t;   GY94[31,47] = t*w;   GY94[31,60] = k*t
    GY94[32,16] = t*w;   GY94[32,20] = t*w; GY94[32,24] = k*t*w;   GY94[32,28] = t*w
      GY94[32,29] = t;   GY94[32,30] = k*t;     GY94[32,31] = t;   GY94[32,48] = t*w
  GY94[32,61] = k*t*w;  GY94[33,1] = k*t*w;   GY94[33,17] = t*w;   GY94[33,34] = t*w
    GY94[33,35] = k*t;   GY94[33,36] = t*w;   GY94[33,37] = t*w; GY94[33,41] = k*t*w
    GY94[33,45] = t*w;  GY94[34,2] = k*t*w;   GY94[34,18] = t*w;   GY94[34,33] = t*w
    GY94[34,35] = t*w;   GY94[34,36] = k*t;   GY94[34,38] = t*w; GY94[34,42] = k*t*w
    GY94[34,46] = t*w;   GY94[34,49] = t*w;  GY94[35,3] = k*t*w;   GY94[35,19] = t*w
    GY94[35,33] = k*t;   GY94[35,34] = t*w;   GY94[35,36] = t*w;   GY94[35,39] = t*w
  GY94[35,43] = k*t*w;   GY94[35,47] = t*w;  GY94[36,4] = k*t*w;   GY94[36,20] = t*w
    GY94[36,33] = t*w;   GY94[36,34] = k*t;   GY94[36,35] = t*w;   GY94[36,40] = t*w
  GY94[36,44] = k*t*w;   GY94[36,48] = t*w;   GY94[36,50] = t*w;  GY94[37,5] = k*t*w
    GY94[37,21] = t*w;   GY94[37,33] = t*w;     GY94[37,38] = t;   GY94[37,39] = k*t
      GY94[37,40] = t;   GY94[37,41] = t*w; GY94[37,45] = k*t*w;   GY94[37,51] = t*w
   GY94[38,6] = k*t*w;   GY94[38,22] = t*w;   GY94[38,34] = t*w;     GY94[38,37] = t
      GY94[38,39] = t;   GY94[38,40] = k*t;   GY94[38,42] = t*w; GY94[38,46] = k*t*w
    GY94[38,52] = t*w;  GY94[39,7] = k*t*w;   GY94[39,23] = t*w;   GY94[39,35] = t*w
    GY94[39,37] = k*t;     GY94[39,38] = t;     GY94[39,40] = t;   GY94[39,43] = t*w
  GY94[39,47] = k*t*w;   GY94[39,53] = t*w;  GY94[40,8] = k*t*w;   GY94[40,24] = t*w
    GY94[40,36] = t*w;     GY94[40,37] = t;   GY94[40,38] = k*t;     GY94[40,39] = t
    GY94[40,44] = t*w; GY94[40,48] = k*t*w;   GY94[40,54] = t*w;  GY94[41,9] = k*t*w
    GY94[41,25] = t*w; GY94[41,33] = k*t*w;   GY94[41,37] = t*w;     GY94[41,42] = t
    GY94[41,43] = k*t;     GY94[41,44] = t;   GY94[41,45] = t*w; GY94[42,10] = k*t*w
    GY94[42,26] = t*w; GY94[42,34] = k*t*w;   GY94[42,38] = t*w;     GY94[42,41] = t
      GY94[42,43] = t;   GY94[42,44] = k*t;   GY94[42,46] = t*w;   GY94[42,55] = t*w
  GY94[43,11] = k*t*w;   GY94[43,27] = t*w; GY94[43,35] = k*t*w;   GY94[43,39] = t*w
    GY94[43,41] = k*t;     GY94[43,42] = t;     GY94[43,44] = t;   GY94[43,47] = t*w
    GY94[43,56] = t*w; GY94[44,12] = k*t*w;   GY94[44,28] = t*w; GY94[44,36] = k*t*w
    GY94[44,40] = t*w;     GY94[44,41] = t;   GY94[44,42] = k*t;     GY94[44,43] = t
    GY94[44,48] = t*w;   GY94[44,57] = t*w; GY94[45,13] = k*t*w;   GY94[45,29] = t*w
    GY94[45,33] = t*w; GY94[45,37] = k*t*w;   GY94[45,41] = t*w;     GY94[45,46] = t
    GY94[45,47] = k*t;     GY94[45,48] = t;   GY94[45,58] = t*w; GY94[46,14] = k*t*w
    GY94[46,30] = t*w;   GY94[46,34] = t*w; GY94[46,38] = k*t*w;   GY94[46,42] = t*w
      GY94[46,45] = t;     GY94[46,47] = t;   GY94[46,48] = k*t;   GY94[46,59] = t*w
  GY94[47,15] = k*t*w;   GY94[47,31] = t*w;   GY94[47,35] = t*w; GY94[47,39] = k*t*w
    GY94[47,43] = t*w;   GY94[47,45] = k*t;     GY94[47,46] = t;     GY94[47,48] = t
    GY94[47,60] = t*w; GY94[48,16] = k*t*w;   GY94[48,32] = t*w;   GY94[48,36] = t*w
  GY94[48,40] = k*t*w;   GY94[48,44] = t*w;     GY94[48,45] = t;   GY94[48,46] = k*t
      GY94[48,47] = t;   GY94[48,61] = t*w;    GY94[49,2] = t*w; GY94[49,18] = k*t*w
    GY94[49,34] = t*w;   GY94[49,50] = k*t;   GY94[49,52] = t*w; GY94[49,55] = k*t*w
    GY94[49,59] = t*w;    GY94[50,4] = t*w; GY94[50,20] = k*t*w;   GY94[50,36] = t*w
    GY94[50,49] = k*t;   GY94[50,54] = t*w; GY94[50,57] = k*t*w;   GY94[50,61] = t*w
     GY94[51,5] = t*w; GY94[51,21] = k*t*w;   GY94[51,37] = t*w;     GY94[51,52] = t
    GY94[51,53] = k*t;     GY94[51,54] = t; GY94[51,58] = k*t*w;    GY94[52,6] = t*w
  GY94[52,22] = k*t*w;   GY94[52,38] = t*w;   GY94[52,49] = t*w;     GY94[52,51] = t
      GY94[52,53] = t;   GY94[52,54] = k*t;   GY94[52,55] = t*w; GY94[52,59] = k*t*w
     GY94[53,7] = t*w; GY94[53,23] = k*t*w;   GY94[53,39] = t*w;   GY94[53,51] = k*t
      GY94[53,52] = t;     GY94[53,54] = t;   GY94[53,56] = t*w; GY94[53,60] = k*t*w
     GY94[54,8] = t*w; GY94[54,24] = k*t*w;   GY94[54,40] = t*w;   GY94[54,50] = t*w
      GY94[54,51] = t;   GY94[54,52] = k*t;     GY94[54,53] = t;   GY94[54,57] = t*w
  GY94[54,61] = k*t*w;   GY94[55,10] = t*w; GY94[55,26] = k*t*w;   GY94[55,42] = t*w
  GY94[55,49] = k*t*w;   GY94[55,52] = t*w;   GY94[55,56] = t*w;   GY94[55,57] = k*t
    GY94[55,59] = t*w;   GY94[56,11] = t*w; GY94[56,27] = k*t*w;   GY94[56,43] = t*w
    GY94[56,53] = t*w;   GY94[56,55] = t*w;   GY94[56,57] = t*w;   GY94[56,60] = t*w
    GY94[57,12] = t*w; GY94[57,28] = k*t*w;   GY94[57,44] = t*w; GY94[57,50] = k*t*w
    GY94[57,54] = t*w;   GY94[57,55] = k*t;   GY94[57,56] = t*w;   GY94[57,61] = t*w
    GY94[58,13] = t*w;   GY94[58,29] = k*t;   GY94[58,45] = t*w; GY94[58,51] = k*t*w
    GY94[58,59] = t*w;   GY94[58,60] = k*t;   GY94[58,61] = t*w;   GY94[59,14] = t*w
  GY94[59,30] = k*t*w;   GY94[59,46] = t*w;   GY94[59,49] = t*w; GY94[59,52] = k*t*w
    GY94[59,55] = t*w;   GY94[59,58] = t*w;   GY94[59,60] = t*w;   GY94[59,61] = k*t
    GY94[60,15] = t*w;   GY94[60,31] = k*t;   GY94[60,47] = t*w; GY94[60,53] = k*t*w
    GY94[60,56] = t*w;   GY94[60,58] = k*t;   GY94[60,59] = t*w;   GY94[60,61] = t*w
    GY94[61,16] = t*w; GY94[61,32] = k*t*w;   GY94[61,48] = t*w;   GY94[61,50] = t*w
  GY94[61,54] = k*t*w;   GY94[61,57] = t*w;   GY94[61,58] = t*w;   GY94[61,59] = k*t
    GY94[61,60] = t*w; diag(GY94) <- rowSums(GY94) * (-1);

# ><><                      Simulation Function Composition                      ><>< #
codonIndex <- 0
codonTriplets <- c()
stopCodons <- c("TAA","TAG","TGA")
nucleotides <- c("A", "C", "G", "T")

for (i in 1:4){   for (j in 1:4){   for (k in 1:4){
      newCodonXter <- paste0(nucleotides[i], nucleotides[j], nucleotides[k])
      if (!(newCodonXter %in% stopCodons)){
        codonIndex <- codonIndex + 1
        codonTriplets[codonIndex] <- newCodonXter
      }
}  } }

litRevSim <- function(nSite, nTaxa=1024, brLength=.05, eqFreqs=rep(1/61,61)){
  progenySize <- log(nTaxa, 2)
  fullLength <- brLength * progenySize
  alignmentMatrix <- matrix(NA, nTaxa, nSite)
  timedGY94 <- GY94 * fullLength; pMatrix <- expm(timedGY94)
  taxaNames <- paste0("Species", sprintf("%04.0f",1:nTaxa))
  rootSeqs <- sapply(1:nSite, function(a) sample(1:61,1,prob=eqFreqs))
  
  for(j in 1:nSite){
    alignmentMatrix[,j] <- sample(codonTriplets, nTaxa, T, pMatrix[rootSeqs[j],]) }
  
  seqFile <- paste0(getwd(), "/cranrSeqs.nex")
  openText <- paste0("  ", nTaxa, "  ", nSite*3, "\n")
  write.table(openText, seqFile, F, F, row.names=F, col.names=F)
  
  for(i in 1:nTaxa){
    newSequence <- paste0(taxaNames[i],"   ",paste0(alignmentMatrix[i,],collapse=""))
    write.table(newSequence, seqFile, T, F, row.names=F, col.names=F)  }
}

# ><><                         Simulation Implementation                         ><>< #
timeFile <- paste0(getwd(), "/cranrTimes.txt")
siteSizes <- c(5e+01, 5e+02, 5e+03, 5e+04, 5e+05)
write.table("", timeFile, F, F, row.names=F, col.names=F)

for(i in siteSizes){
  for(j in 1:51){
    startTime <- proc.time()["elapsed"]
    run <- litRevSim(i)
    cpuTime <- proc.time()["elapsed"] - startTime
    write.table(cpuTime, timeFile, T, F, row.names=F, col.names=F, eol="\t")
  }
  write.table("", timeFile, T, F, row.names=F, col.names=F)
  writeLines(paste0("Completed simulations for size ", i, "!."))
}

# ><>< ========================================================================= ><>< #
# ><><                              CODE ENDS HERE.                              ><>< #
# ><>< ========================================================================= ><>< #
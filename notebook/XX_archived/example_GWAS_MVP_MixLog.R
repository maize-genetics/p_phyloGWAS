#'==================================================================================================
#'
#' Title......: GWAS example using rMVP and mixed logistic regression (for binary traits)
#' Author.....: Germano Costa Neto
#' Date.......: 2022-08-01
#' Modified...: 2022-08-05
#' 
#'==================================================================================================

rm(list = ls())

require(tidyverse)
require(reshape2)
require(plyr)
require(rMVP)      # R version of the MVP package
require(milorGWAS) # mixed logistic regression, see Slack pan_and for reference

home.dir = getwd()
setwd(home.dir)


#'-----------------------------------------------
# #### Data set                            #####
#'-----------------------------------------------


# Data from Hackathon
pheno      = read.csv("Original99_phenomatrix.csv")
K_matrix   = read.csv('Original99_kmatrix.csv'    )
PAV_matrix = read.csv('Original99_pavMatrix.csv'  )


# NOT RUN, go to the line 49  map_panand2 = readRDS('chr_map_panand') #######
map_panand = ape::read.gff('Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3')

head(map_panand)

map_panand <- data.frame(map_panand,colsplit(map_panand$attributes,pattern = ';',names = c('ID','Parent','Protein')))


map_panand$ID=gsub(map_panand$ID,pattern = 'ID=',replacement = '')
map_panand$ID=gsub(map_panand$ID,pattern = 'Parent=',replacement = '')

map_panand <- map_panand %>% filter(type %in% 'mRNA')

map_panand2 <- data.frame(gene=map_panand$ID,chr = map_panand$seqid,pos=map_panand$start)
head(map_panand2)
names(map_panand2 )[3]='POS'


saveRDS(object = map_panand2,file = 'chr_map_panand')
##############

map_panand2 = readRDS('chr_map_panand')

# removing scafolds
levels(map_panand2$chr)
levels(map_panand2$chr)[11:nlevels(map_panand2$chr)] = 'scaf'
map_panand2 <- map_panand2 %>% filter(!chr %in% 'scaf') %>% droplevels()
dim(map_panand2)

map_panand2$gene = as.character(map_panand2$gene)


# organizing the phenotypic inputs

perenniality = data.frame(sample=pheno$assembly_name,perenniality=pheno$perenniality)


perenniality=perenniality[match(perenniality$sample,row.names(PAV_matrix)),]
perenniality[is.na(perenniality$perenniality),]

# I searched in google and it seems these two species are perennial
perenniality$perenniality[is.na(perenniality$perenniality)]   = 1

dim(PAV_matrix)

PAV_matrix_noscaf = PAV_matrix[,which(colnames(PAV_matrix) %in% map_panand2$gene)]
dim(PAV_matrix_noscaf)
dim(map_panand2)

colnames(K_matrix) = rownames(K_matrix)
K_matrix92 = K_matrix[which(rownames(K_matrix) %in% perenniality$sample),which(colnames(K_matrix) %in% perenniality$sample)]

K_matrix92 = data.matrix(K_matrix92)
dim(K_matrix92 )

dim(map_panand2)
dim(PAV_matrix_noscaf)

map_panand2 <- map_panand2 %>% filter(gene %in% colnames(PAV_matrix_noscaf))
dim(map_panand2)


#'-----------------------------------------------
# #### MVP                           #####
#'-----------------------------------------------
#' GCN - Fast and flexible. My plan is to use it for envGWAS not for perenniality
#' It seems GLM and MLM, which are based on a normal distribution, is not suitable for binary traits
#' 
trait = 'perenniality_v7'
path = paste0(home.dir,'/',trait)
dir.create(path)
setwd(path)


# require(AGHmatrix) I did this test to check if the issue was in the Kmatrix. It seems it is not.
#K_matrix92 = Gmatrix(SNPmatrix = data.matrix(PAV_matrix_noscaf),ploidy = 2)

# preparing the inputs (PAV, K matrix and map) according to MVP needs

genotype_PAV     = as.big.matrix(t(PAV_matrix_noscaf)) # dimension: PAV x samples
genotype_Kmatrix = as.big.matrix(K_matrix92)           # dimension: samples x samples

# perenniality  # data.frame, dimension: samples x (sample id + traits)
# map_panand2   # data.frame, dimension: PAVs x 3 (gene, chr and position)

#head(map_panand2)
#map_panand2$chr = as.numeric(gsub(map_panand2$chr,pattern = 'chr',replacement = '')) # ignore it. does not matter

GWAS_perenniality <- rMVP::MVP(
                              # Data Inputs ----
                              phe  = perenniality, 
                              geno = genotype_PAV,
                              K    = genotype_Kmatrix,
                              map  = map_panand2,
                              
                              # Method ------
                              method=c("GLM", "MLM", "FarmCPU"), 
                              
                              # Bin parameters ----
                              # method.bin = 'FaST-LMM', 
                              #  bin.selection = seq(1, 100, 1),
                              
                              # Number of Principal components for each mode ---
                              nPC.FarmCPU = 2,nPC.GLM = 2,nPC.MLM = 2,
                              
                              # Threshold ----
                              permutation.threshold = T, threshold = 0.05, maxLoop=3,
                              
                              # Outputs and computing parameters
                              file.output=T, ncpus=1)

# check your path directory


#'-----------------------------------------------
# #### mixed logistic regression                          #####
#'-----------------------------------------------
#'
#' GCN  - suitable for binary traits.
#' 


trait = 'perenniality_v7_miloreg'
path = paste0(home.dir,'/',trait)
dir.create(path)
setwd(path)


head(map_panand2)

# GCN I am just following the same object naming the authors use

TTN.gen2 = data.matrix(PAV_matrix_noscaf) 
dim(TTN.gen2)

# GCN: they ask for a family pedigree data
TTN.fam2 = data.frame(famid=rownames(TTN.gen2),id=rownames(TTN.gen2),father=0,mother=0,sex=0,pheno=NA)

TTN.bim2 = data.frame(chr=as.character(map_panand2$chr),id =map_panand2$gene,dist=0,pos=map_panand2$POS,A1=NA,A2=NA)
TTN.bim2$chr=as.numeric(gsub(TTN.bim2$chr,pattern = 'chr',replacement = '')) # need to be numeric


# now you put everything together (????)
x <- as.bed.matrix(TTN.gen2, TTN.fam2, TTN.bim2)
x

# and create your matrix of incidence (intercept, in this case)
X = matrix(1,nrow=length(perenniality$perenniality))

# and compute eigenvalues
eigenK = eigen(K_matrix92)


model_1  <- association.test(
                             # Data Inputs ----
                                x = x, # general inputs for map, PAV and pedigree
                                X = X, # indicende matrices for fixed effects. Here is an example of using only the intercept
                                Y = perenniality$perenniality, # vector of phenotypes
                                K = K_matrix92,   # Kinship matrix
                             
                             eigenK = eigenK,p = 2, #number of components
                             
                             # statistical parameters
                             response = 'bin', end = ncol(PAV_matrix_noscaf),# binary
                             test = 'wald',
                       )
head(model_1 )



model_1  %>% ggplot(aes(x=pos,y=freqA2,colour=as.factor(chr)))+
  geom_jitter()+facet_grid(~as.factor(chr),scales = 'free',space = 'free')+
  theme_classic()+theme(axis.text.x = element_text(size=9,angle=90,hjust = 1))


model_1  %>%unique() %>%  ggplot(aes(x=pos,y=-log10(p),colour=as.factor(chr)))+
  geom_jitter()+facet_grid(~as.factor(chr),scales = 'free',space = 'free')+
  theme_classic()+theme(axis.text.x = element_text(size=9,angle=90,hjust = 1))+
  geom_hline(yintercept = -log10(0.05/ncol(PAV_matrix_noscaf)),colour='red',size=1.1)

ggsave(plot = model_1  %>%unique() %>%  ggplot(aes(x=pos,y=-log10(p),colour=as.factor(chr)))+
         geom_jitter()+facet_grid(~as.factor(chr),scales = 'free',space = 'free')+
         theme_classic()+theme(axis.text.x = element_text(size=9,angle=90,hjust = 1))+
         geom_hline(yintercept = -log10(0.05/ncol(PAV_matrix_noscaf)),colour='red',size=1.1),filename = 'test.png')



# GSE155882 validation — data acquisition and sample audit
# Reconstructed from the executed analysis transcript retained with the project.
# Raw GEO files are downloaded directly from NCBI GEO.

library(Seurat)
library(Matrix)

project_dir <- Sys.getenv("GSE155882_PROJECT_DIR", unset = ".")
raw_dir <- file.path(project_dir, "data", "raw")
results_dir <- file.path(project_dir, "results")
dir.create(raw_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

sample_manifest <- data.frame(
  GSM=c("GSM4715045","GSM4715046","GSM4715047","GSM4715048"),
  Sample_ID=c("Sham_Rep1","Sham_Rep2","TAC_Rep1","TAC_Rep2"),
  Condition=c("Sham","Sham","TAC","TAC"),
  Replicate=c(1,2,1,2),
  stringsAsFactors=FALSE
)
write.csv(sample_manifest,file.path(results_dir,"V01_Primary_Validation_Sample_Manifest.csv"),row.names=FALSE)

prefix_map <- c(
 GSM4715045="Sham_scRNA_rep1", GSM4715046="Sham_scRNA_rep2",
 GSM4715047="TAC_scRNA_rep1", GSM4715048="TAC_scRNA_rep2"
)
file_manifest <- expand.grid(
 GSM=sample_manifest$GSM,
 File_Type=c("barcodes","features","matrix"),
 stringsAsFactors=FALSE
)
file_manifest$Sample_ID <- sample_manifest$Sample_ID[match(file_manifest$GSM,sample_manifest$GSM)]
file_manifest$Filename <- mapply(function(gsm,type){
 ext <- if(type=="matrix") "matrix.mtx.gz" else paste0(type,".tsv.gz")
 paste0(gsm,"_",prefix_map[[gsm]],"_",ext)
},file_manifest$GSM,file_manifest$File_Type)
file_manifest$URL <- paste0(
 "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
 substr(file_manifest$GSM,1,7),"nnn/",file_manifest$GSM,"/suppl/",file_manifest$Filename
)
file_manifest$Local_Path <- file.path(raw_dir,file_manifest$Filename)

options(timeout=3600)
for(i in seq_len(nrow(file_manifest))){
 if(!file.exists(file_manifest$Local_Path[i])){
   download.file(file_manifest$URL[i],file_manifest$Local_Path[i],mode="wb")
 }
}
write.csv(file_manifest,file.path(results_dir,"V01_GEO_File_Manifest.csv"),row.names=FALSE)

read_geo_sample <- function(gsm,sample_id,condition,replicate_number){
 sf <- file_manifest[file_manifest$GSM==gsm,]
 counts <- ReadMtx(
  mtx=sf$Local_Path[sf$File_Type=="matrix"],
  cells=sf$Local_Path[sf$File_Type=="barcodes"],
  features=sf$Local_Path[sf$File_Type=="features"],
  feature.column=2, unique.features=TRUE
 )
 obj <- CreateSeuratObject(counts=counts,project="GSE155882",min.cells=0,min.features=0)
 obj$GSE <- "GSE155882"; obj$GSM <- gsm; obj$sample_id <- sample_id
 obj$condition <- condition; obj$replicate <- replicate_number
 obj$validation_role <- "Independent_validation"
 obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern="^mt-")
 obj
}

raw_objects <- list()
for(i in seq_len(nrow(sample_manifest))){
 s <- sample_manifest[i,]
 raw_objects[[s$Sample_ID]] <- read_geo_sample(s$GSM,s$Sample_ID,s$Condition,s$Replicate)
}

shared_genes <- Reduce(intersect,lapply(raw_objects,rownames))
priority_genes <- c("Cyp1b1","Cyp2j6","Cyp2j9","Cyp4f13","Cyp4f16","Cyp4f17","Cyp4f18","Ephx2")
priority_availability <- data.frame(Gene=priority_genes,Present_All_4=priority_genes %in% shared_genes)
write.csv(priority_availability,file.path(results_dir,"V01_Priority_CYP_AA_Target_Availability.csv"),row.names=FALSE)

strict_cyp <- sort(grep("^Cyp[0-9]",shared_genes,value=TRUE))
write.csv(data.frame(Gene=strict_cyp),file.path(results_dir,"V01_Strict_CYP_Genes_Shared.csv"),row.names=FALSE)

saveRDS(raw_objects,file.path(project_dir,"GSE155882_raw_objects.rds"))

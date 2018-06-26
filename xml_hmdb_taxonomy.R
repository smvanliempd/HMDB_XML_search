
# WARNING: if run from rstudio, memory is not released and tends to crash R/Rstudio
# If ran from Rstudio make sure to resart R after script is finished
# Better run as Rscript from command line

# Execute with Rscript on command line (HMDB_selection.txt contains a list of HMDB accession numbers)
# All files (xml, HMDB_selection.txt, script) should be in the same directory
# Make sure that dir of Rscript.exe is in startup path

# To run : Rscript xml_hmdb_taxonomy.R [HMDB accession numebers] [database] [save dir]
# example: Rscript xml_hmdb_taxonomy.R HMDB_selection.txt urine C:/foo/bar


require(XML)
require(data.table)

# the data retrieval function
get_hmdb <- function(file , mets){
  out <- list()
  i   <- 1L
  cnt <- 1L
  get_data <- function( ) {
    metabolite <- function(n) {
      
      # check progress (for % there are a total of 114100 mets in db)
      cnt <<- cnt + 1L
      if (cnt %% 100L == 0L) {
        message(cnt)
      }
      
      # check accession number against selected values, if true get xpath data.
      # gc(reset = T) ##SLOOOWWWWW
      acc <- xpathSApply(n, "//metabolite/accession", xmlValue)
      if (acc %in% mets){
        
        # get values
        name                <- xpathSApply(n, "//metabolite/name", xmlValue) 
        kegg_id             <- xpathSApply(n, "//metabolite/kegg_id", xmlValue)
        kegg_map_id         <- xpathSApply(n, "//metabolite/pathways/pathway/kegg_map_id", xmlValue)
        smpdb_id            <- xpathSApply(n, "//metabolite/pathways/pathway/smpdb_id", xmlValue)
        direct_parent       <- xpathSApply(n, "//metabolite/taxonomy/direct_parent", xmlValue)
        kingdom             <- xpathSApply(n, "//metabolite/taxonomy/kingdom", xmlValue)
        super_class         <- xpathSApply(n, "//metabolite/taxonomy/super_class", xmlValue)
        class               <- xpathSApply(n, "//metabolite/taxonomy/class", xmlValue)
        sub_class           <- xpathSApply(n, "//metabolite/taxonomy/sub_class", xmlValue)
        molecular_framework <- xpathSApply(n, "//metabolite/taxonomy/molecular_framework", xmlValue)
        
        # set to NA when node is empty
        name                <- ifelse(length(name) == 0L,                "", name)
        kegg_id             <- ifelse(length(kegg_id) == 0L,             "", kegg_id)
        direct_parent       <- ifelse(length(direct_parent) == 0L,       "", direct_parent)
        kingdom             <- ifelse(length(kingdom) == 0L,             "", kingdom)
        super_class         <- ifelse(length(super_class) == 0L,         "", super_class)
        class               <- ifelse(length(class) == 0L,               "", class)
        sub_class           <- ifelse(length(sub_class) == 0L,           "", sub_class)
        molecular_framework <- ifelse(length(molecular_framework) == 0L, "", molecular_framework)
        smpdb_id            <- if(length(smpdb_id) == 0L)    "" else smpdb_id[smpdb_id != ""]
        kegg_map_id         <- if(length(kegg_map_id) == 0L) "" else kegg_map_id[kegg_map_id != ""]
        
        # write dataframe to list
        # if more than one value is returned use list() , e.g. kegg_map_id
        out[[i]] <<- data.table(name                = name , 
                                accession_number    = acc , 
                                smpdb_id            = list(smpdb_id),
                                kegg_id             = kegg_id,
                                kegg_map_id         = list(kegg_map_id),
                                direct_parent       = direct_parent,
                                kingdom             = kingdom,
                                super_class         = super_class,
                                class               = class,
                                sub_class           = sub_class,
                                molecular_framework = molecular_framework)
        i <<- i + 1
        rm(name,kegg_id,direct_parent,kingdom,super_class,class,sub_class,molecular_framework,smpdb_id,kegg_map_id)
      }
      
      # free nodes from memory and remove objects
      removeNodes(n, free = T)
      rm(n,acc)
      
    }
    branches <- list(metabolite = metabolite)
    return(branches)
  }
  
  # the SAX parser
  xmlEventParse(file, handlers = list(), branches = get_data())
  
  # make 
  out <- do.call(rbind, out)
  return(out)
}

# arguments for Rscript
args <- commandArgs(trailingOnly = T)

# selected accession numbers
mets_selct <- scan(args[1], what = "character")

# the available databases
if (args[2] == "sweat") xmldoc <- "sweat_metabolites.xml"
if (args[2] == "urine") xmldoc <- "urine_metabolites.xml"
if (args[2] == "all")   xmldoc <- "hmdb_metabolites.xml"

# write results to object
hmdb_tax_data <- get_hmdb(xmldoc,mets_selct)

# save rds file
saveRDS(hmdb_tax_data,paste0(args[3],"/HMDB_out.rds"))



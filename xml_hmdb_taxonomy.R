# TODO set argument for HMDB database directory and database to use

# Execute with Rscript on command line (HMDB_selection.txt contains a list of HMDB accession numbers)
# All files (xml, HMDB_selection.txt, script) should be in the same directory
# Make sure that Rscript.exe is in PATH

# To run:
# Rscript xml_hmdb_taxonomy.R HMDB_selection.txt



require(XML)

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
        direct_parent       <- xpathSApply(n, "//metabolite/taxonomy/direct_parent", xmlValue)
        kingdom             <- xpathSApply(n, "//metabolite/taxonomy/super_class", xmlValue)
        super_class         <- xpathSApply(n, "//metabolite/taxonomy/super_class", xmlValue)
        class               <- xpathSApply(n, "//metabolite/taxonomy/class", xmlValue)
        sub_class           <- xpathSApply(n, "//metabolite/taxonomy/sub_class", xmlValue)
        molecular_framework <- xpathSApply(n, "//metabolite/taxonomy/molecular_framework", xmlValue)
        
        # set to NA when node is empty
        name                <- ifelse(length(name) == 0L,NA,name)
        direct_parent       <- ifelse(length(direct_parent) == 0L,NA,direct_parent)
        kingdom             <- ifelse(length(kingdom) == 0L,NA,kingdom)
        super_class         <- ifelse(length(super_class) == 0L,NA,super_class)
        class               <- ifelse(length(class) == 0L,NA,class)
        sub_class           <- ifelse(length(sub_class) == 0L,NA,sub_class)
        molecular_framework <- ifelse(length(molecular_framework) == 0L,NA,molecular_framework)
        
        # write dataframe to list
        out[[i]] <<- data.frame( accession_number    = acc , 
                                 name                = name , 
                                 direct_parent       = direct_parent,
                                 kingdom             = kingdom,
                                 super_class         = super_class,
                                 class               = class,
                                 sub_class           = sub_class,
                                 molecular_framework = molecular_framework,
                                 stringsAsFactors = F)
        i <<- i + 1
        rm(name,direct_parent,kingdom,super_class,class,sub_class,molecular_framework)
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
  out <- do.call(rbind,lapply(out, unlist))
  return(out)
}

# arguments for Rscript
args <- commandArgs(trailingOnly = T)

# selected accession numbers
mets_selct <- scan(args[1], what = "character")

# the available databases
# xmldoc <- "sweat_metabolites.xml"
# xmldoc <- "urine_metabolites.xml"
xmldoc <- "hmdb_metabolites.xml"

hmdb_tax_data <- get_hmdb(xmldoc,mets_selct)

# separation with semi colon because commas present in text strings
write.table(hmdb_tax_data,"HMDB_out.csv",quote = F,sep=";",row.names = F)




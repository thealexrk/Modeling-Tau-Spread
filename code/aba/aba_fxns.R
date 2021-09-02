ontology.name.children <- function(ontology){ 
  # INPUTS:
  # ontology: nested list from ABA structure ontology
  
  # OUTPUTS:
  # ontology: input ontology list, with children recursively renamed for their acronyms
  
  names(ontology$children) <- sapply(ontology$children, function(x) x$acronym)
  ontology$children <- lapply(ontology$children, function(x) ontology.name.children(x))
  
  return(ontology)
}

ontology.remove.fluff <- function(ontology,unique.identifier = 'ElIcOrNbLaTh'){
  # INPUTS:
  # ontology: nested list from ABA structure ontology
  
  # OUTPUTS:
  # ontology: input ontology list, with all fields removed except for id, children and acronym
  
  ontology[!names(ontology) %in% c('acronym','id','children')] <- NULL
  names(ontology)[names(ontology) == 'id'] <- paste0('id',unique.identifier) # give it a unique identifier so no chance of ambiguity later
  names(ontology)[names(ontology) == 'acronym'] <- paste0('acronym',unique.identifier) # give it a unique identifier so no chance of ambiguity later
  ontology$children <- lapply(ontology$children, function(x) ontology.remove.fluff(x))
  return(ontology)
}

ontology.get.name.id.key <- function(name.vec,df.key){
  # INPUTS
  # name.vec: vector of region names. assumes that each name has i or c as first character, which will get cut off in matching to ABA
  # df.key: data frame with structure_id column containing numeric ABA structure_ids and acronym column containing character structure name
  # OUTPUTS:
  # names.id.key: list whose element names are the names in name.vec, and whose elements contain corresponding structure ids to that name
  
  names.id.key <- as.list(rep(NA,length(name.vec))) # default id is NA
  names(names.id.key) <- name.vec
  for(i.name in name.vec){ # get structure_id corresponding to each name and store in the corresponding list element
    i.name.id <- df.key$structure_id[df.key$acronym == substr(i.name,start=2,stop=nchar(i.name))]
    if(length(i.name.id) > 0){names.id.key[i.name] <- i.name.id} # leave as NA if name doesn't exist in acronym
  }
  return(names.id.key)
}

ontology.assign.expression <- function(names.id.key,goi.exp,expression.metric = 'expression_energy'){
  # INPUTS:
  # names.id.key: list whose element names are the names in name.vec, and whose elements contain corresponding structure ids to that name
  # goi.exp: dataframe containing expression values for each structure id
  # expression.metric: character corresponding to column of goi.exp that names a metric of expression to extract
  # ^ i.e. energy, intensity, density
  
  # OUTPUTS:
  # expression: named vector of numeric gene expression values for given expression metric and each region in names.id.key
  # NA means names.id.key had no match to a structure_id to begin with, or that the region's structure_id was not found
  # in the expression .csv file
  # if a region's structure_id was found multiple times, then average over the expression values
  
  expression <- rep(NA,length(names.id.key))
  names.id.key <- unlist(names.id.key)
  names(expression) <- names(names.id.key)
  for(region in names(names.id.key[!is.na(names.id.key)])){
    region.expression <- goi.exp[goi.exp$structure_id == names.id.key[[region]],expression.metric]
    if(length(region.expression) == 1){expression[region] <- region.expression}
    if(length(region.expression) == 0){expression[region] <- NA}
    if(length(region.expression) > 1){expression[region] <- mean(region.expression)}
  }
  return(expression)
}

ontology.get.expression.df <- function(name.vec,names.id.key,gois,datasets,download=TRUE){
  # INPUTS: 
  # name.vec: vector of region names. assumes that each name has i or c as first character, which will get cut off in matching to ABA
  # names.id.key: list whose element names are the names in name.vec, and whose elements contain corresponding structure ids to that name
  # gois: list of lists for each gene containing:
      # $goi: character vector of names of genes of interest corresponding to ABA gene names
      # $probe: character vector of probe name (pick one)
  # datasets: dataframe containing urls for every gene
  # download: logical indicating whether to download gene expression or not... if already done can skip for speed's sake
  
  # short summary:
  # get gene expression for regions of interest and genes of interest, given mapping between region names and structure_ids, and vector of gene names
  
  # OUTPUTS:
  # df.expression: data frame of gene expression values for each region of interest (row) and each gene of interest (column)
  
  goi.names <- sapply(gois, function(X) paste0(X$goi,'_',X$probe))
  df.expression <- as.data.frame(matrix(NA,nrow=length(name.vec),ncol=length(gois),dimnames = list(name.vec,goi.names)))
  for(goi.info in gois){
    goi <- goi.info$goi
    probe <- goi.info$probe
    goi_mask <- datasets$gene_symbol == goi  # select rows of dataset table that correspond to gene of interest
    if(length(probe)>0){goi_mask <- goi_mask & datasets$probe_name == probe} # if probe input is given, then select specific probe
    goi_urls <- as.character(datasets$structure_unionizes_file_url[goi_mask]) # get urls to experiments measuring expression for gene of interest
    if(length(goi_urls) == 0){
      print('Gene not found in datasets file.')
      proceed <- NA
    } else if(download){
      print(datasets[goi_mask,c('data_set_id','probe_name','gene_symbol')])
      proceed <- tryDownload(url = goi_urls[1], destfile = paste0(abadir,'expression/',goi,probe,'.csv'))
    } else if(!download & length(goi_urls) > 0){
      proceed <- 0 # if you have downloaded gene and it exists, then go ahead
    }
    if(proceed == 0){
        goi.exp <- read.csv(file = paste0(abadir,'expression/',goi,probe,'.csv'))
        df.expression[,paste0(goi,'_',probe)] <- ontology.assign.expression(names.id.key,goi.exp) # default is expression energy
      } else {
        print('bad URL')
    }
  }
  return(df.expression)
}

tryDownload <- function(url,destfile){
  # INPUTS:
  # url: character of url to file
  # destfile: character of path for where to save file
  
  # OUTPUTS:
  # 'good' if it worked, 'bad' if it didn't
  out <- tryCatch({
    download.file(url = url,destfile = destfile)
  }, # try to download expression for gene of interest
  error= function(cond) {return('bad')})
  return(out)
}

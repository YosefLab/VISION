library(jug)
library(jsonlite)

# Can use res$text(body) to set body of response
# Can use res$json(obj) to convert object to json
#  and return

# Options for 'serve_it'
# serve_it<-function(jug,
#    host="127.0.0.1",
#    port=8080,
#    daemonized=FALSE,
#    verbose=FALSE)

# ---------------
# Sample API BELOW
# ----------------

fpout <- readRDS(arg1)


browseURL("http://127.0.0.1:8080/FastProject_Output/html/Results.html")
# Launch the server
jug() %>%
  get("/Signature/Scores/(?<sig_name1>.*)", function(req, res, err) {
    sigMatrix <- fpout@sigMatrix
    name <- URLdecode(req$params$sig_name1)
    out <- "Signature does not exist!"
    if (name %in% rownames(sigMatrix)) {
      out <- sigScoresToJSON(sigMatrix[name,])
    }
    return(out)
  }) %>%
  get("/Signature/ListPrecomputed", function(req, res, err){
    signatures <- fpout@sigList
    keys <- lapply(signatures, function(x) x@name)
    vals <- lapply(signatures, function(x) x@isPrecomputed)
    names(vals) <- keys
    out <- toJSON(vals, auto_unbox=TRUE)
    return(out)
  }) %>%
  get("/Signature/Info/(?<sig_name2>.*)", function(req, res, err){
    signatures <- fpout@sigList
    name <- URLdecode(req$params$sig_name2)
    out <- "Signature does not exist!"
    for (sig in signatures) {
      if (sig@name == name) {
        out <- signatureToJSON(sig)
        break
      }
    }
    return(out)
  }) %>%
  get("/Signature/Ranks/(?<sig_name3>.*)", function(req, res, err) {
    sigMatrix <- fpout@sigMatrix
    name <- URLdecode(req$params$sig_name3)
    out <- "Signature does not exist!"
    if (name %in% rownames(sigMatrix)) {
      out <- sigRanksToJSON(sigMatrix[name,])
    }
    return(out)
  }) %>%
  get("/Signature/Expression/(?<sig_name4>.*)", function(req, res, err) {
    all_names = sapply(fpout@sigList, function(x){return(x@name)})
    name <- URLdecode(req$params$sig_name4)
    index = match(name, all_names)
    if(is.na(index)){
        out <- "Signature does not exist!"
    }
    else{
        sig = fpout@sigList[[index]]
        genes = names(sig@sigDict)
        expMat = fpout@exprData@data
        return(expressionToJSON(expMat, genes))
    }
    return(out)
  }) %>%
  get("/FilterGroup/(?<filter_name1>.*)/(?<proj_name1>.*)/coordinates", function(req, res, err) {
    projData <- fpout@projData
    filter <- URLdecode(req$params$filter_name1)
    proj <- URLdecode(req$params$proj_name1)
    out <- "Filter, Projection pair does not exist!"
    for (pd in projData) {
      if (pd@filter == filter) {
        for (p in pd@projections) {
          if (p@name == proj) {
            print('here')
            out <- coordinatesToJSON(p@pData)
            break
          }
        }
      }
    }
    return(out)
  }) %>%
  get("/FilterGroup/(?<filter_name2>.*)/SigProjMatrix", function(req, res, err) {
    projData <- fpout@projData
    filter <- URLdecode(req$params$filter_name2)
    out <- "Filter, Projection pair does not exist!"
    for (pd in projData) {
      if (pd@filter == filter) {
        out <- sigProjMatrixToJSON(pd@sigProjMatrix)
        break
      }
    }
    return(out)
  }) %>%
  get("/FilterGroup/(?<filter_name3>.*)/SigProjMatrix_P", function(req, res, err) {
    projData <- fpout@projData
    filter <- URLdecode(req$params$filter_name3)
    out <- "Filter, Projection pair does not exist!"
    for (pd in projData) {
      if (pd@filter == filter) {
        out <- sigProjMatrixToJSON(pd@pMatrix)
        break
      }
    }
    return(out)
  }) %>%
  get("/FilterGroup/(?<filter_name4>.*)/(?<proj_name2>.*)/clusters/(?<cluster_procedure>.*)/(?<param>.*)", function(req, res, err) {
    projData <- fpout@projData

    filter <- URLdecode(req$params$filter_name4)
    proj <- URLdecode(req$params$proj_name2)
    method <- URLdecode(req$params$cluster_procedure)
    param <- as.numeric(URLdecode(req$params$param))

    out <- "Filter, Projection pair does not exist!"
    for (pd in projData) {
      if (pd@filter == filter) {
        for (projection in pd@projections) {
          if (projection@name == proj) {
              clust = cluster(projection, method, param)
              out <- clusterToJSON(clust)
              break
          }
        }
      }
    }
    return(out)
  }) %>%
  get("/FilterGroup/(?<filter_name5>.*)/genes", function(req, res, err) {
    projData <- fpout@projData
    filter <- URLdecode(req$params$filter_name5)
    out <- "Filter, Projection pair does not exist!"
    for (pd in projData) {
      if (pd@filter == filter) {
        out <- toJSON(pd@genes)
        break
      }
    }
    return(out)
  }) %>%
  get("/FilterGroup/list", function(req, res, err) {
    projData <- fpout@projData
    filters <- sapply(projData, function(x){
                      return(x@filter);
                      });
    return(toJSON(filters))
  }) %>%
  get("/Expression", function(req, res, err) {
    return(expressionToJSON(fpout@exprData@data, matrix(NA)))
  }) %>%
  get("/path2", function(req, res, err){
    return("Hello Mars!")
  }) %>%
  get("/person/(?<your_name1>.*)/eats/(?<your_food>.*)",
      function(req, res, err){
        name = URLdecode(req$params$your_name1)
        food = URLdecode(req$params$your_food)
        if ('howmany' %in% names(req$params)) {
          howmany = URLdecode(req$params$howmany)
          out = paste0(name, ' would like to eat ',
                       howmany, " ", food, "s!")
        }
        else{
          out = paste(name, 'likes to eat', food)
        }
        return(out)
      }) %>%
  get("/person/(?<your_name2>.*)", function(req, res, err){
    out = paste0("Hello ", 
                 URLdecode(req$params$your_name2),
                 "!")
    return(out)
  }) %>%
  get("/genes", function(req, res, err){
    if ('number' %in% names(req$params)) {
      number = URLdecode(req$params$number)
      number = as.integer(number) # All params are string by default
      
      gene_list = c( "PPP2R1A", "CDK7", "SEM1", 
                     "PSMF1", "PSMB10", "DNM2", "CEP63", "CDC25B", 
                     "CDC25B", "CEP57", "RAB11A", "KHDRBS1", "RCC2", 
                     "ARPP19", "MUS81", "BACH1", "PSMA7", "CDC25C", 
                     "CDC25C", "PSMC5", "PLCB1", "PSMC1", "APP", 
                     "WEE1", "PPP1CB", "CCNB2", "DYNC1H1", "DCTN1", 
                     "PSMA6", "PBX1", "NAE1", "UBC", "UBC", 
                     "UBC", "FOXN3", "CEP41", "TUBG1", "CDK4", 
                     "HAUS6", "CDK6", "PSMA8", "CDK3", "WNT10B", 
                     "WNT10B", "PSMC4", "CLASP1", "PSMD4", "PSMB11", 
                     "KDM8", "ABCB1", "TICRR", "CENPJ", "SFI1", 
                     "DYNLL1", "PPP2R2A", "CHEK2", "YWHAG", "PKIA", 
                     "PSMD1", "PSMD7", "CDKN1A", "PHOX2B", "PSMB7", 
                     "NEK2", "MNAT1", "CCNH", "PHLDA1", "ORAOV1", 
                     "ATM", "SKP2", "CEP72")
      
      res$json(gene_list[1:number])
      
    }
    else{
      out = "ERROR!"
    }
    return(out)
  }) %>%
  serve_static_files("html") %>%
  serve_static_files("../html") %>%
  simple_error_handler_json() %>%
  serve_it()

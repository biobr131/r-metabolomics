# First create a list containing a vector of the compounds to be queried (separated by a semi-colon)  
# and another character vector containing the compound id type.
# The items in the list MUST be queryList and inputType
# Valid input types are: "name", "hmdb", "kegg", "pubchem", "chebi", "metlin"

name.vec <- c("1,3-Diaminopropane;2-Ketobutyric acid;2-Hydroxybutyric acid;2-Methoxyestrone")
toSend = list(queryList = name.vec, inputType = "name")

library(httr)

# The MetaboAnalyst API url
call <- "https://rest.xialab.ca/api/mapcompounds"

# Use httr::POST to send the request to the MetaboAnalyst API
# The response will be saved in query_results
query_results <- httr::POST(call, body = toSend, encode = "json")

# Check if response is ok (TRUE)
# 200 is ok! 401 means an error has occured on the user's end.
query_results$status_code==200

# Parse the response into a table
# Will show mapping to "hmdb_id", "kegg_id", "pubchem_id", "chebi_id", "metlin_id", "smiles" 
query_results_text <- content(query_results, "text", encoding = "UTF-8")
query_results_json <- rjson::fromJSON(query_results_text, flatten = TRUE)
query_results_table <- t(rbind.data.frame(query_results_json))
rownames(query_results_table) <- query_results_table[,1]

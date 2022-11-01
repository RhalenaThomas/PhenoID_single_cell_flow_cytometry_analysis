#### function to test out saving a list as manuliple worksheet in an excel doc


require("XLConnect")
wb <- loadWorkbook("var1.xlsx", create = TRUE)
createSheet(wb, names(var1))
writeWorksheet(wb, var1, names(var1),header=FALSE)
saveWorkbook(wb)



### another an way

write_list <-function(my_list, wb_name = 'var1.xlsx') {    
  wb <- loadWorkbook(wb_name, create = TRUE)
  createSheet(wb, names(my_list))
  writeWorksheet(wb, my_list, names(my_list),header=FALSE)
  saveWorkbook(wb)
}

## should add the name of df in the list (I think)

write_list_name <-function(.list, default = 'var1', path = ''){
  .name <- as.list(match.call())[2]
  if(is.language(.name[[1]])){
    wb_name <- sprintf("%s/%s.xlsx", path, default)
  }
  if(is.symbol(.name[[1]])) {
    wb_name <- sprintf("%s/%s.xlsx", path, as.character(.name))
  }
  wb <- loadWorkbook(wb_name, create = TRUE)
  createSheet(wb, names(.list))
  writeWorksheet(wb,.list, names(.list),header=FALSE)
  saveWorkbook(wb)
}
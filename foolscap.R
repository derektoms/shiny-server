> filter(gpl,gpl=='GPL1261') %>% select(gpl,title)
# Source:   lazy query [?? x 2]
# Database: sqlite 3.19.3 [/Volumes/ULTRA/across_array/GEOmetadb.sqlite]
  gpl     title                                             
  <chr>   <chr>                                             
1 GPL1261 [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
> filter(gpl,gpl=='GPL570') %>% select(gpl,title)
# Source:   lazy query [?? x 2]
# Database: sqlite 3.19.3 [/Volumes/ULTRA/across_array/GEOmetadb.sqlite]
  gpl    title                                                       
  <chr>  <chr>                                                       
1 GPL570 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array


gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL570") %>% 
  select(title,gse) %>%
  collect() %>% ## this was needed here to be able to filter a second time
  filter(str_detect(title, search_string))
 
 
  
  gse %>% 
    left_join(gse_gpl, copy = TRUE) %>% 
    filter(gpl == "GPL570" & str_detect(title, search_string)) %>% 
    select(title,gse) %>%
    collect() %>% ## this was needed here to actually be able to filter a second time
    filter(str_detect(title, search_string))
    
gse %>% 
    left_join(gse_gpl, copy = TRUE) %>% 
  select(title,gse) %>%
  filter(gpl == "GPL570" & str_detect(title, search_string)) %>% 
  
  collect() %>% ## this was needed here to actually be able to filter a second time
  filter(str_detect(title, search_string))
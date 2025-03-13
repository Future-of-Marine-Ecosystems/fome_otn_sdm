# Read in data
otn_read = function(directory, pattern = 'matched_detections'){
  
  # File list
  file_list = list.files(data_dir, 'matched_detections', recursive = TRUE)
  
  # load receiver data (OTN matched detections)
  rec_data = list()
  for(i in 1:length(file_list)){
    
    rec_data[[i]] = read_otn_detections(paste0(data_dir, file_list[i]))
    
  }
  
  rec_data = rbindlist(rec_data)
  
  return(rec_data)
  
} # end otn_read

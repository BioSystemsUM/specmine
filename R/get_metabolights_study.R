is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

get_metabolights_study=function(mtblsID, directory){

  if(!reticulate::py_available()){
    p=reticulate::py_discover_config()
    python_versions=p$python_versions
    python_path=""
    if(is_windows() & length(python_versions)>0) python_path=python_versions[length(python_versions)]
    else python_path=grep(".*python3.*", p$python_versions, value=T)
    if(length(python_path)==0) stop("Python version 3 must be installed")
    reticulate::use_python(python_path)
    reticulate::py_config()
  }
  
  isatools_modified_py_file=system.file("isatools_funcs_modified.py", package="specmine")
  get_data_metadata_py_file=system.file("get_data_metadata.py", package="specmine")
  
  reticulate::py_run_file(isatools_modified_py_file)
  reticulate::py_run_file(get_data_metadata_py_file)
  env=new.env()
  reticulate::source_python(get_data_metadata_py_file,  envir = env)
  
  cat("Downloading files... please wait\n")
  env$load_data_metadata_metabolights(mtblsID, directory)
  cat("\nDone.")
  
}

get_metabolights_study_list=function(){
  
  if(!reticulate::py_available()){
    p=reticulate::py_discover_config()
    python_versions=p$python_versions
    python_path=""
    if(is_windows() & length(python_versions)>0) python_path=python_versions[length(python_versions)]
    else python_path=grep(".*python3.*", p$python_versions, value=T)
    if(length(python_path)==0) stop("Python version 3 must be installed")
    reticulate::use_python(python_path)
    reticulate::py_config()
  }
  
  get_data_metadata_py_file=system.file("get_data_metadata.py", package="specmine")
  env=new.env()
  reticulate::source_python(get_data_metadata_py_file, envir=env)
  
  env$metabolights_studies_list()
}
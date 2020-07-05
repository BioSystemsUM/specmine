
# Internal data:
codes = readRDS('~/Documents/specmine/metabolomicspackage/system_data/codes.rda')
maps_con = readRDS('~/Documents/specmine/metabolomicspackage/system_data/maps_con.rda')
orgs_cpds = readRDS('~/Documents/specmine/metabolomicspackage/system_data/orgs_cpds.rda')
spectra_list = readRDS('~/Documents/specmine/metabolomicspackage/system_data/spectra_list.rda')

usethis::use_data(codes, maps_con, orgs_cpds, spectra_list, internal=T, overwrite = TRUE)


# Extermal data:
#spectra_options = readRDS('./data/spectra_options.rda')

#usethis::use_data(spectra_options, internal=F, overwrite=T)

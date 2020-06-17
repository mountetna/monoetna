packages = read.table('requirements.txt', header=T)

if ("CRAN" %in% packages$source) {
  for (package in packages[ packages$source == "CRAN", ]$package) {
    if (!package %in% rownames(installed.packages())) {
      install.packages(package, repos='http://cran.us.r-project.org')
    }
  }
}

if ("bioconductor" %in% packages$source) {
  source("https://bioconductor.org/biocLite.R")
  for (package in packages[ packages$source == "bioconductor", ]$package) {
    if (!package %in% rownames(installed.packages())) {
      biocLite(package)
    }
  }
}

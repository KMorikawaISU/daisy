usethis::use_description(fields = list(
  Title       = "Data-Adaptive Integration with Summary Data",
  Description = "Estimation and bootstrap tools for generalized entropy balancing that integrate individual-level data with biased summary data. The second-step KL setting follows the original implementation for numerical compatibility.",
  License     = "MIT + file LICENSE",
  URL         = "https://github.com/KMorikawaISU/daisy",
  BugReports  = "https://github.com/KMorikawaISU/daisy/issues"
))
usethis::use_mit_license("Kosuke Morikawa")

# 実行時に必要（Imports）
usethis::use_package("lamW", type = "Imports")
usethis::use_package("gsl",  type = "Imports")
usethis::use_package("stats",type = "Imports")

# あると便利（Suggests：なければフォールバック）
usethis::use_package("MASS",    type = "Suggests")  # rmvnorm 代替
usethis::use_package("mvtnorm", type = "Suggests")

# roxygen2 で Rd/NAMESPACE をMD記法で生成・テスト雛形
usethis::use_roxygen_md()
usethis::use_testthat()
usethis::use_readme_rmd()

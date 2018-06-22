fls <- list.files("R")
fls <- fls[!grepl('sysdata.rda', fls)]
tst <- list.files("tests/manual/test-code")
wd <- getwd()
coverage <- covr::file_coverage(file.path(wd,"R", fls), file.path(wd, 'tests/manual/test-code', tst))

class(coverage)

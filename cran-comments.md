
── R CMD check results ────────────────────── survstan 0.0.7.1 ────
Duration: 11m 49.1s

❯ checking compiled code ... OK
   WARNING
  ‘qpdf’ is needed for checks on size reduction of PDFs

❯ checking installed package size ... NOTE
    installed size is 78.9Mb
    sub-directories of 1Mb or more:
      libs  78.3Mb

❯ checking dependencies in R code ... NOTE
  Namespaces in Imports field not imported from:
    ‘RcppParallel’ ‘rstantools’
    All declared Imports should be used.

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors ✔ | 1 warning ✖ | 3 notes ✖

* Package also checked using GHA workflow.

* This is an update from version 0.0.6.1 to 0.0.7.

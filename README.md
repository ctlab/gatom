[![Travis-CI Build Status](https://travis-ci.org/ctlab/gatom.svg?branch=master)](https://travis-ci.org/ctlab/gatom)
[![codecov](https://codecov.io/gh/ctlab/gatom/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/gatom)


# gatom

An R-package for finding active metabolic modules in atom transition network.

---

* To work with **gatom** you need to download [sgmwcs solver](https://github.com/ctlab/sgmwcs-solver/releases) depending on [CPLEX solver](https://www.ibm.com/uk-en/marketplace/ibm-ilog-cplex) into PATH_TO_SGMWCS and PATH_TO_CPLEX respectively. Make sure you have downloaded cplex.jar and libcplex.so into one folder and both of the same version. Note that CPLEX solver is a proprietary software, still you can find community version on the official site.

* Create the following wrapper for your sgmwcs solver, name it `sgmwcs` for consistency with [gatom-tutorial.Rmd](https://github.com/ctlab/gatom/blob/master/vignettes/gatom-tutorial.Rmd#pre-generated-annotations), fill it out as shown below:
    ```{bash}
    #!/bin/sh
    exec java -Djava.library.path=PATH_TO_CPLEX \
        -Xmx2G \
        -cp PATH_TO_CPLEX/cplex.jar:/PATH_TO_SGMWCS/sgmwcs-solver.jar \
        ru.ifmo.ctddev.gmwcs.Main "$@"
    ```
    and put this file into the directory included in your $PATH.

* Install **gatom** as R-package running the following command inside R:
    ```{R}
    #install.packages("devtools")
    devtools::install_github("ctlab/gatom")
    ```


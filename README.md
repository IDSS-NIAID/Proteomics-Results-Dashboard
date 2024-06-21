# Proteomics-Results-Dashboard

## Conventions

We will generally follow the conventions described in [R Packages](https://r-pkgs.org/) by Hadley Wickham. [JoesFlow](https://github.com/niaid/JoesFlow) is set up using the same conventions we will be using here.

### Directories

* R: R code (i.e. functions)
* inst/extdata: Raw data files that are not generated by the package itself
* shiny: Final Shiny app code (i.e. not for Shiny components, they will reside in `R/`, but for the app that will be installed on the Shiny server)
* tests: Unit testing code
* docs: This is where web documentation (i.e. vignettes) will be generated when we get to that part.

Other directories (i.e. `man`) will be auto-generated during package compilation.

### Documentation

We will use [roxygen2](https://r-pkgs.org/man.html) tags for documentation.

### Building

We will use [devtools](https://bookdown.org/rdpeng/RProgDA/the-devtools-package.html) to help with package building. The following commands will be used:

```{r}
# documentation generation / updates
devtools::document()

# package checking
devtools::check()

# package testing (also run as part of `check`, but if you only want to run tests)
devtools::test()

# package installation
devtools::install()

# Once the package is ready, you can start the Shiny app with:
ProtResDash::run_app()
```
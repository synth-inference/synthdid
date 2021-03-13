### Developing

Contributions are welcome. This repository follows the standard open source protocol and setup with git where there is an abundance of existing resources to get up to speed. Condensed greatly, the workflow is to fork this repository, check out a branch, commit your changes (forming an ideally legible commit history), then submitting a pull request explaining your contribution, ideally referring to the issue you created, or the issue you chose to work on.


### Worfklow

Developing and building can be done through RStudio, or by using the `devtools` package as for example:

1. Edit code
2. Load everything into an R session with `devtools::load_all(".")`
3. Test your changes with `devtools::check()'. 
   This also builds auto-generated documentation, possibly raising warnings that arise from errors in your in-code Markdown. 
   For example: Square brackets are special syntax and must be escaped if meant literally: in Markdown, the unit interval is \[0,1\].
4.Commit the changes (including the generated documentation).

Vignettes for the online documentation are stored in `R/vignettes`. This page can be rendered locally using `pkgdown::build_site()`. The layout for this page is defined in `_pkgdown.yml`.

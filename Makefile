# Makefile that updates README from Rmd vignette
# adapted from https://github.com/daattali/shinyjs/blob/master/Makefile
# based on discussion https://community.rstudio.com/t/readme-md-vs-package-vignette-vs-package-documentation/1359/3

#!usr/bin/make -f
# All commands are run as R functions rather than shell commands so that it will work easily on any Windows machine, even if the Windows machine isn't properly set up with all the right tools

all: README.md

clean:
	Rscript -e 'suppressWarnings(file.remove("README.md", "vignettes/OCMSutility.md"))'

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

README.md : vignettes/OCMSutility.Rmd
#	echo "Rendering the OCMSutility vignette"
	Rscript -e 'rmarkdown::render("vignettes/OCMSutility.Rmd", output_format = "md_document", output_options = list(pandoc_args = c("-t", "markdown")))'
#	echo "Correcting image paths"
#	sed -i -- 's,../inst,inst,g' vignettes/OCMSutility.md
	Rscript -e 'file <- gsub("\\.\\./inst", "inst", readLines("vignettes/OCMSutility.md")); writeLines(file, "vignettes/OCMSutility.md")'
#	echo "Correcting paths to other reference vignettes"
	Rscript -e 'file <- gsub("\\((.*)\\.([rR]md)","(vignettes/\\1.\\2", readLines("vignettes/OCMSutility.md")); writeLines(file, "vignettes/OCMSutility.md")'
#	echo "Copying output to README.md"
#	cp vignettes/OCMSutility.md README.md
	Rscript -e 'file.copy("vignettes/OCMSutility.md", "README.md", overwrite = TRUE)'
	Rscript -e 'suppressWarnings(file.remove("vignettes/OCMSutlity.md"))'

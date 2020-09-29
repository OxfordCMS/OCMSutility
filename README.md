# OCMS_Sandbox
Home for wip pipelines, utility scripts/functions
* all code must have one-liner of what it does and who is working on it
* when done/ready for use, move code to respective utlity folders
* code in utility folders should have preamble documentation and use case example

## OCMS apps

To use scripts that are in the Py_utility directory you should be able to install and run. The commandline interface has been ripped from cgat-apps although is a more rudimentary version.


Install

```
python setup.py install

```

Example script

```
ocms combine_lanes -h
```

To use scripts that are in the R_utilty directory, you should be able to install the R package OCMSutility.

Install in R

```
devtools::install_github("https://github.com/OxfordCMS/OCMS_Sandbox/tree/master/R_utility/OCMSutility")
```

Package documentation

```
library(OCMSutility)
?OCMSutlity
```

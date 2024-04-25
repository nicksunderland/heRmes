
<!-- README.md is generated from README.Rmd. Please edit that file -->

# heRmes <img src="man/figures/hex.png" align="right" width="190"/>

<!-- badges: start -->
<!-- badges: end -->

The goal of **heRmes** is to standardise the heart failure phenotyping
of collections of electronic health records.

## Installation

You can install the latest version of **heRmes** like so:

``` r
# install.packages("devtools")
devtools::install_github("nicksunderland/heRmes")
```

## Phenotypes

The code lists underpinning the various phenotypes are stored in text
files within the package structure at: `inst/extdata/ukhdr_phenotypes`.
The format of the file matches that used by the [UKHDR Phenotype
Library](https://phenotypes.healthdatagateway.org), but the important
columns are: `code`, `description`, `coding_system.name`, `phenotype_id`
and `phenotype_name`. Below is an example of how to view the available
phenotypes and obtain the codes.

### Avaiable phenotypes

For example, view the first 6 phenotypes.

``` r
get_phenotypes()[1:5]
#>                         CCU002_02 Cardiomyopathy 
#>                                         "PH1002" 
#>                Acute Myocardial Infarction (AMI) 
#>                                         "PH1024" 
#>                  Heart Failure (fatal/non-fatal) 
#>                                         "PH1028" 
#> Congestive heart failure - Charlson primary care 
#>                                         "PH1055" 
#>    Myocardial infarction - Charlson primary care 
#>                                         "PH1062"
```

### Codes

View the codes for `PH1645` corresponding to the HERMES Heart Failure
phenotype.

``` r
# top 5 codes
get_codes(pheno_id = "PH1645")[1:5, c("phenotype_id", "phenotype_name", "coding_system.name", "code")]
#>    phenotype_id phenotype_name coding_system.name   code
#>          <char>         <char>             <char> <char>
#> 1:       PH1645  Heart failure         ICD9 codes  40201
#> 2:       PH1645  Heart failure         ICD9 codes  42830
#> 3:       PH1645  Heart failure         ICD9 codes  42832
#> 4:       PH1645  Heart failure         ICD9 codes   4281
#> 5:       PH1645  Heart failure         ICD9 codes  42843
```

### Phenotyping a dataset

Create sample data in long format (multiple ID entries per asscoaited
EHR code).

``` r
set.seed(2020)
n   <- 10
dat <- data.frame(ids   = paste0("ID_", c(1:(n/2), 1:(n/2))), 
                  codes = sample(c("I420", "foo", "bar", "baz"), n, replace = TRUE))
dat
#>     ids codes
#> 1  ID_1   baz
#> 2  ID_2   baz
#> 3  ID_3   bar
#> 4  ID_4   foo
#> 5  ID_5   baz
#> 6  ID_1  I420
#> 7  ID_2  I420
#> 8  ID_3   baz
#> 9  ID_4   foo
#> 10 ID_5   foo
```

Phenotype the individuals with phenotype `PH1645` (heart failure),
excluding phenotype `PH1637` (congenital heart disease). There can be
multiple included or excluded phenotypes given in a list.

``` r
result <- phenotype(dat$ids, dat$codes, 
                    name    = "Heart Failure", 
                    include = list(HF = "PH1645"), 
                    exclude = list(congHD = "PH1637"))
result[]
#>        id include exclude Heart Failure
#>    <char>  <lgcl>  <lgcl>        <lgcl>
#> 1:   ID_1    TRUE   FALSE          TRUE
#> 2:   ID_2    TRUE   FALSE          TRUE
#> 3:   ID_3   FALSE   FALSE         FALSE
#> 4:   ID_4   FALSE   FALSE         FALSE
#> 5:   ID_5   FALSE   FALSE         FALSE
```

### Update library from UKHDR

The package phenotype library can be updated from the [UKHDR Phenotype
Library API](https://phenotypes.healthdatagateway.org/api/v1/) using the
below function. This queries the library for phenotypes matching
enteries in the `search_terms` argument.

``` r
update_library(search_terms = c("heart failure", "cardiomyopathy", "myocardial infarction"))
#> [i] reading phenotype id: PH25 - skipping, already exists
#> [i] reading phenotype id: PH182 - skipping, already exists
#> [i] reading phenotype id: PH530 - skipping, already exists
#> [i] reading phenotype id: PH531 - skipping, already exists
#> [i] reading phenotype id: PH631 - skipping, already exists
#> [i] reading phenotype id: PH687 - skipping, already exists
#> [i] reading phenotype id: PH968 - skipping, already exists
#> [i] reading phenotype id: PH993 - skipping, already exists
#> [i] reading phenotype id: PH1028 - skipping, already exists
#> [i] reading phenotype id: PH1055 - skipping, already exists
#> [i] reading phenotype id: PH1074 - skipping, already exists
#> [i] reading phenotype id: PH1603 - skipping, already exists
#> [i] reading phenotype id: PH129 - skipping, already exists
#> [i] reading phenotype id: PH145 - skipping, already exists
#> [i] reading phenotype id: PH185 - skipping, already exists
#> [i] reading phenotype id: PH961 - skipping, already exists
#> [i] reading phenotype id: PH1002 - skipping, already exists
#> [i] reading phenotype id: PH215 - skipping, already exists
#> [i] reading phenotype id: PH356 - skipping, already exists
#> [i] reading phenotype id: PH481 - skipping, already exists
#> [i] reading phenotype id: PH530 - skipping, already exists
#> [i] reading phenotype id: PH611 - skipping, already exists
#> [i] reading phenotype id: PH612 - skipping, already exists
#> [i] reading phenotype id: PH613 - skipping, already exists
#> [i] reading phenotype id: PH741 - skipping, already exists
#> [i] reading phenotype id: PH886 - skipping, already exists
#> [i] reading phenotype id: PH942 - skipping, already exists
#> [i] reading phenotype id: PH949 - skipping, already exists
#> [i] reading phenotype id: PH988 - skipping, already exists
#> [i] reading phenotype id: PH1024 - skipping, already exists
#> [i] reading phenotype id: PH1062 - skipping, already exists
```

### Plotting phenotype

To see the intersection of the codes in two or more phenotype files use
the `plot_code_overlap()` function.

``` r
plot_code_overlap(pheno_ids = c("PH1645", "PH1028", "PH1055", "PH1074", "PH182", "PH25", "PH530", "PH531", "PH631", "PH687", "PH968", "PH993"), 
                  types = c("ICD10 codes", "ICD9 codes", "OPCS4 codes", "Read codes v2"))
```

<img src="man/figures/README-plot-1.png" width="80%" style="display: block; margin: auto;" />

### Update library from UKHDR (unpublished)

The package phenotype library can be updated with
unpublished/development phenotypes from the [UKHDR Phenotype Library
API](https://phenotypes.healthdatagateway.org/api/v1/) using the below
function. However, since unpublished phenotypes are not searchable by
name, we need to pass the exact ID and also login details for the
website (stored in a local `.Renviron` file in this example.)

``` r
# development phenotypes, ids named for readability only
hermes_phenos <- c(`Congenital heart disease`    = "PH1637", 
                   `Myocardial infarction`       = "PH1636", 
                   `Secondary cardiomyopathies`  = "PH1642", 
                   `Hypertrophic cardiomyopathy` = "PH1640", 
                   `Dilated cardiomyopathy`      = "PH1638", 
                   `Cardiomyopathy`              = "PH1646", 
                   `Heart failure`               = "PH1645", 
                   `Non-ischaemic cardiomyopathy`= "PH1639", 
                   `Heart failure syndrome`      = "PH1643")

# update
update_library(search_terms = c(), 
               ids          = hermes_phenos, 
               UKHDR_UN     = Sys.getenv("UKHDR_UN"), 
               UKHDR_PW     = Sys.getenv("UKHDR_PW"))
#> [i] reading phenotype id: PH1637 - skipping, already exists
#> [i] reading phenotype id: PH1636 - skipping, already exists
#> [i] reading phenotype id: PH1642 - skipping, already exists
#> [i] reading phenotype id: PH1640 - skipping, already exists
#> [i] reading phenotype id: PH1638 - skipping, already exists
#> [i] reading phenotype id: PH1646 - skipping, already exists
#> [i] reading phenotype id: PH1645 - skipping, already exists
#> [i] reading phenotype id: PH1639 - skipping, already exists
#> [i] reading phenotype id: PH1643 - skipping, already exists
```

### Plot the ICD-10 HERMES phenotypes

``` r
plot_code_overlap(pheno_ids = hermes_phenos, types = c("ICD10 codes"))
```

<img src="man/figures/README-plot_hermes-1.png" width="100%" />

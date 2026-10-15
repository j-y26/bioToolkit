# Extracting sample info based on regex keywords

Extracting sample info based on regex keywords found in a string. If
multiple keywords are found, the first keyword found will be used to
extract the info.

## Usage

``` r
keywords_to_info(sample, keywords, info)
```

## Arguments

- sample:

  A string containing the sample info, must be a single string

- keywords:

  A vector of regex keywords to search for

- info:

  A vector of info represented by the keywords, must be the same length
  as keywords

## Value

the info extracted from the sample based on the keywords found

## Examples

``` r
# Extracting sample info based on regex keywords
sample <- "Sample_2wks_1"
keywords <- c("[0-9]wks", "p7")
info <- c("postnatal_2wks", "postnatal_day_7")

# Extracting the stage info
keywords_to_info(sample, keywords, info)
#>    Sample_2wks_1 
#> "postnatal_2wks" 
```

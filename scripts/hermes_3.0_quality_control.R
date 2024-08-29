library(heRmes)

# -------------------
# Step 0 - check column names
# -------------------
example <- data.table::fread("/Users/xx20081/Desktop/hermes_pheno1_eur_harmonised.gz")
example <- example[, list(
  chromosome = chr,
  position   = bp,
  reference  = oa,
  alt        = ea,
  beta       = beta,
  stdErr     = se,
  n          = n,
  eaf        = eaf,
  pValue     = p,
)]
data.table::fwrite(example, "/Users/xx20081/Desktop/hermes_pheno1_eur_harmonised.gz"

# -------------------
# Step 1 - read file
# -------------------
dat <- data.table::fread(file)


# -------------------
# Step 2 - do something
# -------------------
dat <- heRmes::some_function(dat)


# -------------------
# Step 3 - do something else
# -------------------

library(heRmes)

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

# Parse maintained R sources without executing data access or model fitting.

files <- list.files(
    'src', pattern = '[.]R$', recursive = TRUE, full.names = TRUE
)

failures <- character()
for (file in files) {
    error <- tryCatch({
        parse(file = file)
        NULL
    }, error = function(condition) conditionMessage(condition))
    if (!is.null(error)) failures <- c(failures, paste(file, error, sep = ': '))
}

if (length(failures)) {
    cat(paste(failures, collapse = '\n'), '\n')
    quit(status = 1L)
}

cat('Parsed', length(files), 'R source files successfully.\n')

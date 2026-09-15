#!/usr/bin/env Rscript

registry_path <- 'config/experiments.yml'
or_default <- function(x, y) if (is.null(x)) y else x
word_count <- function(x) {
  words <- strsplit(trimws(x), '[[:space:]]+')[[1]]
  if (length(words) == 1L && !nzchar(words)) 0L else length(words)
}

if (!file.exists(registry_path)) {
  stop('Run this validator from the repository root.', call. = FALSE)
}
if (!requireNamespace('yaml', quietly = TRUE)) {
  stop('Package yaml is required.', call. = FALSE)
}

registry <- yaml::read_yaml(registry_path)
experiments <- registry$experiments
word_limit <- as.integer(or_default(registry$finding_word_limit, 45L))
allowed_status <- c(
  'selected', 'retained', 'not_promoted',
  'superseded', 'diagnostic', 'data_pipeline'
)
required <- c(
  'id', 'title', 'track', 'status', 'question', 'finding',
  'report', 'evidence', 'updated'
)
errors <- character()

if (!is.list(experiments) || length(experiments) == 0L) {
  stop('The experiment registry is empty.', call. = FALSE)
}

for (i in seq_along(experiments)) {
  x <- experiments[[i]]
  label <- or_default(x$id, paste0('record_', i))
  missing <- required[vapply(required, function(field) is.null(x[[field]]), logical(1))]
  if (length(missing)) {
    errors <- c(errors, sprintf('%s: missing %s', label, paste(missing, collapse = ', ')))
    next
  }
  if (!x$status %in% allowed_status) {
    errors <- c(errors, sprintf('%s: invalid status %s', label, x$status))
  }
  if (!x$track %in% c('current', 'investigation')) {
    errors <- c(errors, sprintf('%s: invalid track %s', label, x$track))
  }
  if (word_count(x$finding) > word_limit) {
    errors <- c(errors, sprintf('%s: finding exceeds %d words', label, word_limit))
  }
  if (word_count(x$question) > 25L) {
    errors <- c(errors, sprintf('%s: question exceeds 25 words', label))
  }
  if (!file.exists(x$report)) {
    errors <- c(errors, sprintf('%s: report does not exist: %s', label, x$report))
  }
  if (length(unlist(x$evidence, use.names = FALSE)) == 0L) {
    errors <- c(errors, sprintf('%s: evidence list is empty', label))
  }
}

ids <- vapply(experiments, function(x) or_default(x$id, ''), character(1))
reports <- vapply(experiments, function(x) or_default(x$report, ''), character(1))
if (anyDuplicated(ids)) errors <- c(errors, 'Experiment IDs are not unique.')
if (anyDuplicated(reports)) errors <- c(errors, 'Experiment report paths are not unique.')

analysis_reports <- c(
  list.files('analysis/current', pattern = '[.]qmd$', recursive = TRUE, full.names = TRUE),
  list.files('analysis/investigations', pattern = '[.]qmd$', recursive = TRUE, full.names = TRUE)
)
analysis_reports <- chartr('\\', '/', analysis_reports)
registered_reports <- chartr('\\', '/', reports)
unregistered <- setdiff(analysis_reports, registered_reports)
orphaned <- setdiff(registered_reports, analysis_reports)
if (length(unregistered)) {
  errors <- c(errors, paste('Unregistered reports:', paste(unregistered, collapse = ', ')))
}
if (length(orphaned)) {
  errors <- c(errors, paste('Registry reports outside analysis folders:', paste(orphaned, collapse = ', ')))
}

if (length(errors)) {
  stop(paste(c('Experiment registry validation failed:', paste0('- ', errors)), collapse = '\n'),
       call. = FALSE)
}

evidence <- unique(unlist(lapply(experiments, function(x) x$evidence), use.names = FALSE))
missing_evidence <- sum(!file.exists(evidence))
cat(sprintf(
  'Experiment registry valid: %d records, %d current, %d generated evidence paths absent locally.\n',
  length(experiments),
  sum(vapply(experiments, function(x) x$track == 'current', logical(1))),
  missing_evidence
))

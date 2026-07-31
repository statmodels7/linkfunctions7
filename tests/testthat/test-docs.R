# Documentation invariants, and one of them is here because CI caught it first.

test_that("every exported object appears in the pkgdown index", {
  # pkgdown refuses to build a site with a topic missing from the reference
  # index, and it refuses IN CI, minutes after the push, when the fix is one
  # line in a YAML file. Exporting link_bounds_clamp() without adding it there
  # turned a green push red for that reason alone.
  #
  # check_pkgdown() asks exactly this question and runs in a second. Asking it
  # here means the answer arrives while the export is still being written.
  skip_if_not_installed("pkgdown")
  root <- normalizePath("../..", mustWork = FALSE)
  skip_if_not(file.exists(file.path(root, "_pkgdown.yml")),
              "_pkgdown.yml not reachable from here")

  expect_no_error(pkgdown::check_pkgdown(pkg = root))
})


test_that("every exported topic has a value and an executable example", {
  # The two commonest reasons a first CRAN submission comes back, and R CMD
  # check raises neither locally.
  man <- normalizePath("../../man", mustWork = FALSE)
  skip_if_not(dir.exists(man), "man/ not reachable from here")

  exported <- getNamespaceExports(asNamespace("linkfunctions7"))
  info <- lapply(list.files(man, pattern = "[.]Rd$", full.names = TRUE),
                 function(f) {
    rd <- tools::parse_Rd(f)
    tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
    aliases <- unlist(lapply(rd[tags == "\\alias"], function(x)
      trimws(paste(unlist(x), collapse = ""))))
    list(file = basename(f), aliases = aliases,
         value = "\\value" %in% tags, example = "\\examples" %in% tags)
  })

  expect_equal(setdiff(exported, unlist(lapply(info, `[[`, "aliases"))),
               character())

  is_exp <- vapply(info, function(i) any(i$aliases %in% exported), logical(1))
  no_value <- vapply(info, function(i) !i$value, logical(1))
  no_example <- is_exp & vapply(info, function(i) !i$example, logical(1))

  expect_equal(vapply(info[no_value], `[[`, character(1), "file"), character())
  expect_equal(vapply(info[no_example], `[[`, character(1), "file"),
               character())
})

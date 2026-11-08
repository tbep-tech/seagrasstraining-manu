docx_files <- Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES")
docx_files <- strsplit(docx_files, "\n")[[1]]
docx_files <- docx_files[grepl("\\.docx$", docx_files)]

for (f in docx_files) {
  flextable::repair_docx(f, target = f)
}

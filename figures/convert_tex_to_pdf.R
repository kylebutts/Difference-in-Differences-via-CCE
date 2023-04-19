library(here)
library(glue)
library(xfun)

delete_if_exists <- function(filename) {
  if (file.exists(filename)) file.remove(filename)
}

#' @param path Character. Either folder or filename
#' @param redo Logical. If FALSE, doesn't update if pdf date
#'   is newer than the tex date.
#'
#' @return nothing. Used for sideeffect of creating pdf
compile_tikz <- function(path, redo = FALSE) {
  # https://tex.stackexchange.com/questions/138677/why-does-standalone-not-detect-the-tikz-crop-correctly
  # right margin is added cause it usually looks cropped too tight to the right
  WRAPPER_START <- r'(
\documentclass[margin={0mm 0mm 10mm 0mm}]{standalone}
\usepackage{amsfonts}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{a4wide}
\usepackage{lscape}
\usepackage{mathpazo}
\usepackage{rotating}
\usepackage{subcaption}
\usepackage{float}
\usepackage{epstopdf}
\usepackage{rotating}
\usepackage{adjustbox}
\usepackage{setspace}
\usepackage[flushleft]{threeparttable}
\setlength\labelsep{0pt}
\usepackage{booktabs}
\usepackage{adjustbox}
\usepackage{dsfont}
\usepackage{bm}
\usepackage{tikz}
%\numberwithin{equation}{section}
\def\*#1{\mathbf{#1}}
\def\+#1{\boldsymbol{#1}}

% Fix \input with tables
% \input fails when \\ is at end of external .tex file
\makeatletter
\let\input\@@input
\makeatother

\begin{document}

)'

  WRAPPER_END <- r'(

\end{document}
  )'


  is_file <- grepl("\\.tex", path)

  if (is_file == TRUE) {
    tikz_files <- path
    dir <- dirname(path)
  } else {
    tikz_files <- list.files(
      path = here(path), pattern = "\\.tex",
      full.names = TRUE
    )
    dir <- path
  }

  delete_if_exists(here(dir, "temp.tex"))

  for (f in tikz_files) {
    pdf_name <- gsub(".tex", ".pdf", f)

    # Don't recompile new files
    if (file.exists(pdf_name)) {
      tikz_compile_time <- file.info(f)$mtime
      pdf_compile_time <- file.info(pdf_name)$mtime

      if (
        (pdf_compile_time > tikz_compile_time) & redo == FALSE
      ) {
        break
      }
    }

    standalone_str <- paste(
      WRAPPER_START,
      paste(xfun::read_utf8(f), collapse = "\n"),
      WRAPPER_END,
      collapse = "\n"
    )

    xfun::write_utf8(standalone_str, here(dir, "temp.tex"))
    # cat(standalone_str, file = here(dir, "temp.tex"))

    system(glue(r'(
      cd "{here('figures')}" &&
      latexmk temp.tex
    )'))

    file.rename(
      here(dir, "temp.pdf"),
      here(dir, basename(pdf_name))
    )

    cat(paste("\n\nCreated pdf:", pdf_name, "\n\n"))

    delete_if_exists(here(dir, "temp.aux"))
    delete_if_exists(here(dir, "temp.fdb_latexmk"))
    delete_if_exists(here(dir, "temp.flx"))
    delete_if_exists(here(dir, "temp.log"))
    delete_if_exists(here(dir, "temp.xdv"))
    delete_if_exists(here(dir, "temp.fls"))
    delete_if_exists(here(dir, "temp.tex"))
  }
}

# compile_tikz("figures", redo = TRUE)

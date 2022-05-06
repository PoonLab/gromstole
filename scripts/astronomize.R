#' re.findall
#'
#' Emulate behaviour of Python's re.findall() function. Lovingly stolen from \url{https://github.com/PoonLab/gromstole}.
#'
#' @param pat:  regex pattern
#' @param s:  character, a single string
#' @return character, vector of all matching substrings
re.findall <- function(pat, s) {
  if (!is.character(s)) {
    stop("re.findall() requires a character object for input 's'")
  }
  matches <- gregexpr(pat, s)
  index <- as.integer(matches[[1]])
  match.length <- attr(matches[[1]], 'match.length')
  
  sapply(1:length(index), function(i) {
    start <- index[i]
    stop <- start + match.length[i] - 1
    substr(s, start, stop)
  })
}

#' Extract mutation list from a directory of constellation files.
#' 
#' "Constellations" are files produced from \url{https://github.com/cov-lineages/constellations}, which represent the cov-lineage team's best knowledge about which mutations define a lineage. These are updated frequently, please \code{git pull} in your cloned copy accordingly. See details.
#' 
#' Code is adapted from \code{scripts/estimate-freqs.R} in \url{https://github.com/PoonLab/gromstole}.
#' 
#' @param path Path to the constellations folder in the cov-lineages/constellations repository. The default assumes that the current project and the constellations repo are in the same directory.
#' 
#' @return A variant matrix for use with provoc.
#' @export
#' 
#' @details From the repo, a constellation "a collection of mutations which are functionally meaningful, but which may arise independently a number of times".
#' 
#' This function requires a local clone of the constellations repository. By default, I assume that your current project is in the same parent directory as the constellation file (e.g. your working directory is \code{.../parent/current_project} and the constellations repo is cloned as \code{.../parent/constellations}).
#' 
#' If this is not the case, a path to the root folder of the constellations repo is required. The path is to the root folder, not the folder with the constellation files.
#' 
#' There are no options to specify which lineages to include. I am working on a \code{annihilate(coco, varmat)} function to remove mutations that aren't shared between coco and varmat and squash lineages that have too few observed mutations or are too similar to other lineages.
#' 
#' @examples
#' if(dir.exists("../constellations")) {
#'     varmat <- astronomize()
#' }
#' 
astronomize <- function(path = "../constellations") {
    orfs <- list(
        'orf1a'= c(265, 13468),
        'orf1b'= c(13467, 21555),
        'S'= c(21562, 25384),
        'orf3a'= c(25392, 26220),
        'E'= c(26244, 26472),
        'M'= c(26522, 27191),
        'orf6'= c(27201, 27387),
        'orf7a'= c(27393, 27759),
        'orf7b'= c(27755, 27887),
        'orf8'= c(27893, 28259),
        'N'= c(28273, 29533),
        'orf10'= c(29557, 29674),
        'nsp2' = c(806, 2719),
        'nsp3' = c(2720, 8554),
        'nsp4' = c(8555, 10054),
        'nsp5' = c(10055, 10972),
        'nsp6' = c(10973, 11842),
        'nsp7' = c(11263, 11509),
        'nsp12' = c(13442, 16236),
        'nsp13' = c(16237, 18039),
        'nsp15' = c(19621, 20658)
    )

    stelpath <- paste0(path, "/constellations/definitions")

    sitelist <- lapply(list.files(stelpath, full.names = TRUE), function(stelfile){
        stelfile <- stelfile
        lineage <- gsub("^c|\\.json$", "", basename(stelfile))
        constellation <- jsonlite::read_json(stelfile, simplifyVector = TRUE)
        constellation$sites <- unique(constellation$sites)

        len_1a <- (orfs[['orf1a']][2]-orfs[['orf1a']][1])/3 + 1

        # convert constellation to label notation in the mapped files
        sites <- lapply(unique(constellation$sites), function(d) {
            d <- d
            toks <- toupper(strsplit(d, ":")[[1]])

            if (toks[1] != "DEL" && toks[1] != "NUC")
            toks <- c("aa", toks)

            if (toks[2] == "S" || toks[2] == "SPIKE") {
                toks[[2]] <- "S"
            } else if (toks[1] == "DEL") {
                toks[[1]] <- "del"
            } else if (toks[1] == "NUC") {
                toks <- toks[-1]
                toks[[1]] <- substring(toks[1], 2, nchar(toks[1]))
            } else if (toks[2] == "8") {
                toks[[2]] <- "orf8"
            } else if (toks[2] == "ORF1AB" || toks[2] == "1AB") {
                num <- as.numeric(re.findall("\\d+", toks[3]))

                if (num <= len_1a) {
                    toks[[2]] <- "orf1a"
                } else {
                # Determine nucleotide position relative to the start of orf1b
                new_pos <- (((num-1) * 3 + orfs[['orf1a']][1]) - orfs[['orf1b']][1])/3
                toks[[3]] <- gsub(num, floor(new_pos) + 1, toks[[3]])
                toks[[2]] <- "orf1b"
                }

            } else if (nchar(toks[2]) >= 3 && substring(toks[2], 1, 3) == "ORF") {
                toks[[2]] <- tolower(toks[2])
            } else if (substring(toks[2], 1, 3) == "NSP") {
                start_pos <- orfs[[tolower(toks[2])]][1]
                codon <- as.integer(re.findall("\\d+", toks[3]))
                nuc_pos <- start_pos + (codon-1)*3
                if (nuc_pos >= orfs[['orf1a']][1] && nuc_pos <= orfs[['orf1a']][2]) {
                    toks[[2]] <- 'orf1a'
                } else if (nuc_pos >= orfs[['orf1b']][1] && nuc_pos <= orfs[['orf1b']][2]) {
                    toks[[2]] <- 'orf1b'
                } else {
                    stop("Could not convert nsp to orf1a/b")
                }
                new_pos <- ((nuc_pos - orfs[[toks[2]]][1])/3) + 1
                toks[[3]] <- gsub(codon, floor(new_pos), toks[[3]])
            }

            if (grepl("+", toks[1], fixed = TRUE)) {
            ins <- strsplit(toks[[1]], split = "[+]")[[1]]
            toks[[1]] <- gsub(" ", "", paste("+", ins[1], ".", ins[2]))
            }
            toks <- paste(toks, collapse = ":")
        })

        sites <- unlist(sites, recursive = FALSE)
    })

    varmat <- as.matrix(dplyr::bind_rows(sapply(sitelist, function(sites) {
        x <- rep(1, length(sites))
        names(x) <- sites
        x
    })))

    varmat[is.na(varmat)] <- 0
    rownames(varmat) <- gsub(".json", "", list.files(stelpath, full.names = FALSE))
    rownames(varmat) <- gsub("c", "", rownames(varmat))    

    varmat <- varmat[, colSums(varmat) > 0]
    varmat
}













#' Read unisens binary files
#'
#' Function to read unisens binary files into a matrix. Reading happens in
#' blocks of 1 000 000 values at a time, to keep track of progress and to keep
#' the process manageable in size. This function will display a nice
#' progress bar in interactive use.
#'
#' @param folder the folder containing the xml file and the .bin file(s)
#' @param file_id the filename of the binary file. This should also exist in the
#'   xml file as an ID. Upon error, a nice warning will be given.
#'
#' @return A unisens data matrix. The named columns of the matrix are the
#'   different channels, and the rows of the matrix are the samples. The
#'   attributes of the matrix will include further information about the file,
#'   such as the sample rate.
#'
#' @details The reading will happen in blocks of 1 000 000 values at a time, to
#'   keep track of progress and to keep the process manageable in size.
#'
#' @seealso \code{\link{sample_rate}}, \code{\link{sub_sample}}
#'
#' @export
read_unisens <- function(folder, file_id) {
  folder_files <- list.files(folder)
  xml_file     <- folder_files[grepl("\\.xml", folder_files)][1]
  xml <- xml2::xml_ns_strip(xml2::read_xml(file.path(folder, xml_file)))

  xpath <- paste0("//*[@id=\"", file_id, "\"]")
  node  <- xml2::xml_find_first(xml, xpath)

  if (is.na(node)) {
    stop("The file ID could not be found. Available IDs:\n - ",
         paste(xml2::xml_attr(xml2::xml_find_all(xml, "signalEntry"), "id"),
               collapse = "\n - "))
  }

  # parse attributes
  attrs <- as.list(xml2::xml_attrs(node))
  attrs$endianess <- xml2::xml_attr(
    xml2::xml_find_first(node, "binFileFormat"),
    "endianess"
  )
  attrs$channels <-  xml2::xml_attr(xml2::xml_find_all(node, "channel"), "name")

  # actually read in block-wise manner with R
  dat <- .read_bin_unisens(file.path(folder, file_id), attrs = attrs)
  attributes(dat) <- append(attributes(dat), attrs)
  return(dat)
}

#' Read unisens sample properties
#'
#' Parses the xml file into a list to extract the attributes (such as starting
#' timestamp) of the sample.
#'
#' @param folder path to the unisens data folder
#'
#' @return list of attributes associated with the sample
#'
#' @export
read_props <- function(folder) {
  folder_files <- list.files(folder)
  xml_file     <- folder_files[grepl("\\.xml", folder_files)][1]
  xml <- xml2::xml_ns_strip(xml2::read_xml(file.path(folder, xml_file)))
  attr_nodes <- xml2::xml_find_all(xml, ".//customAttribute")
  attrs <- xml2::xml_attr(attr_nodes, "value")
  names(attrs) <- xml2::xml_attr(attr_nodes, "key")

  return(append(as.list(attrs), as.list(xml2::xml_attrs(xml))))
}


#' Binary reading by block
#'
#' @keywords internal
.read_bin_unisens <- function(path, attrs) {
  if (attrs$dataType == "int16") {
    rbSize <- 2
  } else if (attrs$dataType == "int32") {
    rbSize <- 4
  } else {
    stop("Unsupported data type: ", attrs$dataType)
  }

  info <- file.info(path)

  if (is.na(info$size)) {
    stop("File ", path, " not found. Make sure the file is correctly named!",
         call. = FALSE)
  }

  # open a connection to the file
  conn <- file(description = path, open = "rb", raw = TRUE)
  on.exit(close(conn))

  # prepare a vector
  len <- info$size / rbSize
  blk <- 1e6L
  vec <- integer(len)
  if (interactive())
    pb <- progress::progress_bar$new(
      total = ceiling(len/blk),
      format = "[:spin] [:bar] :percent :eta"
    )

  # perform block-wise reading of the data
  for (i in 1:ceiling(len/blk)) {
    if (i < ceiling(len/blk)) {
      idx <- ((i - 1) * blk + 1):(i*blk)
      vec[idx] <- readBin(
        con    = conn,
        what   = "integer",
        n      = blk,
        size   = rbSize,
        signed = TRUE,
        endian = tolower(attrs$endianess)
      )
    } else {
      # final block
      idx <- ((i - 1) * blk + 1):len
      vec[idx] <- readBin(
        con    = conn,
        what   = "integer",
        n      = length(idx),
        size   = rbSize,
        signed = TRUE,
        endian = tolower(attrs$endianess)
      )
    }
    if (interactive()) pb$tick()
  }


  # post-processing
  lsb      <- as.numeric(attrs$lsbValue)
  if (length(lsb) == 0) lsb <- 1
  baseline <- as.numeric(attrs$baseline)
  if (length(baseline) == 0) baseline <- 0

  if (baseline != 0 || lsb != 1) {
    if (interactive()) message("Adjusting baseline and range...")
    vec <- (vec - baseline) * lsb
  }

  dat <- matrix(vec, ncol = length(attrs$channels), byrow = TRUE)
  colnames(dat) <- attrs$channels
  dat
}

#' Check the sample rate of unisens data.
#'
#' @param uni_mat a data matrix as returned by the read_unisens function
#'
#' @export
sample_rate <- function(uni_mat) {
  rate <- as.numeric(attr(uni_mat, "sampleRate"))
  if (length(rate) == 0) {
    stop("Please enter a valid unisens matrix from the read_unisens function.")
  }
  rate
}

#' Subsample unisens data to a target sample rate
#'
#' This function subsamples a unisens data matrix to match a target sample rate.
#' The target sample rate in Hz is approximately matched by selecting the rows
#' closest to the target sample rows.
#'
#' @param uni_mat a data matrix as returned by the read_unisens function
#' @param tgt_rate a numeric value containing the target sample rate in Hz
#'
#' @return a unisens data matrix
#'
#' @export
sub_sample <- function(uni_mat, tgt_rate) {
  rate <- sample_rate(uni_mat)
  if (rate < tgt_rate) stop("Target rate is higher than data sample rate!")
  dat <- uni_mat[round(seq(1L, nrow(uni_mat), rate/tgt_rate)),]
  attr(dat, "sampleRate") <- tgt_rate
  dat
}



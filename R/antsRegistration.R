#' Perform registration between two images.
#'
#' Register a pair of images either through the full or simplified interface
#' to the ANTs registration method.
#'
#' @param fixed fixed image to which we register the moving image.
#' @param moving moving image to be mapped to fixed space.
#' @param typeofTransform A linear or non-linear registration type.  Mutual
#' information metric by default. See \code{Details.}
#' @param initialTransform transforms to prepend
#' @param outprefix output will be named with this prefix.
#' @param mask mask the registration.
#' @param gradStep gradient step size (not for all tx)
#' @param verbose request verbose output (useful for debugging)
#' @param ... additional options see antsRegistration in ANTs
#' @details
#' typeofTransform can be one of:
#' \itemize{
#'   \item{"Rigid": }{Rigid transformation: Only rotation and translation.}
#'   \item{"QuickRigid": }{Rigid transformation: Only rotation and translation.
#'   May be useful for quick visualization fixes.'}
#'   \item{"BOLDRigid": }{Rigid transformation: Parameters typical for BOLD
#'   to BOLD intrasubject registration'.'}
#'   \item{"Affine": }{Affine transformation: Rigid + scaling.}
#'   \item{"AffineFast": }{Fast version of \code{Affine}.}
#'   \item{"BOLDAffine": }{Affine transformation: Parameters typical for BOLD
#'   to BOLD intrasubject registration'.'}
#'   \item{"SyN": }{Symmetric normalization: Affine + deformable transformation,
#'     with mutual information as optimization metric.}
#'   \item{"SyNRA": }{Symmetric normalization: Rigid + Affine + deformable transformation,
#'     with mutual information as optimization metric.}
#'   \item{"SyNOnly": }{Symmetric normalization: no initial transformation,
#'     with mutual information as optimization metric.  Assumes images are
#'     aligned by an inital transformation. Can be useful if you want to run
#'     an unmasked affine followed by masked deformable registration.}
#'   \item{"SyNCC": }{SyN, but with cross-correlation as the metric.}
#'   \item{"SyNBold": }{SyN, but optimized for registrations between
#'     BOLD and T1 images.}
#'   \item{"SyNBoldAff": }{SyN, but optimized for registrations between
#'     BOLD and T1 images, with additional affine step.}
#'   \item{"SyNAggro": }{SyN, but with more aggressive registration
#'     (fine-scale matching and more deformation).  Takes more time than \code{SyN}.}
#'   \item{"TVMSQ": }{time-varying diffeomorphism with mean square metric}
#' }
#' @return outputs a list containing:
#' \itemize{
#'   \item{warpedmovout: }{Moving image warped to space of fixed image.}
#'   \item{warpedfixout: }{Fixed image warped to space of moving image.}
#'   \item{fwdtransforms: }{Transforms to move from moving to fixed image.}
#'   \item{invtransforms: }{Transforms to move from fixed to moving image.}
#' }
#' Ouptut of 1 indicates failure
#' @author Shrinidhi KL, Tustison NJ, Avants BB
#' @examples
#'
#' fi <- antsImageRead(getANTsRData("r16") ,2)
#' mi <- antsImageRead(getANTsRData("r64") ,2)
#' fi<-resampleImage(fi,c(60,60),1,0)
#' mi<-resampleImage(mi,c(60,60),1,0) # speed up
#' mytx <- antsRegistration(fixed=fi, moving=mi, typeofTransform = c('SyN') )
#' mywarpedimage <- antsApplyTransforms( fixed=fi, moving=mi,
#'   transformlist=mytx$fwdtransforms )
#'
#' \dontrun{ # quick visualization fix for images with odd orientation
#' mni = antsImageRead( getANTsRData( "mni" ) )
#' strokt1=antsImageRead('strokt1.nii.gz')
#' strokt1reg=antsRegistration(
#'   fixed=mni,
#'   moving=strokt1,
#'   typeofTransform = "QuickRigid",verbose=TRUE )
#'  plot(  strokt1reg$warpedmovout, axis=3, nslices=20)
#' }
#'
#' @export antsRegistration
antsRegistration <- function( fixed = NA, moving = NA,
  typeofTransform = "SyN", initialTransform = NA,
  outprefix = "", mask = NA, gradStep=NA, verbose=FALSE, ... ) {
  numargs <- nargs()
  if (numargs == 1 & typeof(fixed) == "list") {
    .Call("antsRegistration", .int_antsProcessArguments(c(fixed)), PACKAGE = "ANTsR")
    return(0)
  }
  if (nchar(typeofTransform) == 0)
    typeofTransform = "SyN"
  if (nchar(outprefix) == 0)
    outprefix = tempfile()
  if ( numargs < 1 | missing(fixed) | missing(moving)
        | missing(typeofTransform) | missing(outprefix) )
    {
    cat("for simplified mode: \n")
    cat(" antsRegistration( fixed , moving , typeofTransform = c(\"Rigid\",\"Affine\",\"AffineFast\",\"SyN\",\"SyNCC\"),  outputPrefix=\"./antsRegOut\" \n")
    cat("")
    cat("For full mode: use standard ants call , e.g. : \n")
    cat(" ANTsR::antsRegistration( list( d=2,m=\"mi[r16slice.nii.gz,r64slice.nii.gz,1,20,Regular,0.05]\", t=\"affine[1.0]\", c=\"2100x1200x1200x0\",  s=\"3x2x1x0\", f=\"4x3x2x1\", u=\"1\", o=\"[xtest,xtest.nii.gz,xtest_inv.nii.gz]\" ) )\n")
    cat("full help: \n")
    .Call("antsRegistration", .int_antsProcessArguments(c(list("--help"))), PACKAGE = "ANTsR")
    return(0)
  }
  args <- list(fixed, moving, typeofTransform, outprefix, ...)
  myl=0
  myfAff="6x4x2x1"
  mysAff="3x2x1x0"
  myiterations <- "2100x1200x1200x10"
  if ( typeofTransform == "AffineFast" ) {
    typeofTransform <- "Affine"
    myiterations <- "2100x1200x0x0"
  }
  if ( typeofTransform == "BOLDAffine" ) {
    typeofTransform <- "Affine"
    myfAff="2x1"
    mysAff="1x0"
    myiterations <- "100x20"
    myl=1
  }
  if ( typeofTransform == "QuickRigid" ) {
    typeofTransform <- "Rigid"
    myiterations <- "20x20x0x0"
  }
  if ( typeofTransform == "BOLDRigid" ) {
    typeofTransform <- "Rigid"
    myfAff="2x1"
    mysAff="1x0"
    myiterations <- "100x20"
    myl=1
  }
  if (!is.character(fixed)) {
    if (fixed@class[[1]] == "antsImage" & moving@class[[1]] == "antsImage") {
      inpixeltype <- fixed@pixeltype
      ttexists <- FALSE
      allowableTx <- c("Rigid", "Affine", "SyN","SyNRA","SyNOnly","SyNCC",
        "SyNBold", "SyNBoldAff", "SyNAggro", "SyNLessAggro", "TVMSQ")
      ttexists <- typeofTransform %in% allowableTx
      if (ttexists) {
        initx = initialTransform
        moving <- antsImageClone(moving, "float")
        fixed <- antsImageClone(fixed, "float")
        warpedfixout <- antsImageClone(moving)
        warpedmovout <- antsImageClone(fixed)
        f <- antsrGetPointerName(fixed)
        m <- antsrGetPointerName(moving)
        wfo <- antsrGetPointerName(warpedfixout)
        wmo <- antsrGetPointerName(warpedmovout)
        if (!is.na(mask)) {
          maskScale = mask - min( mask )
          maskScale = maskScale / max( maskScale ) * 255
          charmask <- antsImageClone( maskScale , "unsigned char")
          maskopt <- antsrGetPointerName(charmask)
        } else maskopt=NA
        if (is.na(initx)) {
          initx = paste("[", f, ",", m, ",1]", sep = "")
        }
        if (typeofTransform == "SyNBold") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Rigid[0.25]", "-c", "[1200x1200x100,1e-6,5]", "-s", "2x1x0",
          "-f", "4x2x1", "-m", paste("cc[", f, ",", m, ",1,2]", sep = ""),
          "-t", paste("SyN[0.1,3,0]", sep = ""), "-c", "[200x10,1e-6,5]",
          "-s", "1x0", "-f", "2x1", "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na( maskopt )  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNBoldAff") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Rigid[0.25]", "-c", "[1200x1200x100,1e-6,5]", "-s", "2x1x0",
          "-f", "4x2x1",
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[0.25]", "-c", "[200x20,1e-6,5]", "-s", "1x0",
          "-f", "2x1",
          "-m", paste("cc[", f, ",", m, ",1,2]", sep = ""),
          "-t", paste("SyN[0.1,3,0]", sep = ""), "-c", "[200x10,1e-6,5]",
          "-s", "1x0", "-f", "2x1", "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na( maskopt )  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyN") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[0.25]", "-c", "2100x1200x1200x0", "-s", "3x2x1x0",
          "-f", "4x2x2x1",
          "-m", paste("mattes[", f, ",", m, ",1,32]", sep = ""),
          "-t", paste(typeofTransform, "[0.25,3,0]", sep = ""),
          "-c", "2100x1200x1200x0", "-s", "3x2x1x0", "-f", "4x3x2x1", "-u",
          "1", "-z", "1", "-l", myl, "-o", paste("[", outprefix, ",",
            wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNRA") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Rigid[0.25]", "-c", "2100x1200x1200x0", "-s", "3x2x1x0",
          "-f", "4x2x2x1",
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[0.25]", "-c", "2100x1200x1200x0", "-s", "3x2x1x0",
          "-f", "4x2x2x1",
          "-m", paste("mattes[", f, ",", m, ",1,32]", sep = ""),
          "-t", paste("SyN[0.25,3,0]", sep = ""),
          "-c", "2100x1200x1200x0", "-s", "3x2x1x0", "-f", "4x3x2x1", "-u",
          "1", "-z", "1", "-l", myl, "-o", paste("[", outprefix, ",",
            wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNOnly") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32]", sep = ""),
          "-t", paste("SyN[0.25,3,0]", sep = ""),
          "-c", "2100x1200x1200x0", "-s", "3x2x1x0", "-f", "4x3x2x1", "-u",
          "1", "-z", "1", "-l", myl, "-o", paste("[", outprefix, ",",
            wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNAggro") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[0.25]", "-c", "2100x1200x1200x100", "-s", "3x2x1x0",
          "-f", "4x2x2x1", "-m", paste("meansquares[", f, ",", m, ",1,2]",
            sep = ""), "-t", paste("SyN[0.1,3,0]", sep = ""), "-c", "2100x1200x1200x20",
          "-s", "3x2x1x0", "-f", "4x3x2x1", "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNCC") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Rigid[1]", "-c", "2100x1200x1200x0", "-s", "3x2x1x0",
          "-f", "4x4x2x1",
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[1]", "-c", "1200x1200x100", "-s", "2x1x0",
          "-f", "4x2x1",
          "-m", paste("CC[", f, ",", m, ",1,3]", sep = ""),
          "-t", paste("SyN[0.15,3,0]", sep = ""), "-c", "2100x1200x1200x20",
          "-s", "3x2x1x0", "-f", "4x3x2x1", "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "SyNLessAggro") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", "Affine[0.25]", "-c", "2100x1200x1200x100", "-s", "3x2x1x0",
          "-f", "4x2x2x1", "-m", paste("meansquares[", f, ",", m, ",1,2]",
            sep = ""), "-t", paste("SyN[0.1,3,0.5]", sep = ""), "-c", "2100x1200x1200x20",
          "-s", "3x2x1x0", "-f", "4x3x2x1", "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
        }
        if ( typeofTransform == "TVMSQ" ) {
          if ( is.na(gradStep) ) gradStep=1.0
          tvtx=paste("TimeVaryingVelocityField[",
            gradStep,", 4, 0.0,0.0, 0.5,0 ]",sep='')
          args <- list("-d", as.character(fixed@dimension), # "-r", initx,
            "-m", paste("meansquares[", f, ",", m, ",1,0]", sep = ""),
            "-t", tvtx,
            "-c", "[100,1.e-5,5]",
            "-s", "0",
            "-f", "1",
            "-u", "1", "-z", "1", "-l", myl,
            "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "1Warp.nii.gz", sep = ""),
          paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""),
          paste(outprefix, "1InverseWarp.nii.gz", sep = ""))
          fwdtransforms <- c( paste( outprefix, "0Warp.nii.gz", sep = "") )
          invtransforms <- c( paste(outprefix, "0InverseWarp.nii.gz", sep = ""))
        }
        if (typeofTransform == "Rigid" | typeofTransform == "Affine") {
          args <- list("-d", as.character(fixed@dimension), "-r", initx,
          "-m", paste("mattes[", f, ",", m, ",1,32,regular,0.2]", sep = ""),
          "-t", paste(typeofTransform, "[0.25]", sep = ""), "-c", myiterations,
          "-s", mysAff, "-f", myfAff, "-u", "1", "-z", "1", "-l", myl,
          "-o", paste("[", outprefix, ",", wmo, ",", wfo, "]", sep = ""))
          if ( !is.na(maskopt)  )
            args=lappend( list( "-x", maskopt ), args )
          fwdtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""))
          invtransforms <- c(paste(outprefix, "0GenericAffine.mat", sep = ""))
        }
        args[[ length(args)+1]]="--float"
        args[[ length(args)+1]]="1"
        if ( verbose ) {
          args[[ length(args)+1]]="-v"
          args[[ length(args)+1]]="1"
          }
        args = .int_antsProcessArguments(c(args))
        .Call("antsRegistration", args, PACKAGE = "ANTsR")
        # unlink(ffn) unlink(mfn) outvar<-basename(outprefix) outpath<-dirname(outprefix)
        # txlist<-list.files( path = outpath, pattern = glob2rx( paste(outvar,'*',sep='')
        # ), full.names = TRUE, recursive = FALSE )
        if ( sum(  fixed - warpedmovout ) == 0  ) # FIXME better error catching
          stop( "Registration failed. Use verbose mode to diagnose." )
        return(
          list( warpedmovout = antsImageClone(warpedmovout, inpixeltype),
                warpedfixout = antsImageClone(warpedfixout, inpixeltype),
                fwdtransforms = fwdtransforms,
                invtransforms = invtransforms )
              )
      }
      if (!ttexists) {
        stop("Unrecognized transform type.")
      }
    }
    return(0)
  }
  args[[ length(args)+1]]="--float"
  args[[ length(args)+1]]="1"
  if ( verbose ) {
    args[[ length(args)+1]]="-v"
    args[[ length(args)+1]]="1"
    }
  args = .int_antsProcessArguments(c(args))
  .Call("antsRegistration", args, PACKAGE = "ANTsR")
}

############################################################### .antsrmakeRandomString(n, length) function generates a random string random
############################################################### string of the length (length), made up of numbers, small and capital letters
############################################################### helper function
.antsrmakeRandomString <- function(n = 1, mylength = 12) {
  randomString <- c(1:n)  # initialize vector
  for (i in 1:n) {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS), mylength, replace = TRUE),
      collapse = "")
  }
  return(randomString)
}


#' Return the antsImage pointer string
#'
#' This is a low-level function that may be useful for debugging.  It will
#' return the pointer name itself from an antsImage object.
#'
#' @param img  image whose pointer we want
#' @return string
#' @author Avants BB
#' @examples
#'
#' img = antsImageRead( getANTsRData("r16") )
#' antsrGetPointerName( img )
#' antsrGetPointerName( antsImageClone( img ))
#'
#' @export antsrGetPointerName
antsrGetPointerName <- function(img) {
  # if ( Sys.info()['sysname'] == 'Linux' ) endofpointer<-20 if (
  # Sys.info()['sysname'] == 'Darwin' ) endofpointer<-21 pname<- substr(
  # .int_antsProcessArguments( list( img ) ) , 11 , endofpointer )
  splitptrname <- strsplit(.int_antsProcessArguments(list(img)), " ")[[1]][2]
  pname <- strsplit(splitptrname, ">")
  return(pname[[1]])
}

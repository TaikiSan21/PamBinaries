# Reads the module footer information common to all modules. Differs from
# the legacy code in that it does not read in or skip any information
# specific to a module

readStdModuleFooter <- function(file) {
      footer$length <- pamBinRead(file, 'int32', n=1)
      footer$identifier <- pamBinRead(file, 'int32', n=1)
      footer$binaryLength <- pamBinRead(file, 'int32', n=1)
      return(footer)
}
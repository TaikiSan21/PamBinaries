# Reads module footer informatiopn for the Click Detector module. Note that
# sometimes there is no additional footer information, so check first
# whether or not the binaryLength variable is 0

readClickFooter <- function(file) {
    footer <- readStdModuleFooter(file)
    if(footer$binaryLength != 0) {
        footer$typesCountLength <- pamBinRead(file, 'int16', n=1)
        footer$typesCount <- pamBinRead(file, 'int32', n=footer$typesCountLength)
    }
    return(footer)
}
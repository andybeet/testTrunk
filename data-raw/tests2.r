# Runs a series of tests on Atlantis output to signify "passing"
#
#

test_compare <- function() {
  # read in bio out file for both runs and output the last line
  current <- readLines("example/testFolder/outputSETASBiomIndx.txt")
  # select last line
  print(tail(current, n=1))
  
  release <- readLines("example/testVersionFolder/outputSETASBiomIndx.txt")
  
  print(tail(release, n=1))

}
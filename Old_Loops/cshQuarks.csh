#!/bin/csh

# Set run parameters
set exeFlags = "-np 4 "
set beta = "b460"
set size = "s16t32"
set imp  = "IMP"
set dir  = "su3"$beta$size$imp
set baseConfig = "/orioni/dleinweb/"$dir"/"$dir"c"
set alpha = "0.7"
set nsteps = "4"
set yorn = ".false."
set prefix = "./"

set exeName = "ConstituentQuarks"

foreach cfg ( 072 073 074 )

set thisConfig = $baseConfig$cfg
echo $baseConfig$cfg
echo $alpha
echo $nsteps
echo $yorn

# Run the parallel program
  echo `date`
  echo  "mprun -Mf nodes05 $exeFlags $exeName" 
  mprun -Mf nodes05 $exeFlags $exeName << ....END
$thisConfig
$prefix
3  three-loop improved fMuNu
1  1: action and topological charge, 2: electric and magnetic fields
$yorn
$alpha
$nsteps
....END

end



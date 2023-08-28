#!/usr/bin/env python3

import re
import sys
import os

subroutine = re.compile( r"subroutine\s*(?P<routine>\w+)(?P<def>.*?)end\s*subroutine\s*(?P=routine)", re.M | re.S | re.I )

def main():
  module         = sys.argv[1]
  source         = sys.argv[2]
  prefix         = os.path.splitext( os.path.basename( source ) )[0]
  suffix         = os.path.splitext( os.path.basename( source ) )[1]

  outputDir      = sys.argv[3]
  submodTemplate = sys.argv[4]
  numSubmodFiles = int( sys.argv[5] )

  declSuffix     = "_decl" + suffix
  defSuffix      = "_def"  + ".F"

  fp    = open( source, 'r' )
  lines = fp.read()
  fp.close()

  fp     = open( submodTemplate, "r" )
  submod = fp.read()
  fp.close()

  

  subroutineFormat = "module subroutine {name}{rest}\nend subroutine {name}\n"
  submodDecl = ""
  submodDef  = [""] * numSubmodFiles

  idx      = 0
  for sub in subroutine.finditer( lines ) :
    subDef = sub.group( "def" )

    posLastIntent = subDef.rfind( "INTENT" )
    if posLastIntent == -1 :
      posLastIntent = subDef.rfind( "intent" )

    if posLastIntent == -1 :
      # This subroutine does not have any formal args, just skip all this
      submodDecl += subroutineFormat.format( name=sub.group( "routine" ), rest="" )
    else:
      endPos = subDef.find( "\n", posLastIntent )
      submodDecl += subroutineFormat.format( name=sub.group( "routine" ), rest=subDef[:endPos] )
    
    submodDef[idx] += subroutineFormat.format( name=sub.group( "routine" ), rest=subDef )

    idx = ( idx + 1 ) % numSubmodFiles

  with open( outputDir + "/" + prefix + declSuffix, "w" ) as f :
    f.write( submodDecl )
    f.close()

  # Now write out source files
  for i in range( 0, numSubmodFiles ) :
    submodSource = outputDir + "/" + prefix + "_{id}".format( id=i ) + defSuffix
    with open( submodSource, "w" ) as f :
      f.write( submod.format( id=i, source=submodDef[i], module=module ) )
      f.close()

if __name__ == "__main__" :
  main()
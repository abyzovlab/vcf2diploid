#!/bin/csh

set VERSION    = "v0.2.6"
set distr      = $PWD/"vcf2diploid_"$VERSION".zip"

set tmpdir     = "/tmp"
set vcfdir     = "vcf2diploid"
set maindir    = $tmpdir"/"$vcfdir

rm -rf $maindir
mkdir  $maindir

cp *.java      $maindir
cp *Manifest   $maindir
cp Makefile    $maindir
cp README      $maindir
cp CITATION    $maindir
cp license.rtf $maindir
cp -r example  $maindir

rm -fr $maindir/example/.svn

cd $tmpdir
zip -r $distr  $vcfdir

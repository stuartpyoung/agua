#!/bin/sh

#### GENERATE CHROMOSOME SUB-DIRECTORIES AND COPY FILES INTO SUBDIRS
#echo "directory: " $1
#cd /nethome/bioinfo/data/sequence/chromosomes/mouse/mm9/maq    

cd $1
echo "\$#: " $#

#### CONFIRM ARGUMENTS
echo "Arguments: ";
#for ((index = 0; index < $#; index++)); ## OK
for argument in "$@"
do
    echo "    $argument"
done;

#### CHECK ARGUMENT NUMBER
if [ $# = 0 ];
    then
        echo "Please use two arguments: $0 <directory String> <pattern String>";
    exit;
fi

#if [ $($# > 2) ];
if [ $# != 2 ];
    then
        echo "Please use two arguments:"
        echo "";
        echo "$0 <directory String> <pattern String>";
        echo "";
    exit;
fi

#### 
DIRS="chr*$2"
for DIR in $DIRS;
do
    CHROMO=$(echo $DIR | egrep "^chr[0-9XYM]+$2" ) ### OK!
    if [ "$CHROMO" != "" ]; then
        CHROM=${CHROMO/$2/};
        echo "Doing chromosome: " $CHROM;
        mkdir `pwd`/$CHROM;
        cp -f `pwd`/$CHROM.* `pwd`/$CHROM;        
    fi;
done;


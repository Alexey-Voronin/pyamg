 DIR1=/Users/lexey/pyamg/build/lib.macosx-10.7-x86_64-3.7/pyamg/amg_core
 #build/lib.macosx-10.7-x86_64-3.7/pyamg/amg_core/
 cd $DIR1
 echo `pwd`
 for filename in ./*.so; do
     name=`basename $filename .cpython-37m-darwin.so`
     ext='.so'
     newname="$name$ext"
     echo $filename
     echo $name
     echo $filanme
     echo $newname
     cp $filename ~/pyamg/pyamg/amg_core/$newname
 done

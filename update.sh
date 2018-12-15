#export CXX='clang++'
#export CC='clang'

cd pyamg/amg_core
./generate.sh
rm -f *.so # runtime linked files
cd ../..
rm -rf build/ dist/ pyamg.egg-info/ pyamgc
python setup.py develop

DIR1=/Users/lexey/pyamg/build/lib.macosx-10.7-x86_64-3.7/pyamg/amg_core
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


python -c 'import pyamg'

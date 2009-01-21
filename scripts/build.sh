automake --add-missing
autoreconf
# -pedantic 
# 
export CXXFLAGS="-O6 -Wall -std=c++98 -fomit-frame-pointer -static"
export LDFLAGS=$CXXFLAGS
mkdir -p build
./configure CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS" --prefix=$(pwd)/build
make install

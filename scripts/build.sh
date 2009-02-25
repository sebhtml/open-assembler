aclocal
autoconf
automake --add-missing
autoreconf
# -pedantic 
# 

#export CXXFLAGS="-O6 -Wall -pedantic -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -std=c++98 -fomit-frame-pointer -static -funroll-loops -Weffc++ -Wold-style-cast "
export CXXFLAGS=" -O3 -Wall " # -I /usr/include/mysql -L/usr/lib/mysql -lmysqlclient"
export LDFLAGS=$CXXFLAGS
mkdir -p build
./configure CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS" --prefix=$(pwd)/build
make install

rm -r dist
mkdir -p dist
revision=$(svn info https://denovoassembler.svn.sourceforge.net/svnroot/denovoassembler|grep Revision|awk '{print $2}')
cd dist
svn export https://denovoassembler.svn.sourceforge.net/svnroot/denovoassembler/trunk/src
mv src 454dna-v$revision-source
tar cjf 454dna-v$revision-source.tar.bz2 454dna-v$revision-source
cd 454dna-v$revision-source
echo "Making dist $revision"
make
cd ..


mkdir 454dna-v$revision-x86_64
ls
for i in 454dna LICENSE README 
do
	cp 454dna-v$revision-source/$i 454dna-v$revision-x86_64
done
tar cjf 454dna-v$revision-x86_64.tar.bz2 454dna-v$revision-x86_64

mkdir 454dna-v$revision-x86_64-static
ls
for i in 454dna.static LICENSE README 
do
	cp 454dna-v$revision-source/$i 454dna-v$revision-x86_64-static
done
tar cjf 454dna-v$revision-x86_64-static.tar.bz2 454dna-v$revision-x86_64-static

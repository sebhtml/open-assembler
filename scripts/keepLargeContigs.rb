#!/usr/bin/ruby

if ARGV.size!=3
	puts "usage"
	puts "keepLargeContigs.rb <contigsFile> <minimumContigSize> <largeContigsFile>"
	exit
end

seq=""

contigs=[]
f=File.open ARGV[0]
while l=f.gets
    if l[0..0]=='>'
        contigs<< seq
        seq=""
    else
        seq<< l.strip
    end
end

contigs<< seq
f.close

threshold=ARGV[1].to_i
k=1
out=File.open ARGV[2],"w+"
contigs.each do |i|
    if i.length<threshold
        next
    end
    j=0
    out.puts ">#{k} #{i.length}"
    columns=60
    while j<i.length
        out.puts i[j..(j+columns-1)]
        j+=columns
    end
    k+=1
end

out.close

#!/usr/bin/ruby

if ARGV.size!=2
	puts "usage"
	puts "dumpQuality.rb <contigsFile> <output>"
	exit
end

seq=""
name=""
contigs=[]
f=File.open ARGV[0]
while l=f.gets
    if l[0..0]=='>'
	if name!=""
        	contigs<< [name,seq]
	end
	name=l.split(" ").first.strip
        seq=""
    else
        seq<< l.strip
    end
end

contigs<< [name,seq]
f.close

k=1
out=File.open ARGV[1],"w+"
contigs.each do |i|
    out.puts "#{i[0]} #{i[1].length}"
	i[1].length.times do
	out.print " 40"
    	end
	out.puts ""
end

out.close

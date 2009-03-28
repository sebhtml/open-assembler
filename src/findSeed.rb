#!/usr/bin/ruby


if ARGV.size!=2
	puts "usage findSeed.rb sequenceFile seed"
	exit
end

seq=""
f=File.open ARGV[0]
f.gets

seed=ARGV[1]

#puts "Using seed #{seed}"
while l=f.gets
	seq<< l.strip
end
f.close

#puts "Sequence: #{seq}"
i=0
while i<seq.length
	word=seq[i..(i+seed.length()-1)]
	if word==seed
		puts "Offset: #{i}"
	end
	i+=1
end


#!/usr/bin/ruby

if ARGV.size==0
	puts "You must provide a file"
	exit
end

def revComp a
	b=""
	i=a.length-1
	while i>=0
		s=a[i..i]
		if s=='A'
			b<< 'T'
		elsif s=='T'
			b<< 'A'
		elsif s=='C'
			b<< 'G'
		elsif s=='G'
			b<< 'C'
		end
		i-=1
	end
	b
end

chromosomes=[]
f=File.open  ARGV[0]
seq=""
while l=f.gets
	l=l.upcase
	if l[0..0]=='>'
		if seq!=""
			chromosomes<< seq
		end
		seq=""
	else
		seq<< l.strip
	end
end

chromosomes<< seq
f.close

coverage=20
readLength=250
readID=1
chromosomes.each do |genome|
	gSize=genome.length
	position=0
	while position<gSize
		coverage.times do |t|
			read_length=readLength+rand(100)-50
			start=position+rand(read_length)-read_length/2
			if start<0
				start=0
			end
			if start>=gSize
				next
			end
			sequence=genome[start..(start+read_length)]
			p=rand(sequence.length-30)+30
			if sequence.length < p
				next
			end
			if sequence[p..p]=='T'
				sequence[p..p]='A'
			else
				sequence[p..p]='T'
			end
			if rand(2)==0
				puts ">#{readID}_#{start}_#{read_length}_F"
				puts sequence
			else
				puts ">#{readID}_#{start}_#{read_length}_R"
				puts revComp(sequence)
			end
			readID+=1
		end
		position+=readLength
	end
end

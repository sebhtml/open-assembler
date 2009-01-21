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

coverage=10
readLength=250
readID=1
errors=0
chromosomes.each do |genome|
	gSize=genome.length
	position=0
	while position<gSize
		coverage.times do |t|
			read_length=readLength
			start=position
			if start<0
				start=0
			end
			sequence=genome[start..(start+read_length)]
			if sequence.nil?
				next
			end
			errorsInRead=errors+rand(4)-2
			errorsInRead.times do
				break
				n=rand(4)
				p=rand(sequence.length)
				if n==0
					sequence[p..p]='A'
				elsif n==1
					sequence[p..p]='T'
				elsif n==2
			sequence[p..p]='C'
		elsif n==3
			sequence[p..p]='G'
				end
			end
			if rand(2)==0
		puts ">#{readID}_#{start}_#{read_length}_F_#{errorsInRead}"
		puts sequence
			else
		puts ">#{readID}_#{start}_#{read_length}_R_#{errorsInRead}"
		puts revComp(sequence)
			end
			readID+=1
		end
		position+=readLength
	end
end

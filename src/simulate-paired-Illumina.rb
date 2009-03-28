#!/usr/bin/ruby

if ARGV.size==0
	puts "You must provide a file"
	puts "usage:"
	puts "ruby script fastaFile fragmentSize readLength"
	puts "example: "
	puts "ruby script file.fasta 200 36"
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

prefix="#{ARGV[1]}x#{ARGV[2]}x#{ARGV[2]}-"
leftFile=File.open prefix+ARGV[0]+"_1.fasta","w+"
rightFile=File.open prefix+ARGV[0]+"_2.fasta","w+"

sequencedLength=ARGV[2].to_i
readLength=ARGV[1].to_i
readID=1
errors=0
chromosomes.each do |genome|
	gSize=genome.length
	position=0
	while position<gSize
		if position%1000==0
			puts "#{position} / #{gSize}"
		end
		2.times do |t|
			start=position
			if start<0
				start=0
			end
			if start>=gSize
				next
			end
			sequence=genome[start..(start+readLength)]
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
				sequence=revComp sequence
			end
			leftFile.puts ">r#{readID}_1"
			leftFile.puts sequence[0..(sequencedLength-1)]
			rightFile.puts ">r#{readID}_2"
			rightFile.puts sequence[(sequence.length-sequencedLength)..(sequence.length-1)]
			readID+=1
		end
		position+=sequencedLength/2
	end
end

leftFile.close
rightFile.close

puts ""

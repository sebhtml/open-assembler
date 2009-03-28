#!/usr/bin/ruby

if ARGV.size!=3
	puts "usage"
	puts "script reference.fasta contigs.fasta output.txt"
	puts "You need blat!"
	exit
end

reference=ARGV[0]
contigs=ARGV[1]
output=ARGV[2]

outputStream=File.open output,"w+"
system "blat #{reference} #{contigs} #{output}_blat"
system "cat #{contigs}|grep '>'|awk '{print $1}'|sed 's/>//' >  #{contigs}.names"

contigNames=[]
f=File.open "#{contigs}.names"
while l=f.gets
	if l.strip!=""
		contigNames<< l.strip
	end
end

f.close

blasts={}
f=File.open "#{output}_blat"
5.times do
	f.gets
end

contigSizes={}

while l=f.gets
	tokens=l.split "\t"
	matches=tokens[1-1].to_i
	qName=tokens[10-1].strip
	qSize=tokens[11-1].to_i
	contigSizes[qName]=qSize
	if blasts[qName].nil?
		blasts[qName]=[]
	end
	blasts[qName]<< tokens
end
f.close

misassembled=[]
notFound=[]
assembledSize=0
notFoundSize=0

contigNames.each do |contig|
	outputStream.puts "#{contig} - #{contigSizes[contig]}"
	if contigSizes[contig].nil?
		outputStream.puts "NOT FOUND!"
		notFound<< contig
		next
	end
	assembledSize+=contigSizes[contig]
	outputStream.puts "blocks "
	unless blasts[contig].nil?
		largestBlock=blasts[contig].first[1-1].to_i
		blasts[contig].each do |blast|
			outputStream.puts "  #{blast[1-1]} matches, on target: #{blast[16-1]}-#{blast[17-1]}"
			if blast[1-1].to_i > largestBlock
				largestBlock=blast[1-1].to_i
			end
		end
		diff=largestBlock-contigSizes[contig]
		if diff<0
			diff=-diff
		end
		ratio=diff/contigSizes[contig].to_f
		outputStream.puts "Ratio: #{ratio}"
		if ratio > 0.1
			misassembled<< contig
		end
	end

end

outputStream.puts  "Statistics"
outputStream.puts "total assembled contigs found in  reference: #{assembledSize} bases"
#puts "total assembled  contigs not found in reference: #{notFoundSize} bases"
outputStream.puts "Misassembled contigs: #{misassembled.size}"
outputStream.puts "#{misassembled.join "\n"}"
outputStream.puts "Not found: #{notFound.size}"
outputStream.puts "#{notFound.join "\n"}"

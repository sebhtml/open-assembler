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

system "#blat #{reference} #{contigs} #{output}_blat"
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

contigNames.each do |contig|
	puts "#{contig} - #{contigSizes[contig]}"
	if contigSizes[contig].nil?
		puts "NOT FOUND!"
		next
	end
	puts "blocks "
	unless blasts[contig].nil?
		blasts[contig].each do |blast|
			puts "  #{blast[1-1]} matches, on target: #{blast[16-1]}-#{blast[17-1]}"
		end
	end
end

#!/usr/bin/ruby


puts "Welcome to DNA, the de novo assembler"
puts "DNA.rb is an assembler that calls the c++ modules"
puts "usage:"
puts "DNA.rb [-wordSize 21] [-minimumCoverage 2] -directory output <sequence files (fasta or sff or fastq)"
puts "Note: if your sff files contain  read pairs,  you must utilize dna_ConvertSffToFasta to extract paired information."
puts "DNA.rb will detect fastq or 454 paired information, provided files are named the same except a '_1.*' and '_2.*' suffixes."
puts ""
files=[]
directory="unamed-assembly"
k=21
c=2

i=0

puts "#{ARGV.size} arguments"

while i<ARGV.size
	if ARGV[i]=="-wordSize"&&i<ARGV.size-1
		i+=1
		k=ARGV[i].to_i
	elsif ARGV[i]=="-minimumCoverage"&&i<ARGV.size-1
		i+=1
		c=ARGV[i].to_i
	elsif ARGV[i]=="-directory"&&i<ARGV.size-1
		i+=1
		directory=ARGV[i]
	else
		files<< ARGV[i]
	end
	i+=1
end

if files.size==0
	puts "Error: no files provided."
	exit
end

puts "Parameters"
puts "-wordSize #{k}"
puts "-minimimCoverage #{c}"
puts "-directory #{directory}"
puts "files #{files.join(',')}"

puts "Starting now."


system "mkdir -p #{directory}"

system "dna_GetPairedInformation.rb -directory #{directory} #{files.join ' '} > #{directory}/dna_GetPairedInformation.rb.log"

system "
date > #{directory}/START
dna_BuildGraph -directory #{directory} -minimumCoverage #{c} -wordSize #{k} #{files.join ' '} > #{directory}/dna_BuildGraph.log
dna_ExtractContigs -directory #{directory} > #{directory}/dna_ExtractContigs.log
dna_KeepLargeContigs #{directory}/contigs.fasta 500 #{directory}/2LargeContigs.fasta > #{directory}/dna_KeepLargeContigs.log
dna_MergeContigs  #{directory}/2LargeContigs.fasta  #{directory}/3MergedContigs.fasta > #{directory}/dna_MergeContigs.log
dna_JoinContigs #{k}  #{directory}/3MergedContigs.fasta #{directory}/4JoinedContigs.fasta #{files.join ' '} > #{directory}/dna_JoinContigs.log
date > #{directory}/END
"


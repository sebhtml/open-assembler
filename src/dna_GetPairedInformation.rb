#!/usr/bin/ruby

puts "dna_GetPairedInformation"
puts "This program gets paired information"

puts "usage: "
puts "dna_GetPairedInformation -directory <directory> file1 file2 [....] filen"
puts "it is provided as is and might require changes as it uses default values"
puts "Default insert sizes: "
puts "   Illumina/Solexa: 150"
puts "   Roche/454: 2500"
puts "To edit insert sizes, simply edit <directory>/PairedReads.txt"

files=[]
directory=nil

i=0

while i<ARGV.size
	if ARGV[i]=="-directory"&&i<ARGV.size-1
		i+=1
		directory=ARGV[i]
	else
		files<< ARGV[i]
	end
	i+=1
end

if files.size==0
	puts "No file provided..."
	puts "exiting..."
	exit
end

if directory.nil?
	puts "Directory is nil, exiting..."
	exit
end

system "mkdir -p #{directory}"

pairedFile=File.open "#{directory}/PairedReads.txt","w+"
fileHash={}
files.each do |file|
	fileWithoutExtension=file
	supportedFormat=['_1.fasta', '_2.fasta', '_1.fastq','_2.fastq']
	supportedFormat.each do |aFormat|
		if file.include? aFormat
			fileWithoutExtension=file.gsub aFormat,''
		end
	end
	if fileWithoutExtension!=file
		if fileHash[fileWithoutExtension].nil?
			fileHash[fileWithoutExtension]=[]
		end
		fileHash[fileWithoutExtension]<< file
	end
end

pairedFile.puts "#{fileHash.size}"
fileHash.each do |group,files|
	if files.size==2
		distance=150
		if files[0].include? 'sff'
			distance=2500
		end
		pairedFile.puts "#{files[0]} #{files[1]} #{distance}"
	end
end
pairedFile.close


#!/usr/bin/ruby

if ARGV.size!=1
	puts "Usage: "
	puts "findMutations.rb blast-output-file"
	exit
end

sequence=[]
f=File.open ARGV[0]
while l=f.gets
	if (l.include? 'Query:')||(l.include? 'Sbjct:')
		sequence<< l.strip
	end
end

f.close

i=0

while i<sequence.size
	querySequence=sequence[i]
	i+=1
	subjectSequence=sequence[i]
	i+=1
	referenceBase=subjectSequence.split(" ")[1].to_i
	sampleSequence=querySequence.split(" ")[2].strip	
	referenceSequence=subjectSequence.split(" ")[2].strip
	j=0
	while j<sampleSequence.length
		referenceNucleotide=referenceSequence[j..j].upcase
		sampleNucleotide=sampleSequence[j..j].upcase
		if sampleNucleotide!=referenceNucleotide
			puts "#{referenceBase+j}#{referenceNucleotide}->#{sampleNucleotide}"
		end
		j+=1
	end
end

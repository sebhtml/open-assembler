#!/usr/bin/ruby

require 'set'

if ARGV.size!=2
	puts "usage"
	puts "plotGraph.rb fastaFile graphFile"
	puts "  Use graphviz to plot the graph"
	puts "  dot -T svg graphFile -o graphFile.svg"
	exit
end

f=File.open ARGV[0] 
f.gets
seq=""

while l=f.gets
	seq<< l.strip
end
f.close

#puts seq

graph={}
k=22
i=0
while i<seq.length
	word=seq[i..(i+k-1)]
	i+=1
	if word.length!=k
		next
	end
	prefix=word[0..(k-2)]
	suffix=word[1..(k-1)]
	if graph[prefix].nil?
		graph[prefix]=Set.new
	end
	graph[prefix]<<  suffix
end

f=File.open  ARGV[1],'w+'
f.puts "digraph G{"
graph.each do |prefix,suffixes|
	suffixes.each do |suffix|
		f.puts "#{prefix} -> #{suffix}"
	end
end
f.puts "}"

f.close

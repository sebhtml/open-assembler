#!/usr/bin/ruby

require 'set'

if ARGV.size!=3
	puts "usage"
	puts "plotGraph.rb wordSize fastaFile graphFile"
	puts "  Use graphviz to plot the graph"
	puts "  dot -T svg graphFile -o graphFile.svg"
	exit
end

f=File.open ARGV[1] 
f.gets
seq=""

while l=f.gets
	seq<< l.strip.gsub(' ','')
end
f.close

#puts seq

graph={}
k=ARGV[0].to_i
i=0
while i<seq.length
	word=seq[i..(i+k)]
	i+=1
	if word.length!=(k+1)
		next
	end
	prefix=word[0..(k-1)]
	suffix=word[1..(k-0)]
	if graph[prefix].nil?
		graph[prefix]=Set.new
	end
	graph[prefix]<<  suffix
end

f=File.open  ARGV[2],'w+'
f.puts "digraph G{"
graph.each do |prefix,suffixes|
	suffixes.each do |suffix|
		f.puts "#{prefix} -> #{suffix}"
	end
end
f.puts "}"

f.close

puts "Type this now:"
puts "dot -T svg #{ARGV[2]} -o #{ARGV[2]}.svg"



while l=STDIN.gets
	if l.strip=="Writing contigs"
		STDIN.gets
		contigs=STDIN.gets.split(" ").first.to_i
		puts contigs.to_s
		contigs.times do
			puts STDIN.gets
		end
	end
end

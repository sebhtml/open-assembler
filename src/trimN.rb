while l=STDIN.gets
	puts ">"+l.split(' ').first.gsub('>','')+"__no_N"
	seq=STDIN.gets
	j=0
	while seq[j..j]!='N'&&j<seq.length
		j+=1
	end
	j-=1
	puts seq[0..j]
end

#!/usr/bin/ruby

puts ">"
10000.times do 
	i=rand(4)
	if i==0
		print 'A'
	elsif i==1
		print 'T'
	elsif i==2
		print 'C'
	elsif i==3
		print 'G'
	end
end

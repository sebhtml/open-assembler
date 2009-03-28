#!/usr/bin/ruby

theLength.times do 
	i=rand(4)
	if i==0
		g<< 'A'
	elsif i==1
		g<< 'T'
	elsif i==2
		g<< 'C'
	elsif i==3
		g<< 'G'
	end
end

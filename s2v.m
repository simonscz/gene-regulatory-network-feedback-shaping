function r = s2v(a) 
	if a ~= 0 && a ~= 1
		error('The input argument must be 0 or 1');
	end
    r = [1-a, a]';
end
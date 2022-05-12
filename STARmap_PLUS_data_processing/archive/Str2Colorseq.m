function cs = Str2Colorseq( s )
%STR2COLORSEQ convert string color sequence back to num sequence 
cs = [];
for i=1:numel(s)
    cs = [cs str2num(s(i))];
end

end


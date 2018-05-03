clear
m=1;
n=25;
a=1;
for ii=m:n
    a=lcm(a,ii);
end;
b=1;
for jj=21:25
    b=lcm(b,jj);
end
a=a+111;
expand(lcm(sym(a),sym(b)))
save a.txt a -ascii -double
26771144400

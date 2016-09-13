%function out=circularize(in)
%
%Maps the input to the interval [0,2pi], such that negative values are
%mapped as 2pi-in

function out=circularize(in)

out=in;

for k=1:length(out)
    while out(k)<0 || out(k)>2*pi
        if out(k)>2*pi
            out(k)=out(k)-2*pi;
        elseif out(k)<0
            out(k)=2*pi-abs(out(k));
        end
    end
end

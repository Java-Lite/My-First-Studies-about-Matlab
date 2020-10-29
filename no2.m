n=[1 2 3 4 5];

filter(n,1,2)


function fl = filter(matrix,var1,var2)
for i = 1:5
  fil = matrix(1,i);
  
        if fil == var1 || fil == var2
%              disp(matrix(1,var1))
        else
          out(fil) = fil;
        end
    end
    out
end
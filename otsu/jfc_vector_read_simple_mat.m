function data=jfc_vector_read_simple_mat(filename)

fid = fopen(filename,'r');
m=fscanf(fid,'%d\n',1);
n=fscanf(fid,'%d\n',1);
data=zeros(m,n);
fprintf('size: m= %d, n= %d\n',m,n);
for i=1:m
    data(i,:)=fscanf(fid,'%g\n',n);
end
fclose(fid);

end
function jfc_vector_save_simple_mat(data,filename)

[m,n]=size(data);
fid = fopen(filename,'w');
fprintf(fid,'%d\n',m);
fprintf(fid,'%d\n',n);
for i=1:m
    fprintf(fid,'%24.18g\n',data(i,:));
end
fclose(fid);

end
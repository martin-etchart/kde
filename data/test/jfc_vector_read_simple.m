function data=jfc_vector_read_simple(filename)

	fid = fopen(filename,'r');
	n=fscanf(fid,'%d\n',1);
    fprintf('n: %d\n',n);
	data=fscanf(fid,'%g\n',n);
	fclose(fid);

end
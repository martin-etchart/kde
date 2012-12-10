function jfc_vector_save_simple(data,filename)

	fid = fopen(filename,'w');
	fprintf(fid,'%d\n',length(data));
	fprintf(fid,'%24.18g\n',data);
	fclose(fid);

end
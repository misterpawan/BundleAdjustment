[i,j,val] = find(JTJ);     % find to get index & value vectors
data_dump = [i,j,val];     % save this in data_dump
fid = fopen('/home/siddhant/Desktop/JTJ.txt','w');    % open file in write mode
fprintf( fid,'%d %d %f\n', transpose(data_dump));     % print the JTJ to the file in COO 
fclose(fid);     % close file 

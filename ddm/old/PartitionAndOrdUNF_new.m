function [perm, sizes] = PartitionAndOrdUNF_new(matrix, nc)
   
input_parameters_ndp; 
metis_aggregation_ndp; 

perm = load('perm_array.txt');
sizes = load('size_array.txt');

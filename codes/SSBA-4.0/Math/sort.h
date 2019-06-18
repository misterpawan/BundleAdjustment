#ifndef V3D_SORT_H
#define V3D_SORT_H

// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(int arr[], double val[],int l, int m, int r) 
{ 
	int i, j, k; 
	int n1 = m - l + 1; 
	int n2 = r - m; 

	/* create temp arrays */
	int L[n1], R[n2]; 
	double Lval[n1], Rval[n2];

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++) 
	{
		L[i] = arr[l + i]; 
		Lval[i] = val[l+i];
	}
	for (j = 0; j < n2; j++)
	{ 
		R[j] = arr[m + 1+ j]; 
		Rval[j] = val[m +1 + j];
	}
	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray 
	j = 0; // Initial index of second subarray 
	k = l; // Initial index of merged subarray 
	while (i < n1 && j < n2) 
	{ 
		if (L[i] <= R[j]) 
		{ 
			arr[k] = L[i]; 
			val[k] = Lval[i];
			i++; 
		} 
		else
		{ 
			arr[k] = R[j]; 
			val[k] = Rval[j];
			j++; 
		} 
		k++; 
	} 

	/* Copy the remaining elements of L[], if there 
	are any */
	while (i < n1) 
	{ 
		arr[k] = L[i]; 
		val[k] = Lval[i];
		i++; 
		k++; 
	} 

	/* Copy the remaining elements of R[], if there 
	are any */
	while (j < n2) 
	{ 
		arr[k] = R[j]; 
		val[k] = Rval[j];
		j++; 
		k++; 
	} 
} 

/* l is for left index and r is right index of the 
sub-array of arr to be sorted */
void mergeSort(int arr[], double val[],int l, int r) 
{ 
	if (l < r) 
	{ 
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		int m = l+(r-l)/2; 

		// Sort first and second halves 
		mergeSort(arr, val,l, m); 
		mergeSort(arr, val, m+1, r); 

		merge(arr, val, l, m, r); 
	} 
} 

/* UTILITY FUNCTIONS */
/* Function to print an array */
void printArray(int A[], int size) 
{ 
	int i; 
	for (i=0; i < size; i++) 
		printf("%d ", A[i]); 
	printf("\n"); 
} 

/* Driver program to test above functions */
void sort_rows_val(int* i, double* x, int nzmax) 
{ 
	//int arr[] = {12, 11, 13, 5, 6, 7}; 
	//int arr_size = sizeof(arr)/sizeof(arr[0]); 

	//printf("Given array is \n"); 
	//printArray(arr, arr_size); 

	mergeSort(i, x, 0, nzmax-1); 

	//printf("\nSorted array is \n"); 
	//printArray(i, nzmax); 
	//return 0; 
} 
#endif
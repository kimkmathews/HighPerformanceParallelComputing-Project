#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <string.h>
#include <math.h>

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

// Function to find the maximum element of the array
int max_element(int array[], int size) 
{  
    // Initializing max variable to minimum value so that it can be updated
    // when we encounter any element which is greater than it.
    int max = INT_MIN;  
    for (int i = 0; i < size; i++)
    {
        //Updating max when array[i] is greater than max
        if (array[i] > max)  
        max = array[i];  
    }
    //return the max element
    return max;  
}

int min_element(int array[], int size) {
    // Initialize min variable to the first element of the array
    int min = array[0];
    for (int i = 1; i < size; i++) {
        // Update min if the current element is smaller
        if (array[i] < min) {
            min = array[i];
        }
    }
    // Return the minimum element
    return min;
}
void printArray(int array[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}

int sort_check(int array[], int size) {
  // Iterate through the array from the second element onwards
  for (int i = 1; i < size; i++) {
    // If the current element is less than the previous element, the array is not sorted
    if (array[i] < array[i - 1]) {
      printf("Array is not sorted.\n");
      return 0; // Indicate unsuccessful sort check (not sorted)
    }
  }

  return 1; // Indicate successful sort check (sorted)
}
// Function to sort a bucket using bubble sort
void bubbleSort(int* bucket, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (bucket[j] > bucket[j + 1]) {
                int temp = bucket[j];
                bucket[j] = bucket[j + 1];
                bucket[j + 1] = temp;
            }
        }
    }
}

void mergeSort(int arr[], int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;

        // Sort first and second halves
        mergeSort(arr, left, mid);
        mergeSort(arr, mid + 1, right);

        // Merge the sorted halves
        int i, j, k;
        int n1 = mid - left + 1;
        int n2 = right - mid;

        int* L = (int*)malloc(n1 * sizeof(int));
        int* R = (int*)malloc(n2 * sizeof(int));

        if (L == NULL || R == NULL) {
            printf("Memory allocation failed. Exiting...\n");
            return;
        }

        for (i = 0; i < n1; i++)
            L[i] = arr[left + i];
        for (j = 0; j < n2; j++)
            R[j] = arr[mid + 1 + j];

        i = 0;
        j = 0;
        k = left;
        while (i < n1 && j < n2) {
            if (L[i] <= R[j]) {
                arr[k] = L[i];
                i++;
            } else {
                arr[k] = R[j];
                j++;
            }
            k++;
        }

        while (i < n1) {
            arr[k] = L[i];
            i++;
            k++;
        }

        while (j < n2) {
            arr[k] = R[j];
            j++;
            k++;
        }
    }
}
void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
int partition(int array[], int low, int high) {
  
  // select the rightmost element as pivot
  int pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swap(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSort(int array[], int low, int high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partition(array, low, high);
    
    // recursive call on the left of pivot
    quickSort(array, low, pi - 1);
    
    // recursive call on the right of pivot
    quickSort(array, pi + 1, high);
  }
}
// Implementing bucket sort using OpenMP 
void Bucket_Sort_1(int array[], int size) 
{  
    //Finding max element of array which we will use to create buckets
    int max = max_element(array, size); 
 
    int* bucket = (int*)malloc((max + 1) * sizeof(int));
    if (bucket == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
    //Initializing buckets to zero
    for (int i = 0; i <= max; i++)  
    bucket[i] = 0;  
 
    // Pushing elements in their corresponding buckets
    for (int i = 0; i < size; i++)  
        bucket[array[i]]++;
 
    // Merging buckets effectively
    int j=0;   // j is a variable which points at the index we are updating
    for (int i = 0; i <= max; i++)  
    { 
        // Running while loop until there is an element in the bucket
        while (bucket[i] > 0)  
        {  
            // Updating array and increment j          
            array[j++] = i;  
 
            // Decreasing count of bucket element
            bucket[i]--;   
        }  
    }  
}  

// Implementing bucket sort using OpenMP 
void Bucket_Sort_2(int array[], int size, int nThreads) {
    // Finding the max element of the array which we will use to create buckets
    int max = max_element(array, size);
    int threadNum , threadStart ,threadEnd ; 

    int* bucket = (int*)malloc((max + 1) * sizeof(int));
    if (bucket == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
    
    // Initializing buckets to zero
    #pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < max + 1; i++)
        bucket[i] = 0;


    // Pushing elements into their corresponding buckets
/*   #pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < size; i++) {
        int index = (array[i] * nBuckets) / (max + 1);
        #pragma omp atomic
        bucket[index]++;
    }
*/
    #pragma omp parallel for num_threads(nThreads) 
    for (int i = 0; i < size; i++)
        bucket[array[i]]++;
/*
for (int i = 0; i < max+1; i++)
        printf("%d ", bucket[i]);   
    printf("\n"); 
*/
   
   // Merging buckets effectively
    int j=0;   // j is a variable which points at the index we are updating
    //#pragma omp parallel for num_threads(nThreads) 
    for (int i = 0; i <= max; i++)  
    { 
        // Running while loop until there is an element in the bucket
        while (bucket[i] > 0)  
        {  
            // Updating array and increment j          
            // #pragma omp critical
            array[j++] = i;  
 
            // Decreasing count of bucket element
            // #pragma omp critical
            bucket[i]--;   
        }  
    }  
   free(bucket);
  
}
void Bucket_Sort_3(int array[], int size, int numBuckets) {
    
    int** buckets = (int**)malloc(numBuckets * sizeof(int *));

    if (buckets == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
     // Initialize each bucket pointer to NULL 
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = NULL;
    }
    
    int max = max_element(array, size);

    int min = min_element(array,size);

    // Calculate bucket size and elements per bucket
    int bucket_size = (max - min) / numBuckets;

    int elements_per_bucket = (size / numBuckets) + 1;
    
    //Allocate memory for bucket size limit and filled counts
    int bucket_size_limit[numBuckets];

    int bucket_size_fill[numBuckets];

    // Step 1: Initialize an array of empty buckets
    for (int i = 0; i < numBuckets; i++) {
        //buckets[i] = 0;
        bucket_size_limit[i] = elements_per_bucket;
        bucket_size_fill[i] =   0;
        buckets[i] = (int *)malloc(elements_per_bucket * sizeof(int));
        }
    // Step 2: Bucket assignment
    for (int i = 0; i < size; i++) {
        int bucket_id = array[i]/bucket_size;
        
        if (bucket_id >= numBuckets) {
            bucket_id = numBuckets - 1;
        }
        if (bucket_size_fill[bucket_id] >= bucket_size_limit[bucket_id]) {
                bucket_size_limit[bucket_id] *= 1.1; // Double capacity if needed
                buckets[bucket_id] = (int *)realloc(buckets[bucket_id], bucket_size_limit[bucket_id] * sizeof(int));
                if (buckets[bucket_id] == NULL) {
                    printf("Memory allocation failed for bucket reallocation.\n");
                    // Handle error appropriately (e.g., exit or use alternative strategy)
                }
            }
        
         // Assign element to the bucket
        buckets[bucket_id][bucket_size_fill[bucket_id]] = array[i];
        bucket_size_fill[bucket_id]++;
    
        
    } 
    // Step 3: Sort individual buckets
    for (int i = 0; i < numBuckets; i++){
        int current_bucket_size = bucket_size_fill[i];
        //int* current_bucket_elements = buckets[i];
        
        quickSort(buckets[i],0,current_bucket_size-1);

    }  

    //Step 4: Merge  sorted buckets into one array
    int  index_array = 0;
    for (int i = 0; i < numBuckets; i++){
        for (int j = 0;  j < bucket_size_fill[i]; j++ ){
            array[index_array++]=buckets[i][j];
        }
        free(buckets[i]);
    }
    
}

void Bucket_Sort_4(int array[], int size, int numBuckets, int numThreads) {
    
    int** buckets = (int**)malloc(numBuckets * sizeof(int *));

    if (buckets == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
     // Initialize each bucket pointer to NULL (optional for clarity)
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = NULL;
    }
    
    int max = max_element(array, size);

    int min = min_element(array,size);

    // Calculate bucket size and elements per bucket
    
    int bucket_size = max / numBuckets;

    int elements_per_bucket = (size / numBuckets) + 1;
 
    //Allocate memory for bucket size limit and filled counts
    int bucket_size_limit[numBuckets];

    int bucket_size_fill[numBuckets];

    // Step 1: Initialize an array of empty buckets
    for (int i = 0; i < numBuckets; i++) {
        //buckets[i] = 0;
        bucket_size_limit[i] = elements_per_bucket;
        bucket_size_fill[i] =   0;
        buckets[i] = (int *)malloc(elements_per_bucket * sizeof(int));
        }
    // Step 2: Bucket assignment
    for (int i = 0; i < size; i++) {
        int bucket_id = array[i]/bucket_size;
        
        if (bucket_id >= numBuckets) {
            bucket_id = numBuckets - 1;
        }
        {
            if (bucket_size_fill[bucket_id] >= bucket_size_limit[bucket_id]) {
                bucket_size_limit[bucket_id] *= 1.1; // Double capacity if needed
                buckets[bucket_id] = (int *)realloc(buckets[bucket_id], bucket_size_limit[bucket_id] * sizeof(int));
                if (buckets[bucket_id] == NULL) {
                    printf("Memory allocation failed for bucket reallocation.\n");
                    // Handle error appropriately (e.g., exit or use alternative strategy)
                }
            }
        
         // Assign element to the bucket
        buckets[bucket_id][bucket_size_fill[bucket_id]] = array[i];
        bucket_size_fill[bucket_id]++;
        }
        
    }
    /*for (int i = 0; i < numBuckets; i++){
        int current_bucket_size = bucket_size_fill[i];
        printf(" Bucket %d : ", i);
        for ( int j = 0; j < current_bucket_size; j++){
            printf( " %d ", buckets[i][j]);
        }
        printf("\n");
    }
    */
    // Step 3: Sort individual buckets
    #pragma omp parallel for schedule(dynamic, 1) num_threads(numBuckets)
    for (int i = 0; i < numBuckets; i++){
        int current_bucket_size = bucket_size_fill[i];
        //int* current_bucket_elements = buckets[i];
        quickSort(buckets[i],0,current_bucket_size-1);
        //mergeSort(buckets[i],0,current_bucket_size-1);
        //bubbleSort(buckets[i],current_bucket_size);

    }  
    /*for (int i = 0; i < numBuckets; i++){
        int current_bucket_size = bucket_size_fill[i];
        printf(" Bucket %d : ", i);
        for ( int j = 0; j < current_bucket_size; j++){
            printf( " %d ", buckets[i][j]);
        }
        printf("\n");
    }
    */
    
    // Calculate  start and end index of elements in all buckets and allocate
    int start_index[numBuckets], end_index[numBuckets];
    for (int i = 0; i < numBuckets; i++) {
          start_index[i] = (i == 0) ? 0 : end_index[i - 1] + 1;
          end_index[i] = start_index[i] + bucket_size_fill[i] - 1; 
    }

    
    //Step 4: Merge  sorted buckets into one array
    #pragma omp parallel for schedule(dynamic, 1) num_threads(numBuckets)
    for (int i = 0;  i < numBuckets; i++) {
        for (int j = start_index[i]; j <= end_index[i]; j++) {
            array[j] = buckets[i][j - start_index[i]];
        }
    }
    free(buckets);

}
/*
// Below function are from the Individual project from VT 23
void Bucket_Sort_3(int array[], int size, int numBuckets) {
    int* buckets = (int*)malloc(numBuckets * sizeof(int));

    if (buckets == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
    int* copyArray = (int*)malloc(size * sizeof(int));
    memcpy(copyArray, array, size * sizeof(int));

    // Step 1: Initialize an array of empty buckets
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = 0;
    }

    // Step 2: Partition the original array into k sets
    int elementsPerThread = size / numBuckets;

    // Step 3: Put elements in buckets
    int max = max_element(array, size);
    int bucketRange = max / numBuckets;
    for (int i = 0; i < size; i++) {
        int findIndex = array[i] / bucketRange;
        if (findIndex > (numBuckets - 1))
            findIndex = numBuckets - 1;
        int bucketIndex = findIndex;
        buckets[bucketIndex]++;
    }
    //printf("Step 3 completed\n");
    // Step 4: Sort non-empty buckets
    int copyIndex = 0;
    for (int i = 0; i < numBuckets; i++) {
        if (buckets[i] > 0) {
            int bucketSize = buckets[i];
            int* bucketElements = (int*)malloc(bucketSize * sizeof(int));
            // Copy elements from the original array to the bucket
            int count = 0;

            for (int j = 0; j < size; j++) {
                int findIndex = array[j] / bucketRange;
                if (findIndex > (numBuckets - 1))
                    findIndex = numBuckets - 1;
                int bucketIndex = findIndex;
                if (bucketIndex == i) {
                    bucketElements[count++] = array[j];
                }
            }

            // Sort the bucket
            //bubbleSort(bucketElements, count);
            //mergeSort(bucketElements,0,count-1);
            quick_sort_old(bucketElements,0,count-1);
            // Copy the sorted bucket back to the copyArray
            for (int j = 0; j < bucketSize; j++) {
                copyArray[copyIndex++] = bucketElements[j];
            }

            free(bucketElements);
        }
    }

    free(buckets);

    // Replace the contents of the original array with the sorted copyArray
    for (int i = 0; i < size; i++) {
        array[i] = copyArray[i];
    }

    free(copyArray);
}
// Function to perform parallel bucket sort
void Bucket_Sort_4(int array[], int size, int numBuckets, int numThreads) {
    int* buckets = (int*)malloc(numBuckets * sizeof(int));

    if (buckets == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return;
    }
    int* copyArray = (int*)malloc(size * sizeof(int));
    memcpy(copyArray, array, size * sizeof(int));
    
    // Step 1: Initialize an array of empty buckets
    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = 0;
    } 
    
    // Step 2: Partition the original array into k sets
    int elementsPerThread = size / numThreads;
     
    // Step 3: In parallel, put elements in buckets
    #pragma omp parallel num_threads(numThreads)
    {
        int threadId = omp_get_thread_num();
        int startIndex = threadId * elementsPerThread;
        int endIndex = (threadId == numThreads - 1) ? size - 1 : startIndex + elementsPerThread - 1;
        int max = max_element(array, size);
        // Put elements in buckets
        for (int i = startIndex; i <= endIndex; i++) {
            int findIndex = array[i] /  (max/numBuckets ) ;
            if (findIndex > (numBuckets - 1))
                findIndex = numBuckets-1;
            int bucketIndex = findIndex;
            // Increment the count of elements in the bucket
            #pragma omp atomic
            buckets[bucketIndex]++;
        }
    }
    
    // Step 4: Partition buckets into k sets
    int bucketsPerThread = numBuckets / numThreads;
    
    // Step 5: In parallel, sort non-empty buckets
    #pragma omp parallel num_threads(numThreads) 
    {
        int threadId = omp_get_thread_num();
        int startIndex = threadId * bucketsPerThread;
        int endIndex = (threadId == numThreads - 1) ? numBuckets - 1 : startIndex + bucketsPerThread - 1;
        int max = max_element(array, size);
        // Sort non-empty buckets
        for (int i = startIndex; i <= endIndex; i++) {
            if (buckets[i] > 0) {
                int bucketSize = buckets[i];
                int* bucketElements = (int*)malloc(bucketSize * sizeof(int));
                // Copy elements from the original array to the bucket
                int count = 0;
                
                for (int j = 0; j < size; j++) {
                    int findIndex = array[j] /  (max/numBuckets ) ;
                    if (findIndex > (numBuckets - 1))
                        findIndex = numBuckets-1;
                    int bucketIndex = findIndex;
                    if (bucketIndex == i) {
                        bucketElements[count++] = array[j];
                        
                        
                    }
                }
               
                
                // Sort the bucket
                
                //bubbleSort(bucketElements, count);
                //mergeSort(bucketElements,0,count-1);
                quickSort(bucketElements,0,count-1);

                // Copy the sorted bucket back to the original array
                int arrayIndex = 0;
                for (int j = 0; j<i; j++) {
                    arrayIndex+=buckets[j];
                }
                
                for (int j = 0; j < bucketSize; j++) {
                    copyArray[arrayIndex++] = bucketElements[j];
                }

                free(bucketElements);
                
            }
        }
    }
    
    free(buckets);
    for (int i = 0; i < size; i++) {
        array[i] = copyArray[i];
    }
    free(copyArray);
}
*/
/* The main() begins */
int main(int argc, char *argv[]) {
    
    if (argc != 3)
    {
        printf("There are some missing parameters! Exiting the code...\n ArraySize nThreads");
        return -1;
    }

    int arraySize = atoi(argv[1]);

    int nThreads = atoi(argv[2]);


    // Dynamically allocate memory for the array
    int* uniform_array = (int*)malloc(arraySize * sizeof(int));
    int* expo_array = (int*)malloc(arraySize * sizeof(int));
    int* normal_array = (int*)malloc(arraySize * sizeof(int));

    int* uniform_array_cp = (int*)malloc(arraySize * sizeof(int));
    int* expo_array_cp = (int*)malloc(arraySize * sizeof(int));
    int* normal_array_cp = (int*)malloc(arraySize * sizeof(int));


    //int* array2 = (int*)malloc(arraySize * sizeof(int));
    if (uniform_array == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return 1;
    }
    int MAXLIMIT = 100000000;
    // Seed the random number generator
    srand(time(0));

    int flag = 1; //flag = 0 for bucketsize = arraysize , flag = 1 bucketsize= nbuckets

    // Generate random values for the array elements
    for (int i = 0; i < arraySize; i++) {
        
        uniform_array[i] = rand() % MAXLIMIT; // Adjust the range as needed
    }
    // Generate random values with exponential distribution for the array elements
    for (int i = 0; i < arraySize; i++) {
        double u = (double)rand() / RAND_MAX;
        double lambda = 0.0001; // Adjust the lambda parameter as needed
        expo_array[i] = fmod((-log(1 - u) / lambda), MAXLIMIT);
    }

    // Generate random values with normal distribution for the array elements
    for (int i = 0; i < arraySize; i++) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        normal_array[i] = fmod(abs(z * MAXLIMIT),MAXLIMIT); // Adjust the range as needed
}
    memcpy(uniform_array_cp, uniform_array, arraySize * sizeof(int));
    memcpy(expo_array_cp, expo_array, arraySize * sizeof(int));
    memcpy(normal_array_cp, normal_array, arraySize * sizeof(int));
    
    double start_time,timeTaken;
    
    if (flag == 0 ) {
        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_1(uniform_array, arraySize);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Uniform Distribution Serial : %f\n", timeTaken);


        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_1(expo_array, arraySize);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Exponential Distribution Serial : %f\n", timeTaken);



        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_1(normal_array, arraySize);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Normal Distribution Serial : %f\n", timeTaken);

        //printf("---------------------------------------------\n");

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_2(uniform_array_cp, arraySize,nThreads);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Uniform Distribution Parallel : %f\n", timeTaken);

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_2(expo_array_cp, arraySize,nThreads);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Exponential Distribution Parallel : %f\n", timeTaken);

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_2(normal_array_cp, arraySize,nThreads);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Normal Distribution Parallel : %f\n", timeTaken);

        //printf("---------------------------------------------\n");
    }
    else {
        int nBuckets = nThreads;
        //printArray(uniform_array,arraySize);   
        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_3(uniform_array, arraySize,nBuckets);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Uniform Distribution Serial with Buckets : %f\n", timeTaken);
        //printArray(uniform_array,arraySize);
        sort_check(uniform_array,arraySize);

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_3(expo_array, arraySize,nBuckets);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Exponential Distribution Serial with Buckets : %f\n", timeTaken);
        sort_check(expo_array,arraySize);
        
        start_time = omp_get_wtime();
        // Calling the bucket sort function
        Bucket_Sort_3(normal_array, arraySize,nBuckets);

        timeTaken = omp_get_wtime() - start_time;
        printf("Time Normal Distribution Serial with Buckets : %f\n", timeTaken);
        sort_check(normal_array,arraySize);

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        //printArray(uniform_array_cp,arraySize);
        Bucket_Sort_4(uniform_array_cp, arraySize,nBuckets,nThreads);
        //printArray(uniform_array_cp,arraySize);
        timeTaken = omp_get_wtime() - start_time;
        printf("Time Uniform Distribution Parallel with Buckets : %f\n", timeTaken);
        sort_check(uniform_array_cp,arraySize);
        
        start_time = omp_get_wtime();
        // Calling the bucket sort function
        //printArray(expo_array_cp,arraySize);
        Bucket_Sort_4(expo_array_cp, arraySize,nBuckets,nThreads);
        //printArray(expo_array_cp,arraySize);
        timeTaken = omp_get_wtime() - start_time;
        printf("Time Exponential Distribution Parallel with Buckets : %f\n", timeTaken);
        sort_check(expo_array_cp,arraySize);

        start_time = omp_get_wtime();
        // Calling the bucket sort function
        //printArray(normal_array_cp,arraySize);
        Bucket_Sort_4(normal_array_cp, arraySize,nBuckets,nThreads);
        //printArray(normal_array_cp,arraySize);
        timeTaken = omp_get_wtime() - start_time;
        printf("Time Normal Distribution Parallel with Buckets : %f\n", timeTaken);
        sort_check(normal_array_cp,arraySize);
        
    }
    // Free the dynamically allocated memory
    free(uniform_array);
    free(normal_array);
    free(expo_array);
    free(uniform_array_cp);
    free(normal_array_cp);
    free(expo_array_cp);

    printf("\n");
    return 0;
}

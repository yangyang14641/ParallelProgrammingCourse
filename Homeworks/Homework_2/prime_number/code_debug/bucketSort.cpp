// This program is designed to sort the bucket list for bucket sort algorithm

//head files
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> 
//#include <math.h>
#include <stdlib.h>




// function sortBucket
//--------------------------------------------------------------------------------------
/*long int *sortBucket(long int *bucket,long int *nElementInBucket,int nBucket
,int *indexBucket,int *nIndexBucket)*/

int *sortBucket(long int *bucket,long int *nElementInBucket,int nBucket
,int *nIndexBucket)
{  
   int i, j;
   long int bucketHead;
   long int temp_new, temp_old;
   int index_new, index_old;
   int *temp_int;
   long int *temp_longint;
   long int tempShift, tempShift_old, tempShift_new ;
   long int *listHead;
   int  *listHeadIndex;
   
   
       
   
/*   printf("Function sort Bucket:\n");
   tempShift = 0;
   for (i = 0; i < nBucket; i++)
   {  
      printf("Bucket id: %d,   Elements:\n",i);
       
      for (j = 0; j < *(nElementInBucket+i); j++)
        {                   
            printf("%ld  ", *(bucket + tempShift + j));
        }
           
     tempShift += *(nElementInBucket + i);
     printf("\n");
   }  
*/



// list head record 
   listHead = (long int *)malloc(sizeof(long int));
   listHeadIndex = (int *)malloc(sizeof(int));
   *nIndexBucket = 0;
   tempShift = 0;
   for (i = 0; i < nBucket; i++)
   {  

      bucketHead = *(bucket + tempShift);               // find the first element in every bucket.
//      printf("Pointer shift: %ld\nCurrent bucketHead: %ld\n",tempShift,bucketHead);
      if (bucketHead == 0 ) 
        {
            tempShift += *(nElementInBucket + i);
            continue;
        }
      
      *nIndexBucket += 1;
       
      temp_int = listHeadIndex;
      temp_longint = listHead;
      
      listHeadIndex = (int *)malloc(*nIndexBucket*sizeof(int));
      listHead = (long int *)malloc(*nIndexBucket*sizeof(long int));
      
      memcpy(listHeadIndex, temp_int, ((*nIndexBucket) - 1 )*sizeof(int));
      memcpy(listHead, temp_longint, ((*nIndexBucket) - 1 )*sizeof(long int));
     
      *(listHeadIndex + (*nIndexBucket) - 1 ) = i;
      *(listHead + (*nIndexBucket) - 1 ) = bucketHead;

      free(temp_int);
      free(temp_longint);
      
      tempShift += *(nElementInBucket + i);

   }
   
   printf("\n\n");
   printf("List Head has %d elements\n",*nIndexBucket);
   for (i = 0; i< *nIndexBucket;i++)
  {
     printf("Element value %ld and id is %d\n",*(listHead+i), *(listHeadIndex+i));
  }


//bubble sort algorithm O(n^2)
 
  for (i = 0; i < *nIndexBucket;i++)
  {

     temp_old = *(listHead+i);
     index_old = *(listHeadIndex+i); 
    
     for(j = i + 1; j < *nIndexBucket;j++)
     {
        temp_new = *(listHead+j);
        index_new = *(listHeadIndex+j); 
        if(temp_new < temp_old)     // swap operator
        {
           *(listHead+j) = temp_old;            // value old to new
           *(listHeadIndex+j) = index_old;      // index old to new
           *(listHead+i) = temp_new;            // value new to old
           *(listHeadIndex+i) = index_new;      // index new to old
        }
        temp_old = *(listHead+i);
        index_old = *(listHeadIndex+i);  
     }
  }
  
   printf("\n\n");
   printf("Sorted list Head has %d elements\n",*nIndexBucket);
   for (i = 0; i< *nIndexBucket;i++)
  {
     printf("Element value %ld and id is %d\n",*(listHead+i), *(listHeadIndex+i));
  }
  
      
  free(listHead);
  
  return listHeadIndex;
}
//--------------------------------------------------------------------------------------




// function sortAllElements
//--------------------------------------------------------------------------------------

long int *sortAllElements(long int *bucket,long int *nElementInBucket,int nBucket
,int *indexBucket,int *nIndexBucket)
{  
   int i,j;
   long int *temp;
   long int tempShift, bucketShift;
   int index_temp;
   
   long int nNonZero;
   

   nNonZero = 0;
   for (i=0; i < *nIndexBucket; i++)
   {
      nNonZero += *(nElementInBucket+*(indexBucket + i));
   }



   temp = bucket;
   bucket = (long int*)malloc((nNonZero) * sizeof(long int));
   
   bucketShift = 0;
// Bomb found in here 
   for (i=0; i < *nIndexBucket; i++)
   {
       index_temp = *(indexBucket + i);
       
       tempShift = 0;        
       for (j = 0; j < index_temp; j++)
         {
            tempShift += *(nElementInBucket + j);
         }

       memcpy((bucket + bucketShift), (temp + tempShift), (*(nElementInBucket + index_temp)) * sizeof(long int));
       bucketShift += *(nElementInBucket + index_temp);

      
   }



   printf("\n");
   printf("We have %d non-empty buckets.\n",*nIndexBucket);
   printf("Bucket sort index are showed below:\n");
   nNonZero = 0;
   for (i=0; i < *nIndexBucket; i++)
   {
      printf("%d\t", *(indexBucket + i));
      nNonZero += *(nElementInBucket+*(indexBucket + i));
   }
   printf("\n");



   printf("\n\nAfter sort, elements in buckets:\n");
   for (i=0; i < nNonZero; i++) printf("%ld  ",*(bucket+i));
     printf("\n\n\n\n");

  
  free(temp);

  return bucket;
}
//--------------------------------------------------------------------------------------







// function main
//--------------------------------------------------------------------------------------
int main(int argc, char* argv[]) 
{  
   int i,j;

   int nBucket;
   long int *nElementInBucket,nTotalElements;                            // This pointer is to show how many 
   long int *bucket;                                                     // This pointer point to  all elements in the bucket

   int *indexBucket;                                                     // this pointer point to buckets' index
   int *nIndexBucket;                                                    // counter to record how many non-empty buckets 
   
   long int temp,nNonZero;
       
   printf("Please input number of Bucket: \n"); 
   scanf("%d",&nBucket);                                                 // input number of bucket
   nElementInBucket = (long int*)malloc(nBucket*sizeof(long int));       // alloc memory for number of elements in each bucket
   printf("memory alloc complected for nBucket\n");
   
   printf("Please input number of elements in each bucket:\n");
   for (i = 0; i < nBucket; i++)
   {  
      printf("Bucket id: %d: ",i);
      scanf("%ld",nElementInBucket+i);                                   //scan for number of elements in each bucket
   }
   
   
   nTotalElements = 0;
   for (i = 0; i < nBucket; i++)
   {  
      printf("Bucket id: %d, number of elements: %ld \n",i,*(nElementInBucket+i));
      nTotalElements  +=  *(nElementInBucket+i);
   }
   printf("number of total elements: %ld \n",nTotalElements);
   
   
   bucket = (long int*)malloc(nTotalElements*sizeof(long int));
   printf("memory alloc complected for bucket\n");
    



// scan elements in each bucket

   temp = 0;
   for (i = 0; i < nBucket; i++)
   {  
      printf("please enter elements for bucket %d\n",i);
  
      for (j = 0; j < *(nElementInBucket+i); j++)
         {   
            //printf("temp = %ld , j = %d \n",temp,j);
            scanf("%ld", (bucket + temp + j));
            //printf("elements: %ld\n",*(bucket + temp + j));
         }
      printf("scan elements for bucket %d finished! \n",i);
      temp += *(nElementInBucket + i);              // pointer shift for the list head
   }
   
   printf("\n\n");

// print element in each bucket
   temp = 0;
   for (i = 0; i < nBucket; i++)
   {  
      printf("Bucket id: %d,   Elements:\n",i);
       
      for (j = 0; j < *(nElementInBucket+i); j++)
        {                   
            printf("%ld  ", *(bucket + temp + j));
        }
           
     temp += *(nElementInBucket + i);
     printf("\n");
   }   



printf("\n\n elements in buckets:\n");
for (i=0;i<nTotalElements;i++) printf("%ld  ",*(bucket+i));
   printf("\n\n\n\n");



// call function sort buckets 
   nIndexBucket = (int *)malloc(sizeof(int));
//   indexBucket = (int *)malloc(sizeof(int));

//   sortBucket(bucket,nElementInBucket,nBucket,indexBucket,nIndexBucket);
   indexBucket = sortBucket(bucket,nElementInBucket,nBucket,nIndexBucket);
   
   printf("\n");
   printf("We have %d non-empty buckets.\n",*nIndexBucket);
   printf("Bucket sort index are showed below:\n");
   nNonZero = 0;
   for (i=0; i < *nIndexBucket; i++)
   {
      printf("%d\t", *(indexBucket + i));
      nNonZero += *(nElementInBucket+*(indexBucket + i));
   }
   printf("\n");


// call function to sort all element in the bucket
bucket = sortAllElements(bucket,nElementInBucket,nBucket,indexBucket,nIndexBucket);


   
   printf("\n\nAfter sort, elements in buckets:\n");
   for (i=0;i<nNonZero;i++) printf("%ld  ",*(bucket+i));
     printf("\n\n\n\n");


}




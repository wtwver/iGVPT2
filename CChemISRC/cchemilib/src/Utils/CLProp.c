/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2010 Abdulrahman Allouche (University Lyon 1)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
********************************************************************************/

/* CLProp.c */
#include <stdio.h>
#include <stdlib.h>
#include "../Utils/Types.h"

#ifdef ENABLE_CL
#include <CL/cl.h>
static CLProp clProp={0,0};

/*****************************************************************************************/
CLProp getCLProp()
{
	return clProp;
}
/*****************************************************************************************/
void initCLProp()
{
	cl_context_properties properties[3];
	cl_int err;
	cl_uint num_of_platforms=0;
	cl_platform_id platform_id;
	cl_platform_id pid[2];
	cl_device_id device_id;
	cl_uint num_of_devices=0;
	cl_ulong long_entries;
        size_t p_size;
        cl_uint entries;
	char dname[500];
        cl_uint k;
	if(clProp.context!=0) return; /* already initialized */

	// retreive a list of platforms avaible
	if (clGetPlatformIDs(2, pid, &num_of_platforms)!= CL_SUCCESS)
	{
		printf("Unable to get CL platform_id\n");
		exit(1);
	}
	printf("num_of_platforms=%d\n",num_of_platforms);
#ifdef CPU_CL
	platform_id = pid[1];
	// try to get a supported CPU device
	if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &num_of_devices) != CL_SUCCESS)
	{
		printf("Unable to get CL device_id\n");
		exit(1);
	}
#else
	platform_id = pid[0];
	// try to get a supported GPU device
	if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_of_devices) != CL_SUCCESS)
	{
		printf("Unable to get CL device_id\n");
		exit(1);
	}
#endif

	/* printf("num_of_devices=%d\n",num_of_devices);*/
	clGetPlatformInfo(platform_id,CL_PLATFORM_NAME,500,dname,NULL);
	printf("CL_PLATFORM_NAME = %s\n", dname);
	clGetPlatformInfo(platform_id,CL_PLATFORM_VERSION,500,dname,NULL);
	printf("CL_PLATFORM_VERSION = %s\n", dname);
	{
		clGetDeviceInfo(device_id, CL_DEVICE_NAME, 500, dname,NULL);
		printf("Device name = %s\n", dname);
		clGetDeviceInfo(device_id,CL_DRIVER_VERSION, 500, dname,NULL);
		printf("\tDriver version = %s\n", dname);
		clGetDeviceInfo(device_id,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
		printf("\tGlobal Memory (MB):\t%llu\n",long_entries/1024/1024);
		clGetDeviceInfo(device_id,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
		printf("\tGlobal Memory Cache (MB):\t%llu\n",long_entries/1024/1024);
		clGetDeviceInfo(device_id,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
		printf("\tLocal Memory (KB):\t%llu\n",long_entries/1024);
		clGetDeviceInfo(device_id,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
		printf("\tMax clock (MHz) :\t%llu\n",long_entries);
		clGetDeviceInfo(device_id,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
		printf("\tMax Work Group Size:\t%d\n",p_size);
		clGetDeviceInfo(device_id,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&entries,NULL);
		printf("\tNumber of parallel compute cores:\t%d\n",entries);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR : \t%d\n",k);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT : \t%d\n",k);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_INT : \t%d\n",k);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG : \t%d\n",k);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT : \t%d\n",k);
		clGetDeviceInfo(device_id,CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),&k,NULL);
		printf("\tCL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE : \t%d\n",k);
	}



	// context properties list - must be terminated with 0
	properties[0]= CL_CONTEXT_PLATFORM;
	properties[1]= (cl_context_properties) platform_id;
	properties[2]= 0;

	// create a context with the GPU device
	clProp.context = clCreateContext(properties,1,&device_id,NULL,NULL,&err);

	// create command queue using the context and device
	clProp.command_queue = clCreateCommandQueue(clProp.context, device_id, 0, &err);
	clProp.device_id = device_id;
}
/*****************************************************************************************/
void freeCLProp()
{
	clReleaseCommandQueue(clProp.command_queue);
	clReleaseContext(clProp.context);
}
void printErrorCLRun(cl_int err)
{
switch(err)
{
case CL_INVALID_PROGRAM_EXECUTABLE : printf("there is no successfully built program executable available for device associated with command_queue\n"); break;
case CL_INVALID_COMMAND_QUEUE :  printf("command_queue is not a valid command-queue.\n");  break;
case CL_INVALID_KERNEL :  printf("kernel is not a valid kernel object.\n");  break;
case CL_INVALID_CONTEXT :  printf("context associated with command_queue and kernel is not the same or if the context associated with command_queue and events in event_wait_list are not the same.\n");  break;
case CL_INVALID_KERNEL_ARGS :  printf("the kernel argument values have not been specified.\n");  break;
case CL_INVALID_WORK_DIMENSION :  printf("work_dim is not a valid value (i.e. a value between 1 and 3).\n");  break;
case CL_INVALID_WORK_GROUP_SIZE :  printf("local_work_size is specified and number of work-items specified by global_work_size is not evenly divisable by size of work-group given by local_work_size or does not match the work-group size specified for kernel using the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier in program source.\n");  break;
printf("local_work_size is specified and the total number of work-items in the work-group computed as local_work_size[0] *... local_work_size[work_dim - 1] is greater than the value specified by CL_DEVICE_MAX_WORK_GROUP_SIZE in the table of OpenCL Device Queries for clGetDeviceInfo.\n");  break;
 printf("local_work_size is NULL and the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier is used to declare the work-group size for kernel in the program source.\n");  break;
case CL_INVALID_WORK_ITEM_SIZE :  printf("the number of work-items specified in any of local_work_size[0], ... local_work_size[work_dim - 1] is greater than the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0], .... CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim - 1].\n");  break;
case CL_INVALID_GLOBAL_OFFSET :  printf("global_work_offset is not NULL.\n");  break;
case CL_OUT_OF_RESOURCES :  printf("there is a failure to queue the execution instance of kernel on the command-queue because of insufficient resources needed to execute the kernel. For example, the explicitly specified local_work_size causes a failure to execute the kernel because of insufficient resources such as registers or local memory. Another example would be the number of read-only image args used in kernel exceed the CL_DEVICE_MAX_READ_IMAGE_ARGS value for device or the number of write-only image args used in kernel exceed the CL_DEVICE_MAX_WRITE_IMAGE_ARGS value for device or the number of samplers used in kernel exceed CL_DEVICE_MAX_SAMPLERS for device.\n");  break;
case CL_MEM_OBJECT_ALLOCATION_FAILURE :  printf("there is a failure to allocate memory for data store associated with image or buffer objects specified as arguments to kernel.\n");  break;
case CL_INVALID_EVENT_WAIT_LIST :  printf("event_wait_list is NULL and num_events_in_wait_list > 0, or event_wait_list is not NULL and num_events_in_wait_list is 0, or if event objects in event_wait_list are not valid events.\n");  break;
case CL_OUT_OF_HOST_MEMORY :  printf("there is a failure to allocate resources required by the OpenCL implementation on the host.\n"); break;
default:;
}
}
#endif




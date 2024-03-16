/********************************************************************************
 cchemi is an interface to ab initio computational chemistry programs 
 designed for add them many functionalities non available in these packages.
 Copyright (C) 2020 Abdulrahman Allouche (University Lyon 1)

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

#ifndef __CCHEMILIB_SLIST_H__
#define __CCHEMILIB_SLIST_H__

typedef struct _LNode  LNode;
typedef struct _SListClass  SListClass;
typedef struct _SList  SList;


struct _LNode {
        void* data;
        LNode* next;
};
struct _SList {
        LNode* head;
        LNode* tail;
        LNode* current;
        int size;
	SListClass* klass;
};
struct _SListClass
{
	void (*destroy)(SList* slist);
	int (*pushFront)(SList* slist, void* data);
	int (*pushBack)(SList* slist, void* data);
	int (*append)(SList* slist, void* data);
	void* (*popFront)(SList* slist);
	void* (*popBack)(SList* slist);
	int (*remove)(SList* slist, void* data);
	int (*step)(SList* slist);
	LNode* (*last)(SList* slist);
	void* (*readIndex)(SList* slist, int index);
	int (*insertAfter)(SList* slist, int index, void *data);
	void* (*extractAfter)(SList* slist, int index);
	
};
SList* newSList();

#endif

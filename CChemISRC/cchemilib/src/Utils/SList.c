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

#include <stdlib.h>
#include "SList.h"

static void slistDestroy(SList *SList);
static int slistPushFront(SList *SList, void *data);
static int slistPushBack(SList *SList, void *data);
static void* slistPopFront(SList *SList);
static void* slistPopBack(SList *SList);
static int slistStep(SList *SList);
static void* slistReadIndex(SList *SList, int index);
static int slistInsertAfter(SList *SList, int index, void *data);
static void* slistExtractAfter(SList *SList, int index);
static LNode* slistLast (SList *slist);
static int slistRemove(SList *slist, void* data);
/*************************************************************************************************************/
SList* newSList()
{
        SList *slist;
        slist = malloc(sizeof(SList));
        if (slist != NULL) {
                slist->head = NULL;
                slist->tail = NULL;
                slist->current = NULL;
                slist->size = 0;

                slist->klass = malloc(sizeof(SListClass));
                slist->klass->destroy = slistDestroy;
                slist->klass->pushFront =slistPushFront;
                slist->klass->pushBack = slistPushBack;
                slist->klass->append = slistPushBack;
                slist->klass->popFront = slistPopFront;
                slist->klass->popBack = slistPopBack;
                slist->klass->remove = slistRemove;
                slist->klass->step = slistStep;
                slist->klass->last = slistLast;
                slist->klass->readIndex = slistReadIndex;
                slist->klass->insertAfter = slistInsertAfter;
                slist->klass->extractAfter = slistExtractAfter;
        }
        return slist;
}
/*************************************************************************************************************/
static void slistDestroy(SList *slist)
{
        LNode *save_next;
        slist->current = slist->head;
        while(slist->current != NULL) {
                save_next = slist->current->next;
                free(slist->current->data);
                free(slist->current);
                slist->current = save_next;
        }
        free(slist);
}
/*************************************************************************************************************/
static int slistPushFront(SList *slist, void *data)
{
	LNode *node;
	node = malloc(sizeof(LNode));
	if (node == NULL) return -1;
	node->data = data;
	node->next = slist->head;
	slist->head = node;
	if (slist->size == 0) slist->tail = node;
	slist->size++;
	return 0;
}
/*************************************************************************************************************/
static int slistPushBack(SList *slist, void *data)
{
	LNode *node;
	node = malloc(sizeof(LNode));
	if (node ==NULL) return -1;
	node->data = data;
	node->next = NULL;
	if (slist->size == 0) { slist->head = node; slist->tail = node; } 
	else { slist->tail->next = node; slist->tail = node; }
	slist->size++;
	return 0;
}
/*************************************************************************************************************/
static void* slistPopFront(SList *slist)
{
	if (slist->size == 0) return NULL;
	void* data = slist->head->data;
	LNode *save_head = slist->head;
	if (slist->size == 1)
	{
		slist->head = NULL;
		slist->tail = NULL;
		slist->current = NULL;
	}
	else
	{
		slist->head = slist->head->next;
	}
	free(save_head);
	slist->size--;
	return data;
}
/*************************************************************************************************************/
static void* slistPopBack(SList *slist)
{
	if (slist->size == 0) return NULL;
	void *data = slist->tail->data;
	LNode *save_tail = slist->tail;
	if (slist->size == 1)
	{
		slist->head = NULL;
		slist->tail = NULL;
		slist->current = NULL;
	}
	else
	{
		LNode *new_tail = slist->head;
		while(new_tail->next->next != NULL) new_tail = new_tail->next;
		slist->tail = new_tail;
		slist->tail->next = NULL;
	}
	free(save_tail);
	slist->size--;
	return data;
}
/*************************************************************************************************************/
static int slistStep(SList *slist)
{
	if (slist->current == NULL) return 1;
        else
	{
		slist->current = slist->current->next;
		return 0;
	}
}
/*************************************************************************************************************/
static LNode* slistLast (SList *slist)
{
	if(!slist) return NULL;
	return slist->tail;
}
/*************************************************************************************************************/
static void* slistReadIndex(SList *slist, int index)
{
	LNode *target;
	int i;
	if ( ((slist->size - index - 1) < 0 ) || (index < 0) ) return NULL;
	target = slist->head;
	for(i = 0; i < index; i++) target = target->next;
	return (target->data);
}
/*************************************************************************************************************/
static int slistRemove(SList *slist, void* data)
{
	LNode *current;
	LNode *prev;
	if(!slist) return 1;
	if(!data) return 2;
	if(data==slist->head->data)
	{ 
		if (slist->size == 1)
		{
			slist->head = NULL;
			slist->tail = NULL;
			slist->current = NULL;
		}
		else
		{
			slist->head = slist->head->next;
		}
		slist->size--;
		return 0;
	}
	if(data==slist->tail->data)
	{ 
		if (slist->size == 1)
		{
			slist->head = NULL;
			slist->tail = NULL;
			slist->current = NULL;
		}
		else
		{
			LNode *new_tail = slist->head;
			while(new_tail->next->next != NULL) new_tail = new_tail->next;
			slist->tail = new_tail;
			slist->tail->next = NULL;
		}
		slist->size--;
		return 0;
	}
	
	prev = slist->head;
	current = prev->next;
        while(current)
        {
                void* cdata = current->data;
		if(cdata==data)
		{
			prev->next = current->next;
			slist->size--;
			return 0;
		}
		prev = current;
		current = prev->next;
	}
	return 1;
}
/*************************************************************************************************************/
static int slistInsertAfter(SList *slist, int index, void *data)
{
	LNode *target;
	LNode *node;
	int i;
	if ( ((slist->size - index - 1) < 0 ) || (index < 0) ) return 1;
	target = slist->head;
	for(i = 0; i < index; i++) target = target->next;
	node = malloc(sizeof(LNode));
	if (node == NULL) return -1;
	node->data = data;
	node->next = target->next;
	target->next = node;
	if (index == slist->size - 1) slist->tail = node;
	slist->size++;
	return 0;
}
/*************************************************************************************************************/
static void* slistExtractAfter(SList *slist, int index)
{
	LNode *target;
	void *data;
	LNode *saveObsolete;
	int i;
	if ( ((slist->size - index - 2) < 0 ) || (index < 0) ) return NULL;
	target = slist->head;
	for(i = 0; i < index; i++) target = target->next;
	if (index == slist->size - 1) slist->tail = target;
	data = target->next->data;
	saveObsolete = target->next;
	target->next = target->next->next;
	free(saveObsolete);
	slist->size--;
	return data;
}
/*************************************************************************************************************/

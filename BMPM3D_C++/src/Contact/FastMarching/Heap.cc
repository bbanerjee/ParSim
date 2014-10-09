/*
 * Heap.cc
 *
 *  Created on: 4/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/heap.cpp
 *
 *  The heap index (returned by the push() method) is a map from the grid
 *  address to the index of heap's internal data representation.
 *
 * d_distance is the unsigned distance from the zero level set to each
 * element in the heap.
 *
 * d_address is the (original) grid address of each element in the heap.
 *
 * d_backPointer is a map from the index of d_distance (or d_address) to the
 * current location of the element in the heap.
 *
 * d_heap is an array of integer indices into the d_distance (or d_address)
 * array the heap invariant is maintained by moving elements in this list.
 *
 * if d_selfTest is true a consistency check is done after each operation.
 *
 * Currently, the heap constructor needs to know the number of elements
 * that will enter the heap during the calculation.
 */

#include <Contact/FastMarching/Heap.h>

#include <iostream>
#include <limits>
#include <Exception.h>

using namespace BrMPM;

Heap::Heap(int maxLength, bool selfTest)
  : d_maxLength(maxLength), d_heapLength(0), d_listLength(0),
    d_distance(maxLength, 0.0), d_heap(maxLength, 0), d_address(maxLength, 0),
    d_backPointer(maxLength, 0), d_selfTest(selfTest)
{
}

Heap::~Heap()
{
}

int
Heap::push(int address, double value)
{
  if (d_heapLength == d_maxLength) {
    throw Exception("FastMarching contact algorithm::Heap push error::Heap is full, increase max_length.",
                     __FILE__, __LINE__);
  }
  d_heap[d_heapLength] = d_listLength;
  d_address[d_listLength] = address;
  d_distance[d_listLength] = value;
  d_backPointer[d_listLength] = d_heapLength;
  d_listLength++;
  d_heapLength++;

  siftDown(0, d_heapLength-1);

  if (d_selfTest) test();
  return d_listLength-1;
}

void
Heap::pop(int& address, double& value)
{
  if (d_heapLength == 0) {
    throw Exception("FastMarching contact algorithm::Heap pop error::Heap is empty.",
                    __FILE__, __LINE__);
  }
  int loc = d_heap[0];
  value = d_distance[loc];
  address = d_address[loc];
  d_heap[0] = d_heap[d_heapLength-1];
  d_backPointer[d_heap[0]] = 0;
  d_heapLength--;

  siftUp(0);

  if (d_selfTest) test();
}

void Heap::set(int index, double newDistance)
{
  double oldDistance = d_distance[index];
  int pos = d_backPointer[index];
  d_distance[index]=newDistance;
  if (newDistance > oldDistance)
  {
    siftUp(pos);
  }
  if (d_distance[d_heap[pos]] != newDistance)
  {
    if (d_selfTest) test();
    return;
  }
  siftDown(0,pos);
  if (d_selfTest) test();
}

bool Heap::empty() const
{
  if (d_heapLength == 0) return true;
  return false;
}

void Heap::test() const
{
  for (int i=0; i < d_heapLength; i++)
   {
     int c[2];
     c[0]=2*i+1;
     c[1]=c[0]+1;
     for (int j=0; j<2; j++)
     {
       if (c[j] < d_heapLength-1)
       {
         double dp = d_distance[d_heap[i]];
         double dc = d_distance[d_heap[c[j]]];
         if (! (dp<=dc))
         {
           throw Exception("heap invariant violation", __FILE__, __LINE__);
         }
       }
     }
   }
   for (int i=0; i<d_heapLength; i++)
   {
     if (! d_backPointer[d_heap[i]]==i)
     {
       std::cerr << "error " << i << std::endl;
       throw Exception("heap backpointer inconsistancy", __FILE__, __LINE__);
     }
   }

}

void
Heap::siftUp(int pos)
{
  int endPos = d_heapLength;
  int startPos = pos;
  int rightPos;
  int newItem = d_heap[pos];
  int childPos = 2*pos + 1;
  while (childPos < endPos)
  {
    rightPos = childPos + 1;
    if ((rightPos < endPos) && !
        (d_distance[d_heap[childPos]] < d_distance[d_heap[rightPos]]))
    {
      childPos = rightPos;
    }
    d_heap[pos]=d_heap[childPos];
    d_backPointer[d_heap[childPos]]=pos;
    pos = childPos;
    childPos = 2*pos + 1;
  }
  d_heap[pos] = newItem;
  siftDown(startPos, pos);
}

void
Heap::siftDown(int startPos, int pos)
{
  int parent;
  int parentPos;
  int newItem = d_heap[pos];
  while (pos > startPos)
  {
    parentPos = (pos-1)>>1;
    parent = d_heap[parentPos];
    if (d_distance[newItem] < d_distance[parent])
    {
      d_heap[pos]=parent;
      d_backPointer[parent]=pos;
      pos=parentPos;
      continue;
    }
    break;
  }
  d_heap[pos]=newItem;
  d_backPointer[newItem]=pos;
}

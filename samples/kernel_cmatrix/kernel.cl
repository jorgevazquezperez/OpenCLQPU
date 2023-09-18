/**
 * clIndexedQueue.cl
 * Indexed queue implementation in OpenCL
 * Copyright (C) 2013-2014 Roy Spliet, Delft University of Technology
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
 * USA
 */
#include "clIndexedQueue.h"

/* Lines get too long... */
#define uptr uintptr_t
#define POISON 1

struct mem_item {
	clIndexedQueue_item q;
	uint32_t filler;
	uintptr_t thisobj;
};

/**
 * clIndexedQueue_ptr2idx() - Convert a pointer and tag to an index value
 * @q: Indexed queue object
 * @ptr: Pointer to convert
 * This helper function produces a uint32_t variable with the following bit-layout:
 * 0     Poison
 * 21:1  Index (base + stride * (n-1)
 * 32:22 Tag
 *
 * Returns 0 on failure.
 */
uint32_t
clIndexedQueue_ptr2idx(clIndexedQueue __global *q, void __global *ptr)
{
	uptr idx = (uptr) ptr;

	/* Before the base isn't possible */
	if(idx < (uptr) q->base)
		return 0;

	idx -= (uptr) q->base;
	/* Does it align to a stride boundary nicely? */
	if(idx & ((1 << q->stride) - 1))
		return 0;

	idx >>= q->stride;
	idx++;
	/* Does the index still fit? */
	if(idx > ((1 << 20) - 1))
		return 0;

	return (uint32_t) idx << 1;
}

/**
 * clIndexedQueue_idx2ptr() - Convert a pointer and tag to an index value
 * @q: Indexed queue object
 * @idx: Index to convert
 */
inline void __global*
clIndexedQueue_idx2ptr(clIndexedQueue __global *q, uint32_t idx)
{
	size_t i = idx;
	idx >>= 1;
	idx &= 0xfffff;
	if(idx == 0)
		return 0;

	idx--;

	i = idx;
	i <<= q->stride;

	return &q->base[i];
}

#define PTR(i,t) ((i & 0x1fffe) | ((((t)+1) & 0x7ff) << 22))
#define TAG(i) ((i >> 22) & 0x7ff)
#define IDX(i) ((i & 0x1fffe) >> 1)

__kernel void
clIndexedQueue_init(void __global *queue, void __global *base,
		uint32_t stride_l2, void __global *i)
{
	clIndexedQueue __global *q = (clIndexedQueue __global *) queue;
	uint32_t idx, tag;
	clIndexedQueue_item __global *item = (clIndexedQueue_item __global *)i;

	q->base = base;
	q->stride = stride_l2;

	tag = TAG(item->next);
	item->next = PTR(0, tag);
	idx = clIndexedQueue_ptr2idx(q, item);

	q->head = PTR(idx, tag);
	q->tail = PTR(idx, tag);
	mem_fence(CLK_GLOBAL_MEM_FENCE);
}

/**
 * idxd_enqueue() - Add an item to the indexed queue
 * @q: Queue to add the item to
 * @item: Item to add to the queue
 * @return 1 iff enqueuing succeeded, 0 otherwise
 */
int
idxd_enqueue(clIndexedQueue __global *q, clIndexedQueue_item __global *item)
{
	clIndexedQueue_item __global *tail;
	unsigned int ret = 0;
	uint32_t idx, tag,
	tailidx,
	nextidx;

	if(item == NULL)
		return 0;

	tag = TAG(item->next);
	item->next = PTR(0, tag-1);
	tag++;
	idx = clIndexedQueue_ptr2idx(q, item);

	while(1) {
		tailidx = q->tail;
		mem_fence(CLK_GLOBAL_MEM_FENCE);

		tail = (clIndexedQueue_item __global *)
				clIndexedQueue_idx2ptr(q, tailidx);
		nextidx = tail->next;
		mem_fence(CLK_GLOBAL_MEM_FENCE);

		/* Did I read a consistent state? */
		if(q->tail == tailidx) {
			if(IDX(nextidx) == 0) {
				tag = TAG(nextidx);
				if(atom_cmpxchg(&tail->next, nextidx, PTR(idx, tag)) == nextidx) {
					mem_fence(CLK_GLOBAL_MEM_FENCE);
					ret = 1;
					break;
				}
			} else {
				tag = TAG(tailidx);
				atom_cmpxchg(&q->tail, tailidx, PTR(nextidx, tag));
				mem_fence(CLK_GLOBAL_MEM_FENCE);
			}
		}
	}
	tag = TAG(tailidx);
	atom_cmpxchg(&q->tail, tailidx, PTR(idx, tag));
	mem_fence(CLK_GLOBAL_MEM_FENCE);

	return ret;
}

/**
 * idxd_dequeue() - Remove and return the next item in the queue
 * @q: Queue to get the item from
 * @return The next queue item, NULL on failure.
 */
clIndexedQueue_item __global *
idxd_dequeue(clIndexedQueue __global *q)
{
	clIndexedQueue_item __global *head;
	uint32_t 	tag,
			nextidx,
			tailidx,
			headidx;

	while(1) {
		headidx = q->head;
		mem_fence(CLK_GLOBAL_MEM_FENCE);
		head = (clIndexedQueue_item __global *)
				clIndexedQueue_idx2ptr(q, headidx);
		tailidx = q->tail;
		nextidx = head->next;
		mem_fence(CLK_GLOBAL_MEM_FENCE);

		if(headidx == q->head) {
			if(IDX(headidx) == IDX(tailidx)) {
				if(IDX(nextidx) == 0) {
					return NULL;
				}
				tag = TAG(tailidx);
				atom_cmpxchg(&q->tail, tailidx, PTR(nextidx, tag));
				mem_fence(CLK_GLOBAL_MEM_FENCE);
			} else {
				tag = TAG(headidx);
				if(atom_cmpxchg(&q->head, headidx, PTR(nextidx, tag)) == headidx)
					break;
			}
		}

	}
	mem_fence(CLK_GLOBAL_MEM_FENCE);
	return head;
}

/* Test enqueueing. */
__kernel void
clIndexedQueue_test_enqueue(void __global *queue, unsigned int __global *mem)
{
	clIndexedQueue __global *q = (clIndexedQueue __global *) queue;
	size_t pid = 0, j;
	unsigned int i;
	struct mem_item __global *item, *prev;

	/* First find global unique ID */
	for(i = 0, j = 1; i < get_work_dim(); i++) {
		pid += j * get_global_id(i);
		j *= get_global_size(i);
	}
	pid++;

	item = (struct mem_item __global *) &mem[pid * 4];

	idxd_enqueue(q, &item->q);

	if(pid == 1) {
		prev = (struct mem_item __global *) &mem[0];
		prev->thisobj = (uintptr_t) &mem[0];
	}
	item->thisobj = (uintptr_t) item;
}

/* Test enqueueing, dequeueing. */
__kernel void
clIndexedQueue_test_dequeue(void __global *queue, unsigned int __global *mem)
{
	clIndexedQueue __global *q = (clIndexedQueue __global *) queue;
	size_t pid = 0, j;
	unsigned int i;
	struct mem_item __global *item, *prev;

	/* First find global unique ID */
	for(i = 0, j = 1; i < get_work_dim(); i++) {
		pid += j * get_global_id(i);
		j *= get_global_size(i);
	}
	pid++;

	item = (struct mem_item __global *) &mem[pid * 4];
	item->thisobj = (uintptr_t) item;
	if(pid == 1) {
		prev = (struct mem_item __global *) &mem[0];
		prev->thisobj = (uintptr_t) &mem[0];
	}
	idxd_enqueue(q, (clIndexedQueue_item __global *) item);

	/* Do the shuffle */
	for(i = 0; i < 10; i++) {
		item = (struct mem_item __global *) idxd_dequeue(q);
		if(item != NULL) {
			if(!idxd_enqueue(q, (clIndexedQueue_item __global *) item)) {
			}
		} else {
		}
	}
}


/**
 * kma.cl
 * Main Kernel Memory Allocator OpenCL implementation
 * Copyright (C) 2013-2014 Roy Spliet, Delft University of Technology
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
 * USA
 */
#include "kma.h"

/** Initialise the heap
 * @param heap The heap to initialise
 * Sets the pointers to null, initialises the free list with all pages
 */
__kernel void
clheap_init(void __global *hp)
{
	__global struct clheap *heap = (__global struct clheap *)hp;
	struct kma_sb __global *sb;
	char __global *ptr;
	unsigned int pages;
	unsigned int i;

	/* Empty the superblock hashtable */
	for(i = 0; i < KMA_SB_SIZE_BUCKETS; i++) {
		heap->sb[i] = NULL;
	}

	/* Add all pages to the free list and initialise them */
	pages = (heap->bytes >> KMA_SB_SIZE_LOG2) - 1;
	ptr = (char __global *)heap;
	ptr += sizeof(struct clheap);
	sb = (struct clSuperBlock __global *)ptr;

	/* Initialise the free-list */
	clIndexedQueue_init(&heap->free, &sb[0], KMA_SB_SIZE_LOG2, &sb[0]);
	for(i = 1; i < pages; i++) {
		idxd_enqueue(&heap->free, &sb[i].q);
	}
}

/* For given sbid, return the size of a block in bytes */
size_t
_kma_size_by_sbid(int block)
{
	size_t size, total, bytes;
	unsigned int items, alloc_bytes;

	if(block > KMA_SB_SIZE_BUCKETS || block < 0)
		return 0;

	/* Blocks of 1 and 2 bytes don't exist for this malloc */
	block += 2;

	/* A bit of heuristic: For smaller blocks (sqrt(KMA_SB_SIZE))
	 * sacrifice one or two at the end of the superblock. */
	if(block <= (KMA_SB_SIZE_LOG2 >> 1)) {
		return 1 << (block);
	}

	/* For bigger blocks, adjust size to waste as little as possible
	 * 1. Calculate number of items when size is 2^(block+2)
	 * 2. Find total amount of space left when subtracting "alloc bits"
	 * 3. Divide this by the number of blocks*/
	bytes = 1 << block;
	items = KMA_SB_SIZE / bytes;

	/* Bytes rounded up to the nearest doubleword */
	alloc_bytes = items >> 3;
	alloc_bytes += (items & 0x7) ? 1 : 0;
	alloc_bytes += (alloc_bytes & 0x3) ? (4 - (alloc_bytes & 0x3)) : 0;

	/* Total space */
	total = KMA_SB_SIZE - (12 + alloc_bytes);
	size = total / items;

	return size & ~0x3;
}

/** For given size, return the superblock ID in the array
 * @param size size of the desired block
 * Returns the desired superblock id, -1 if no such superblock exists */
int
_kma_slots_by_size(size_t size)
{
	int slots;
	size_t space;

	space = KMA_SB_SIZE - 12;

	/* First approx */
	slots = space / size;
	/* Find out how many doublewords we need at the end of the block */
	if(slots & 0x1f)
		slots += 32;
	slots >>= 3;
	slots &= ~0x3;
	space -= slots;

	/* Recalculate with new space */
	return space / size;
}

/** For given size, return the superblock ID in the array
 * @param size size of the desired block
 * Returns the desired superblock id, -1 if no such superblock exists */
int
_kma_sbid_by_size(size_t size)
{
	int sbid, i;
	int test;
	size_t ssize;

	if(size > KMA_SB_SIZE - 16)
		return -1;

	ssize = size >> 2;

	/* First find a good approximation */
	for(i = 0; i < KMA_SB_SIZE_BUCKETS; i++, ssize >>= 1) {
		if(ssize & 1)
			sbid = i;
	}

	/* When block size increases, this no longer holds */
	while(true) {
		test = _kma_size_by_sbid(sbid);
		if(size > test && test > 0) {
			sbid++;
		} else if (size < _kma_size_by_sbid(sbid - 1)) {
			sbid--;
		} else {
			return sbid;
		}
	}
}

struct clSuperBlock __global *
_kma_reserve_block(__global struct clheap *heap, int block,
		unsigned int *slot)
{
	volatile struct kma_sb __global *cursor;
	unsigned int state, state_old, slots, i;
	volatile unsigned int __global *abits_ptr;

	if(block < 0 || block >= KMA_SB_SIZE_BUCKETS)
		return NULL;

	while(1) {
		/* Is there a superblock available */
		cursor = (struct kma_sb __global *)
				atom_cmpxchg((volatile uintptr_t __global *) &heap->sb[block], 0, POISON);
		mem_fence(CLK_GLOBAL_MEM_FENCE);
		if(cursor == 0) {
			/* No, let's reserve one */
			cursor = (struct kma_sb __global *)
					idxd_dequeue(&heap->free);
			if(!cursor) {
				/* No free pages left, return NULL */
				atom_xchg((volatile uintptr_t __global *) &heap->sb[block], 0);
				mem_fence(CLK_GLOBAL_MEM_FENCE);
				return NULL;
			}

			/* Reserve one for me */
			cursor->size = _kma_size_by_sbid(block);
			slots = _kma_slots_by_size(cursor->size);
			state = (slots << 16) | (slots - 1);
			cursor->state = state;

			/* Set all allocation bits to 0 (unallocated) */
			abits_ptr = (unsigned int __global *)cursor + (KMA_SB_SIZE >> 2);
			for(i = 0; i < slots; i += 32) {
				abits_ptr--;
				*abits_ptr = 0;
			}
			mem_fence(CLK_GLOBAL_MEM_FENCE);

			atom_cmpxchg((volatile uintptr_t __global *) &heap->sb[block], POISON, (uintptr_t) cursor);
			*slot = slots - 1;
			mem_fence(CLK_GLOBAL_MEM_FENCE);
			return cursor;
		}

		if((uintptr_t) cursor > POISON) {
			/* First reserve a slot */
			state_old = atom_add(&cursor->state, 0);
			mem_fence(CLK_GLOBAL_MEM_FENCE);
			state = state_old & 0xffff;
			slots = (state_old & 0xffff0000) >> 16;
			if(state == 0)
				continue;

			/* Decrease counter by 1 */
			state--;
			*slot = state;
			state |= (state_old & 0xffff0000);

			if(atom_cmpxchg(&cursor->state, state_old, state) != state_old)
				continue;
			mem_fence(CLK_GLOBAL_MEM_FENCE);

			/* If this was the last block in the SB, unlink */
			if((state & 0xffff) == 0) {
				atom_xchg((__global volatile uintptr_t *)&heap->sb[block], 0);
				mem_fence(CLK_GLOBAL_MEM_FENCE);
			}
			//heap->sb[8] += 1;
			return cursor;

		}
	}
}

/*
 * Return a pointer to a free block
 * @param heap Heap to allocate from
 * @pre Block has been reserved in state
 */
void __global *
_kma_get_block(struct kma_sb __global *sb, unsigned int slot)
{
	unsigned int abits = 0;
	volatile unsigned int __global *abits_ptr;
	unsigned int slots;
	uintptr_t ptr;
	unsigned int slot_orig = slot;

	slots = ((sb->state & 0xffff0000) >> 16);

	while(true) {
		abits_ptr = (volatile unsigned int __global *)sb;
		abits_ptr += (KMA_SB_SIZE >> 2);
		abits_ptr -= (slot >> 5);
		if(slot & 0x1f)
			abits_ptr--;
		if(slot)
			abits = *abits_ptr;
		abits >>= (slot & 0x1f);
		for(; slot < slots; slot++, abits >>= 1) {
			if((slot & 0x1f) == 0) {
				abits_ptr--;
				mem_fence(CLK_GLOBAL_MEM_FENCE);
				abits = *abits_ptr;
			}

			if((abits & 0x1) == 0) {
				/* Try setting the bit */
				if((atom_or(abits_ptr, (1 << (slot & 0x1f))) & (1 << (slot & 0x1f))) == 0) {
					mem_fence(CLK_GLOBAL_MEM_FENCE);
					/* Gotcha, I have block i */
					ptr = (uintptr_t) &sb->data;
					ptr += (slot * sb->size);
					*(unsigned int __global *)ptr = slot_orig;
					return (void __global *)ptr;
				}
			}
		}
		slot = 0;
	}
}

/** Allocate memory
 * @param heap Heap
 * @param size Size of the desired block
 */
void __global *
malloc(__global struct clheap *heap, size_t size)
{
	int block;
	unsigned int slot,i;
	struct kma_sb __global *sb;
	/* Sizes all come in log2 for now. This means a lot of wastage for
	 * medium-sized memory blocks.
	 * Earlier experiments showed that traversing a linked list could lead
	 * to a corrupted cursor, with unpredictable behaviour. We can improve
	 * by increasing the granularity and adding more size buckets, at the
	 * cost of possibly more internal fragmentation
	 *
	 * Let's find a suitable superblock */
	block = _kma_sbid_by_size(size);
	if(block < 0)
		return NULL;

	sb = _kma_reserve_block(heap, block, &slot);
	if(!sb) {
		return NULL;
	}

	return _kma_get_block(sb, slot);
}

void
free(__global struct clheap *heap, uintptr_t block)
{
	unsigned int size, mask;
	volatile struct kma_sb __global *sb;
	uintptr_t first_sb, off;
	unsigned int state_old, state, slots, sbid;
	bool enq;
	volatile unsigned int __global *abits_ptr;
	volatile unsigned int __global *b = (volatile unsigned int __global *) block;

	if(block == NULL)
		return;

	/* Find superblock */
	first_sb = ((uintptr_t) heap + sizeof(struct clheap));
	off = block - first_sb;
	mask = (1 << KMA_SB_SIZE_LOG2) - 1;

	/* Find size of block */
	sb = (volatile struct clSuperBlock __global *)(first_sb + (off & ~mask));
	size = sb->size;

	/* Index of this block */
	block -= (uintptr_t) sb;
	block -= (8 + sizeof(clIndexedQueue_item));
	block /= size;

	/* Update the "taken" bit
	 * XXX: If you try to free a block that isn't taken, "free slots"
	 * does get incremented. Corrupting the state! */
	abits_ptr = (volatile unsigned int __global *)(sb+1);
	abits_ptr -= ((block >> 5) + 1);
	*b = atom_and(abits_ptr, ~(1 << (block & 0x1f)));

	/* Update free slots */
	do {
		mem_fence(CLK_GLOBAL_MEM_FENCE);
		state_old = atom_add(&sb->state, 0);

		state = state_old & 0xffff;
		slots = (state_old & 0xffff0000) >> 16;
		state++;

		/* Enqueue this superblock and "unlink" */
		if(state == slots) {
			enq = 1;
			state = 0;
		} else {
			enq = 0;
		}

		state |= (slots << 16);
		mem_fence(CLK_GLOBAL_MEM_FENCE);
	} while (atom_cmpxchg(&sb->state, state_old, state) != state_old);
	mem_fence(CLK_GLOBAL_MEM_FENCE);

	/* find the right sbid and enqueue superblock if required */
	if(enq) {
		sbid = _kma_sbid_by_size(sb->size);
		atom_cmpxchg((volatile uintptr_t __global *)&heap->sb[sbid], (uintptr_t) sb, 0);
		mem_fence(CLK_GLOBAL_MEM_FENCE);
		idxd_enqueue(&heap->free, &sb->q);
	} else {
		/* Try to re-attach to avoid wasting too much mem */
		//atom_cmpxchg((volatile uintptr_t __global *)&heap->sb[sbid], NULL, (uintptr_t) sb);
		//mem_fence(CLK_GLOBAL_MEM_FENCE);
	}
}

//------------------------------------------------------------------------------------------------

typedef struct {
   double real;
   double imag;
} Complex;

Complex double_to_complex(double real, double imag)
{
    Complex c;

    c.real = real;
    c.imag = imag;

    return c;
}

double real(Complex c) 
{
   return c.real;
}

double imag(Complex c) 
{
   return c.imag;
}

Complex conj(Complex c)
{
    
    c.imag = -c.imag;
    
    return c;
}

Complex add(Complex c1, Complex c2)
{
    Complex c3;

    c3.real = c1.real + c2.real;
    c3.imag = c1.imag + c2.imag;

    return c3;
}

Complex sub(Complex c1, Complex c2)
{
    Complex c3;

    c3.real = c1.real - c2.real;
    c3.imag = c1.imag - c2.imag;

    return c3;
}
Complex mul(Complex c1, Complex c2)
{
    Complex c3;

    c3.real = (c1.real * c2.real) - (c1.imag * c2.imag);
    c3.imag = (c1.real * c2.imag) + (c1.imag * c2.real);

    return c3;
}

Complex div(Complex c1, Complex c2)
{
    Complex c3;

    c3.real = ((c1.real * c2.real) + (c1.imag * c2.imag)) / (pow(c2.real, 2) + pow(c2.imag, 2));
    c3.imag = -1 * (((c1.real * c2.imag) - (c1.imag * c2.real)) / (pow(c2.real, 2) + pow(c2.imag, 2)));

    if(c3.imag == 0)
        c3.imag = 0;

    return c3;
}

double abs(Complex c)
{
    return sqrt(c.real * c.real + c.imag * c.imag);
}

void printComplex(Complex complex)
{
    printf("Complejo: (%f, %f) \n", complex.real, complex.imag);
}

//--------------------------------------------------------------------------------------------------
//------------------------------------ Matrix4cd y Vector4cd----------------------------------------
typedef Complex Matrix4cd[4][4];
typedef Complex Vector4cd[4];

void printMatrix4cd(Matrix4cd matrix)
{
    int i, j;
    printf("MATRIZ 4x4: \n");
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            printf("(%f, %f) ", matrix[i][j].real, matrix[i][j].imag);
        }
        printf("\n");
    }
}

void printVector4cd(Vector4cd vector)
{
    int i;
    printf("Vector: \n");
    for (i = 0; i < 4; i++) {
        printf("(%f, %f) ", vector[i].real, vector[i].imag);
    }
    printf("\n");
}

void initialise_vector4cd(Vector4cd vector, size_t length)
{
	int i;
    for (i = 0; i < length; i++) {
		vector[i] = double_to_complex(0, 0);
    }
}

void mul4cd(Matrix4cd matrix, Vector4cd vector, Vector4cd result)
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            result[i] = add(result[i], mul(matrix[i][j], vector[j]));
        }
    }
}

//----------------------------------------------------------------------------------------
//------------------------------MatrixXcd y VectorXcd-------------------------------------
typedef Complex __global *MatrixXcd;
typedef Complex __global *VectorXcd;

void printMatrixXcd(MatrixXcd matrix, size_t rows, size_t cols)
{
    int i, j;
    printf("MATRIZ XxX: \n");
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("(%f, %f) ", matrix[i*cols + j].real, matrix[i*cols + j].imag);
        }
        printf("\n");
    }
}

void printVectorXcd(VectorXcd vector, size_t length)
{
    int i;
    printf("Vector X: \n");
    for (i = 0; i < length; i++) {
        printf("(%f, %f) ", vector[i].real, vector[i].imag);
    }
	printf("\n");
}

MatrixXcd initialise_matrixXcd(struct clheap __global *heap, size_t rows, size_t cols)
{
    MatrixXcd matrix = (MatrixXcd)malloc(heap, sizeof(Complex)*rows*cols);
	int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            matrix[i*cols + j] = double_to_complex(1, 0);
        }
    }

	return matrix;
}

VectorXcd initialise_vectorXcd(struct clheap __global *heap, size_t length)
{
    VectorXcd vector = (VectorXcd)malloc(heap, sizeof(Complex)*length);
	int i;
    for (i = 0; i < length; i++) {
		vector[i] = double_to_complex(0, 0);
    }

	return vector;
}

VectorXcd mulXcd(struct clheap __global *heap, MatrixXcd matrix, VectorXcd vector, size_t rows, size_t cols)
{
	VectorXcd result = initialise_vectorXcd(heap, rows);

	int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            result[i] = add(result[i], mul(matrix[i*cols + j], vector[j]));
        }
    }
	
	return result;
}

//----------------------------------------------------------------------------------------
//----------------------------------Sparse matrix-----------------------------------------
//FALTA

__kernel void ComplexMatrix(struct clheap __global *heap) 
{
    /* TESTS Matrix4cd y Vector4cd
    Matrix4cd matriz = {
        {{1, 2}, {2, 2.3}, {3, 4.5}, {0, 1}},
        {{1, 2}, {2, 2.3}, {3, 4.5}, {0, 1}},
        {{1, 2}, {2, 2.3}, {3, 4.5}, {0, 2}},
        {{1, 2}, {2, 2.3}, {3, 4.5}, {0, 1}}
    };

    Vector4cd vector = {{1, 0}, {1, 0}, {1, 0}, {1, 0}};
    Vector4cd result = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
    
    mul4cd(matriz, vector, result);
    printVector4cd(result);
    */

	/* TESTS MatrixXcd y VectorXcd
	size_t rows = 3, cols = 3;
    MatrixXcd matrix = initialise_matrixXcd(heap, rows, cols);
	printMatrixXcd(matrix, rows, cols);
	
	size_t length = 3;
	VectorXcd vector = initialise_vectorXcd(heap, length);
	vector[0] = double_to_complex(1, 0);

	VectorXcd result = mulXcd(heap, matrix, vector, rows, cols);
	printVectorXcd(result, rows);
	*/
}
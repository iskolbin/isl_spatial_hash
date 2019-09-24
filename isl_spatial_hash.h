#ifndef ISL_SPATIAL_HASH_H_
#define ISL_SPATIAL_HASH_H_
/* 
 isl_spatial_hash.h - v0.1 - public domain spatial hash
  
 This is as single-header-file library that provides ease-to-use
 spatial hash structure for fast broadphase collision detection

 To use this library, do this in *one* C or C++ file:
   #define ISL_SPATIAL_HASH_IMPLEMENTATAION	 
   #include "isl_spatial_hash.h"

 Dynamic array implementation taken from stb: https://github.com/nothings/stb/blob/master/stb_ds.h

 author: Ilya Kolbin (iskolbin@gmail.com)
 url: github.com/iskolbin/isl_spatial_hash 

 LICENSE
 
 See end of file for license information.
*/


#ifndef ISLSH_FREE
#define ISLSH_FREE(ctx,ptr) free(ptr)
#endif

#ifndef ISLSH_REALLOC
#define ISLSH_REALLOC(ctx,ptr,size) realloc((ptr),(size))
#endif

#ifndef ISLSH_STATIC
#define ISLSH_API extern
#else
#define ISLSH_API static
#endif

struct islsh_bucket {
	int x;
	int y;
	int *ids;
};

struct islsh_xywh {
	float x;
	float y;
	float w;
	float h;
};

struct islsh_spatial_hash {
	float cell_size;
	struct islsh_bucket **buckets;
	int spatial_buckets;
	struct islsh_xywh *xywhs;
	int **intersections;
	int id_counter;
	int *released_ids;
};

ISLSH_API struct islsh_spatial_hash *islsh_create(float cell_size, int spatial_buckets);
ISLSH_API void islsh_destroy(struct islsh_spatial_hash *self);
ISLSH_API int islsh_put(struct islsh_spatial_hash *self, float x, float y, float w, float h);
ISLSH_API int islsh_del(struct islsh_spatial_hash *self, int id);
#define islsh_contains(self,id) (((self)->id_counter>(id))&&((self)->xywhs[(id)].x==(self)->xywhs[(id)].x))
#define islsh_update(self,id,x,y,w,h) {if (islsh_del((self),(id)) == 0) islsh_put((self),(x),(y),(w),(h));}

ISLSH_API int islsh_intersections_len(struct islsh_spatial_hash *self, int id);
ISLSH_API int islsh_bucket_len(struct islsh_spatial_hash *self, int x, int y);

#endif

#ifdef ISL_SPATIAL_HASH_IMPLEMENTATION
#ifndef ISL_SPATIAL_HASH_IMPLEMENTATION_ONCE
#define ISL_SPATIAL_HASH_IMPLEMENTATION_ONCE

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

// Based on stb_ds.h version 0.62 (see https://github.com/nothings/stb/blob/master/stb_ds.h)
struct islsh_stbds_array_header {
  size_t length;
  size_t capacity;
};

#define islsh_stbds_header(t)  ((struct islsh_stbds_array_header *) (t) - 1)
#define islsh_stbds_arrgrow(a,b,c)   ((a) = islsh_stbds_arrgrowf((a), sizeof *(a), (b), (c)))
#define islsh_stbds_arrmaybegrow(a,n)  ((!(a) || islsh_stbds_header(a)->length + (n) > islsh_stbds_header(a)->capacity) ? (islsh_stbds_arrgrow(a,n,0),0) : 0)
#define islsh_stbds_arrput(a,v)     (islsh_stbds_arrmaybegrow(a,1), (a)[islsh_stbds_header(a)->length++] = (v))
#define islsh_stbds_arrpop(a)       (islsh_stbds_header(a)->length--, (a)[islsh_stbds_header(a)->length])
#define islsh_stbds_arrfree(a)      ((void) ((a) ? ISLSH_FREE(NULL,islsh_stbds_header(a)) : (void)0), (a)=NULL)
#define islsh_stbds_arrlen(a)       ((a) ? (ptrdiff_t) islsh_stbds_header(a)->length : 0)
#define islsh_stbds_arrcap(a)       ((a) ? islsh_stbds_header(a)->capacity : 0)
#define islsh_stbds_arrdeln(a,i,n)  (memmove(&(a)[i], &(a)[(i)+(n)], sizeof *(a) * (islsh_stbds_header(a)->length-(n)-(i))), islsh_stbds_header(a)->length -= (n))
#define islsh_stbds_arrdel(a,i)     islsh_stbds_arrdeln(a,i,1)

static void *islsh_stbds_arrgrowf(void *a, size_t elemsize, size_t addlen, size_t min_cap) {
  void *b;
  size_t min_len = islsh_stbds_arrlen(a) + addlen;

  // compute the minimum capacity needed
  if (min_len > min_cap)
    min_cap = min_len;

  if (min_cap <= islsh_stbds_arrcap(a))
    return a;

  // increase needed capacity to guarantee O(1) amortized
  if (min_cap < 2 * islsh_stbds_arrcap(a))
    min_cap = 2 * islsh_stbds_arrcap(a);
  else if (min_cap < 4)
    min_cap = 4;

  b = ISLSH_REALLOC(NULL, (a) ? islsh_stbds_header(a) : 0, elemsize * min_cap + sizeof(struct islsh_stbds_array_header));
  b = (char *) b + sizeof(struct islsh_stbds_array_header);
  if (a == NULL) {
    islsh_stbds_header(b)->length = 0;
  } 
  islsh_stbds_header(b)->capacity = min_cap;
  return b;
}

static int islsh_array_index_of(int *ids, int id) {
	for (int i = 0; i < islsh_stbds_arrlen(ids); i++) {
		if (ids[i] == id) {
			return i;
		}
	}
	return -1;
}

static int islsh_array_remove(int *ids, int id) {
	for (int i = 0; i < islsh_stbds_arrlen(ids); i++) {
		if (ids[i] == id) {
			islsh_stbds_arrdel(ids, i);
			return i;
		}
	}
	return -1;
}

static int islsh_intersects(struct islsh_xywh rec1, struct islsh_xywh rec2) {
	return rec1.x <= rec2.x + rec2.w && rec1.x + rec1.w >= rec2.x &&
		rec1.y <= rec2.y + rec2.h && rec1.y + rec1.h >= rec2.y;
}

#define ISLSH_BUCKET_INDEX(self,x,y) ((size_t)((x)*65479+(y)))%((self)->spatial_buckets)

static int *islsh_get_ids(struct islsh_spatial_hash *self, int x, int y) {
	struct islsh_bucket *buckets_list = self->buckets[ISLSH_BUCKET_INDEX(self,x,y)];
	for (int i = 0; i < islsh_stbds_arrlen(buckets_list); i++) {
		if (buckets_list[i].x == x && buckets_list[i].y == y) {
			return buckets_list[i].ids;
		}
	}
	return NULL;
}

static void islsh_set_ids(struct islsh_spatial_hash *self, int x, int y, int *ids) {
	size_t index = ISLSH_BUCKET_INDEX(self,x,y);
	struct islsh_bucket *buckets_list = self->buckets[index];
	for (int i = 0; i < islsh_stbds_arrlen(buckets_list); i++) {
		if (buckets_list[i].x == x && buckets_list[i].y == y) {
			buckets_list[i].ids = ids;
			return;
		}
	}
	struct islsh_bucket new_bucket = {x, y, ids};
	islsh_stbds_arrput(buckets_list, new_bucket);
	self->buckets[index] = buckets_list;
}

struct islsh_spatial_hash *islsh_create(float cell_size, int spatial_buckets) {
	struct islsh_spatial_hash *self = ISLSH_REALLOC(NULL, NULL, sizeof *self);
	memset(self, 0, sizeof *self);
	self->spatial_buckets = spatial_buckets;
	self->cell_size = cell_size;
	self->buckets = malloc(sizeof(struct islsh_bucket *) * spatial_buckets);
	for (int i = 0; i < spatial_buckets; i++) {
		self->buckets[i] = NULL;
	}
	return self;
}

int islsh_put(struct islsh_spatial_hash *self, float x, float y, float w, float h) {
	struct islsh_xywh xywh = {x, y, w, h};
	int id;
	if (islsh_stbds_arrlen(self->released_ids) == 0) {
		id = self->id_counter++;
		islsh_stbds_arrput(self->xywhs, xywh);
		islsh_stbds_arrput(self->intersections, NULL);
	} else {
		id = islsh_stbds_arrpop(self->released_ids);
		self->xywhs[id] = xywh;
	}
	int *intersections = NULL;
	int x0 = floor(x/self->cell_size), y0 = floor(y/self->cell_size),
			x1 = ceil((x+w)/self->cell_size), y1 = ceil((y+h)/self->cell_size);
	for (int x_ = x0; x_ < x1; x_++) {
		for (int y_ = y0; y_ < y1; y_++) {
			int *ids = islsh_get_ids(self, x_, y_);
			for (int i = 0; i < islsh_stbds_arrlen(ids); i++) {
				int id_ = ids[i];
				if (id_ != id) {
					struct islsh_xywh xywh_ = self->xywhs[id_];
					if (islsh_intersects(xywh, xywh_)) {
						if (islsh_array_index_of(intersections, id_) < 0) {
							islsh_stbds_arrput(intersections, id_);
						}
						int *intersections_ = self->intersections[id_];
						if (islsh_array_index_of(intersections_, id) < 0) {
							islsh_stbds_arrput(intersections_, id_);
							self->intersections[id_] = intersections_;
						}
					}
				}
			}
			islsh_stbds_arrput(ids, id);
			islsh_set_ids(self, x_, y_, ids);
		}
	}
	self->intersections[id] = intersections;
	return id;
}

int islsh_del(struct islsh_spatial_hash *self, int id) {
	if (islsh_contains(self,id)) {
		struct islsh_xywh xywh = self->xywhs[id];
		float x = xywh.x, y = xywh.y, w = xywh.w, h = xywh.h;
		int x0 = floor(x/self->cell_size), y0 = floor(y/self->cell_size),
				x1 = ceil((x+w)/self->cell_size), y1 = ceil((y+h)/self->cell_size);
		int *intersections = self->intersections[id];
		for (int i = 0; i < islsh_stbds_arrlen(intersections); i++) {
			int id_ = intersections[i];
			islsh_array_remove(self->intersections[id_], id);
		}
		islsh_stbds_arrfree(self->intersections[id]);
		self->intersections[id] = NULL;
		for (int x_ = x0; x_ < x1; x_++) {
			for (int y_ = y0; y_ < y1; y_++) {
				int *ids = islsh_get_ids(self, x_, y_);
				islsh_array_remove(ids, id);
				islsh_set_ids(self, x_, y_, ids);
			}
		}
		self->xywhs[id].x = (0.0/0.0); // NAN
		islsh_stbds_arrput(self->released_ids, id);
		return 0;
	} else {
		return 1;
	}	
}

void islsh_destroy(struct islsh_spatial_hash *self) {
	if (self != NULL) {
		islsh_stbds_arrfree(self->xywhs);
		self->xywhs = NULL;
		islsh_stbds_arrfree(self->released_ids);
		self->released_ids = NULL;
		for (int i = 0; i < islsh_stbds_arrlen(self->intersections); i++) {
			islsh_stbds_arrfree(self->intersections[i]);
		}
		islsh_stbds_arrfree(self->intersections);
		self->intersections = NULL;
		self->id_counter = 0;
		for (int i = 0; i < self->spatial_buckets; i++) {
			islsh_stbds_arrfree(self->buckets[i]);
			self->buckets[i] = NULL;
		}
		free(self->buckets);
		self->buckets = NULL;
		self->spatial_buckets = 0;
		free(self);
	}	
}

int islsh_bucket_len(struct islsh_spatial_hash *self, int x, int y) {
	return islsh_stbds_arrlen(islsh_get_ids(self, x, y));
}

int islsh_intersections_len(struct islsh_spatial_hash *self, int id) {
	return islsh_contains(self,id) ? islsh_stbds_arrlen(self->intersections[id]) : 0;
}

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2019 Ilya Kolbin
Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this 
software, either in source code form or as a compiled binary, for any purpose, 
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this 
software dedicate any and all copyright interest in the software to the public 
domain. We make this dedication for the benefit of the public at large and to 
the detriment of our heirs and successors. We intend this dedication to be an 
overt act of relinquishment in perpetuity of all present and future rights to 
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/
#endif
#endif

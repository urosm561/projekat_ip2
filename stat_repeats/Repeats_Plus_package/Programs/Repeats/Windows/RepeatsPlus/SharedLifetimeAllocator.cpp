#include <cstdlib>
#include "SharedLifetimeAllocator.h"

shared_lifetime_allocator_base::shared_lifetime_memory::shared_lifetime_memory() {
	for (int x = 0; x < sizeof (free_blocks_) / sizeof (free_blocks_[0]); ++x) {
		free_blocks_[x] = nullptr;
	}

	nextFree_ = nullptr;
	nextFreeSize_ = 0;
	fast_free_on_ = false;

	ref_count = 0;
}

shared_lifetime_allocator_base::shared_lifetime_memory::~shared_lifetime_memory() {
	for (auto ptr : allocated_blocks_) {
		std::free(ptr);
	}
}


void* shared_lifetime_allocator_base::shared_lifetime_memory::allocate(size_t bytes) {
	bytes = (bytes + (sizeof(void*)-1)) & ~(sizeof(void*)-1);

	if (bytes >= use_default_allocator_treshold_) {
		return std::malloc(bytes);
	}
	else {
		auto head = free_blocks_[bytes / sizeof(void*)];
		if (head == nullptr) {
			if (nextFreeSize_ < bytes) {
				nextFree_ = std::malloc(block_size_);
				allocated_blocks_.push_back(nextFree_);
				nextFreeSize_ = block_size_;
			}

			auto ret = nextFree_;
			nextFree_ = static_cast<char*> (nextFree_)+bytes;
			nextFreeSize_ -= bytes;

			return ret;
		}
		else {
			void* nextFree = *static_cast<void**> (head);
			free_blocks_[bytes / sizeof (void*)] = nextFree;
			return head;
		}
	}
}

void shared_lifetime_allocator_base::shared_lifetime_memory::deallocate(void* p, size_t bytes) {
	if (p == nullptr) return;

	bytes = (bytes + (sizeof(void*)-1)) & ~(sizeof(void*)-1);

	if (bytes >= use_default_allocator_treshold_) {
		std::free(p);
	}
	else {
		if (fast_free_on_) return;

		auto prevHead = free_blocks_[bytes / sizeof (void*)];
		*static_cast<void**> (p) = prevHead;
		free_blocks_[bytes / sizeof (void*)] = p;
	}
}

void shared_lifetime_allocator_base::set_scope(shared_lifetime_memory* scope) {
	if (scope == scope_) return;

	if (scope_ != nullptr) {
		if (--scope_->ref_count == 0) delete scope_;
	}

	scope_ = scope;

	if (scope_ != nullptr) ++scope_->ref_count;
}

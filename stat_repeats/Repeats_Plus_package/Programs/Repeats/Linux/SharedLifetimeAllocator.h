#pragma once

#ifndef SHARED_LIFETIME_ALLOCATOR_H
#define SHARED_LIFETIME_ALLOCATOR_H

#include <vector>
#include <cstddef>


class shared_lifetime_allocator_base {
	class shared_lifetime_memory {
		enum {
			use_default_allocator_treshold_ = 1024,
			block_size_ = 16384 - 50
		};

		std::vector<void*> allocated_blocks_;

		void* free_blocks_[use_default_allocator_treshold_ / sizeof (void*)];
		void* nextFree_;
		size_t nextFreeSize_;

		shared_lifetime_memory(const shared_lifetime_memory& other) = delete;
		void operator= (const shared_lifetime_memory& other) = delete;

		bool fast_free_on_;

	public:
		int ref_count;

		shared_lifetime_memory();

		~shared_lifetime_memory();

		void fast_free_on() {
			fast_free_on_ = true;
		}

		void* allocate(size_t bytes);

		void deallocate(void* p, size_t bytes);
	};

	template<typename T>
	friend class shared_lifetime_allocator;

	shared_lifetime_memory* scope_;

	shared_lifetime_allocator_base() {
		scope_ = nullptr;
	}

	~shared_lifetime_allocator_base() {
		set_scope(nullptr);
	}

	shared_lifetime_memory* get_scope() const {
		return scope_;
	}

	void set_scope(shared_lifetime_memory* scope);
};

template <class T>
class shared_lifetime_allocator : public shared_lifetime_allocator_base {
public:
	typedef size_t    size_type;
	typedef std::ptrdiff_t difference_type;
	typedef T*        pointer;
	typedef const T*  const_pointer;
	typedef T&        reference;
	typedef const T&  const_reference;
	typedef T         value_type;

	shared_lifetime_allocator(bool create) {
		if (create) {
			set_scope(new shared_lifetime_memory());
		}
	}

	shared_lifetime_allocator(const shared_lifetime_allocator& other) {
		set_scope(other.get_scope());
	}

	shared_lifetime_allocator(shared_lifetime_allocator&& other) {
		scope_ = other.scope_;
		other.scope_ = nullptr;
	}

	template <class U>
	shared_lifetime_allocator(const shared_lifetime_allocator<U>& other) {
		set_scope(other.get_scope());
	}

	template <class U>
	shared_lifetime_allocator(shared_lifetime_allocator<U>&& other) {
		scope_ = other.scope_;
		other.scope_ = nullptr;
	}

	pointer allocate(size_type n, const void * = 0) {
		T* t = static_cast<T*> (get_scope()->allocate(n * sizeof(T)));
		return t;
	}

	void deallocate(void* p, size_type count) {
		get_scope()->deallocate(p, sizeof(T)* count);
	}

	pointer address(reference x) const { return &x; }
	const_pointer address(const_reference x) const { return &x; }
	shared_lifetime_allocator<T>& operator=(const shared_lifetime_allocator&) { return *this; }

	void construct(pointer p, const T& val) {
		new ((T*)p) T(val);
	}

	void destroy(pointer p) { p->~T(); }

	size_type max_size() const { return size_t(-1); }

	template <class U>
	struct rebind { typedef shared_lifetime_allocator<U> other; };

	template <class U>
	shared_lifetime_allocator& operator=(const shared_lifetime_allocator<U>& other) {
		set_scope(other.get_scope());
		return *this;
	}

	void fast_free_on() {
		get_scope()->fast_free_on();
	}
};

#endif
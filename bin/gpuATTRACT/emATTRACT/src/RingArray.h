/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef RINGARRAY_H_
#define RINGARRAY_H_

#include <list>

namespace ema {

template <typename T>
class RingArray {

public:
	using value_type = T;
private:

	using internal_iterator = typename std::list<T>::iterator;

public:

	template <typename _Tp>
	class ring_iterator {
	public:
		using _Self = ring_iterator<_Tp>;
		using iterator_category = std::bidirectional_iterator_tag;
        using value_type = _Tp ;
        using pointer    = _Tp*;
        using reference  = _Tp&;

		ring_iterator() : _iter(), _begin(), _end() {}
		ring_iterator(const internal_iterator& iter, const internal_iterator& begin, const internal_iterator& end) :
			_iter(iter), _begin(begin), _end(end) {}

		ring_iterator(const _Self& iter) = default;
		_Self& operator= (const _Self&) = default;
		ring_iterator (_Self&&) = default;
		_Self& operator= (_Self&&) = default;
//		ring_iterator (_Self&&) = delete;
//		_Self& operator= (_Self&&) = delete;


		reference operator*() const {
			return *_iter;
		}

		pointer operator->() const {
			return static_cast<pointer>(_iter);
		}

		_Self& operator++() { // Prefix
			if (++_iter == _end) {
				_iter = _begin;
			}
			return *this;
		}

		_Self operator++(int) { // Postfix
			_Self tmp = *this;
			++*this;
			return tmp;
		}

		_Self& operator--() { // Prefix
			if (_iter == _begin) {
				_iter = _end;
				--_iter;
			} else {
				--_iter;
			}
			return *this;
		}

		_Self operator--(int) { // Postfix
			_Self tmp = *this;
			--*this;
			return tmp;
		}

		bool operator==(const _Self& x) const {
//			assert(_begin == x._begin);
//			assert(_end == x._end);
			return _iter == x._iter;
		}

		bool operator!=(const _Self& x) const {
//			assert(_begin == x._begin);
//			assert(_end == x._end);
			return _iter != x._iter;
		}

	private:

		internal_iterator _iter;
		internal_iterator _begin;
		internal_iterator _end;

		template <typename _Tf> friend class RingArray;

	};

	using iterator = ring_iterator<T>;


	/****************************
	 * public member functions
	 ****************************/
	iterator begin() {
		return iterator(_ring.begin(), _ring.begin(), _ring.end());
	}

	iterator end() {
		return iterator(_ring.end(), _ring.begin(), _ring.end());
	}

	unsigned size() {
		return _ring.size();
	}


	/*
	 ** @brief: The container is extended by inserting new elements before the element at the specified position.
	 */
	iterator insert(const iterator& iter, const T& value) {
		internal_iterator tmp = _ring.insert(iter._iter, value);
		return iterator(tmp, iter._begin, iter._end);
	}

	iterator insert(const iterator& iter, T&& value) {
		internal_iterator tmp = _ring.insert(iter._iter, value);
		return iterator(tmp, iter._begin, iter._end);
	}

	iterator erase(const iterator& iter) {
		internal_iterator tmp = _ring.erase(iter._iter);
		return iterator(tmp, iter._begin, iter._end);
	}


	void push_back(const T& value) {
		_ring.push_back(value);
	}

	void push_back(T&& value) {
		_ring.push_back(value);
	}

	void push_front(const T& value) {
		_ring.push_front(value);
	}

	void push_front(T&& value) {
		_ring.push_front(value);
	}

	void pop_front() {
		_ring.pop_front();
	}

	void pop_back() {
		_ring.pop_back();
	}



private:
	/****************************
	 * private member functions
	 ****************************/



	/****************************
	 * private member variables
	 ****************************/
	std::list<T> _ring;

};

}  // namespace


#endif /* RINGARRAY_H_ */

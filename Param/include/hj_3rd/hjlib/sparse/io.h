#ifndef _HJ_SPARSE_IO_H_
#define _HJ_SPARSE_IO_H_

#include "format.h"

namespace hj { namespace sparse {

template <typename OS, typename VAL_TYPE, typename INT_TYPE,
		  template <typename K, typename T> class MAP_TYPE>
OS &operator << (OS &os, const map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE> &v) {
	typedef typename map_vec<VAL_TYPE, INT_TYPE, MAP_TYPE>::container container;
	os << "[ " << v.size() << ':' << nnz(v) << '\n';
	for(typename container::const_iterator i = v.data_.begin(); i != v.data_.end(); ++i) {
		os << "  (" << i->first << ':' << i->second << ")\n";
	}
	os << "]\n";
	return os;
}

}}

#endif

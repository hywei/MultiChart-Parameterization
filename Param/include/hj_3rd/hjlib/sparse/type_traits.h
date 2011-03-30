#ifndef HJ_TYPE_TRAITS_H_
#define HJ_TYPE_TRAITS_H_

namespace hj { namespace sparse {

template <typename T>
class value_type
{
public:
	static const unsigned char type = 'r';
	static const unsigned char size = sizeof(T);
};

template <typename T>
class value_type<std::complex<T> >
{
public:
	static const unsigned char type = 'c';
	static const unsigned char size = sizeof(T);
};

}}

#endif

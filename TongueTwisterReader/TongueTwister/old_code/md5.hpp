#ifndef MD5_H
#define MD5_H

#include <cstring>
#include <iostream>


class MD5Hash
	{
	public:
		MD5Hash();
		MD5Hash(const std::string& text);
		inline int Compare(MD5Hash& m) { return memcmp(digest, m.digest, sizeof(digest)); }

	private:
		unsigned char digest[16];
    };

#endif
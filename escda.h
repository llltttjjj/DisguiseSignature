#ifndef ESCDA_H_
#define ESCDA_H_

#include<iostream>
#include<NTL/ZZ.h>
using namespace NTL;
using namespace std;

static ZZ mod, n;
static ZZ a, b;     //Elliptic Curve Function: y ^ 2 = x ^ 3 + a * x + b
const int SIZE = 256;
static const ZZ pow(uint32_t n) {
	if (n == 0)
		return (ZZ)1;
	return 2 * pow(n - 1);
}
static const ZZ modPow(const ZZ& num, ZZ t, const ZZ& n) {
	ZZ temp, ans;
	ans = 1;
	temp = num;
	while (t > 0) {
		if (t % 2 == 1)
			ans = ans * temp % n;
		temp *= temp;
		temp %= n;
		t /= 2;
	}
	return ans;
}
const uint32_t cut(const ZZ& num) {     //Take Low 32 Bits Of ZZ, Transfer Into Uint32_t
	uint32_t n = 0;
	uint32_t temp = 1;
	ZZ x = num;
	for (int i = 0; i < 32; i++) {
		if (x % 2 == 1)
			n += temp;
		temp += temp;
		x /= 2;
	}
	return n;
}
static const ZZ reverse(const ZZ& num, const ZZ& n) { return modPow(num, n - 2, n); }

class ellPoint {
public:
	ZZ x, y;
	ellPoint(ZZ _x = (ZZ)0, ZZ _y = (ZZ)0) { x = _x; y = _y; }
	const bool setEllPoint(const ZZ& _x) {
		x = _x;
		ZZ num = (x * x * x + a * x + b) % mod;
		if (modPow(num, (mod - 1) / 2, mod) == 1) {
			y = modPow(num, (mod + 1) / 4, mod);
			return 1;
		}
		return 0;
	}
	void pointHash(uint32_t* hash) {     //x,y are 256 bits
		uint32_t XY[16];
		for (uint32_t i = 0; i < 8; i++) {
			XY[i] = cut(x / pow(7 - i));
			XY[i + 8] = cut(y / pow(7 - i));
		}
	}
	const bool OnCurve() { return y * y % mod == (x * x * x + a * x + b) % mod; }
	const ellPoint operator+(const ellPoint& n)const {
		ellPoint p;
		ZZ lambda;
		if (x == n.x && y == n.y)
			lambda = (3 * x * x + a) * reverse(2 * y, mod) % mod;
		else
			lambda = (n.y - y) * reverse(n.x - x, mod) % mod;
		p.x = (lambda * lambda - x - n.x) % mod;
		p.y = (lambda * (x - p.x) - y) % mod;
		return p;
	}
	const ellPoint& operator+=(const ellPoint& n) { *this = *this + n; return *this; }
	const ellPoint operator*(const ZZ& n)const {
		ellPoint p;
		ellPoint temp = *this;
		ZZ num = n;
		while (num % 2 == 0) {
			num /= (ZZ)2;
			temp += temp;
		}
		p = temp;
		while (num > 0) {
			num /= 2;
			temp += temp;
			if (num % 2 == 1)
				p += temp;
		}
		return p;
	}
	friend const ellPoint operator*(const ZZ& n, const ellPoint& p) { return p * n; }
	const ellPoint& operator*=(const ZZ& n) { *this = *this * n; return *this; }
	friend ostream& operator<<(ostream& os, const ellPoint& n) {
		os << "(" << n.x << ", " << n.y << ")";
		return os;
	}
	const bool operator!=(const ellPoint& p) const { if (x == p.x && y == p.y) return 0; return 1; }
	const bool operator==(const ellPoint& p) const { if (*this != p) return 0; return 1; }
	void place(uint32_t* pos)const {
		ellPoint p = *this;
		uint32_t* mem = pos;
		for (int i = SIZE / 32 - 1; i >= 0; i--) {
			mem[i] = cut(p.x);
			p.x /= pow((uint32_t)32);
		}
		mem = &pos[SIZE / 32];
		for (int i = SIZE / 32 - 1; i >= 0; i--) {
			mem[i] = cut(p.y);
			p.y /= pow((uint32_t)32);
		}
	}
};
const ZZ transfer(const uint32_t* ptr) {
	ZZ ans = (ZZ)0;
	ZZ temp;
	for (int i = 0; i < 8; i++) {
		if (ptr[i] >= 2147483648) {
			temp = ptr[i];
			temp += pow((uint32_t)32);
		}
		else
			temp = ptr[i];
		ans += temp * pow((uint32_t)(7 - i) * 32);
	}
	return ans;
}

void GenerateSignature() {
	const uint32_t arrmod[8] = { 0x8542D69E, 0x4C044F18 ,0xE8B92435 ,0xBF6FF7DE ,0x45728391 ,0x5C45517D ,0x722EDB8B ,0x08F1DFC3 };
	mod = transfer(arrmod);
	const uint32_t arra[8] = { 0x787968B4 ,0xFA32C3FD ,0x2417842E ,0x73BBFEFF ,0x2F3C848B ,0x6831D7E0 ,0xEC65228B ,0x3937E498 };
	const uint32_t arrb[8] = { 0x63E4C6D3 ,0xB23B0C84 ,0x9CF84241 ,0x484BFE48 ,0xF61D59A5 ,0xB16BA06E ,0x6E12D1DA ,0x27C5249A };
	a = transfer(arra); b = transfer(arrb);
	const uint32_t arrn[8] = { 0x8542D69E, 0x4C044F18, 0xE8B92435, 0xBF6FF7DD, 0x29772063, 0x0485628D, 0x5AE74EE7, 0xC32E79B7 };
	ZZ n = transfer(arrn);
	ellPoint G, P;
	const uint32_t arrX_G[8] = { 0x421DEBD6, 0x1B62EAB6, 0x746434EB, 0xC3CC315E, 0x32220B3B, 0xADD50BDC, 0x4C4E6C14, 0x7FEDD43D };
	const uint32_t arrY_G[8] = { 0x0680512B, 0xCBB42C07, 0xD47349D2, 0x153B70C4, 0xE5D7FDFC, 0xBFA36EA1, 0xA85841B9, 0xE46E09A2 };
	G.x = transfer(arrX_G); G.y = transfer(arrY_G);
	ZZ d = RandomBits_ZZ(SIZE) % n;     //have trouble getting Satoshi's secret key, take a random one instead :(
	P = d * G;
	ZZ u, v, r, s, e;     //basic numbers in signature
	u = RandomBits_ZZ(100);
	v = RandomBits_ZZ(100);
	ellPoint R = u * G + v * P;
	r = R.x;
	e = r * u * reverse(v, n) % n;
	s = r * reverse(v, n) % n;
	if (e * reverse(s, n) * G + r * reverse(s, n) * P == R)
		cout << "r = " << r << endl << "s = " << s << endl 
		<< "is a valid signature of e = " << e << endl << "with secret key d." << endl;
	else
		cout << "Invalid key!" << endl;
}
#endif

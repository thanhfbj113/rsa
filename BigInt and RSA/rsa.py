import numpy as np
import time
import random


maxInt = np.iinfo(np.uint64).max + 1

class BigInt :

	def __init__(self, inp):
		self.dwords = np.zeros(16,dtype = np.uint64)
		if inp == '':
			return
		inpBin = bin(inp)[2:]
		for i in reversed(range(16)) :
			left  = max(0,len(inpBin) - (16-i) * 64  ) 
			right = len(inpBin) - (15-i)*64 
			tmp = int(inpBin[left: right],2) 

			self.dwords[i] = tmp
			if left == 0 :
				break
	def Is_Negative(self):
		if self.dwords[0] & 0XA000000000000000:
			return True
		else: 
			return False
	def __gt__(self,opt): # So sánh lớn hơn
		for i in range(16):
			if self.dwords[i] > opt.dwords[i]:
				return True
			elif self.dwords[i] < opt.dwords[i]:
				return False
	def __eq__(self,opt): # So sánh bằng
		for i in range(16):
			if self.dwords[i] != opt.dwords[i]:
				return False
		return True
	def __ge__(self,opt): # So sánh lớn hơn hoặc bằng
		if self > opt or self == opt:
			return True
		return False
	# def __eq__(self,i):
	# 	if i < 0:
	# 		return False
	# 	tmp = BigInt(i)
	# 	return self == tmp
	def __neg__(self):
		pat = 1
		for i in reversed(range(16)):
			self.dwords[i] = self.dwords[i] ^ 0XFFFFFFFFFFFFFFFF
			if pat == 1 :
				if self.dwords[i] == 0XFFFFFFFFFFFFFFFF :
					print("i= ",i,":",self.dwords[i])
					self.dwords[i] = 0
				else:
					self.dwords[i] += 1
					pat = 0
		return self
	def __add__(self, other): # Cộng bigInt
		result = BigInt(0)

		surplus = 0
		for i in reversed(range(16)):
			res = int(self.dwords[i]) + int(other.dwords[i]) + surplus
			surplus = res // maxInt
			result.dwords[i] = res % maxInt

		return result
	def __rshift__(self, other): # Dịch phải other bit
		result = BigInt(0)
		k = other // 64
		if k > 0:
			for i in reversed(range(16 - k)):
				result.dwords[i + k] = self.dwords[i]
		else:
			result.dwords = self.dwords.copy()
		k = other % 64
		surplus = 0
		for i in range(16):
			temp = int(result.dwords[i]) << 64 - k
			result.dwords[i] = (int(result.dwords[i]) >> k) | surplus
			surplus = temp % maxInt
		return result

	def __lshift__(self, other): # Dịch trái other bit
		result = BigInt(0)
		k = other // 64
		if k > 0:
			for i in range(k, 16):
				result.dwords[i - k] = self.dwords[i]
		else:
			result.dwords = self.dwords.copy()
		k = other % 64
		surplus = 0
		for i in reversed(range(16)):
			temp = int(result.dwords[i]) >> 64 - k
			result.dwords[i] = ((int(result.dwords[i]) << k) % maxInt) | surplus
			surplus = temp
		return result
	def __sub__(self, other): # Trừ bigInt
		result = BigInt(0)

		sursub = 0
		for i in reversed(range(16)):
			res = int(self.dwords[i]) - int(other.dwords[i]) - sursub
			sursub = 1 if res < 0 else 0
			result.dwords[i] = res if res>=0 else maxInt + res
		return result

	def __str__(self): # to String
		res = 0
		for i in range(16):
			res = (res << 64) + int(self.dwords[i])
		return str(res)

	# def __div__(self, other):
	# 	if other > self:
	# 		return BigInt(0)

	# 	x = self
	# 	y = other

	# 	diff = Countbit(x) - Countbit(y)

	# 	y = y << diff
	# 	if y > x:
	# 		y = y >> 1
	# 		diff -= 1
	# 	result = BigInt(2**diff)
		
	# 	while x > other:

	# 	x = x - (y << diff)




	def __mod__(self, other): # modulo cơ số other
		if other == BigInt(0) :
			return -1 
		if other > self :
			return self
		if other == self:
			return 0
		x = self
		y = other 

		diff = Countbit(x) - Countbit(y) # Độ chênh lệch bits giữa x và y
		# Dịch phải cho đến khi y gần bằng x (y = k*other)
		if diff > 0:
			y = y << diff
		if (y > x):
			y = y >> 1

		while x >= other: # Kiểm tra đến khi x nhỏ hơn cơ số modulo thì ngừng
			x = x - y
			if y > x: 
				diff = Countbit(y) - Countbit(x) # Độ chênh lệch bits giữa y và x
				# đưa y về dạng k*other nhưng nhỏ hơn x
				if diff > 0:
					y = y >> diff 
				if y > x:
					y = y >> 1
		return x

	def Getbit(self, pos): # Kiểm tra bits là 0 hay 1 tại vị trí cho trước
		div = pos // 64 
		sur = pos % 64
		return (int(self.dwords[15-div]) >> sur) & 1
def AddMod(a,b,n): # Cộng modulo
	res = a + b
	res = res % n
	return res
def Countbit(bInt): # Đếm số lượng bit của một bigInt
	result = 0
	count = 16
	for i in range(16):
		if int(bInt.dwords[i]) == 0:
			count -= 1
		else:
			result = len(bin(int(bInt.dwords[i])))-2
			break
	if count == 0:
		return 0
	result += (count-1)*64
	return result

def MulMod(a, b, m): # Nhân modulo
	result = BigInt(0)
	if a > m:
		a = a % m
	if b > m:
		b = b % m
	countB = Countbit(b)
	# Nhân từng bit của b với 2**i và a
	for i in range(countB): 
		if b.Getbit(i) == 1:
			temp = (a << i) % m
			result = (result + temp) % m
	return result


def PowMod(a,b,m): # Lũy thừa modulo
	if b == ZERO:
		return BigInt(1)
	res = BigInt(1)
	if a > m:
		a = a % m
	if b.Getbit(0) == 1 :
		res = a
	tmp = a

	countB = Countbit(b)
	for i in range(1,countB): 
		tmp = MulMod(tmp,tmp,m) # a**(2i)
		if b.Getbit(i) == 1 :
			res = MulMod(res,tmp,m)
	return res
def GCD(a,b):
	zero = BigInt(0)
	if a == zero or b == zero:
		return a+b

	g = BigInt(1)
	# kiem tra a, b chan
	while a.Getbit(0) == 0 and b.Getbit(0) == 0  :
		a = a >> 1	# a /=2
		b = b >> 1	# b/=2
		g = g << 1  #g *=2
	while a > zero : 
		#dua a ve so len
		while a.Getbit(0) == 0:
			a = a >> 1
		#dua b ve so le
		while b.Getbit(0) == 0:
			b = b >> 1
		# t = |a-b|/2
		if a >= b :
			a = a-b
			a = a >> 1
		else:
			b = b - a
			b = b >> 1

	return MulMod(g,b,BigInt(2**1024 -1 ))

def GenKey(p,q):
	#n
	ONE = BigInt(1)

	n = MulMod(p,q,BigInt(2**1024-1))

	q= q - ONE
	p=p-ONE
	
	print(p,q)
	phi = MulMod(p,q,n)	
	#find e
	e = BigInt(3)
	TWO = BigInt(2)
	while GCD(e,phi) > ONE :
		e+= TWO
	# tmp = random.randrange(1,int(str(phi))/2-1)
	# e = BigInt(tmp*2 + 1)
	# while GCD(e,phi) > ONE :
	# 	tmp = random.randrange(1,int(str(phi))/2-1)
	# 	e = BigInt(tmp*2 + 1)
	tmp = pow(int(str(e)),-1,int(str(phi)))
	d = BigInt(tmp)
	d = d % phi
	return ((e,n),(d,n))
def Encrypt(pubKey , plaintext):
	e,n = pubKey
	cipher = []
	for c in plaintext:
		crs = int(str(PowMod(BigInt(ord(c)),e,n)))
		char = ''
		while crs > 0 :
			char += chr(crs % 2**7 + 65)
			crs //= 2**7
		cipher.append(char)
	return cipher
def FastDecrypt(prvKey,p,q, ciphers):
	TWO = BigInt(2)
	'''
	prvKey la khoa bao mat
	p,q la 2 so tu nhien tao ra n va thuoc kieu BIGINT
	cipher la 
	p < q
	'''
	d, n = prvKey

	# lay message
	cip = []
	for cipher in ciphers:
		chars = 0 
		t = 1
		for char in cipher:
			chars += (ord(char) - 65) * t
			t*= 128
		cip.append(BigInt(chars))

	

	# Dp = d mod phi(p-1)  rut gon d
	Dp = d % (p-ONE)
	# Dq = d mod phi(q-1)  rut gon d
	Dq = d % (q-ONE)
	ret = ''
	# ret = ''
	for c in cip :
		# Rut gon lai c
		Cp = c % p
		Cq = c % q

		# Tinh char theo p
		c1 = PowMod(Cp,Dp,p)
		# Tinh char theo q
		c2 = PowMod(Cq,Dq,q)

		#Ap dung m = c1 mod p 
		#		 m = c2 mod q
		#      =================
		#       c2 +q*x = c1 mod p          vi q > p

		#rut gon p theo q  va tinh q^-1 trong p
		q1 = PowMod((q%p),p-TWO,p)		# q^p-2 ap dung dinh ly EUler
		# char = c1 - c2       x = (c2-c1)*q1 mod p
		char = c1 -c2
		x = MulMod(char,q1,p)
		# thế x nguoc vào m = c2 + q*x
		ret += chr(int(str(c2 + MulMod(x,q,n))))


		#Cach khac ap dung dong du trung hoa:
	return ret





ZERO = BigInt(0)
ONE = BigInt(1)
		#Big int 1
#p < q
p = BigInt(57787)
q = BigInt(60521)
t0 = time.time()
pubKey  , prvKey = GenKey(p,q)

print("public Key", pubKey[0],pubKey[1])
print("Private Key", prvKey[0],prvKey[1])


text = input('Chuoi can ma hoa: ')
print('Day chua ma hoa: ', text)
cipher = Encrypt(pubKey, text)
t1 = time.time()
print("da ma hoa",[x for x in cipher])
t2 = time.time()
mes = FastDecrypt(prvKey,p,q,cipher)
t3 = time.time()

print('processing time: ',t3 - t2)
print(mes)


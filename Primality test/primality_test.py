from math import sqrt
import random

count_assign = 0
count_compare = 0


# from var_global import count_assign, count_compare
def miller_rabin_base_2(n):
    global count_assign, count_compare
    """Perform the Miller Rabin primality test base 2"""
    d = n-1
    s = 0
    count_assign += 2

    count_compare += 1
    while not d & 1: # Check for divisibility by 2
        count_compare += 1
        d = d >> 1 # Divide by 2 using a binary right shift
        s += 1
        count_assign += 2

    x = pow(2, d, n)
    count_assign += 1

    count_compare += 2
    if x == 1 or x == n-1:
        return True

    count_assign += 1
    count_compare += 1
    for i in range(s-1):
        count_assign += 1
        count_compare += 1
        x = pow(x, 2, n)
        count_assign += 1


        count_compare += 2
        if x == 1:
            count_compare -= 1
            return False
        elif x == n - 1:
            return True
    return False

def jacobi_symbol(a, n):
    global count_assign, count_compare
    """Calculate the Jacobi symbol (a/n)"""
    

    count_compare += 5
    if n == 1:
        count_compare -= 4
        return 1
    elif a == 0:
        count_compare -= 3
        return 0
    elif a == 1:
        count_compare -= 2
        return 1
    elif a == 2:
        count_compare -= 1
        if n % 8 in [3, 5]:
            count_compare += 1
            return -1
        elif n % 8 in [1, 7]:
            count_compare += 2
            return 1
    elif a < 0:
        return (-1)**((n-1)/2) * jacobi_symbol(-1*a, n)

    if a % 2 == 0:
        count_compare += 1
        return jacobi_symbol(2, n) * jacobi_symbol(a / 2, n)
    elif a % n != a:
        count_compare += 2
        return jacobi_symbol(a % n, n)
    else:
        count_compare += 2

        count_compare += 2
        if a % 4 == n % 4 == 3:
            return -1 * jacobi_symbol(n, a)
        else:
            return jacobi_symbol(n, a)


def U_V_subscript(k, n, U, V, P, Q, D):
    global count_assign, count_compare
    k, n, U, V, P, Q, D = map(int, (k, n, U, V, P, Q, D))
    digits = list(map(int, str(bin(k))[2:]))
    subscript = 1
    count_assign += 9

    count_assign += 1
    count_compare += 1
    # print(len(digits))
    for digit in digits[1:]:
        count_assign += 1
        count_compare += 1
        U, V = U*V % n, (pow(V, 2, n) - 2*pow(Q, subscript, n)) % n
        subscript *= 2
        count_assign += 3


        count_compare += 1
        if digit == 1:
            count_compare += 2
            if not (P*U + V) & 1:
                if not (D*U + P*V) & 1:
                    U, V = (P*U + V) >> 1, (D*U + P*V) >> 1
                    count_assign += 2
                else:
                    U, V = (P*U + V) >> 1, (D*U + P*V + n) >> 1
                    count_assign += 2
            elif not (D*U + P*V) & 1:
                U, V = (P*U + V + n) >> 1, (D*U + P*V) >> 1
                count_assign += 2
            else:
                U, V = (P*U + V + n) >> 1, (D*U + P*V + n) >> 1
                count_assign += 2
            subscript += 1
            U, V = U % n, V % n
            count_assign += 3
    return U, V

def lucas_pp(n, D, P, Q):                                                                                                                                                                                                                         
    global count_assign, count_compare
    """Perform the Lucas probable prime test"""
    U, V = U_V_subscript(n+1, n, 1, P, P, Q, D)
    count_assign += 2

    count_compare += 1
    if U != 0:
        return False

    d = n + 1
    s = 0
    count_assign += 2

    count_compare += 1
    while not d & 1:
        count_compare += 1
        d = d >> 1
        s += 1
        count_assign += 2

    U, V = U_V_subscript(n+1, n, 1, P, P, Q, D)
    count_assign += 2

    count_compare += 1
    if U == 0:
        return True

    count_assign += 1
    count_compare += 1
    for r in xrange(s):
        count_assign += 1
        count_compare += 1

        U, V = (U*V) % n, (pow(V, 2, n) - 2*pow(Q, d*(2**r), n)) % n
        count_assign += 2
        count_compare += 1
        if V == 0:
            return True

    return False


def D_chooser(candidate):
    global count_assign, count_compare
    """Choose a D value suitable for the Baillie-PSW test"""
    D = 5
    count_assign += 1

    count_compare += 1
    while jacobi_symbol(D, candidate) != -1:
        count_compare += 1
        D += 2 if D > 0 else -2
        D *= -1
        count_assign += 2
    return D

def baillie_psw(candidate):
    global count_assign, count_compare
    """Perform the Baillie-PSW probabilistic primality test on candidate"""
    # Check divisibility by a short list of primes less than 50
    count_assign += 1
    count_compare += 1
    for known_prime in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        count_assign += 1
        count_compare += 1

        count_compare += 2
        if candidate == known_prime:
            count_compare -= 1
            return True
        elif candidate % known_prime == 0:
            return False

    # Now perform the Miller-Rabin primality test base 2
    count_compare += 1
    if not miller_rabin_base_2(candidate):
        return False
    
    # Check that the number isn't a square number, as this will throw out 
    # calculating the correct value of D later on (and means we have a
    # composite number)
    # the slight ugliness is from having to deal with floating point numbers
    count_compare += 1
    if int(sqrt(candidate) + 0.5) ** 2 == candidate:
        return False

    # Finally perform the Lucas primality test
    count_assign += 1
    D = D_chooser(candidate)

    count_compare += 1
    if not lucas_pp(candidate, D, 1, (1-D)/4):
        return False

    # You've probably got a prime!

    return True

def IsprimeSimple(candidate):
    global count_assign, count_compare

    count_compare += 3
    if candidate < 3:
        count_compare -= 2
        return candidate > 1
    elif candidate % 2 == 0 or candidate % 3 == 0:
        return False

    count_compare += 1
    count_assign += 1
    for i in range(5, int(sqrt(candidate)) + 1, 6):
        count_compare += 1
        count_assign += 1

        count_compare += 2
        if candidate % i == 0 or candidate % (i + 2) == 0:
            return False
    return True


def power(x, y, p):
    global count_assign
    global count_compare
    res = 1;
    x = x % p;
    count_assign+=2
    while (y > 0):
        count_compare+=2
        if (y & 1):
            res = (res * x) % p;
            count_assign+=1
        y = y>>1;
        x = (x * x) % p;
        count_assign+=2
    count_compare+=1
    return res;

def miillerTest(d, n):
    global count_assign
    global count_compare
    a = 2 + random.randint(1, n - 4);
    x = power(a, d, n);
    count_assign+=2
    count_compare+1
    if (x == 1 or x == n - 1):
        return True;
    while (d != n - 1):
        x = (x * x) % n;
        d *= 2;
        count_assign+=2
        count_compare+=1
        if (x == 1):
            return False;
        count_compare+=1
        if (x == n - 1):
            return True;
    count_compare+=1
    return False;


def isPrime(n, k):
    global count_assign
    global count_compare
    count_compare+=1
    if (n <= 1 or n == 4):
        return False;
    count_compare+=1
    if (n <= 3):
        return True;
    d = n - 1;
    count_assign+=1
    while (d % 2 == 0):
        d //= 2;
        count_compare+=1
        count_assign+=1
    count_compare+=1
    for i in range(k):
        count_compare+=2
        if (miillerTest(d, n) == False):
            return False;
    count_compare+=1
    return True;


def readFile(nameFi):
    f = open(nameFi, 'r')
    buf = f.read()
    f.close()
    tmp= buf.split()
    result = []
    for i in tmp:
        result.append(int(i))
    return result

def main():
    global count_assign, count_compare
    print('===   PRIMALITY TEST   ===')
    print('1. Simple')
    print('2. Miller Rabin')
    print('3. Baillie PSW')
    print('4. Simple with file test')
    print('5. Miller Rabin with file test')
    print('6. Baillie PSW with file test')
    try:
        choice = int(input('Your choice:'))
    except:
        print('ERRO!')
        return

    if choice not in range(1,7):
        print('ERRO!')
        return
    if choice == 4:
        try:
            fi = open('test_simple.txt')
            fo = open('test_simple_out.txt', 'w')
        except:
            print("Can not open file.")
            return
        line = fi.readline()
        while line:
            fo.write(str(int(line)) + ' ' + str(IsprimeSimple(int(line))) + '\n')
            line = fi.readline()
        fi.close()
        fo.close()
        return
    if choice == 5:
        try:
            fi = open('test.txt')
            fo = open('output_millerRabin.txt', 'w')
            k = int(input('Input k:'))
        except:
            print("Can not open file.")
            return
        line = fi.readline()
        while line:
            fo.write(str(int(line)) + ' ' + str(isPrime(int(line), k)) + '\n')
            line = fi.readline()
        fi.close()
        fo.close()
        return
    if choice == 6:
        try:
            fi = open('test.txt')
            fo = open('output_bailliePSW.txt', 'w')
        except:
            print("Can not open file.")
            return
        line = fi.readline()
        while line:
            fo.write(str(int(line)) + ' ' + str(baillie_psw(int(line))) + '\n')
            line = fi.readline()
        fi.close()
        fo.close()
        return

    try:
        inp = int(input('Input candidate:'))
    except:
        print('Input ERRO!')
        return

    if inp < 1:
        print('ERRO!')
        return
    if choice == 1:
        if IsprimeSimple(inp) == True:
            print('It\'s a prime number.')
        else:
            print('It\'s not a prime number.')
        print('Count compare:',count_compare)
        print('Count assign:',count_assign)
    elif choice == 2:
        try:
            k = int(input('Input k:'))
        except:
            print('Input ERRO!')
            return
        if isPrime(inp, k) == True:
            print('It\'s a prime number.')
            print('Reliability:', (1 - 1 / (4 ** k)) * 100, '%')
        else:
            print('It\'s not a prime number.')
        print('Count compare:',count_compare)
        print('Count assign:',count_assign)
    elif choice == 3:
        if baillie_psw(inp) == True:
            print('It\'s a prime number.')
            k = len(bin(inp)) - 2
            print('Reliability:', (1 - (4 / 15) ** k) * 100, '%')
        else:
            print('It\'s not a prime number.')
        print('Count compare:',count_compare)
        print('Count assign:',count_assign)


main()
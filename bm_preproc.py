#!/usr/bin/env python

import unittest
from collections import defaultdict

# Boyer-Moore preprocessing functions
def z_array(s):
    """ Z algorithm for preprocessing string s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                z[k] = zkp
            else:
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z

def n_array(s):
    """ Return the N array from the Z array (reversed) """
    return z_array(s[::-1])[::-1]

def big_l_prime_array(p, n):
    """ Compute L' array using pattern and N array """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp

def big_l_array(p, lp):
    """ Compute L array from L' array """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l

def small_l_prime_array(n):
    """ Compute lp' array from N array """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp

def good_suffix_table(p):
    """ Return good suffix tables (L', L, lp') """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)

def dense_bad_char_tab(p, amap):
    """ Create a dense bad character table """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

class BoyerMoore(object):
    """ Boyer-Moore pattern matching class """
    
    def __init__(self, p, alphabet='ACGT'):
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        self.bad_char = dense_bad_char_tab(p, self.amap)
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Return shift as determined by good suffix rule """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return shift if full match occurs """
        return len(self.small_l_prime) - self.small_l_prime[1]

# Reading the FASTA file and applying Boyer-Moore search
def search_pattern_in_fasta(fasta_file, pattern):
    # Read the sequence from the FASTA file
    with open(fasta_file, 'r') as f:
        sequence = ''.join([line.strip() for line in f.readlines() if not line.startswith('>')])

    # Create the Boyer-Moore object with the given pattern
    bm = BoyerMoore(pattern)

    # Start Boyer-Moore matching process
    matches = []
    i = 0
    while i <= len(sequence) - len(pattern):
        j = len(pattern) - 1
        while j >= 0 and pattern[j] == sequence[i + j]:
            j -= 1
        if j < 0:
            matches.append(i)
            i += bm.match_skip()
        else:
            shift_bad = bm.bad_character_rule(j, sequence[i + j])
            shift_good = bm.good_suffix_rule(j)
            i += max(shift_bad, shift_good)

    return matches

# Example of using the search
if __name__ == '__main__':
    fasta_file = 'chr1.GRCh38.excerpt.fasta'  # path to your FASTA file
    pattern = 'ATCG'  # change to your desired pattern
    
    matches = search_pattern_in_fasta(fasta_file, pattern)
    
    print(f"Pattern {pattern} found at positions: {matches}")

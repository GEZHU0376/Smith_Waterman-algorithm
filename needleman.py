#!/usr/bin/env python3
#This program include multiple functions and needleman-Wunsch method

def create_submat (match, mismatch, alphabet):
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if (c1==c2):
                sm[c1+c2] = match
            else:
                sm[c1+c2] = mismatch
    return sm

def score_pos(c1,c2,sm,g):
    if c1 == "-" or c2 == "-":
        return g
    else:
        return sm[c1+c2]
    
def print_mat(mat):
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print('{0:4d}'.format(mat[i][j]),end = '')
        print()
    print()

def max3t (v1, v2, v3):
    if max(v1, v2, v3) == v1:
        return 1
    if max(v1, v2, v3) == v2:
        return 2
    if max(v1, v2, v3) == v3:
        return 3
    
def recover_alignment (T, seq1, seq2):
    res = ["",""]
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if T[i][j] == 1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0]
            res[1] = seq2[j-1] + res[1]
            j -= 1
        else:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res

# Needleman_Wunsch
def needleman_Wunsch (seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    for j in range(1, len(seq2) + 1):
        S[0].append(g*j)
        T[0].append(3)
    
    for i in range(1, len(seq1) + 1):
        S.append([g*i])
        T.append([2])
    
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            S[i+1].append(max(s1,s2,s3))
            T[i+1].append(max3t(s1,s2,s3))
    return (S,T)

def test_global():
   sm = create_submat(1,-1,"ACGT")
   S,T = needleman_Wunsch("ACGT","ATT",sm,-2)
   res = recover_alignment(T,"ACGT","ATT")
   print("The needleman wunsch global alignment: ")
   print (res)

if  __name__ == "__main__":
   test_global()

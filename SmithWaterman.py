#!/usr/bin/env python3
#This program import the functions from the needleman.py
#This program also contains smith waterman method

import needleman as needle

def max_mat(mat):
    maxval = mat[0][0]
    maxrow = 0
    maxcol = 0
    for i in range(0,len(mat)):
        for j in range(0, len(mat[i])):
            if mat[i][j] > maxval:
                 maxval = mat[i][j]
                 maxrow = i
                 maxcol = j
    return (maxrow,maxcol)

def recover_align_local(S,T,seq1,seq2):
    res = ["",""]
    i,j = max_mat(S)
    while T[i][j] > 0:
        if T[i][j] == 1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0]
            res[1] = seq2[j-1] + res[1]
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res

def smith_Waterman (seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    maxscore = 0
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + needle.score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1,s2,s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(needle.max3t(s1,s2,s3))
                if b > maxscore:
                    maxscore = b
    return (S, T, maxscore)

def test_local():
   sm = needle.create_submat(1,-1,"ACGT")
   S,T,maxscore = smith_Waterman("ACGT","ATT",sm,-2)
   res = recover_align_local(S,T,"ACGT","ATT")
   print("The smith waterman local alignment result: ")
   print (res)
   print("The maxscore is: ")
   print(maxscore)

if  __name__ == "__main__":
   test_local()

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//  Each Table Entry consists of a value and a pointer to another
//      entry for the traceback
//  Also has three strings for the values of the sequences. Used
//      to compute the scores
struct TableEntry {
    int val;
    string pt;
    string s1;
    string s2;
    string s3;
};

//  Blast Scores
const int match = 5;
const int mismatch2 = -4;
const int indel = -8;

//  -----------------------------GETSCORE WORKS-------------------

// Gets the Scoring Matrix from user & converts all scores into ints
// aa, cc, gg, and tt are the match scores
// ac, ag, at, cg, ct, and gt are the mismatch scores
// ab, cb, gb, tb, and bb are the indel scores
void GetScore(string fileName, int &aa, int &ac, int &ag, int &at, 
    int &cc, int &cg, int &ct, int &gg, int &gt, int &tt, int &ab, 
    int &cb, int &gb, int &tb, int &bb) {
    
    int n,i,j;
    ifstream fin(fileName.c_str());
    vector< vector<int> > vect;
    vect.resize(5);
    for (i = 0; i < 5; ++i) {
        vect.at(i).resize(5);
    }
    i = 0, j = 0;
    while (fin >> n) {
        vect.at(i).at(j) = n;
        if (j == 4) {
            j = 0;
            ++i;
        }
        else
            ++j;
    }
    aa = vect[0][0], ac = vect[0][1], ag = vect[0][2];
    at = vect[0][3], ab = vect[0][4];
    cc = vect[1][1], cg = vect[1][2], ct = vect[1][3], cb = vect[1][4];
    gg = vect[2][2], gt = vect[2][3], gb = vect[2][4];
    tt = vect[3][3], tb = vect[3][4];
    bb = vect[4][4];
}

//  Input:  Seven numbers
//  Output: Maximum integer
//  Used in the scoring function
int max(int a, int b, int c, int d, int e, int f, int g) {
    if (a >= b && a >= c && a >= d && a >= e && a >= f && a >= g) {
        return a;
    }
    else if (b >= a && b >= c && b >= d && b >= e && b >= f && b >= g) {
        return b;
    }
    else if (c >= a && c >= b && c >= d && c >= e && c >= f && c >= g) {
        return c;
    }
    else if (d >= a && d >= b && d >= c && d >= e && d >= f && d >= g) {
        return d;
    }
    else if (e >= a && e >= b && e >= c && e >= d && e >= f && e >= g) {
        return e;
    }
    else if (f >= a && f >= b && f >= c && f >= d && f >= e && f >= g) {
        return f;
    }
    else {
        return g;
    }
}

//  Input: File Name
//  Output: String of nucleotides in the file
string fileSequences(string fileName) {
    string input, seq = "";
    ifstream fin(fileName.c_str());
    while(fin >> input) {
        seq += input;
    }
    fin.close();
    return seq;
}

//  Input:  Table Entries: diag, x, y, z, xy, xz, yz
// //  Output: String of direction of pointer
// string maxPointer(TableEntry xyz, TableEntry x, TableEntry y,
//     TableEntry z, TableEntry xy, TableEntry xz, TableEntry yz) {
string maxPointer(int xyz, int x, int y, int z, int xy, int xz, int yz){
    //int max2 = max(xyz.val, x.val, y.val, z.val, xy.val, xz.val, yz.val);
    //if (max2 == xyz.val)
    int max2 = max(xyz,x,y,z,xy,xz,yz);
    if (max2 == xyz)
        return "xyz";
    //else if (max2 == x.val)
    else if (max2 == x)
        return "x";
    //else if (max2 == y.val)
    else if (max2 == y)
        return "y";
    //else if (max2 == z.val)
    else if (max2 == z)
        return "z";
    //else if (max2 == xy.val)
    else if (max2 == xy)
        return "xy";
    //else if (max2 == xz.val)
    else if (max2 == xz)
        return "xz";
    //else {
    else {
        return "yz";
    }
}

void initTable(vector< vector< vector<TableEntry> > > &scores,
    string seq1, string seq2, string seq3) {
    
    unsigned i, j, k;
    //  Add 1 to sizes to include dummy row/column
    unsigned a = seq1.size()+1, b = seq2.size()+1, c = seq3.size()+1;
    
    //  Create the 3D matrix
    scores.resize(a);
    for (i = 0; i < a; ++i) {
        scores.at(i).resize(b);
    }
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            scores.at(i).at(j).resize(c);
        }
    }
    
    //  Initialize dummy rows and columns
    int dummyScores = 0;
    scores.at(0).at(0).at(0).val = 0;
    scores.at(0).at(0).at(0).pt = "";
    for (i = 0; i < a; ++i) {
        scores.at(i).at(0).at(0).val = 0 + dummyScores;
        scores.at(i).at(0).at(0).pt = "z";
        dummyScores += 2*indel;
        for (j = 1; j < b; ++j) {
            scores.at(i).at(j).at(0).val = 
                scores.at(i).at(j-1).at(0).val + 2 * indel;
            scores.at(i).at(j).at(0).pt = "x";
        }
        for (k = 1; k < c; ++k) {
            scores.at(i).at(0).at(k).val = 
                scores.at(i).at(0).at(k-1).val + 2*indel;
            scores.at(i).at(0).at(k).pt = "y";
        }
    }
    
    //  These next 3 outer loops were originally just adding 
    //      indel instead of 2*indel. Checking if it gives 
    //      correct answer
    for (j = 1; j < b; ++j) {
        for (k = 1; k < c; ++k) {
            scores.at(0).at(j).at(k).val = 
                scores.at(0).at(j-1).at(k-1).val + 2*indel;
            scores.at(0).at(j).at(k).pt = "xy";
        }
    }
    
    for (i = 1; i < a; ++i) {
        for (j = 1; j < b; ++j) {
            scores.at(i).at(j).at(0).val =
                scores.at(i-1).at(j-1).at(0).val + 2*indel;
            scores.at(i).at(j).at(0).pt = "xz";
        }
    }
    
    for (i = 1; i < a; ++i) {
        for (k = 1; k < c; ++k) {
            scores.at(i).at(0).at(k).val =
                scores.at(i-1).at(0).at(k-1).val + 2*indel;
            scores.at(i).at(0).at(k).pt = "xy";
        }
    } 
    
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            for (k = 0; k < c; ++k) {
                if (i == 0) {
                    scores.at(i).at(j).at(k).s1 = " ";
                }
                else {
                    scores.at(i).at(j).at(k).s1 = seq1.at(i-1);
                }
                if (j == 0) {
                    scores.at(i).at(j).at(k).s2 = " ";
                }
                else {
                    scores.at(i).at(j).at(k).s2 = seq2.at(j-1);
                }
                if (k == 0) {
                    scores.at(i).at(j).at(k).s3 = " ";
                }
                else {
                    scores.at(i).at(j).at(k).s3 = seq3.at(k-1);
                }
            }
        }
    }
    // cout << "-----------------" << endl;
    // cout << "000: " << scores.at(0).at(0).at(0).s1 << 
    //     scores.at(0).at(0).at(0).s2 << scores.at(0).at(0).at(0).s3
    //     << endl;
    // cout << "001: " << scores.at(0).at(0).at(1).s1 << 
    //     scores.at(0).at(0).at(1).s2 << scores.at(0).at(0).at(1).s3
    //     << endl;
    // cout << "002: " << scores.at(0).at(0).at(2).s1 << 
    //     scores.at(0).at(0).at(2).s2 << scores.at(0).at(0).at(2).s3
    //     << endl;
    // cout << "010: " << scores.at(0).at(1).at(0).s1 << 
    //     scores.at(0).at(1).at(0).s2 << scores.at(0).at(1).at(0).s3
    //     << endl;
    // cout << "011: " << scores.at(0).at(1).at(1).s1 << 
    //     scores.at(0).at(1).at(1).s2 << scores.at(0).at(1).at(1).s3
    //     << endl;
    // cout << "012: " << scores.at(0).at(1).at(2).s1 << 
    //     scores.at(0).at(1).at(2).s2 << scores.at(0).at(1).at(2).s3
    //     << endl;
    // cout << "020: " << scores.at(0).at(2).at(0).s1 << 
    //     scores.at(0).at(2).at(0).s2 << scores.at(0).at(2).at(0).s3
    //     << endl;
    // cout << "021: " << scores.at(0).at(2).at(1).s1 << 
    //     scores.at(0).at(2).at(1).s2 << scores.at(0).at(2).at(1).s3
    //     << endl;
    // cout << "022: " << scores.at(0).at(2).at(2).s1 << 
    //     scores.at(0).at(2).at(2).s2 << scores.at(0).at(2).at(2).s3
    //     << endl;
    // cout << endl;
    // cout << "100: " << scores.at(1).at(0).at(0).s1 << 
    //     scores.at(1).at(0).at(0).s2 << scores.at(1).at(0).at(0).s3
    //     << endl;
    // cout << "101: " << scores.at(1).at(0).at(1).s1 << 
    //     scores.at(1).at(0).at(1).s2 << scores.at(1).at(0).at(1).s3
    //     << endl;
    // cout << "102: " << scores.at(1).at(0).at(2).s1 << 
    //     scores.at(1).at(0).at(2).s2 << scores.at(1).at(0).at(2).s3
    //     << endl;
    // cout << "110: " << scores.at(1).at(1).at(0).s1 << 
    //     scores.at(1).at(1).at(0).s2 << scores.at(1).at(1).at(0).s3
    //     << endl;
    // cout << "111: " << scores.at(1).at(1).at(1).s1 << 
    //     scores.at(1).at(1).at(1).s2 << scores.at(1).at(1).at(1).s3
    //     << endl;
    // cout << "112: " << scores.at(1).at(1).at(2).s1 << 
    //     scores.at(1).at(1).at(2).s2 << scores.at(1).at(1).at(2).s3
    //     << endl;
    // cout << "120: " << scores.at(1).at(2).at(0).s1 << 
    //     scores.at(1).at(2).at(0).s2 << scores.at(1).at(2).at(0).s3
    //     << endl;
    // cout << "121: " << scores.at(1).at(2).at(1).s1 << 
    //     scores.at(1).at(2).at(1).s2 << scores.at(1).at(2).at(1).s3
    //     << endl;
    // cout << "122: " << scores.at(1).at(2).at(2).s1 << 
    //     scores.at(1).at(2).at(2).s2 << scores.at(1).at(2).at(2).s3
    //     << endl;
    // cout << "-------------------------" << endl;
    
    
    
}

void spaceSavingDP(vector< vector< vector<int> > > arr, 
    string seq1, string seq2, string seq3) {
    
    
    
    
    
}

//  Inputs: 3D array of scores, 3 sequences, 3 strings for results
void DP(vector< vector< vector<TableEntry> > > &scores, string seq1, 
    string seq2, string seq3, string &res1, string &res2, string &res3) {
    
    unsigned i, j, k;
    //  Add 1 to sizes to include dummy row/column
    unsigned a = seq1.size()+1, b = seq2.size()+1, c = seq3.size()+1;
    //  All possibilities of reccurence relation
    int xyz, x, y, z, xy, xz, yz;    
    
    //  Recurrence Relation and filling the Table
    for (i = 1; i < a; ++i) {
        for (j = 1; j < b; ++j) {
            for (k = 1; k < c; ++k) {
                
                //  Diagonal: Add scores of all possible seq pairs (3)
                xyz = scores.at(i-1).at(j-1).at(k-1).val;
                if (scores.at(i).at(j).at(k).s1 == 
                    scores.at(i).at(j).at(k).s2) {
                    xyz += match;
                }
                else {
                    xyz += mismatch2;
                }
                if (scores.at(i).at(j).at(k).s1 ==
                    scores.at(i).at(j).at(k).s3) {
                    xyz += match;
                }
                else {
                    xyz += mismatch2;
                }
                if (scores.at(i).at(j).at(k).s2 == 
                    scores.at(i).at(j).at(k).s3) {
                    xyz += match;
                }
                else {
                    xyz += mismatch2;
                }
                //  ----------------------------------------
                
                x = scores.at(i-1).at(j).at(k).val + 2*indel;
                y = scores.at(i).at(j-1).at(k).val + 2*indel;
                z = scores.at(i).at(j).at(k-1).val + 2*indel;
                
                //  Changed the next three from indel to 2*indel
                xy = scores.at(i-1).at(j-1).at(k).val + 2*indel;
                //  Compare seq1 with seq2
                if (scores.at(i).at(j).at(k).s1 == 
                    scores.at(i).at(j).at(k).s2) {
                    xy += match;
                }
                else {
                    xy += mismatch2;
                }
                //  Compare seq1 with seq3
                xz = scores.at(i-1).at(j).at(k-1).val + 2*indel;
                if (scores.at(i).at(j).at(k).s1 == 
                    scores.at(i).at(j).at(k).s3) {
                    xz += match;
                }
                else {
                    xz += mismatch2;
                }
                //  Compare seq2 with seq3
                yz = scores.at(i).at(j-1).at(k-1).val + 2*indel;
                if (scores.at(i).at(j).at(k).s2 == 
                    scores.at(i).at(j).at(k).s3) {
                    yz += match;
                }
                else {
                    yz += mismatch2;
                }
                scores.at(i).at(j).at(k).val = max(xyz, x, y, z, xy,
                    xz, yz);
                // scores.at(i).at(j).at(k).pt = maxPointer(
                //     scores.at(i-1).at(j-1).at(k-1), 
                //     scores.at(i-1).at(j).at(k), scores.at(i).at(j-1).at(k),
                //     scores.at(i).at(j).at(k-1),
                //     scores.at(i-1).at(j-1).at(k),
                //     scores.at(i-1).at(j).at(k-1),
                //     scores.at(i).at(j-1).at(k-1));
                scores.at(i).at(j).at(k).pt = maxPointer(
                    xyz, 
                    x, y,
                    z,
                    xy,
                    xz,
                    yz);
            }        
        }
        
    }
    
    for (i = 0; i < a; ++i) {
        for (j = 0; j < b; ++j) {
            for (k = 0; k < c; ++k) {
                cout << scores.at(i).at(j).at(k).val << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    
    //  Traceback
    i = a-1, j = b-1, k = c-1;
    while (i > 0 || j > 0 || k > 0) {
        if (i == 0 && j == 0) {
            res1 = "-" + res1;
            res2 = "-" + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --k;
        }
        else if (i == 0 && k == 0) {
            res1 = "-" + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = "-" + res3;
            --j;
        }
        else if (j == 0 && k == 0) {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = "-" + res2;
            res3 = "-" + res3;
            --i;
        }
        else if (i == 0) {
            res1 = "-" + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --j; --k;
        }
        else if (j == 0) {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = "-" + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --i; --k;
        }
        else if (k == 0) {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = "-" + res3;
            --i; --j;
        }
        else if (scores.at(i).at(j).at(k).pt == "xyz") {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --i; --j; --k;
        }
        else if (scores.at(i).at(j).at(k).pt == "xy") {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = "-" + res3;
            --i; --j;
        }
        else if (scores.at(i).at(j).at(k).pt == "xz") {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = "-" + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --i; --k;
        }
        else if (scores.at(i).at(j).at(k).pt == "yz") {
            res1 = "-" + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --j; --k;
        }
        else if (scores.at(i).at(j).at(k).pt == "x") {
            res1 = scores.at(i).at(j).at(k).s1 + res1;
            res2 = "-" + res2;
            res3 = "-" + res3;
            --i;
        }
        else if (scores.at(i).at(j).at(k).pt == "y") {
            res1 = "-" + res1;
            res2 = scores.at(i).at(j).at(k).s2 + res2;
            res3 = "-" + res3;
            --j;
        }
        else if (scores.at(i).at(j).at(k).pt == "z") {
            res1 = "-" + res1;
            res2 = "-" + res2;
            res3 = scores.at(i).at(j).at(k).s3 + res3;
            --k;
        }
    }
}

int main() {
    //  Testing max  ---- WORKS
    // cout << max(41,60,4, 5, 10, 30, 40) << endl;
    
    int aa, ac, ag, at, cc, cg, ct, gg, gt, tt, ab, cb, gb, tb,bb;
    
    //  Testing fileSequences --- WORKS
    string seq1 = "", seq2 = "", seq3 = "";
    seq1 = fileSequences("seq1.txt");
    cout << seq1 << "\tz-axis" << endl;
    seq2 = fileSequences("seq2.txt");
    cout << seq2 << "\tx-axis" << endl;
    seq3 = fileSequences("seq3.txt");
    cout << seq3 << "\ty-axis" << endl;
    
    GetScore("scoreMat.txt", aa, ac, ag, at, cc, cg, ct, gg, gt,
        tt, ab, cb, gb, tb, bb);
    
    vector< vector< vector<TableEntry> > > A;
    string result1 = "", result2 = "", result3 = "";    
    
    initTable(A, seq1, seq2, seq3);
    // cout << "--------------" << endl;
    // for (unsigned i = 0; i <= seq1.size(); ++i) {
    //     for (unsigned j = 0; j <= seq2.size(); ++j) {
    //         for (unsigned k = 0; k <= seq3.size(); ++k) {
    //             cout << A.at(i).at(j).at(k).s1 <<
    //                 A.at(i).at(j).at(k).s2 << A.at(i).at(j).at(k).s3
    //                 << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    // cout << "-----------------" << endl;
    
    DP(A, seq1, seq2, seq3, result1, result2, result3);
    
    cout << result1 << endl << result2 << endl << result3 << endl;
    
    cout << "\nDon't forget to add functionality for applying the "
        << "scores\n";
    
    return 0;
}

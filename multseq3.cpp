#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

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

string ReverseString(string abc) {
    string str = "";
    for (int i = abc.size() - 1; i >= 0; --i) {
        str += abc.at(i);
    }
    return str;
}

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
//  Output: String of direction of pointer
string maxPointer(int xyz, int x, int y, int z, int xy, int xz, int yz){
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

//	Used for normal dynamic programming DP
void initTable(vector< vector< vector<TableEntry> > > &scores,
    string seq1, string seq2, string seq3) {
    
    unsigned i, j, k;
    //  Add 1 to sizes to include dummy row/column
    unsigned a = seq1.size()+1, b = seq2.size()+1, c = seq3.size()+1;
    
	//	Error Check
	//if (a == 1 || b == 1 || c == 1) {
	//	cout << "One sequence is empty" << endl;
	//	cout << "Seq1 size: " << a << " - 1 = " << a - 1 << endl;
	//	cout << "Seq2 size: " << b << " - 1 = " << b - 1 << endl;
	//	cout << "Seq3 size: " << c << " - 1 = " << c - 1 << endl;
	//	//return;
	//}

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
        scores.at(i).at(0).at(0).pt = "x";
        dummyScores += 2*indel;
        for (j = 1; j < b; ++j) {
            scores.at(i).at(j).at(0).val = 
                scores.at(i).at(j-1).at(0).val + 2 * indel;
            scores.at(i).at(j).at(0).pt = "y";
        }
        for (k = 1; k < c; ++k) {
            scores.at(i).at(0).at(k).val = 
                scores.at(i).at(0).at(k-1).val + 2*indel;
            scores.at(i).at(0).at(k).pt = "z";
        }
    }
	//	Error Check
	//if (a == 1 || b == 1 || c == 1)
	//	cout << "InitTable: Survived the dummy column initializations" << endl;
 
    //  These next 3 outer loops were originally just adding 
    //      indel instead of 2*indel. Checking if it gives 
    //      correct answer
    for (j = 1; j < b; ++j) {
        for (k = 1; k < c; ++k) {
            scores.at(0).at(j).at(k).val = 
                scores.at(0).at(j-1).at(k-1).val + 2*indel;
            scores.at(0).at(j).at(k).pt = "yz";
        }
    }
    
    for (i = 1; i < a; ++i) {
        for (j = 1; j < b; ++j) {
            scores.at(i).at(j).at(0).val =
                scores.at(i-1).at(j-1).at(0).val + 2*indel;
            scores.at(i).at(j).at(0).pt = "xy";
        }
    }
    
    for (i = 1; i < a; ++i) {
        for (k = 1; k < c; ++k) {
            scores.at(i).at(0).at(k).val =
                scores.at(i-1).at(0).at(k-1).val + 2*indel;
            scores.at(i).at(0).at(k).pt = "xz";
        }
    } 
    

	//	Error Checks
//	if (a == 1 || b == 1 || c == 1)
//		cout << "InitTable: Survived iniatilizeing mutliple dummies at once" << endl;

    //  Apply strings to the table entries
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

	//	Error Check
//	if (a == 1 || b == 1 || c == 1) {
//		cout << "reached the end of the initTable function and nothing crashed" << endl;
//	}
}

//  Inputs: 3D array of scores, 3 sequences, 3 strings for results
void DP(vector< vector< vector<TableEntry> > > scores, 
    string &res1, string &res2, string &res3) {
    
    unsigned i, j, k;
    //  Add 1 to sizes to include dummy row/column
    	unsigned a = scores.size();
	if (a == 0)
		cout << "Beginning of DP: Something happened in a\n";
	unsigned b = scores.at(0).size();
	if (b == 0)
		cout << "Beginning of DP: Something happened in b\n";
        unsigned c = scores.at(0).at(0).size();
	if (c == 0)
		cout << "Beginning of DP: Something happened in c\n";

	//	Error Check
//	if (a == 1 || b == 1 || c == 1) {
//		cout << "Beginning DP with a degenerate string" << endl;
//	}	

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
    
//	if (a == 1 || b == 1 || c == 1)
//		cout << "Finished the recurrence relation and now beginning traceback with a degenerate string" << endl;
    //  Just outputting to test
    // for (i = 0; i < a; ++i) {
    //     for (j = 0; j < b; ++j) {
    //         for (k = 0; k < c; ++k) {
    //             cout << scores.at(i).at(j).at(k).val << "\t";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    
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
    //cout << res1 << endl << res2 << endl << res3 << endl;
    //cout << "Leaving DP" << endl;
	//	Error Check
//	if (a == 1 || b == 1 || c == 1) {
//		cout << "Made it to the end of the DP function with no errors" << endl;
//	}
}

//  Create the 2D matrices for space saving. Returns the final matrix
//      that will be used to find the midpoint of the traceback line
void Create2DMat(vector< vector<TableEntry> > &vect, 
    vector< vector<TableEntry> >&vect2, int rowsize,
    int colsize, string seq1, string seq2, string seq3) {
    
    unsigned i,j,k;
    int m = seq1.size() + 1;
    
    vect.resize(rowsize);
    vect2.resize(rowsize);
    for (i = 0; i < rowsize; ++i) {
        vect.at(i).resize(colsize);
        vect2.at(i).resize(colsize);
    }
    vect.at(0).at(0).val = 0;
    vect.at(0).at(0).pt = ""; vect.at(0).at(0).s1 = " ";
    vect.at(0).at(0).s2 = " "; vect.at(0).at(0).s3 = " ";
    int dummyScore = 2 * indel;
    for (j = 1; j < rowsize; ++j) {
        vect.at(j).at(0).val = dummyScore;
        vect.at(j).at(0).pt = "x";
        dummyScore += 2 * indel;
    }
    dummyScore = 2*indel;
    for (k = 1; k < colsize; ++k) {
        vect.at(0).at(k).val = dummyScore;
        vect.at(0).at(k).pt = "y";
        dummyScore += 2 * indel;
    }
    for (j = 1; j < rowsize; ++j) {
        for (k = 1; k < colsize; ++k) {
            vect.at(j).at(k).val = vect.at(j-1).at(k-1).val + 
                2 * indel;
            vect.at(j).at(k).pt = "xy";
        }
    }
    
    int xyz, x, y, z, xy, xz, yz;
    for (i = 1; i < m; ++i) {
        vect2.at(0).at(0).val = 2*indel;
        vect2.at(0).at(0).pt = "x";
        for (j = 1; j < rowsize; ++j) {
            vect2.at(j).at(0).val = vect.at(j-1).at(0).val + 2*indel;
        }
        for (k = 1; k < colsize; ++k) {
            vect2.at(0).at(k).val = vect.at(0).at(k-1).val + 2*indel;
        }
        for (j = 1; j < rowsize; ++j) {
            for (k = 1; k < colsize; ++k) {
                xyz = vect.at(j-1).at(k-1).val;
                if (seq1.at(i-1) == seq2.at(j-1))
                    xyz += match;
                else
                    xyz += mismatch2;
                if (seq1.at(i-1) == seq3.at(k-1))
                    xyz += match;
                else
                    xyz += mismatch2;
                if (seq2.at(j-1) == seq3.at(k-1))
                    xyz += match;
                else
                    xyz += mismatch2;
                
                x = vect.at(j).at(k).val + 2*indel;
                y = vect2.at(j-1).at(k).val + 2*indel;
                z = vect2.at(j).at(k-1).val + 2*indel;
                
                //  Compare seq1 with seq2
                xy = vect.at(j-1).at(k).val + 2*indel;
                if (seq1.at(i-1) == seq2.at(j-1))
                    xy += match;
                else
                    xy += mismatch2;
                //  Compare seq1 with seq3
                xz = vect.at(j).at(k-1).val + 2*indel;
                if (seq1.at(i-1) == seq3.at(k-1))
                    xz += match;
                else
                    xz += mismatch2;
                //  Compare seq2 with seq3
                yz = vect2.at(j-1).at(k-1).val + 2*indel;
                if (seq2.at(j-1) == seq3.at(k-1))
                    yz += match;
                else
                    yz += mismatch2;
                vect2.at(j).at(k).val = max(xyz, x, y, z, xy, xz, yz);
                vect2.at(j).at(k).pt = maxPointer(xyz, x, y, z, xy,xz,yz);
            }
        }
        vect = vect2;
    }
}

//  We shall slice the z-axis
void spaceSavingDP( string seq1, string seq2, string seq3,
    string &res1, string &res2, string &res3) {

    unsigned i, j, k,m;
    unsigned a = seq1.size() + 1, b = seq2.size() + 1;
    unsigned c = seq3.size() + 1;
    
    //cout << seq1 << endl << seq2 << endl << seq3 << endl;
    
    m = a/2 + 1;
    //cout << "m = " << m << endl;
    
    //  Base Case:
    //      When a sequence has 0 or 1 characters
    //      Just compute normal Dynamic Programming
    if (m == 0 || m == 1 || a <= 3 || b <= 3 || c <= 3) {
        //cout << "Something degenerated" << endl;
        string r1, r2, r3;
        vector< vector< vector<TableEntry> > > scores;
        initTable(scores, seq1, seq2, seq3);
        DP(scores, r1, r2, r3);
        res1 = res1 + r1;
        res2 = r2 + res2;
        res3 = r3 + res3;
        return;
    }
    
    //  Recursive Step
    //      oh is the vector for the "other half"
    //      fhseq1 is the first half sequence
    string fhseq1 = "", ohseq1 = "", ohseq2 = "", ohseq3 = "";
    ohseq2 = seq2;
    ohseq3 = seq3;
    for (i = 1; i < m; ++i) {
        fhseq1 += seq1.at(i-1);
    }
    for (i = a - 1; i >= m; --i) {
        ohseq1 += seq1.at(i-1);
    }
    ohseq2 = ReverseString(seq2);
    ohseq3 = ReverseString(seq3);
    //  Create the tables of the two halves
    // cout << fhseq1 << endl;
    // cout << ohseq1 << endl << ohseq2 << endl << ohseq3 << endl;
    
    vector< vector<TableEntry> > fh1, fh2;
    vector< vector<TableEntry> > oh1, oh2;
    
    Create2DMat(fh1, fh2, b, c,fhseq1, seq2, seq3);
    Create2DMat(oh1, oh2, b, c,ohseq1, ohseq2, ohseq3);
    
    int max = -1000000;
    int value;
    unsigned maxj, maxk;
    for (j = 0; j < b; ++j) {
        for (k = 0; k < c; ++k) {
            value = fh2.at(j).at(k).val + oh2.at(b-1-j).at(c-1-k).val;
            if (value > max) {
                max = value;
                maxj = j, maxk = k;
            }
        }
    }
    
    // cout << "Final Max: " << max << endl;
    // cout << "Maxj: " << maxj << endl;
    // cout << "Maxk: " << maxk << endl;
    
    //  FH
    //      Cut Half the matrix along i away (fhseq1)
    //      Cut 1 to maxj along j
    //      Cut 1 to maxk along k
    string fhseq2 = "", fhseq3 = "";
    for (i = 0; i < maxj; ++i) {
        fhseq2 += seq2.at(i);
    }
    for (i = 0; i < maxk; ++i) {
        fhseq3 += seq3.at(i);
    }
    //cout << "fhseq2: " << fhseq2 << endl << "fhseq3: " << fhseq3 << endl;
    
    //  OH
    //      seq1 m:a
    //      seq2 maxj:b-2
    //      seq3 maxk:c-2
    string ohseqalt2 = "", ohseqalt3 = "";
    // -----------------------------------------------------------
    // cout << "b: " << b << endl;
    // cout << "maxj: " << maxj << endl;
    // cout << seq1 << endl;
    // cout << seq2 << endl;
    for (int d = b - 2; d >= maxj; --d) {
        ohseqalt2 += seq2.at(d);
    }
    
    for (int d = c - 2; d >= maxk; --d) {
        ohseqalt3 += seq3.at(d);
    }
    // cout << "ohseq2: " << ohseqalt2 << endl;
    // cout << "ohseq3: " << ohseqalt3 << endl;
    
    //cout << "-------------RECURSION 1--------------" << endl;
    string str1 = "", str2 = "", str3 = "", str4 = "", str5 = "", str6= "";
    
    spaceSavingDP(fhseq1, fhseq2, fhseq3, str1, str2, str3);
    
    //cout << "-------------RECURSION 2--------------" << endl;
    spaceSavingDP(ohseq1, ohseqalt2, ohseqalt3, str4, str5, str6);
    
    string str7, str8, str9;
    str7 = ReverseString(str4), str8 = ReverseString(str5);
    str9 = ReverseString(str6);
    
    res1 = str1 + res1 + str7;
    res2 = str2 + res2 + str8;
    res3 = str3 + res3 + str9;
    
}

int main() {
    
    //  For the Scores of the Scoring Matrix
    int aa, ac, ag, at, cc, cg, ct, gg, gt, tt, ab, cb, gb, tb,bb;
    int start_time = time(0);
    unsigned i;

	int setnum;
	cout << "Enter set number: ";   
 	cin >> setnum;

    //  Testing fileSequences --- WORKS
    string seq1 = "", seq2 = "", seq3 = "";
	while (setnum < 1 || setnum > 4) {
		cout << "Invalid Number. Try again: ";
		cin >> setnum;
	}
	if (setnum == 1) {
		cout << "Using set 1" << endl;
		seq1 = fileSequences("OneOne.txt");
		seq2 = fileSequences("OneTwo.txt");
		seq3 = fileSequences("OneThree.txt");
	}
	else if (setnum == 2) {
		cout << "Using set 2" << endl;
		seq1 = fileSequences("TwoOne.txt");
		seq2 = fileSequences("TwoTwo.txt");
		seq3 = fileSequences("TwoThree.txt");
	}
	else if (setnum == 3) {
		cout << "Using set 3" << endl;
		seq1 = fileSequences("ThreeOne.txt");
		seq2 = fileSequences("ThreeTwo.txt");
		seq3 = fileSequences("ThreeThree.txt");
	}
	else if (setnum == 4) {
		cout << "Using set 4" << endl;
		seq1 = fileSequences("FourOne.txt");
		seq2 = fileSequences("FourTwo.txt");
		seq3 = fileSequences("FourThree.txt");
	}

	/*
	seq1 = fileSequences("seq1.txt");
	cout << seq1 << endl;
	seq2 = fileSequences("seq2.txt");
	cout << seq2 << endl;
	seq3 = fileSequences("seq3.txt");
	cout << seq3 << endl;
	*/
    
    GetScore("scoreMat.txt", aa, ac, ag, at, cc, cg, ct, gg, gt,
        tt, ab, cb, gb, tb, bb);
    
    vector< vector< vector<TableEntry> > > A;
    string result1 = "", result2 = "", result3 = "";    
    
    //initTable(A, seq1, seq2, seq3);
    //DP(A, result1, result2, result3);
    
    spaceSavingDP(seq1, seq2, seq3, result1, result2, result3);
    
    cout << "\n\n\nSeqeuence1\n\n\n" << result1 << endl << "\n\n\nSequence2\n\n\n" << result2 << endl << "\n\n\nSequence3\n\n\n" << result3 << endl;
    int length = result1.size();
    int counter = 0, score = 0;
    cout << "Length of alignment: " << length << endl;
   	cout << result2.size() << endl;
	cout << result3.size() << endl; 
    
    for (i = 0; i < length; ++i) {
        //	AAA
	if (result1.at(i) == result2.at(i) && 
            result2.at(i) == result3.at(i)) {
            ++counter;
            score += 3*match;
        }
	//	AAG
        else if (result1.at(i) == result2.at(i) && 
            result1.at(i) != result3.at(i) && result1.at(i) != '-' &&
            result3.at(i) != '-') {
            score += match + 2*mismatch2; 
        }
	//	AGA
        else if (result1.at(i) == result3.at(i) && 
            result1.at(i) != result2.at(i) && result1.at(i) != '-' &&
            result2.at(i) != '-') {
            score += match + 2*mismatch2; 
        }
	//	GAA
        else if (result2.at(i) == result3.at(i) && 
            result1.at(i) != result3.at(i) && result1.at(i) != '-' &&
            result3.at(i) != '-') {
            score += match + 2*mismatch2; 
        }
        
        //	AGC
        else if (result1.at(i) != result2.at(i) && result3.at(i) != '-' &&
            result2.at(i) != result3.at(i) && result2.at(i) != '-' &&
            result3.at(i) != result1.at(i) && result1.at(i) != '-') {
            score += 3*mismatch2;
        }
        
	//	AA-
        else if (result1.at(i) == result2.at(i) && 
            result2.at(i) != result3.at(i) && result2.at(i) != '-' &&
            result3.at(i) == '-') {
            score += match + 2*indel;
        }
	//	A-A
        else if (result1.at(i) == result3.at(i) && 
            result2.at(i) != result3.at(i) && result3.at(i) != '-' &&
            result2.at(i) == '-') {
            score += match + 2*indel;
        }
	//	-AA
        else if (result2.at(i) == result3.at(i) && 
            result1.at(i) != result3.at(i) && result1.at(i) == '-' &&
            result3.at(i) != '-') {
            score += match + 2*indel;
        }
       	
	//	--A
        else if (result1.at(i) == result2.at(i) && 
            result1.at(i) != result3.at(i) && result1.at(i) == '-') {
            score += 2*indel;
        }
	//	-A-
        else if (result1.at(i) == result3.at(i) && 
            result1.at(i) != result2.at(i) && result1.at(i) == '-') {
            score += 2*indel;
        }
	//	A--
        else if (result2.at(i) == result3.at(i) && 
            result1.at(i) != result3.at(i) && result3.at(i) == '-') {
            score += 2*indel;
        }

	//	-AC
        else if (result1.at(i) == '-' && result2.at(i) != result3.at(i) && result2.at(i) != result1.at(i) && result1.at(i) != result3.at(i)) {
		score += 2*indel + mismatch2;
	}
	//	A-C
	else if (result2.at(i) == '-' && result1.at(i) != result2.at(i) && result1.at(i) != result3.at(i) && result2.at(i) != result3.at(i)) {
		score += 2*indel + mismatch2;
	}
	//	AC-
	else if (result3.at(i) == '-' && result1.at(i) != result2.at(i) && result1.at(i) != result3.at(i) && result2.at(i) != result3.at(i)) {
		score += 2*indel + mismatch2;
	}

        else {
		cout << "------\n";
            cout << "I forgot something" << endl;
		cout << result1.at(i) << endl << result2.at(i) << endl << result3.at(i) << endl;
		cout << "-----\n";
        }
    }
    cout << "Score: " << score << endl;
    cout << "Number of columns perfectly matched: " << counter << endl;
    
	int finish_time = time(0);
	int total_time = finish_time - start_time;
	int numHours = total_time / 3600;
	total_time = total_time % 3600;
	int numMinutes = total_time / 60;
	total_time = total_time % 60;
	int numSeconds = total_time;
	cout << "Runtime\n\t" << numHours << " hour(s), " << numMinutes << " minute(s), " << numSeconds << " second(s)" << endl;


    cout << "\nDon't forget to add functionality for applying the "
        << "scores\n";
    
    return 0;
}

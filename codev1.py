#function def to remove header
import time
def fasta_reader(inputpath):
    with open(inputpath, "r") as fileinput:
        line = fileinput.readline()

        header = ''
        sequence = ''
        while line:
            line = line.rstrip('\n')
            if '>' in line:
                header = line
            else:
                sequence = sequence + line
            line = fileinput.readline()
        return sequence



#function defination for nucleotide physiocal properties
def nucleotide(sequence):
    n_s_t = time.time()
    print("\t\t\t\t\tPhysical Properties")
    seq = str.upper(sequence)
    length = len(seq)
    A = seq.count('A')
    T = seq.count('T')
    G = seq.count('G')
    C = seq.count('C')
    extinction_C = ((A * 15.4) + (C * 7.2) + (G * 11.5) + (T * 8.7)) * 0.9 * 1000
    m_weight = ((A * 313.21) + (T * 304.2) + (G * 329.21) + (C * 289.18)) - 61.96
    if length < 14:
        temp1 = 4 * (G + C) + 2 * (A + T)
        print('\tTheoratical Basic temp. is:%.1f C' % temp1)
    elif length >= 14:
        temp2 = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)
        print('\tTheoratical temp. is:%.1f C' % temp2)

    print('\tLenngth of sequence is %d bP:' % length)
    print('\tAnhydrous Molecular weight:%.1f daltons' % m_weight)
    print('\tNo. of Gaunine in sequence are:%d' % G)
    print('\tNo. of Cytocine in sequence are:%d' % C)
    print('\tNo. of Adnine in sequence are:%d' % A)
    print('\tNo. of Thymine in sequence are:%d' % T)
    print('\tGC content of seq. is:%d' % ((G + C) / length * 100), '%')
    print('\tExtenction Coefficient of sequence is: %d /M /cm' % extinction_C)

    '''variable for each IUPAC aminoacid codon'''
    f1 = 0;
    f2 = 0;
    l1 = 0;
    l2 = 0;
    l3 = 0;
    l4 = 0;
    l5 = 0;
    l6 = 0;
    s1 = 0;
    s2 = 0;
    s3 = 0
    s4 = 0;
    s5 = 0;
    s6 = 0;
    r1 = 0;
    r2 = 0;
    r3 = 0;
    r4 = 0;
    r5 = 0;
    r6 = 0;
    p1 = 0;
    p2 = 0
    p3 = 0;
    p4 = 0;
    v1 = 0;
    v2 = 0;
    v3 = 0;
    v4 = 0;
    t1 = 0;
    t2 = 0;
    t3 = 0;
    t4 = 0;
    a1 = 0
    a2 = 0;
    a3 = 0;
    a4 = 0;
    g1 = 0;
    g2 = 0;
    g3 = 0;
    g4 = 0;
    i1 = 0;
    i2 = 0;
    i3 = 0
    stop1 = 0;
    stop2 = 0;
    stop3 = 0;
    y1 = 0;
    y2 = 0;
    h1 = 0;
    h2 = 0;
    q1 = 0;
    q2 = 0;
    n1 = 0
    n2 = 0;
    k1 = 0;
    k2 = 0;
    d1 = 0;
    d2 = 0;
    e1 = 0;
    e2 = 0;
    c1 = 0;
    c2 = 0;
    w1 = 0;
    m = 0

    x = 0
    size = len(seq)
    while (x < size):
        if ((x + 3) <= size):
            if (seq[x:x + 3] == "TTT"):
                f1 = f1 + 1
            elif (seq[x:x + 3] == "TTC"):
                f2 = f2 + 1
            elif (seq[x:x + 3] == "TTA"):
                l1 = l1 + 1
            elif (seq[x:x + 3] == "TTG"):
                l2 = l2 + 1
            elif (seq[x:x + 3] == "CTT"):
                l3 = l3 + 1
            elif (seq[x:x + 3] == "CTC"):
                l4 = l4 + 1
            elif (seq[x:x + 3] == "CTA"):
                l5 = l5 + 1
            elif (seq[x:x + 3] == "CTG"):
                l6 = l6 + 1
            elif (seq[x:x + 3] == "ATT"):
                i1 = i1 + 1
            elif (seq[x:x + 3] == "ATC"):
                i2 = i2 + 1
            elif (seq[x:x + 3] == "ATA"):
                i3 = i3 + 1
            elif (seq[x:x + 3] == "ATG"):
                m = m + 1
            elif (seq[x:x + 3] == "GTT"):
                v1 = v1 + 1
            elif (seq[x:x + 3] == "GTC"):
                v2 = v2 + 1
            elif (seq[x:x + 3] == "GTA"):
                v3 = v3 + 1
            elif (seq[x:x + 3] == "GTG"):
                v4 = v4 + 1
            elif (seq[x:x + 3] == "TCT"):
                s1 = s1 + 1
            elif (seq[x:x + 3] == "TCC"):
                s2 = s2 + 1
            elif (seq[x:x + 3] == "TCA"):
                s3 = s3 + 1
            elif (seq[x:x + 3] == "TCG"):
                s4 = s4 + 1
            elif (seq[x:x + 3] == "CCT"):
                p1 = p1 + 1
            elif (seq[x:x + 3] == "CCC"):
                p2 = p2 + 1
            elif (seq[x:x + 3] == "CCA"):
                p3 = p3 + 1
            elif (seq[x:x + 3] == "CCG"):
                p4 = p4 + 1
            elif (seq[x:x + 3] == "ACT"):
                t1 = t1 + 1
            elif (seq[x:x + 3] == "ACC"):
                t2 = t2 + 1
            elif (seq[x:x + 3] == "ACA"):
                t3 = t3 + 1
            elif (seq[x:x + 3] == "ACG"):
                t4 = t4 + 1
            elif (seq[x:x + 3] == "GCT"):
                a1 = a1 + 1
            elif (seq[x:x + 3] == "GCC"):
                a2 = a2 + 1
            elif (seq[x:x + 3] == "GCA"):
                a3 = a3 + 1
            elif (seq[x:x + 3] == "GCG"):
                a4 += 1
            elif (seq[x:x + 3] == "TAT"):
                y1 = y1 + 1
            elif (seq[x:x + 3] == "TAC"):
                y2 = y2 + 1
            elif (seq[x:x + 3] == "TAA"):
                stop1 = stop1 + 1
            elif (seq[x:x + 3] == "TAG"):
                stop2 = stop2 + 1
            elif (seq[x:x + 3] == "CAT"):
                h1 = h1 + 1
            elif (seq[x:x + 3] == "CAC"):
                h2 = h2 + 1
            elif (seq[x:x + 3] == "CAA"):
                q1 = q1 + 1
            elif (seq[x:x + 3] == "CAG"):
                q2 = q2 + 1
            elif (seq[x:x + 3] == "AAT"):
                n1 = n1 + 1
            elif (seq[x:x + 3] == "AAC"):
                n2 = n2 + 1
            elif (seq[x:x + 3] == "AAA"):
                k1 = k1 + 1
            elif (seq[x:x + 3] == "AAG"):
                k2 = k2 + 1
            elif (seq[x:x + 3] == "GAT"):
                d1 = d1 + 1
            elif (seq[x:x + 3] == "GAC"):
                d2 = d2 + 1
            elif (seq[x:x + 3] == "GAA"):
                e1 = e1 + 1
            elif (seq[x:x + 3] == "GAG"):
                e2 = e2 + 1
            elif (seq[x:x + 3] == "TGT"):
                c1 = c1 + 1
            elif (seq[x:x + 3] == "TGC"):
                c2 = c2 + 1
            elif (seq[x:x + 3] == "TGA"):
                stop3 = stop3 + 1
            elif (seq[x:x + 3] == "TGG"):
                w1 = w1 + 1
            elif (seq[x:x + 3] == "CGT"):
                r1 = r1 + 1
            elif (seq[x:x + 3] == "CGC"):
                r2 = r2 + 1
            elif (seq[x:x + 3] == "CGA"):
                r3 = r3 + 1
            elif (seq[x:x + 3] == "CGG"):
                r4 = r4 + 1
            elif (seq[x:x + 3] == "AGT"):
                s5 = s5 + 1
            elif (seq[x:x + 3] == "AGC"):
                s6 = s6 + 1
            elif (seq[x:x + 3] == "AGA"):
                r5 = r5 + 1
            elif (seq[x:x + 3] == "AGG"):
                r6 = r6 + 1
            elif (seq[x:x + 3] == "GGT"):
                g1 = g1 + 1
            elif (seq[x:x + 3] == "GGC"):
                g2 = g2 + 1
            elif (seq[x:x + 3] == "GGA"):
                g3 = g3 + 1
            elif (seq[x:x + 3] == "GGG"):
                g4 = g4 + 1
        x = x + 3
    print("\n\t\t\t\t\tGenome Biasness")
    print("\tAla \t GCG \t %d" % a4)
    print("\tAla \t GCA \t %d" % a3)
    print("\tAla \t GCT \t %d" % a1)
    print("\tAla \t GCC \t %d" % a2)
    print("\tCys \t TGT \t %d" % c1)
    print("\tCys \t TGC \t %d" % c2)
    print("\tAsp \t GAT \t %d" % d1)
    print("\tAsp \t GAC \t %d" % d2)
    print("\tGlu \t GAG \t %d" % e2)
    print("\tGlu \t GAA \t %d" % e1)
    print("\tPhe \t TTT \t %d" % f1)
    print("\tPhe \t TTC \t %d" % f2)
    print("\tGly \t GGG \t %d" % g4)
    print("\tGly \t GGA \t %d" % g3)
    print("\tGly \t GGT \t %d" % g1)
    print("\tGly \t GGC \t %d" % g2)
    print("\tHis \t CAT \t %d" % h1)
    print("\tHis \t CAC \t %d" % h2)
    print("\tIle \t ATA \t %d" % i3)
    print("\tIle \t ATT \t %d" % i1)
    print("\tIle \t ATC \t %d" % i2)
    print("\tLys \t AAG \t %d" % k2)
    print("\tLys \t AAA \t %d" % k1)
    print("\tLeu \t TTG \t %d" % l2)
    print("\tLeu \t TTA \t %d" % l1)
    print("\tLeu \t CTG \t %d" % l6)
    print("\tLeu \t CTA \t %d" % l5)
    print("\tLeu \t CTT \t %d" % l3)
    print("\tLeu \t CTC \t %d" % l4)
    print("\tMet \t ATG \t %d" % m)
    print("\tAsn \t AAT \t %d" % n1)
    print("\tAsn \t AAC \t %d" % n2)
    print("\tPro \t CCG \t %d" % p4)
    print("\tPro \t CCA \t %d" % p3)
    print("\tPro \t CCT \t %d" % p1)
    print("\tPro \t CCC \t %d" % p2)
    print("\tGlu \t CAG \t %d" % q2)
    print("\tGlu \t CAA \t %d" % q1)
    print("\tArg \t AGG \t %d" % r6)
    print("\tArg \t AGA \t %d" % r5)
    print("\tArg \t CGG \t %d" % r4)
    print("\tArg \t CGA \t %d" % r3)
    print("\tArg \t CGT \t %d" % r1)
    print("\tArg \t CGC \t %d" % r2)
    print("\tSer \t AGT \t %d" % s5)
    print("\tSer \t AGC \t %d" % s6)
    print("\tSer \t TCG \t %d" % s4)
    print("\tSer \t TCA \t %d" % s3)
    print("\tSer \t TCT \t %d" % s1)
    print("\tSer \t TCC \t %d" % s2)
    print("\tThr \t ACG \t %d" % t4)
    print("\tThr \t ACA \t %d" % t3)
    print("\tThr \t ACT \t %d" % t1)
    print("\tThr \t ACC \t %d" % t2)
    print("\tVal \t GTG \t %d" % v4)
    print("\tVal \t GTA \t %d" % v3)
    print("\tVal \t GTT \t %d" % v1)
    print("\tVal \t GTC \t %d" % v2)
    print("\tTry \t TGG \t %d" % w1)
    print("\tTyr \t TAT \t %d" % y1)
    print("\tTyr \t TAC \t %d" % y2)
    print("\tEnd \t TGA \t %d" % stop3)
    print("\tEnd \t TAG \t %d" % stop2)
    print("\tEnd \t TAA \t %d" % stop1)
    n_e_t = time.time()
    n_t = n_e_t - n_s_t
    return n_t

#get all the files present in that directory
import glob
list_of_files = glob.glob('./*.txt')           # create the list of file
for file_name in list_of_files:
    seq=file_name
    sequence=fasta_reader(seq)    #function call to remove header
    nucleotide(sequence)
file=open(r"C:\Users\aadit\OneDrive\Desktop\Thermotoga petrophila genome.txt", "r")
genome=file.read()

def rev_com(s):
    c=''
    for i in range(0,len(s)):
        if(s[i]=='A'):
            c+="T"
        if(s[i]=='T'):
            c+="A"
        if(s[i]=='C'):
            c+="G"
        if(s[i]=='G'):
            c+="C"
    return (c[::-1])

patterns=set()     # final list of strings returned as answer
freqMap={}         # empty map 
nucleotides={'A',"T",'G','C'}
def suffix(pattern):
    return pattern[1:]

d=1             # hamming distance
n=len(genome)
k=9            # k value of kmer

def hamming_dist(p,q):
    count=0
    for i in range(0,len(p)):
        if p[i]!=q[i]:count+=1
    return count

def neighbors(pattern, d):
    if(d==0):
        return pattern
    if(len(pattern)==1):
        return {'A','C','G','T'}
    neighbourhood=set()
    suffix_neighbor=neighbors(suffix(pattern),d)
    for ele in suffix_neighbor:
        if(hamming_dist(suffix(pattern),ele)<d):
            for x in nucleotides:
                neighbourhood.add(x+ele)
        else:
            neighbourhood.add(pattern[0]+ele)

    return (neighbourhood)

for i in range(0,n-k+1):
    pattern=genome[i:i+k]
    rev_pat=rev_com(pattern)
    nayrc=neighbors(rev_pat,d)
    neighborhood=neighbors(pattern,d)
    for naybur in neighborhood:
        if naybur not in freqMap:
            freqMap[naybur]=1
        else:
            freqMap[naybur]+=1
    
    for revnay in nayrc:
        if revnay not in freqMap:
            freqMap[revnay]=1
        else:
            freqMap[revnay]+=1
    
mf=max(freqMap.values())

for k, v in freqMap.items():
    if(v==mf):
        patterns.add(k)

skew=[0]*len(genome)
for i in range(1,len(genome)):
    ele=genome[i]
    if(ele=='C'):
        skew[i]=skew[i-1]-1
    elif(ele=='G'):
        skew[i]=skew[i-1]+1
    else:
        skew[i]=skew[i-1]

ori=set()

min_val=min(skew)
for i in range(0,len(skew)):
    if(skew[i]==min_val):
       ele=genome[i:i+k]
       ele_neighborhood=neighbors(ele,d)
       for element in ele_neighborhood:
           for pat in patterns:
               if(pat==ele):
                   ori.add(pat)
                   break
            
print(ori)

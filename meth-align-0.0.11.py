# coding=utf-8

# ori_seq = original sequence
# con_seq = converted sequence(all c change to t)
# seq_data = sequencing data

# 输入序列
ori_seq = str(input('Please input the original sequence:\n'))
seq_data = str(input('Please input the sequencing result:\n'))



# 将输入的序列统一转换为大写字母，方便对比
u_ori_seq = ori_seq.upper()
u_seq_data = seq_data.upper()

# 过滤序列
get_ori_seq = []
get_seq_data = []
for ori in u_ori_seq:
    if ori == 'A' or ori == 'C' or ori == 'T' or ori == 'G':
        get_ori_seq .append(ori)
# print (get_ori_seq)
u_ori = ''.join(get_ori_seq[:])
for seq in u_seq_data:
    if seq == 'A' or seq == 'C' or seq == 'T' or seq == 'G':
        get_seq_data .append(seq)
# print (get_ori_seq)
u_seq = ''.join(get_seq_data[:])
# print (u_seq)
  
### 计算“C”的含量
# 原始序列
numc = []
numg = []
numt = []
numa = []
numcg = []
for i in u_ori:
    if i == "C":
        numc.append(i)
    elif i == "G":
        numg.append(i)
    elif i == 'T':
        numt.append(i)
    else:
        numa.append(i)
lenc = len(numc)
lenseq = len(u_ori)
crace = lenc/lenseq *100

#测序序列
s_numc = []
s_numg = []
s_numt = []
s_numa = []
s_numcg = []
for i in u_seq:
    if i == "C":
        s_numc.append(i)
    elif i == "G":
        s_numg.append(i)
    elif i == 'T':
        s_numt.append(i)
    else:
        s_numa.append(i)
s_lenc = len(s_numc)
s_lenseq = len(u_seq)
s_crace = s_lenc/s_lenseq *100

### 计算CpG个数及含量
n = 0
while n <= lenseq:
    cpg = u_ori[n:n+2]
    if cpg == "CG":
        numcg.append(cpg)
    n += 1
# print (numcg)

m = 0
while m <= s_lenseq:
    s_cpg = u_seq[m:m+2]
    if s_cpg == 'CG':
        s_numcg.append(s_cpg)
    m += 1

### 将原始序列中的“C”转换为“T”,CpG位点的C不变

base_of_con = []
for base_con in u_ori:
    base_of_con.append(base_con)
    
base_of_con_seq = []
a = 0
while a < len(u_ori):
    if a == len(u_ori)-1:
        if base_of_con[a] == 'C':
            base_of_con[a] = 'T'
            base_of_con_seq.append(base_of_con[a])
        else:
            base_of_con_seq.append(base_of_con[a])
        
    else:
        if base_of_con[a] == 'C' and base_of_con[a+1] != 'G':
            base_of_con[a] = 'T'
            base_of_con_seq.append(base_of_con[a])
        else:
            base_of_con_seq.append(base_of_con[a])
    a += 1
con_seq = "".join(base_of_con_seq[:])
#print (con_seq)

# 将测序序列单独取出再加入到一个元组中，方便统计碱基数
base_of_seq = []
for base_seq in u_seq:
    base_of_seq.append(base_seq)
#print(base_of_seq)


### simth-waterman alignment
n = len(con_seq)+1
m = len(u_seq)+1
# 建立矩阵，计算分值
# match = 1, mismatch = -1, gap = -2
    
score_matrix = [[0 for a in range(n)] for b in range(m)]
max_score = 0
max_pos = None
for i in range(1, m):
    for j in range(1, n):
        score_matrix[i][j] = max(
            score_matrix[i-1][j-1] + (2 if con_seq[j-1] == u_seq[i-1] else -2),
            score_matrix[i-1][j] - 1,
            score_matrix[i][j-1] - 1,
            0)
        # 记录最大分值位置
        score = score_matrix[i][j]
        if score > max_score:
            max_score = score
            max_pos = (i, j)
# 回溯

aligned_con = []
aligned_seq = []

#x, y = len(u_seq), len(con_seq)
x, y = max_pos 
while x > 0 and y > 0:
        diag = score_matrix[x-1][y-1]
        up = score_matrix[x-1][y]
        left = score_matrix[x][y-1]

        if diag >= up and diag >= left:
            
            if x == 1 and y > 1:
                y -=1
                aligned_con.append(con_seq[y-1])
                aligned_seq.append('-')
                
            elif x > 1 and y == 1:
                x -= 1
                aligned_con.append('-')
                aligned_seq.append(u_seq[x-1])
            else:
                x -= 1
                y -= 1
                aligned_con.append(con_seq[y-1])
                aligned_seq.append(u_seq[x-1])

        
        elif up >= left and up > diag:
                aligned_con.append('-')
                aligned_seq.append(u_seq[x-1])
                x -= 1
        else:
                aligned_con.append(con_seq[y-1])
                aligned_seq.append('-')
                y -= 1


get_aligned_con = ''.join(reversed(aligned_con))
get_aligned_seq = ''.join(reversed(aligned_seq))
 


### 加入匹配符号

match, mismatch, gap = 0, 0, 0
align_str = []
for base1, base2 in zip(aligned_con, aligned_seq):
    if base1 == base2:
        align_str.append('|')
        match += 1
    elif '-' in (base1,base2):
        align_str.append(' ')
        gap += 1
    else:
        align_str.append(':')
        mismatch += 1
get_align_str = ''.join(reversed(align_str[:]))

### 计算转化率,CpG中的C不计算
ori_c = []
seq_c = []
k = 0
l =0

while k < len(u_ori):
    if k  == len(u_ori)-1:
        if u_ori[k] == 'C':
            ori_c.append(u_ori[k])
    else:
        if u_ori[k] == 'C' and u_ori[k+1] != 'G':
            ori_c.append(u_ori[k])
    k += 1
    
while l < len(get_aligned_seq):
    if l == len(get_aligned_seq)-1:
        if get_aligned_seq[l] == 'C':
            seq_c.append(l)
    else:
        if get_aligned_seq[l] == 'C' and get_aligned_seq[l+1] != 'G':
            seq_c.append(l)
    l +=1

if len(ori_c) == 0:
    con_rate = 0
else:
    con_rate = (1-len(seq_c)/len(ori_c))*100

### 查找比对序列中CpG的位置，计算甲基化概率

i = 0
n = 1
c_of_seq = []
t_of_seq = []
m_of_seq = []
cpg_pos = []
meth_data = {}
while i < len(get_aligned_con):
    if get_aligned_con[i] == 'C':
        if get_aligned_seq[i] == 'C':
            c_of_seq.append(get_aligned_seq[i])
        elif get_aligned_seq[i] == 'T':
            t_of_seq.append(get_aligned_seq[i])
        else:
            m_of_seq.append(get_aligned_seq[i])
        cpg_pos.append(i)
        meth_data[n] = get_aligned_seq[i]
        n += 1
    i += 1


### 输出结果
print()
print ('Original sequence length:%d'%len(u_ori))
print ('Sequencing result length:%d'%len(u_seq))
print ('Identitied sequence length:%d'%len(get_aligned_con))
print ()

print ('Number of original sequence CpG:%d '%len(numcg))
print ()

print ('conversion rate:%.2f%%'%con_rate)
print()

print ('Methylatied CpG:%d (%.2f%%)'%(len(c_of_seq),len(c_of_seq)/len(numcg)*100))
print ('Unmethylated CpG:%d (%.2f%%)'%(len(t_of_seq),len(t_of_seq)/len(numcg)*100))
print ('Mutation of CpG:%d'%len(m_of_seq))
print()
print ('Details of methylation site\n ',meth_data)
print()
print ('The position of CpG in sequence\n',cpg_pos)
print()
print ('Matchs:%d/%d (%.2f%%), '%(match,len(get_aligned_con),match/len(get_aligned_con)*100),
       'Mismatchs:%d/%d (%.2f%%), '%(mismatch,len(get_aligned_con),mismatch/len(get_aligned_con)*100),
       'Gaps:%d/%d (%.2f%%)'%(gap,len(get_aligned_con),gap/len(get_aligned_con)*100))
print()
length = len(get_aligned_con)
for i in range(0,length,50):
    seq1_len = len(get_aligned_con[i:i+50])
    print ('Original    %d\t'%(i+1),get_aligned_con[i:i+50],'\t%d'%(i+seq1_len))
    print ('                ',get_align_str[i:i+50])
    seq2_len = len(get_aligned_seq[i:i+50])
    print ('Converted   %d\t'%(i+1),get_aligned_seq[i:i+50],'\t%d'%(i+seq2_len))
    print ()
    









    

import matplotlib.pyplot as plt
import numpy as np

# Read numbers from the file
with open("/Users/dnikzad/Documents/code/cpp/results/gb_0_s.txt", "r") as file:
    data = [float(line.strip()) for line in file]  # Convert lines to float

# Plot the data
plt.plot(data[300:])

# Add labels and title
plt.xlabel("Sample")
plt.ylabel("Score")
plt.title("Traceplot")


# Show the plot
plt.show()

def read_matrix_from_txt(file_path):
    matrix = []
    with open(file_path, 'r') as f:
        # Read each line in the file
        for line in f:
            # Split the line by spaces and convert to float
            row = list(map(float, line.split()))
            matrix.append(row)
    return np.array(matrix)

ref= read_matrix_from_txt("/Users/dnikzad/Documents/code/datasets/asia100sparse3ref.txt")
out= read_matrix_from_txt("/Users/dnikzad/Documents/code/cpp/results/gb_0_m.txt")


print('Max loss: ',np.max(np.abs(ref-out)))
print('Norm-L1 loss: ',np.sum(np.abs(ref-out))/np.sum(ref))
# The following code is commented out as it was part of the original request but not needed for the current task.
# import numpy as np

# # Load the adjacency matrix from the space-separated file
# matrix = np.loadtxt("/Users/dnikzad/Documents/code/relevant/DDAG/matrix.txt")

# # Print to verify
# print(matrix)

# import re
# import numpy as np
# import math
# from scipy.special import logsumexp


# def cache_f (path,max_parents,n_nodes):
# #with open('./data/scores/asia.jkl','r') as file:
#     with open(path,'r') as file:
#         filedata = file.readlines()
   
       
#     #line = filedata[1]
#     #print(line)
#     #print(re.findall(r'\d+', line))

#     n_parents_set=[]
#     cache ={}
   
   
   
#     t=1
#     vertices=[]
#     for i in range (n_nodes):
#         line = filedata[t]
#         vertices.append(int(re.findall(r'\d+', line)[0]))
#         t+=int(re.findall(r'\d+', line)[1])+1
#     #print(vertices)
   
#     t=1
#     for i in vertices :
#         line = filedata[t]
#         now_n_parents_set =int(re.findall(r'\d+', line)[1])
#         #print(now_n_parents_set)
#         n_parents_set.append(now_n_parents_set)
#         now_dict= {}
#         #print(filedata[t])
#         for j in range(now_n_parents_set):
#             now_parents=()
#             for k in range (len(re.findall(r'\d+', filedata[t+j+1]))-3):
#                 #print(re.findall(r'\d+', filedata[t+j+1]))
#                 now_parents= now_parents + (int(re.findall(r'\d+', filedata[t+j+1])[k+3]),)
#             #now_dict[float(filedata[t+j+1].split(' ')[0])]=now_parents
#             now_parents=tuple(sorted(now_parents))
#             #if now_parents==(12,18):
#                 #print(i)
#             now_dict[now_parents]=float(filedata[t+j+1].split(' ')[0])
#         cache[i]=now_dict
#         t=t+now_n_parents_set+1
#     return cache


# bida=cache_f('/Users/dnikzad/Documents/code/datasets/asia100sparse3bida.out',3,8)

# n=8
# post_prob=np.zeros((n,n))

# for i in range (n):
#     for j in range (n):
#         if i!=j:
#             x=[]
#             y=[]
#             keys=bida[i].keys()
#             for k in keys:
#                 #print(j,k)
#                 if j in k:
#                     #print(math.comb(len(cache)-1,len(keys)))
#                     x.append(bida[i][k])
#                 y.append(bida[i][k])
#             post_prob[j,i]=np.exp(logsumexp(x)-logsumexp(y))
# np.round(post_prob,2)

# print(post_prob)

# print(np.max(np.abs(matrix-post_prob)))
# print(np.sum(np.abs(matrix-post_prob))/np.sum(post_prob))
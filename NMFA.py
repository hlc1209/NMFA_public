'''
Deciphering the Generating Rules and Functionalities of Complex Networks
Author X. Xiao, H. Chen, P. Bogdan
Nov. 25 2021
https://www.nature.com/articles/s41598-021-02203-4
'''


from scipy.stats import linregress
import numpy as np
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool
from tqdm import tqdm
import time


    

def is_isomorphic(G1,G2,Q = [q/100 for q in range(-2000,2001,10)]):
    ntauls1 = np.array(nfd_nk(G1,Q))
    ntauls2 = np.array(nfd_nk(G2,Q))
    if np.array_equal(ntauls1, ntauls2):
        return True
    else:
        return False

def node_dimension_single(G,node,weight=None,fig=True):
    grow = []
    r_g = []
    num_g = []
    num_nodes = 0
    if weight == None:
        spl = nx.single_source_shortest_path_length(G,node)
    else:
        spl = nx.single_source_dijkstra_path_length(G,node)
    for s in spl.values():
        if s>0:
            grow.append(s)
    grow.sort()
    num = Counter(grow)
    for i,j in num.items():
        num_nodes += j
        if i>0:
            r_g.append(i)
            num_g.append(num_nodes)
    x = np.log(r_g)
    y = np.log(num_g)
    slope, intercept, r_value, _, _  = linregress(x, y)
    dimension = slope
    if fig:
        plt.plot(x,y,'o',label='fitted data')
        plt.plot(x,intercept + slope*x,'r',label='fitted line')
        plt.title("NFD of "+str(node)+'='+str(slope))
    return dimension,r_value**2
    



def nfd(G,Q = [q/100 for q in range(-2000,2001,10)]):
    
    N_list = []
    rmax_list = []
    for node in tqdm(G.nodes(), total=G.number_of_nodes()):
        grow = []
        r_g_all = []
        num_g_all = []
        num_nodes = 0
        num_nodes_all = 0
        spl = nx.single_source_shortest_path_length(G,node)
        for s in spl.values():
            if s>0:
                grow.append(s)
        grow.sort()
        num = Counter(grow)
        rmax_list.append(grow[-1])
        for k,h in num.items():
            num_nodes_all += h
        for i,j in num.items():
            num_nodes += j
            if i>0:
                r_g_all.append(i)
                num_g_all.append(num_nodes/num_nodes_all)
        N_list.append(num_g_all)
    diameter = np.max(rmax_list)
    
    
    r_max = diameter
    r_min=0
    tau_list=[]
    cnt = 0
    Zq_list_q=[]
    for q in Q:
        Zq_list_q.append([])
    for j in range(r_min,r_max):
        Zq_q = [0]*len(Q)
        for k in range(len(N_list)):
            if j < len(N_list[k]):
                for idx, q in enumerate(Q):
                    Zq_q[idx] += (N_list[k][j])**q
            else:
                for idx, q in enumerate(Q):
                    Zq_q[idx] += 1
        for idx, q in enumerate(Q):
            Zq_list_q[idx].append(Zq_q[idx])
    
    for idx, q in enumerate(Q):
        Zq_list = Zq_list_q[idx]
        x = np.log(range(r_min+1,r_max+1)/diameter)
        y = np.log(Zq_list)
        q = format (q, '.0f' ) 
        slope, intercept, _, _, _  = linregress(x, y)
        tau_list.append(slope)
        cnt += 1
    return tau_list

def wnfd(G,Q = [q/100 for q in range(-2000,2001,10)],fdigi=3):
    N_list = []
    r_g_all = []
    num_nodes_all = nx.number_of_nodes(G)
    for node in G.nodes():
        grow = []    
        spl = nx.single_source_dijkstra_path_length(G,node)
        for s in spl.values():
            if s>0:
                s = round(s,fdigi)
                grow.append(s)
                if s not in r_g_all:
                    r_g_all.append(s)
        grow.sort()
        grow = [round(d,fdigi) for d in grow]
        num = Counter(grow)
        N_list.append(num)
    r_g_all.sort()
    
    r_num = len(r_g_all)
    Nw_mat = np.zeros((num_nodes_all,r_num))
    for column,r in enumerate(r_g_all):
        for row,n_ls in enumerate(N_list):        
            num_r = 1
            for key in n_ls:
                if key <= r:
                    num_r += n_ls[key]
            Nw_mat[row,column] = num_r

    diameter = r_g_all[-1]
    
    Zq_list = [[] for i in range(len(Q))]
    for j in range(len(r_g_all)):
        Zq_q = [0 for i in range(len(Q))]
        for i in range(num_nodes_all):
            for qdx, q in enumerate(Q):
                Zq_q[qdx] += (Nw_mat[i][j]/num_nodes_all)**q 
        for idx,q in enumerate(Q):
            Zq_list[idx].append(Zq_q[idx])
            
    tau_list = []
    for idx, q in enumerate(Q):
        r_g_all_np = np.array(r_g_all)
        x = np.log(r_g_all_np/diameter)
        y = np.log(Zq_list[idx])
        slope, _, _, _, _  = linregress(x, y)
        tau_list.append(slope)
           
    return tau_list
    
def nfd_nk(G,Q = [q/100 for q in range(-2000,2001,10)]):
    import networkit as nk
    N_list = []
    rmax_list = []
    G = nx.convert_node_labels_to_integers(G)
    G_nk = nk.nxadapter.nx2nk(G)

    for node in tqdm(G.nodes(), total=G.number_of_nodes()):
        r_g_all = []
        num_g_all = []
        num_nodes = 0
        num_nodes_all = 0
        grow = nk.distance.BFS(G_nk, node, storePaths=False).run().getDistances()
        grow = list(map(int,grow))
        grow.sort()
        grow = grow[1:]
        num = Counter(grow)
        rmax_list.append(grow[-1])
        for k,h in num.items():
            num_nodes_all += h
        for i,j in num.items():
            num_nodes += j
            if i>0:
                r_g_all.append(i)
                num_g_all.append(num_nodes/num_nodes_all)
        N_list.append(num_g_all)
    diameter = np.max(rmax_list)

    r_max = diameter
    r_min=0
    tau_list=[]
    cnt = 0
    Zq_list_q=[]
    for q in Q:
        Zq_list_q.append([])
    for j in range(r_min,r_max):
        Zq_q = [0]*len(Q)
        for k in range(len(N_list)):
            if j < len(N_list[k]):
                for idx, q in enumerate(Q):
                    Zq_q[idx] += (N_list[k][j])**q
            else:
                for idx, q in enumerate(Q):
                    Zq_q[idx] += 1
        for idx, q in enumerate(Q):
            Zq_list_q[idx].append(Zq_q[idx])
    
    for idx, q in enumerate(Q):
        Zq_list = Zq_list_q[idx]
        x = np.log(range(r_min+1,r_max+1)/diameter)
        y = np.log(Zq_list)
        q = format (q, '.0f' ) 
        slope, intercept, _, _, _  = linregress(x, y)
        tau_list.append(slope)
        cnt += 1
    return tau_list

def wnfd_nk(G,Q = [q/100 for q in range(-2000,2001,10)],fdigi=3):
    import networkit as nk
    N_list = []
    r_g_all_set = set()
    num_nodes_all = nx.number_of_nodes(G)
    G = nx.convert_node_labels_to_integers(G)
    G_nk = nk.nxadapter.nx2nk(G,weightAttr='weight')
    for node in tqdm(G.nodes(), total=num_nodes_all):
        grow = nk.distance.Dijkstra(G_nk, node, storePaths=False).run().getDistances()
        grow.sort()
        grow = [round(d,fdigi) for d in grow]
        grow = grow[1:]
        num = Counter(grow)
        r_g_all_set.update(num.keys())
        N_list.append(num)
    r_g_all = sorted(list(r_g_all_set))
    
    r_num = len(r_g_all)
    Nw_mat = np.zeros((num_nodes_all,r_num))
    for column,r in enumerate(r_g_all):
        for row,n_ls in enumerate(N_list):        
            num_r = 1
            for key in n_ls:
                if key <= r:
                    num_r += n_ls[key]
            Nw_mat[row,column] = num_r

    diameter = r_g_all[-1]
    
    Zq_list = [[] for i in range(len(Q))]
    for j in range(len(r_g_all)):
        Zq_q = [0 for i in range(len(Q))]
        for i in range(num_nodes_all):
            for qdx, q in enumerate(Q):
                Zq_q[qdx] += (Nw_mat[i][j]/num_nodes_all)**q 
        for idx,q in enumerate(Q):
            Zq_list[idx].append(Zq_q[idx])
            
    tau_list = []
    for idx, q in tqdm(enumerate(Q),total=len(Q)):
        r_g_all_np = np.array(r_g_all)
        x = np.log(r_g_all_np/diameter)
        y = np.log(Zq_list[idx])
        slope, _, _, _, _  = linregress(x, y)
        tau_list.append(slope)
           
    return tau_list

def _calcSPL_nx(G,node_list):
    result = {}
    for node in node_list:
        result[node] = nx.single_source_shortest_path_length(G,node)
    return result

def nfd_multiprocess(G,Q = [q/100 for q in range(-2000,2001,10)],cpu_num=os.cpu_count()-1):
    node_list_list = []
    node_list = list(G.nodes())
    num_per_core = G.number_of_nodes() // cpu_num
    num_left = G.number_of_nodes() % cpu_num
    num_temp = 0
    for idx in range(cpu_num):
        if num_left:
            node_list_list.append(node_list[idx*num_per_core+num_temp:(idx+1)*num_per_core+1+num_temp])
            num_left -= 1
            num_temp += 1
        else:
            node_list_list.append(node_list[idx*num_per_core+num_temp:(idx+1)*num_per_core+num_temp])

    pool = Pool(processes=cpu_num)
    multi_result = [pool.apply_async(_calcSPL_nx,(G, node_list)) for node_list in node_list_list]
    pool.close()
    pool.join()

    spl_dict = {}
    for res in multi_result:
        spl_dict.update(res.get())

    N_list = []
    rmax_list = []
    for node in tqdm(G.nodes(), total=G.number_of_nodes()):
        grow = []
        r_g_all = []
        num_g_all = []
        num_nodes = 0
        num_nodes_all = 0
        spl = spl_dict[node]
        for s in spl.values():
            if s>0:
                grow.append(s)
        grow.sort()
        num = Counter(grow)
        rmax_list.append(grow[-1])
        for k,h in num.items():
            num_nodes_all += h
        for i,j in num.items():
            num_nodes += j
            if i>0:
                r_g_all.append(i)
                num_g_all.append(num_nodes/num_nodes_all)
        N_list.append(num_g_all)
    diameter = np.max(rmax_list)


    r_max = diameter
    r_min=0
    tau_list=[]
    cnt = 0
    Zq_list_q=[]
    for q in Q:
        Zq_list_q.append([])
    for j in range(r_min,r_max):
        Zq_q = [0]*len(Q)
        for k in range(len(N_list)):
            if j < len(N_list[k]):
                for idx, q in enumerate(Q):
                    Zq_q[idx] += (N_list[k][j])**q
            else:
                for idx, q in enumerate(Q):
                    Zq_q[idx] += 1
        for idx, q in enumerate(Q):
            Zq_list_q[idx].append(Zq_q[idx])
    
    for idx, q in enumerate(Q):
        Zq_list = Zq_list_q[idx]
        x = np.log(range(r_min+1,r_max+1)/diameter)
        y = np.log(Zq_list)
        q = format (q, '.0f' ) 
        slope, _, _, _, _  = linregress(x, y)
        tau_list.append(slope)
        cnt += 1
    return tau_list

def nspectrum(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],name=None,linewidth=5):
    tau_list = tau_list[170:231]
    q_list = q_list[170:231]
    al_list = []
    fal_list = []
    for i in range(1,len(q_list)):
        al=(tau_list[i]-tau_list[i-1])/(q_list[i]-q_list[i-1])
        al_list.append(al)
    for j in range(len(q_list)-1):
        fal=q_list[j]*al_list[j]-tau_list[j]
        fal_list.append(fal)
    if name is not None:
        plt.plot(al_list,fal_list,label=name,linewidth=linewidth)
        plt.legend()
    else:
        plt.plot(al_list,fal_list,linewidth=linewidth)
    plt.xlabel('Lipschitz-Holder exponent, 'r'$\alpha$')
    plt.ylabel('Multi-fractal spectrum, 'r'$f(\alpha)$')
    alpha_0 = np.argmax(fal_list)
    asymmetry_ratio = np.log((al_list[-1] - al_list[alpha_0])/(al_list[alpha_0]-al_list[0]))
    return al_list,fal_list, asymmetry_ratio


def ndimension(tau_list,q_list=[q/100 for q in range(-2000,2001,10)],linewidth=5,name=None):
    dim_list = []
    qd_list = []
    for i in range(len(q_list)):
        if q_list[i] != 0:
            dim = tau_list[i]/q_list[i]
            dim_list.append(dim)
            qd_list.append(q_list[i])
    if name is not None:
        plt.plot(qd_list,dim_list,label=name,linewidth=linewidth)
        plt.legend()
    else:
        plt.plot(qd_list,dim_list,linewidth=linewidth)
    plt.xlabel('Distortion exponent, 'r'$q$')
    plt.ylabel('Generalized fractal dimension, 'r'$D(q)$')
    return dim_list, qd_list

def distance(ndim_list):
    dn = len(ndim_list)
    dist = np.zeros((dn,dn))
    for i in range(dn):
        for j in range(dn):
            d1 = np.array(ndim_list[i])
            d2 = np.array(ndim_list[j])
            dist[i][j] = np.sqrt(np.sum((d1-d2)**2)/len(ndim_list[i]))
    return dist


def nheat(tau_list, q_list=[q/100 for q in range(-2000,2001,10)]):
    heat_list = []
    for i in range(2,len(tau_list)):
        heat = -100 * (tau_list[i] - 2*tau_list[i-1] + tau_list[i-2])
        heat_list.append(heat)
    return heat_list, q_list[2:]
